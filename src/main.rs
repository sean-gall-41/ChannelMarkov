extern crate rand;

use std::io::Error;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
};

const DEFAULT_NUM_CHANNELS: u32 = 1;
const DEFAULT_NUM_TS: u32       = 1000; // ms
const DEFAULT_STIM_MV: f32      = 10.0;
const DEFAULT_OPEN_RATE: f32    = 1.22;
const DEFAULT_CLOSE_RATE: f32   = 0.056;

/*
 * Here's the basic state diagram that we're assuming our Markov Model follows:
 *
 *  close            close            close            close            open
 * +-----+          +-----+          +-----+          +-----+          +-----+
 * |     |    4a    |     |    3a    |     |    2a    |     |    1a    |     |
 * |  1  |  <---->  |  2  |  <---->  |  3  |  <---->  |  4  |  <---->  |  5  |
 * |     |    1b    |     |    2b    |     |    3b    |     |    4b    |     |
 * +-----+          +-----+          +-----+          +-----+          +-----+
 *
 * where here a == OPEN_RATE and b == CLOSE_RATE
 *
 * one thing to change, after checking with real data, is which transitions are allowed:
 * this linear transition diagram is a simplication, and there may be other diagrams
 * that better fit the data
 */

#[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
#[repr(i32)]
pub enum ChannelState {
    Closed1 = 0,
    Closed2,
    Closed3,
    Closed4,
    Open,
}

impl Distribution<ChannelState> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ChannelState {
        unsafe { ::std::mem::transmute(rng.gen_range(0..=4)) }
    }
}

#[derive(Debug, Clone)]
pub struct Channel {
    pub state: ChannelState,
}

impl Channel {
    fn new() -> Self {
        Channel {
            state: rand::random(),
        }
    }

    fn make_transition(&mut self, transition_matrix: &[[f32; 5]; 5]) {
        let transition_row = transition_matrix[self.state as usize];
        let dist = WeightedIndex::new(&transition_row).unwrap();
        let transition_prob = transition_row[dist.sample(&mut rand::thread_rng())];
        self.state = unsafe { ::std::mem::transmute(transition_row.iter()
                                                                  .position(|&r| r == transition_prob)
                                                                  .unwrap() as i32)
        };
    }
}

/* Takes in a vector of tuples that describe a step-wise voltage stimulus
 * e.g., suppose you have this vector:
 * vec![(20, -100), (20, 10), (20, -100)]
 * this says:
 *
 * set the stimulus to -100mV in interval [0, 20),
 * set the stimulus to 10mV in interval [20, 40),
 * set the stimulus to -100mV in interval [40, 60)
 *
 * for safety, asserts that the sum of the p.0s == num_ts
 */
fn generate_stimulus(stim_intervals: &Vec<(u32, f32)>, num_ts: u32) -> Vec<f32> {
    assert_eq!(num_ts, stim_intervals.iter().map(|&p| p.0).sum());
    let mut stim: Vec<f32> = vec![];
    for (t, v) in stim_intervals.iter() {
        for _ in 0..(*t) {
            stim.push(*v);
        }
    }
    stim
}

pub struct Simulation {
    pub num_channels: u32,
    pub num_ts: u32,
    pub open_rate: f32,
    pub close_rate: f32,
    pub channels: Vec<Channel>,
    pub stimulus: Vec<f32>,
    pub transition_matrix: [[f32; 5]; 5],
}

impl Simulation {
    fn from(num_channels: u32,
            num_ts: u32,
            init_open_rate: f32,
            init_close_rate: f32,
            stim_intervals: &Vec<(u32, f32)>) -> Self {
        Simulation {
            num_channels: num_channels,
            num_ts: num_ts,
            open_rate: init_open_rate,
            close_rate: init_close_rate,
            channels: vec![Channel::new(); num_channels as usize],
            stimulus: generate_stimulus(stim_intervals, num_ts),
            transition_matrix: [
    [            0.0,  4.0 * init_open_rate,                   0.0,                   0.0,            0.0],
    [init_close_rate,                   0.0,  3.0 * init_open_rate,                   0.0,            0.0],
    [            0.0, 2.0 * init_close_rate,                   0.0,  2.0 * init_open_rate,            0.0],
    [            0.0,                   0.0, 3.0 * init_close_rate,                   0.0, init_open_rate],
    [            0.0,                   0.0,                   0.0, 4.0 * init_close_rate,            0.0],
            ],
        }
    }

    fn update_rates(&mut self, ts: u32) {
        // for now: hard-coding in values for delayed K+ rectifier conductance
        self.open_rate = DEFAULT_OPEN_RATE * (0.04 * self.stimulus[ts as usize]).exp();
        self.close_rate = DEFAULT_CLOSE_RATE * (-0.0125 * self.stimulus[ts as usize]).exp();
    }

    fn update_transition_matrix(&mut self) {
            self.transition_matrix = [
    [            0.0,  4.0 * self.open_rate,                   0.0,                   0.0,            0.0],
    [self.close_rate,                   0.0,  3.0 * self.open_rate,                   0.0,            0.0],
    [            0.0, 2.0 * self.close_rate,                   0.0,  2.0 * self.open_rate,            0.0],
    [            0.0,                   0.0, 3.0 * self.close_rate,                   0.0, self.open_rate],
    [            0.0,                   0.0,                   0.0, 4.0 * self.close_rate,            0.0],
            ];

    }

    fn run(&mut self) {
        for ts in 0..self.num_ts {
            self.update_rates(ts);
            self.update_transition_matrix();
            for channel in &mut self.channels {
                channel.make_transition(&self.transition_matrix);
                println!("ts: {ts}, state: {:?}", channel.state);
                //if channel.state == ChannelState::Open {
                //}
            }
        }
    }
}

/*
 * Couple of fixes:
 * c. need to create function to generate stimulus
 * d. even later: put all relevant params in json and read from
 */
fn main() -> Result<(), Error> {
    // initial values for delayed K+ rectifier conductance
    let stimulus: Vec<(u32, f32)> = vec![(15, -100.0), (20, 10.0), (15, -100.0)];
    let mut model_sim = Simulation::from(1,
                                         50,
                                         DEFAULT_OPEN_RATE,
                                         DEFAULT_CLOSE_RATE,
                                         &stimulus);
    model_sim.run();

    Ok(())
}

