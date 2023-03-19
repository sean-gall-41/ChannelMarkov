extern crate rand;

use std::fs;
use std::io::Error;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
};
use rand_distr::{Normal};

const PARAM_FILE_NAME: &str = "params.json";

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
    pub emis: f32
}

impl Channel {
    fn new() -> Self {
        Channel {
            state: rand::random(),
            emis: 0.0f32
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

    fn sample_emission(&mut self, emis_dist: &EmisDist) {
        match self.state {
            ChannelState::Closed1 => {
                self.emis = emis_dist.closed_1_dist.sample(&mut rand::thread_rng());
            },
            ChannelState::Closed2 => {
                self.emis = emis_dist.closed_2_dist.sample(&mut rand::thread_rng());
            },
            ChannelState::Closed3 => {
                self.emis = emis_dist.closed_3_dist.sample(&mut rand::thread_rng());
            },
            ChannelState::Closed4 => {
                self.emis = emis_dist.closed_4_dist.sample(&mut rand::thread_rng());
            },
            ChannelState::Open => {
                self.emis = emis_dist.open_dist.sample(&mut rand::thread_rng());
            },
        }
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
fn generate_stimulus(stim_intervals: &json::JsonValue, num_ts: u32) -> Vec<f32> {
    assert!(stim_intervals.is_array());
    let mut stim: Vec<f32> = vec![];
    for pair in stim_intervals.members() {
        for _ in 0..(pair["step_dt"].as_u32().unwrap()) {
            stim.push(pair["step_v"].as_f32().unwrap());
        }
    }
    stim
}

#[derive(Debug)]
pub struct ChannelRateParams {
    pub open_rate: f32,
    pub open_exp_const: f32,
    pub close_rate: f32,
    pub close_exp_const: f32,
}

impl ChannelRateParams {
    fn from(init_open_rate: f32,
            open_exp_const: f32,
            init_close_rate: f32,
            close_exp_const: f32) -> Self {
       ChannelRateParams {
        open_rate: init_open_rate,
        open_exp_const: open_exp_const,
        close_rate: init_close_rate,
        close_exp_const: close_exp_const
       }
    }
}

fn generate_emis_params(emis_params: &json::JsonValue) -> Vec<(f32, f32)> {
    assert!(emis_params.is_array());
    let emis: Vec<(f32, f32)> = emis_params.members().map(|p| {
        (p["mu"].as_f32().unwrap(), p["sigma"].as_f32().unwrap())
    }).collect();
    emis
}

#[derive(Debug)]
pub struct EmisDist {
    pub closed_1_dist: Normal<f32>,
    pub closed_2_dist: Normal<f32>,
    pub closed_3_dist: Normal<f32>,
    pub closed_4_dist: Normal<f32>,
    pub open_dist: Normal<f32>,
}

impl EmisDist {
    fn from(emis_params: &json::JsonValue) -> Self {
        let emis_params_vec: Vec<(f32, f32)> = generate_emis_params(emis_params);
        assert_eq!(emis_params_vec.len(), 5); // the number of different distributions we need
        EmisDist {
            closed_1_dist: Normal::new(emis_params_vec[0].0, emis_params_vec[0].1).unwrap(),
            closed_2_dist: Normal::new(emis_params_vec[1].0, emis_params_vec[1].1).unwrap(),
            closed_3_dist: Normal::new(emis_params_vec[2].0, emis_params_vec[2].1).unwrap(),
            closed_4_dist: Normal::new(emis_params_vec[3].0, emis_params_vec[3].1).unwrap(),
            open_dist: Normal::new(emis_params_vec[4].0, emis_params_vec[4].1).unwrap()
        }
    }
}

#[derive(Debug)]
pub struct Simulation {
    pub num_channels: u32,
    pub num_ts: u32,
    pub channels: Vec<Channel>,
    pub stimulus: Vec<f32>,
    pub emis_dists: EmisDist,
    pub transition_matrix: [[f32; 5]; 5],
    pub channel_rate_params: ChannelRateParams,
}

impl Simulation {
    fn from(num_channels: u32,
            num_ts: u32,
            channel_rate_params: ChannelRateParams,
            stim_intervals: &json::JsonValue,
            emis_params: &json::JsonValue) -> Self {
        Simulation {
            num_channels: num_channels,
            num_ts: num_ts,
            channels: vec![Channel::new(); num_channels as usize],
            stimulus: generate_stimulus(stim_intervals, num_ts),
            emis_dists: EmisDist::from(emis_params),
            transition_matrix: [
                [0.0, 4.0 * channel_rate_params.open_rate, 0.0, 0.0, 0.0],
                [channel_rate_params.close_rate, 0.0, 3.0 * channel_rate_params.open_rate, 0.0, 0.0],
                [0.0, 2.0 * channel_rate_params.close_rate, 0.0, 2.0 * channel_rate_params.open_rate, 0.0],
                [0.0, 0.0, 3.0 * channel_rate_params.close_rate, 0.0, channel_rate_params.open_rate],
                [0.0, 0.0, 0.0, 4.0 * channel_rate_params.close_rate, 0.0],
            ],
            channel_rate_params: channel_rate_params,
        }
    }

    fn update_rates(&mut self, ts: u32) {
        // for now: hard-coding in values for delayed K+ rectifier conductance
        self.channel_rate_params.open_rate =
            DEFAULT_OPEN_RATE *
            (self.channel_rate_params.open_exp_const * self.stimulus[ts as usize]).exp();
        self.channel_rate_params.close_rate =
            DEFAULT_CLOSE_RATE *
            (self.channel_rate_params.close_exp_const * self.stimulus[ts as usize]).exp();
    }

    fn update_transition_matrix(&mut self) {
        self.transition_matrix = [
            [0.0, 4.0 * self.channel_rate_params.open_rate, 0.0, 0.0, 0.0],
            [self.channel_rate_params.close_rate, 0.0, 3.0 * self.channel_rate_params.open_rate, 0.0, 0.0],
            [0.0, 2.0 * self.channel_rate_params.close_rate, 0.0, 2.0 * self.channel_rate_params.open_rate, 0.0],
            [0.0, 0.0, 3.0 * self.channel_rate_params.close_rate, 0.0, self.channel_rate_params.open_rate],
            [0.0, 0.0, 0.0, 4.0 * self.channel_rate_params.close_rate, 0.0]
        ];
    }

    fn run(&mut self) {
        for ts in 0..self.num_ts {
            self.update_rates(ts);
            self.update_transition_matrix();
            for channel in &mut self.channels {
                channel.make_transition(&self.transition_matrix);
                channel.sample_emission(&self.emis_dists);
                println!("ts: {ts}, state: {:?}, emis: {}", channel.state, channel.emis);
            }
        }
    }
}

// FIXME: transition matrix requires diagonal elements to satisfy prob dist req
// (see https://link.springer.com/referenceworkentry/10.1007/978-1-4614-6675-8_131)
//
// FIXME: change make_transition function: use transition rates and simulate channel
// probabilities over time
fn main() -> Result<(), Error> {
    // initial values for delayed K+ rectifier conductance
    let param_str: String = fs::read_to_string(PARAM_FILE_NAME)?;
    let params_json = json::parse(param_str.as_str()).unwrap();
    let mut model_sim = Simulation::from(params_json["num_channels"].as_u32().unwrap(),
                                         params_json["num_ts"].as_u32().unwrap(),
                                         ChannelRateParams::from(
                                            params_json["init_open_rate"].as_f32().unwrap(),
                                            params_json["open_exp_const"].as_f32().unwrap(),
                                            params_json["init_close_rate"].as_f32().unwrap(),
                                            params_json["close_exp_const"].as_f32().unwrap(),
                                        ),
                                         &params_json["stimulus"],
                                         &params_json["emissions"]);

    model_sim.run();

    Ok(())
}

