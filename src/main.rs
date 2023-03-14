extern crate rand;

use std::io::Error;
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

const DEFAULT_NUM_CHANNELS: u32 = 1;
const DEFAULT_NUM_TS: u32 = 1000; // ms
const DEFAULT_STIM_MV: f32 = 10.0;
const DEFAULT_CLOSE_RATE: f32 = 0.01;
const DEFAULT_OPEN_RATE: f32 = 0.01;

/* think about it from matrix mult lens: the col index is the state we come
 * from, and the row index is the state we go to
 *
 *  close            close            close            close            open
 * +-----+          +-----+          +-----+          +-----+          +-----+
 * |     |    4a    |     |    3a    |     |    2a    |     |    1a    |     |
 * |  1  |  <---->  |  2  |  <---->  |  3  |  <---->  |  4  |  <---->  |  5  |
 * |     |    1b    |     |    2b    |     |    3b    |     |    4b    |     |
 * +-----+          +-----+          +-----+          +-----+          +-----+
 *
 * where here a == OPEN_RATE and b == CLOSE_RATE
 */


#[derive(Debug, PartialEq, PartialOrd, Clone)]
pub enum ChannelState {
    Closed1,
    Closed2,
    Closed3,
    Closed4,
    Open,
}

impl Distribution<ChannelState> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ChannelState {
        match rng.gen_range(0..=4) {
            0 => ChannelState::Closed1,
            1 => ChannelState::Closed2,
            2 => ChannelState::Closed3,
            3 => ChannelState::Closed4,
            _ => ChannelState::Open,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Channel {
    pub state: ChannelState,
    pub open_rate: f32,
    pub close_rate: f32,
}

impl Channel {
    fn new(default: bool) -> Self {
        if default {
            Channel {
                state: ChannelState::Open,
                open_rate: DEFAULT_OPEN_RATE,
                close_rate: DEFAULT_CLOSE_RATE,
            }
        } else {
            Channel {
                state: rand::random(),
                open_rate: rand::random(),
                close_rate: rand::random(),
            }
        }
    }

    //fn from_state(in_state: ChannelState) -> Self {
    //    Channel {
    //        state: in_state,
    //        open_rate: DEFAULT_OPEN_RATE,
    //        close_rate: DEFAULT_CLOSE_RATE,
    //    }
    //}

    fn from_rates(open_rate: f32, close_rate: f32) -> Self {
        Channel {
            state: rand::random(),
            open_rate: open_rate,
            close_rate: close_rate,
        }
    }

    fn update_rates(&mut self, stim: &Vec<f32>, ts: u32) {
        // for now: hard-coding in values for delayed K+ rectifier conductance
        self.open_rate = 1.22 * (0.04 * stim[ts as usize]).exp();
        self.close_rate = 0.056 * (-0.0125 * stim[ts as usize]).exp();
    }

    fn make_transition(&mut self) {
        let choice = rand::thread_rng().gen_range(0..2);
        match self.state {
            ChannelState::Closed1 => {
                if rand::thread_rng().gen::<f32>() < 4.0 * DEFAULT_OPEN_RATE {
                    self.state = ChannelState::Closed2;
                    return;
                }
            },
            ChannelState::Closed2 => {
                match choice {
                    0 => {
                        if rand::thread_rng().gen::<f32>() < DEFAULT_CLOSE_RATE {
                            self.state = ChannelState::Closed1;
                        }
                    },
                    1 => {
                        if rand::thread_rng().gen::<f32>() < 3.0 * DEFAULT_OPEN_RATE {
                            self.state = ChannelState::Closed3;
                        }
                    },
                    _ => (),
                }
                return;
            },
            ChannelState::Closed3 => {
                match choice {
                    0 => {
                        if rand::thread_rng().gen::<f32>() < 2.0 * DEFAULT_CLOSE_RATE {
                            self.state = ChannelState::Closed2;
                        }
                    },
                    1 => {
                        if rand::thread_rng().gen::<f32>() < 2.0 * DEFAULT_OPEN_RATE {
                            self.state = ChannelState::Closed4;
                        }
                    },
                    _ => (),
                }
                return;
            },
            ChannelState::Closed4 => {
                match choice {
                    0 => {
                        if rand::thread_rng().gen::<f32>() < 3.0 * DEFAULT_CLOSE_RATE {
                            self.state = ChannelState::Closed3;
                        }
                    },
                    1 => {
                        if rand::thread_rng().gen::<f32>() < DEFAULT_OPEN_RATE {
                            self.state = ChannelState::Open;
                        }
                    },
                    _ => (),
                }
                return;
            },
            ChannelState::Open => {
                if rand::thread_rng().gen::<f32>() < 4.0 * DEFAULT_CLOSE_RATE {
                    self.state = ChannelState::Closed4;
                    return;
                }
            },
        }
    }
}

pub struct Simulation {
    pub num_channels: u32,
    pub num_ts: u32,
    pub channels: Vec<Channel>,
    pub stimulus: Vec<f32>,
}

impl Simulation {
    fn new() -> Self {
        Simulation {
            num_channels: DEFAULT_NUM_CHANNELS,
            num_ts: DEFAULT_NUM_TS,
            channels: vec![Channel::new(true); DEFAULT_NUM_CHANNELS as usize],
            stimulus: vec![DEFAULT_STIM_MV; DEFAULT_NUM_TS as usize],
        }
    }

    fn from(num_channels: u32,
            num_ts: u32,
            init_open_rate: f32,
            init_close_rate: f32,
            stimulus: Vec<f32>) -> Self {
        Simulation {
            num_channels: num_channels,
            num_ts: num_ts,
            channels: vec![Channel::from_rates(init_open_rate, init_close_rate); num_channels as usize],
            stimulus: stimulus,
        }
    }

    fn run(&mut self) {
        for ts in 0..self.num_ts {
            for channel in &mut self.channels {
                channel.update_rates(&self.stimulus, ts);
                channel.make_transition();
                if channel.state == ChannelState::Open {
                    println!("ts: {ts}, Open");
                }
            }
        }
    }
}

/*
 * Couple of fixes:
 *
 * a. need to scale the transition rates so lie in range [0, 1)
 * b. need to initialize rates to *correct* values wrt input current at t == 0
 * c. need to create function to generate stimulus
 * d. even later: put all relevant params in json and read from
 */
fn main() -> Result<(), Error> {
    // initial values for delayed K+ rectifier conductance
    let init_open_rate: f32 = 1.22;
    let init_close_rate: f32 = 0.056;
    let mut model_sim = Simulation::from(1,
                                         50,
                                         init_open_rate,
                                         init_close_rate,
                                         vec![DEFAULT_STIM_MV; 50]);
    //let mut model_sim = Simulation::new();
    model_sim.run();
    Ok(())
}

