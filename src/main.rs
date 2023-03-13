extern crate rand;

use std::io::Error;
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

const DEFAULT_NUM_CHANNELS: u32 = 1;
const DEFAULT_NUM_TS: u32 = 1000; // ms
const CLOSE_RATE: f32 = 0.01;
const OPEN_RATE: f32 = 0.01;

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
const TRANSITION_MATRIX: [[f32; 5]; 5] = [
    [0.0, 4.0 * CLOSE_RATE,              0.0,              0.0,            0.0],
    [OPEN_RATE,        0.0, 3.0 * CLOSE_RATE,              0.0,            0.0],
    [0.0,  2.0 * OPEN_RATE,              0.0, 2.0 * CLOSE_RATE,            0.0],
    [0.0,              0.0,  3.0 * OPEN_RATE,              0.0,     CLOSE_RATE],
    [0.0,              0.0,              0.0,  4.0 * OPEN_RATE,            0.0],
];

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
}

impl Channel {
    fn new(default: bool) -> Self {
        if default {
            Channel {
                state: ChannelState::Open
            }
        } else {
            Channel {
                state: rand::random()
            }
        }
    }

    fn from(in_state: ChannelState) -> Self {
        Channel {
            state: in_state
        }
    }
}

pub struct Simulation<'a> {
    pub num_channels: u32,
    pub num_ts: u32,
    pub channels: Vec<Channel>,
    pub t_mat: &'a[[f32; 5]; 5],
}

impl Simulation<'static> {
    fn new() -> Self {
        Simulation {
            num_channels: DEFAULT_NUM_CHANNELS,
            num_ts: DEFAULT_NUM_TS,
            channels: vec![Channel::new(true); DEFAULT_NUM_CHANNELS as usize],
            t_mat: &TRANSITION_MATRIX,
        }
    }

    fn run(&self) {
        todo!["pls impl me :-)"];
    }
}

fn main() -> Result<(), Error> {
    let model_sim = Simulation::new();
    println!("Hello, world!");
    Ok(())
}

