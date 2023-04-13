extern crate rand;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
};
use rand_distr::{Normal};

use crate::util;
/* ====================================================================================
 * |                                                                                  |
 * |                           Potassium Channel Modeling                             |
 * |                                                                                  |
 * ====================================================================================
 */

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
pub enum KChannelState {
    Closed1 = 0,
    Closed2,
    Closed3,
    Closed4,
    Open,
}

impl Distribution<KChannelState> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> KChannelState {
        unsafe { ::std::mem::transmute(rng.gen_range(0..=4)) }
    }
}

#[derive(Debug, Clone)]
pub struct KChannel {
    pub state: KChannelState,
    pub emis: f32,
    pub state_hist: Vec<u32>,
    pub emis_hist: Vec<f32>
}

impl KChannel {
    fn new() -> Self {
        KChannel {
            state: rand::random(),
            emis: 0.0f32,
            state_hist: vec![], // look, I'm not going to say this is the fastest solution...
            emis_hist: vec![],
        }
    }

    fn make_transition(&mut self, trans_matrix: &[[f32; 5]; 5]) {
        let transition_row = trans_matrix[self.state as usize];
        let dist = WeightedIndex::new(&transition_row).unwrap();
        let transition_prob = transition_row[dist.sample(&mut rand::thread_rng())];
        self.state = unsafe { ::std::mem::transmute(transition_row.iter()
                                                                  .position(|&r| r == transition_prob)
                                                                  .unwrap() as i32)
        };
        self.state_hist.push(self.state as u32);
    }

    fn sample_emission(&mut self, emis_dist: &KEmisDist) {
        match self.state {
            KChannelState::Closed1 => {
                self.emis = emis_dist.closed_1_dist.sample(&mut rand::thread_rng());
            },
            KChannelState::Closed2 => {
                self.emis = emis_dist.closed_2_dist.sample(&mut rand::thread_rng());
            },
            KChannelState::Closed3 => {
                self.emis = emis_dist.closed_3_dist.sample(&mut rand::thread_rng());
            },
            KChannelState::Closed4 => {
                self.emis = emis_dist.closed_4_dist.sample(&mut rand::thread_rng());
            },
            KChannelState::Open => {
                self.emis = emis_dist.open_dist.sample(&mut rand::thread_rng());
            },
        }
        self.emis_hist.push(self.emis);
    }
}

#[derive(Debug)]
pub struct KChannelParams {
    pub e_rev: f32,
    pub init_open_prob: f32,
    pub init_close_prob: f32,
    pub open_prob: f32,
    pub open_exp_const: f32,
    pub close_prob: f32,
    pub close_exp_const: f32,
}

impl KChannelParams {
    pub fn from(e_rev: f32,
            init_open_prob: f32,
            open_exp_const: f32,
            init_close_prob: f32,
            close_exp_const: f32) -> Self {
       KChannelParams {
        e_rev: e_rev,
        init_open_prob: init_open_prob,
        init_close_prob: init_close_prob,
        open_prob: init_open_prob,
        open_exp_const: open_exp_const,
        close_prob: init_close_prob,
        close_exp_const: close_exp_const
       }
    }
}

#[derive(Debug)]
pub struct KEmisDist {
    pub closed_1_dist: Normal<f32>,
    pub closed_2_dist: Normal<f32>,
    pub closed_3_dist: Normal<f32>,
    pub closed_4_dist: Normal<f32>,
    pub open_dist: Normal<f32>,
}

impl KEmisDist {
    fn from(emis_params: &json::JsonValue) -> Self {
        let emis_params_vec: Vec<(f32, f32)> = util::gen_emis_params(emis_params);
        assert_eq!(emis_params_vec.len(), 5); // the number of different distributions we need
        KEmisDist {
            closed_1_dist: Normal::new(emis_params_vec[0].0, emis_params_vec[0].1).unwrap(),
            closed_2_dist: Normal::new(emis_params_vec[1].0, emis_params_vec[1].1).unwrap(),
            closed_3_dist: Normal::new(emis_params_vec[2].0, emis_params_vec[2].1).unwrap(),
            closed_4_dist: Normal::new(emis_params_vec[3].0, emis_params_vec[3].1).unwrap(),
            open_dist:     Normal::new(emis_params_vec[4].0, emis_params_vec[4].1).unwrap()
        }
    }
}

fn gen_trans_matrix(c_p: &KChannelParams) -> [[f32; 5]; 5] {
    [
        [1.0 - 4.0 * c_p.open_prob, 4.0 * c_p.open_prob, 0.0, 0.0, 0.0],
        [c_p.close_prob, 1.0 - c_p.close_prob - 3.0 * c_p.open_prob, 3.0 * c_p.open_prob, 0.0, 0.0],
        [0.0, 2.0 * c_p.close_prob, 1.0 - 2.0 * (c_p.close_prob + c_p.open_prob), 2.0 * c_p.open_prob, 0.0],
        [0.0, 0.0, 3.0 * c_p.close_prob, 1.0 - 3.0 * c_p.close_prob - c_p.open_prob, c_p.open_prob],
        [0.0, 0.0, 0.0, 4.0 * c_p.close_prob, 1.0 - 4.0 * c_p.close_prob],
    ]
}

#[derive(Debug)]
pub struct KHMMSim {
    pub num_channels: u32,
    pub total_time: u32,
    pub dt: f32,
    pub channels: Vec<KChannel>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_dists: KEmisDist,
    pub cum_emis: Vec<f32>,
    pub trans_matrix: [[f32; 5]; 5],
    pub c_p: KChannelParams,
}

impl KHMMSim {
    pub fn from(num_channels: u32,
            total_time: u32,
            dt: f32,
            c_p: KChannelParams,
            stim_intervals: &json::JsonValue,
            emis_params: &json::JsonValue) -> Self {
        KHMMSim {
            num_channels: num_channels,
            total_time: total_time,
            dt: dt,
            channels: vec![KChannel::new(); num_channels as usize],
            stimulus: util::gen_stimulus(stim_intervals),
            emis_dists: KEmisDist::from(emis_params),
            cum_emis: vec![0.0; (total_time as f32 / dt) as usize],
            trans_matrix: gen_trans_matrix(&c_p),
            c_p: c_p,
        }
    }

    fn update_probs(&mut self, voltage: f32) {
        // NOTE: since multiplying by init probs, already taken extra mult factor of dt
        // into acct
        self.c_p.open_prob =
            self.c_p.init_open_prob *
            (self.c_p.open_exp_const * voltage).exp();
        self.c_p.close_prob =
            self.c_p.init_close_prob *
            (self.c_p.close_exp_const * voltage).exp();
    }

    fn update_trans_matrix(&mut self) {
        self.trans_matrix = gen_trans_matrix(&self.c_p);
    }

    pub fn run(&mut self) {
        let mut stim_id: usize = 0;
        for ts in 0..((self.total_time as f32 / self.dt) as u32) {
            if stim_id < self.stimulus.len() && ts == (self.stimulus[stim_id].0 / self.dt) as u32 {
                self.update_probs(self.stimulus[stim_id].1);
                self.update_trans_matrix();
                stim_id += 1;
            }
            let mut cum_emis_ts: f32 = 0.0;
            for channel in &mut self.channels {
                channel.make_transition(&self.trans_matrix);
                channel.sample_emission(&self.emis_dists);
                cum_emis_ts += channel.emis * (self.stimulus[stim_id-1].1 - self.c_p.e_rev);
            }
            self.cum_emis[ts as usize] = cum_emis_ts;
        }
    }
}

