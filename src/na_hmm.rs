extern crate rand;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
};
use rand_distr::{Normal};

use crate::util;

/* ====================================================================================
 * |                                                                                  |
 * |                             Sodium Channel Modeling                              |
 * |                                                                                  |
 * ====================================================================================
 */

/*
 * Here's the basic state diagram that we're assuming our Markov Model follows:
 *
 *                                            k1
 *                      ------------------------------------------------
 *                     /                                ah              \/
 *  close            close            close            open             inact
 * +-----+          +-----+          +-----+          +-----+          +-----+
 * |     |    3am   |     |    2am   |     |    am    |     |    k3    |     |
 * |  1  |  <---->  |  2  |  <---->  |  3  |  <---->  |  4  |   ---->  |  5  |
 * |     |    1bm   |     |    2bm   |     |    3bm   |     |          |     |
 * +-----+          +-----+          +-----+          +-----+          +-----+
 *                                     /\               k2             -/
 *                                       -------------------------------
 *                                                      ah
 *
 * where here am == activation rate and bm == deactivation rate,
 *            k_1 == inactivation rate 1, k_2 == inactivation rate 2,
 *            k_3 inactivation rate 3, ah == deinactivation rate
 *
 */

#[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
#[repr(i32)]
pub enum NaChannelState {
    Closed1 = 0,
    Closed2,
    Closed3,
    Open,
    Inactive,
}

impl Distribution<NaChannelState> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> NaChannelState {
        unsafe { ::std::mem::transmute(rng.gen_range(0..=4)) }
    }
}

#[derive(Debug, Clone)]
pub struct NaChannel {
    pub state: NaChannelState,
    pub emis: f32,
    pub state_hist: Vec<u32>,
    pub emis_hist: Vec<f32>
}

impl NaChannel {
    fn new() -> Self {
        NaChannel {
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

    fn sample_emission(&mut self, emis_dist: &NaEmisDist) {
        match self.state {
            NaChannelState::Closed1 => {
                self.emis = emis_dist.closed_1_dist.sample(&mut rand::thread_rng());
            },
            NaChannelState::Closed2 => {
                self.emis = emis_dist.closed_2_dist.sample(&mut rand::thread_rng());
            },
            NaChannelState::Closed3 => {
                self.emis = emis_dist.closed_3_dist.sample(&mut rand::thread_rng());
            },
            NaChannelState::Open => {
                self.emis = emis_dist.open_dist.sample(&mut rand::thread_rng());
            },
            NaChannelState::Inactive => {
                self.emis = emis_dist.inactive_dist.sample(&mut rand::thread_rng());
            },
        }
        self.emis_hist.push(self.emis);
    }
}

// TODO: update with Na params
#[derive(Debug)]
pub struct NaChannelParams {
    pub k_1: f32,
    pub k_2: f32,
    pub k_3: f32,
    pub m_open_prob: f32,
    pub init_m_open_prob: f32,
    pub m_pre_v_fact_open: f32,
    pub m_v_offset_open: f32,
    pub m_close_prob: f32,
    pub init_m_close_prob: f32,
    pub close_m_exp_const: f32,
    pub m_v_offset_close: f32,
    pub h_open_prob: f32,
    pub init_h_open_prob: f32,
    pub h_v_offset_open: f32,
    pub open_h_exp_const: f32,
    pub h_close_prob: f32,
    pub init_h_close_prob: f32,
    pub h_pre_v_fact_close: f32,
    pub h_v_offset_close: f32,
    pub e_rev: f32,
}

impl NaChannelParams {
    pub fn from(k_1: f32,
            k_2: f32,
            k_3: f32,
            init_m_open_prob: f32,
            m_pre_v_fact_open: f32,
            m_v_offset_open: f32,
            init_m_close_prob: f32,
            close_m_exp_const: f32,
            m_v_offset_close: f32,
            init_h_open_prob: f32,
            h_v_offset_open: f32,
            open_h_exp_const: f32,
            init_h_close_prob: f32,
            h_pre_v_fact_close: f32,
            h_v_offset_close: f32,
            e_rev: f32) -> Self {
       NaChannelParams {
            k_1: k_1,
            k_2: k_2,
            k_3: k_3,
            m_open_prob: init_m_open_prob,
            init_m_open_prob: init_m_open_prob,
            m_pre_v_fact_open: m_pre_v_fact_open,
            m_v_offset_open: m_v_offset_open,
            m_close_prob: init_m_close_prob,
            init_m_close_prob: init_m_close_prob,
            close_m_exp_const: close_m_exp_const,
            m_v_offset_close: m_v_offset_close,
            h_open_prob: init_h_open_prob,
            init_h_open_prob: init_h_open_prob,
            h_v_offset_open: h_v_offset_open,
            open_h_exp_const: open_h_exp_const,
            h_close_prob: init_h_close_prob,
            init_h_close_prob: init_h_close_prob,
            h_pre_v_fact_close: h_pre_v_fact_close,
            h_v_offset_close: h_v_offset_close,
            e_rev: e_rev
       }
    }
}

#[derive(Debug)]
pub struct NaEmisDist {
    pub closed_1_dist: Normal<f32>,
    pub closed_2_dist: Normal<f32>,
    pub closed_3_dist: Normal<f32>,
    pub open_dist: Normal<f32>,
    pub inactive_dist: Normal<f32>,
}

impl NaEmisDist {
    fn from(emis_params: &json::JsonValue) -> Self {
        let emis_params_vec: Vec<(f32, f32)> = util::gen_emis_params(emis_params);
        assert_eq!(emis_params_vec.len(), 5); // the number of different distributions we need
        NaEmisDist {
            closed_1_dist: Normal::new(emis_params_vec[0].0, emis_params_vec[0].1).unwrap(),
            closed_2_dist: Normal::new(emis_params_vec[1].0, emis_params_vec[1].1).unwrap(),
            closed_3_dist: Normal::new(emis_params_vec[2].0, emis_params_vec[2].1).unwrap(),
            open_dist:     Normal::new(emis_params_vec[3].0, emis_params_vec[3].1).unwrap(),
            inactive_dist: Normal::new(emis_params_vec[4].0, emis_params_vec[4].1).unwrap(),
        }
    }
}

fn gen_na_trans_matrix(c_p: &NaChannelParams) -> [[f32; 5]; 5] {
    [
        [1.0 - 3.0 * c_p.m_open_prob, 3.0 * c_p.m_open_prob, 0.0, 0.0, 0.0],
        [c_p.m_close_prob, 1.0 - c_p.m_close_prob - 2.0 * c_p.m_open_prob - c_p.k_1, 2.0 * c_p.m_open_prob, 0.0, c_p.k_1],
        [0.0, 2.0 * c_p.m_close_prob, 1.0 - 2.0 * c_p.m_close_prob - c_p.m_open_prob - c_p.k_2, c_p.m_open_prob, c_p.k_2],
        [0.0, 0.0, 3.0 * c_p.m_close_prob, 1.0 - 3.0 * c_p.m_close_prob - c_p.k_3, c_p.k_3],
        [0.0, 0.0, c_p.h_open_prob, 0.0, 1.0 - c_p.h_open_prob],
    ]
}

#[derive(Debug)]
pub struct NaHMMSim {
    pub num_channels: u32,
    pub total_time: u32,
    pub dt: f32,
    pub channels: Vec<NaChannel>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_dists: NaEmisDist,
    pub cum_emis: Vec<f32>,
    pub trans_matrix: [[f32; 5]; 5],
    pub c_p: NaChannelParams,
}

impl NaHMMSim {
    pub fn from(num_channels: u32,
            total_time: u32,
            dt: f32,
            c_p: NaChannelParams,
            stim_intervals: &json::JsonValue,
            emis_params: &json::JsonValue) -> Self {
        NaHMMSim {
            num_channels: num_channels,
            total_time: total_time,
            dt: dt,
            channels: vec![NaChannel::new(); num_channels as usize],
            stimulus: util::gen_stimulus(stim_intervals),
            emis_dists: NaEmisDist::from(emis_params),
            cum_emis: vec![0.0; (total_time as f32 / dt) as usize],
            trans_matrix: gen_na_trans_matrix(&c_p),
            c_p: c_p,
        }
    }

    fn update_probs(&mut self, voltage: f32) {
        // NOTE: since multiplying by init probs, already taken extra mult factor of dt
        // into acct
        self.c_p.m_open_prob = (self.c_p.init_m_open_prob * (voltage - self.c_p.m_v_offset_open))
                             / (1.0 - (-self.c_p.m_pre_v_fact_open * (voltage - self.c_p.m_v_offset_open)).exp());
        self.c_p.m_close_prob = self.c_p.init_m_close_prob * (self.c_p.close_m_exp_const * (voltage - self.c_p.m_v_offset_close)).exp();
        self.c_p.h_open_prob = self.c_p.init_h_open_prob * (self.c_p.open_h_exp_const * (voltage - self.c_p.h_v_offset_open)).exp();
        self.c_p.h_close_prob = self.c_p.init_h_close_prob / (1.0 + (self.c_p.h_pre_v_fact_close * (voltage - self.c_p.h_v_offset_close)).exp());
    }

    fn update_trans_matrix(&mut self) {
        self.trans_matrix = gen_na_trans_matrix(&self.c_p);
    }

    pub fn run(&mut self) {
        let mut stim_id: usize = 0;
        for ts in 0..((self.total_time as f32 / self.dt) as u32) {
            if stim_id < self.stimulus.len() && ts == (self.stimulus[stim_id].0 / self.dt) as u32 {
                self.update_probs(self.stimulus[stim_id].1);
                self.update_trans_matrix();
                //if ts == 0 { println!("{:#?}", self.trans_matrix); }
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

