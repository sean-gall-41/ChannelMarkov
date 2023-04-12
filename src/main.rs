extern crate rand;

use std::fs;
use std::io::Error;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
};
use rand_distr::{Normal};
use plotters::prelude::*;

const PARAM_FILE_NAME: &'static str = "params.json";
//const OUT_IMG_FILE_NAME: &'static str = "example_trace.png";
//const OUT_IMG_FILE_NAME: &'static str = "example_cum_emis_trace.png";
const OUT_IMG_FILE_NAME: &'static str = "example_hhk_emis_trace.png";

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
    pub emis: f32,
    pub state_hist: Vec<u32>,
    pub emis_hist: Vec<f32>
}

impl Channel {
    fn new() -> Self {
        Channel {
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
        self.emis_hist.push(self.emis);
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
fn gen_stimulus(stim_intervals: &json::JsonValue) -> Vec<(f32, f32)> {
    assert!(stim_intervals.is_array());
    let stim: Vec<(f32, f32)> = stim_intervals.members().map(|p| {
        (p["step_t"].as_f32().unwrap(), p["step_v"].as_f32().unwrap())
    }).collect();
    stim
}

#[derive(Debug)]
pub struct ChannelParams {
    pub e_rev: f32,
    pub init_open_prob: f32,
    pub init_close_prob: f32,
    pub open_prob: f32,
    pub open_exp_const: f32,
    pub close_prob: f32,
    pub close_exp_const: f32,
}

impl ChannelParams {
    fn from(e_rev: f32,
            init_open_prob: f32,
            open_exp_const: f32,
            init_close_prob: f32,
            close_exp_const: f32) -> Self {
       ChannelParams {
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

fn gen_emis_params(emis_params: &json::JsonValue) -> Vec<(f32, f32)> {
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
        let emis_params_vec: Vec<(f32, f32)> = gen_emis_params(emis_params);
        assert_eq!(emis_params_vec.len(), 5); // the number of different distributions we need
        EmisDist {
            closed_1_dist: Normal::new(emis_params_vec[0].0, emis_params_vec[0].1).unwrap(),
            closed_2_dist: Normal::new(emis_params_vec[1].0, emis_params_vec[1].1).unwrap(),
            closed_3_dist: Normal::new(emis_params_vec[2].0, emis_params_vec[2].1).unwrap(),
            closed_4_dist: Normal::new(emis_params_vec[3].0, emis_params_vec[3].1).unwrap(),
            open_dist:     Normal::new(emis_params_vec[4].0, emis_params_vec[4].1).unwrap()
        }
    }
}

fn gen_trans_matrix(c_p: &ChannelParams) -> [[f32; 5]; 5] {
    [
        [1.0 - 4.0 * c_p.open_prob, 4.0 * c_p.open_prob, 0.0, 0.0, 0.0],
        [c_p.close_prob, 1.0 - c_p.close_prob - 3.0 * c_p.open_prob, 3.0 * c_p.open_prob, 0.0, 0.0],
        [0.0, 2.0 * c_p.close_prob, 1.0 - 2.0 * (c_p.close_prob + c_p.open_prob), 2.0 * c_p.open_prob, 0.0],
        [0.0, 0.0, 3.0 * c_p.close_prob, 1.0 - 3.0 * c_p.close_prob - c_p.open_prob, c_p.open_prob],
        [0.0, 0.0, 0.0, 4.0 * c_p.close_prob, 1.0 - 4.0 * c_p.close_prob],
    ]
}

#[derive(Debug)]
pub struct HMMSim {
    pub num_channels: u32,
    pub total_time: u32,
    pub dt: f32,
    pub channels: Vec<Channel>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_dists: EmisDist,
    pub cum_emis: Vec<f32>,
    pub trans_matrix: [[f32; 5]; 5],
    pub c_p: ChannelParams,
}

impl HMMSim {
    fn from(num_channels: u32,
            total_time: u32,
            dt: f32,
            c_p: ChannelParams,
            stim_intervals: &json::JsonValue,
            emis_params: &json::JsonValue) -> Self {
        HMMSim {
            num_channels: num_channels,
            total_time: total_time,
            dt: dt,
            channels: vec![Channel::new(); num_channels as usize],
            stimulus: gen_stimulus(stim_intervals),
            emis_dists: EmisDist::from(emis_params),
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

    fn run(&mut self) {
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

/*
 * struct representing the simulation/ updating of the dynamic var n
 * in the HH scheme
 */
pub struct HHKSim {
    pub total_time: u32,
    pub dt: f32,
    pub k_g_max: f32,
    pub k_probs: Vec<f32>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_hist: Vec<f32>,
    pub c_p: ChannelParams
}

impl HHKSim {
    fn from(total_time: u32,
            dt: f32,
            k_g_max: f32,
            c_p: ChannelParams,
            stim_intervals: &json::JsonValue) -> Self {
        let mut k_probs_init = vec![0.0; (total_time as f32 / dt) as usize];
        k_probs_init[0] = 0.5; // initial prob at ts == 0
        HHKSim {
            total_time: total_time,
            dt: dt,
            k_g_max: k_g_max,
            k_probs: k_probs_init,
            stimulus: gen_stimulus(stim_intervals),
            emis_hist: vec![0.0; (total_time as f32 / dt) as usize],
            c_p: c_p,
        }
    }

    // here we interpret these guys differently: as rates instead of probs.
    // the name remains as a prob (my bad, needs refactor) but they are
    // to be understood as rates
    fn update_rates(&mut self, voltage: f32) {
        self.c_p.open_prob =
            self.c_p.init_open_prob *
            (self.c_p.open_exp_const * voltage).exp();
        self.c_p.close_prob =
            self.c_p.init_close_prob *
            (self.c_p.close_exp_const * voltage).exp();
    }

    fn run(&mut self) {
        let mut stim_id: usize = 0;
        for ts in 1..((self.total_time as f32 / self.dt) as u32) {
            if stim_id < self.stimulus.len() && (ts-1) == (self.stimulus[stim_id].0 / self.dt) as u32 {
                self.update_rates(self.stimulus[stim_id].1);
                stim_id += 1;
            }
            let ts: usize = ts as usize;
            // no need to use time step, considered in probs update
            self.k_probs[ts] = self.k_probs[ts-1]
                        + (1.0 - self.k_probs[ts-1]) * self.c_p.open_prob
                        - self.k_probs[ts-1] * self.c_p.close_prob;

            self.emis_hist[ts] = self.k_g_max
                               * self.k_probs[ts].powf(4.0)
                               * (self.stimulus[stim_id-1].1 - self.c_p.e_rev);
        }
    }

}

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
    pub e_rev: f32,
    pub init_open_prob: f32,
    pub init_close_prob: f32,
    pub open_prob: f32,
    pub open_exp_const: f32,
    pub close_prob: f32,
    pub close_exp_const: f32,
}

impl NaChannelParams {
    fn from(e_rev: f32,
            init_open_prob: f32,
            open_exp_const: f32,
            init_close_prob: f32,
            close_exp_const: f32) -> Self {
       NaChannelParams {
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
pub struct NaEmisDist {
    pub closed_1_dist: Normal<f32>,
    pub closed_2_dist: Normal<f32>,
    pub closed_3_dist: Normal<f32>,
    pub open_dist: Normal<f32>,
    pub inactive_dist: Normal<f32>,
}

impl NaEmisDist {
    fn from(emis_params: &json::JsonValue) -> Self {
        let emis_params_vec: Vec<(f32, f32)> = gen_emis_params(emis_params);
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

// TODO: update to reflect above Na scheme
fn gen_na_trans_matrix(c_p: &NaChannelParams) -> [[f32; 5]; 5] {
    [
        [1.0 - 4.0 * c_p.open_prob, 4.0 * c_p.open_prob, 0.0, 0.0, 0.0],
        [c_p.close_prob, 1.0 - c_p.close_prob - 3.0 * c_p.open_prob, 3.0 * c_p.open_prob, 0.0, 0.0],
        [0.0, 2.0 * c_p.close_prob, 1.0 - 2.0 * (c_p.close_prob + c_p.open_prob), 2.0 * c_p.open_prob, 0.0],
        [0.0, 0.0, 3.0 * c_p.close_prob, 1.0 - 3.0 * c_p.close_prob - c_p.open_prob, c_p.open_prob],
        [0.0, 0.0, 0.0, 4.0 * c_p.close_prob, 1.0 - 4.0 * c_p.close_prob],
    ]
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // initial values for delayed K+ rectifier conductance
    let params_json = json::parse(fs::read_to_string(PARAM_FILE_NAME)?.as_str()).unwrap();
    let k_model_params = &params_json["potassium"];
    let mut model_hmm_sim = HMMSim::from(k_model_params["num_channels"].as_u32().unwrap(),
                                         k_model_params["total_time"].as_u32().unwrap(),
                                         k_model_params["dt"].as_f32().unwrap(),
                                         ChannelParams::from(
                                            k_model_params["e_reverse"].as_f32().unwrap(),
                                            k_model_params["init_open_rate"].as_f32().unwrap()
                                            * k_model_params["dt"].as_f32().unwrap(),
                                            k_model_params["open_exp_const"].as_f32().unwrap(),
                                            k_model_params["init_close_rate"].as_f32().unwrap()
                                            * k_model_params["dt"].as_f32().unwrap(),
                                            k_model_params["close_exp_const"].as_f32().unwrap(),
                                        ),
                                         &k_model_params["stimulus"],
                                         &k_model_params["emissions"]);

    model_hmm_sim.run();

    let mut model_hhk_sim = HHKSim::from(k_model_params["total_time"].as_u32().unwrap(),
                                         k_model_params["dt"].as_f32().unwrap(),
                                         k_model_params["emissions"][4]["mu"].as_f32().unwrap(),
                                         ChannelParams::from(
                                            k_model_params["e_reverse"].as_f32().unwrap(),
                                            k_model_params["init_open_rate"].as_f32().unwrap()
                                            * k_model_params["dt"].as_f32().unwrap(),
                                            k_model_params["open_exp_const"].as_f32().unwrap(),
                                            k_model_params["init_close_rate"].as_f32().unwrap()
                                            * k_model_params["dt"].as_f32().unwrap(),
                                            k_model_params["close_exp_const"].as_f32().unwrap(),
                                        ),
                                         &k_model_params["stimulus"]);

    model_hhk_sim.run();

    // it's plotting time
    let root = BitMapBackend::new(OUT_IMG_FILE_NAME, (1024, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.titled(format!("HMM vs HH emissions num_channels={}", model_hmm_sim.num_channels).as_str(),
                           ("sans-serif", 40))?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 45)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(-50f32..((model_hmm_sim.total_time as f32) / model_hmm_sim.dt + 50f32),
                                                (-15.0 * model_hmm_sim.num_channels as f32)..(40.0 * model_hmm_sim.num_channels as f32))?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_labels(6)
        .x_desc("time (microsecs)")
        .y_desc("Current (pA)")
        .draw()?;

    chart.draw_series(LineSeries::new(
        model_hmm_sim.cum_emis.iter().enumerate().map(|(ts, e)| {
            (ts as f32, *e)
        }),
        &BLUE,
    ))?;
    chart.draw_series(LineSeries::new(
        model_hhk_sim.emis_hist.iter().enumerate().map(|(ts, e)| {
            (ts as f32, model_hmm_sim.num_channels as f32 * *e)
        }),
        &BLACK,
    ))?;

    // uncomment for individual channel emission time series
    //for channel in model_sim.channels {
    //    chart.draw_series(LineSeries::new(
    //        channel.state_hist.iter().enumerate().map(|(ts, s)| {
    //            (ts as f32, 0.65f32 * (*s as f32) + 10f32)
    //        }),
    //        &full_palette::ORANGE,
    //    ))?;
    //    chart.draw_series(LineSeries::new(
    //        channel.emis_hist.iter().enumerate().map(|(ts, e)| {
    //            (ts as f32, *e)
    //        }),
    //        &BLUE,
    //    ))?;
    //}

    root.present().expect("Unable to write result to file.");
    println!("Result have been saved to {}", OUT_IMG_FILE_NAME);

    Ok(())
}

