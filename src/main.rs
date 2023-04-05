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
const OUT_IMG_FILE_NAME: &'static str = "example_trace.png";

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
    pub occ_time: f32,
    pub state: ChannelState,
    pub emis: f32,
    pub state_hist: Vec<(f32, u32)>,
    //pub emis_hist: Vec<f32>
}

impl Channel {
    fn new() -> Self {
        Channel {
            occ_time: 0.0,
            state: rand::random(),
            emis: 0.0f32,
            state_hist: vec![], // look, I'm not going to say this is the fastest solution...
            //emis_hist: vec![],
        }
    }

    fn compute_occupancy_time(&mut self, occ_rates: &[f32; 5]) {
        let occ_rate = occ_rates[self.state as usize];
        let r: f32 = rand::random::<f32>();
        // need to tweak this due to log(u) for u -> 0.0
        self.occ_time = - (1.0 / occ_rate) * r.ln();
        self.state_hist.push((self.occ_time, self.state as u32));
    }

    fn make_transition(&mut self, cond_trans_matrix: &[[f32; 5]; 5]) {
        let transition_row = cond_trans_matrix[self.state as usize];
        let dist = WeightedIndex::new(&transition_row).unwrap();
        let transition_prob = transition_row[dist.sample(&mut rand::thread_rng())];
        self.state = unsafe { ::std::mem::transmute(transition_row.iter()
                                                                  .position(|&r| r == transition_prob)
                                                                  .unwrap() as i32)
        };
        //self.state_hist.push(self.state as u32);
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
        //self.emis_hist.push(self.emis);
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
fn generate_stimulus(stim_intervals: &json::JsonValue) -> Vec<(f32, f32)> {
    assert!(stim_intervals.is_array());
    let stim: Vec<(f32, f32)> = stim_intervals.members().map(|p| {
        (p["step_dt"].as_f32().unwrap(), p["step_v"].as_f32().unwrap())
    }).collect();
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
    pub total_time: u32,
    pub dt: f32,
    pub channels: Vec<Channel>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_dists: EmisDist,
    pub cum_emis_hist: Vec<f32>,
    pub occ_rates: [f32; 5],
    pub cond_trans_matrix: [[f32; 5]; 5],
    pub c_r_p: ChannelRateParams,
}


fn generate_occ_rates(c_r_p: &ChannelRateParams) -> [f32; 5] {
    [4.0 * c_r_p.open_rate,
     3.0 * c_r_p.open_rate + c_r_p.close_rate,
     2.0 * (c_r_p.open_rate + c_r_p.close_rate),
     c_r_p.open_rate + 3.0 * c_r_p.close_rate,
     4.0 * c_r_p.close_rate]
}

fn generate_cond_trans_matrix(c_r_p: &ChannelRateParams) -> [[f32; 5]; 5] {
    [
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [c_r_p.close_rate / (3.0 * c_r_p.open_rate + c_r_p.close_rate), 0.0, (3.0 * c_r_p.open_rate) / (3.0 * c_r_p.open_rate + c_r_p.close_rate), 0.0, 0.0],
        [0.0, c_r_p.close_rate / (c_r_p.open_rate + c_r_p.close_rate), 0.0, c_r_p.open_rate / (c_r_p.open_rate + c_r_p.close_rate), 0.0],
        [0.0, 0.0, (3.0 * c_r_p.close_rate) / (c_r_p.open_rate + 3.0 * c_r_p.close_rate), 0.0, c_r_p.open_rate / (c_r_p.open_rate + 3.0 * c_r_p.close_rate)],
        [0.0, 0.0, 0.0, 1.0, 0.0],
    ]
}

fn update_rates(c_r_p: &mut ChannelRateParams, stim: &Vec<(f32, f32)>, index: usize) {
    // for now: hard-coding in values for delayed K+ rectifier conductance
    c_r_p.open_rate =
        DEFAULT_OPEN_RATE *
        (c_r_p.open_exp_const * stim[index].1).exp();
    c_r_p.close_rate =
        DEFAULT_CLOSE_RATE *
        (c_r_p.close_exp_const * stim[index].1).exp();
}

//fn update_occ_rates(c_r_p: &ChannelRateParams, occ_rates: &mut [f32; 5]) {
//    occ_rates = generate_occ_rates(c_r_p);
//}

//fn update_cond_trans_matrix(&mut self) {
//    self.cond_trans_matrix =  generate_cond_trans_matrix(&self.c_r_p);
//}



impl Simulation {
    fn from(num_channels: u32,
            total_time: u32,
            dt: f32,
            c_r_p: ChannelRateParams,
            stim_intervals: &json::JsonValue,
            emis_params: &json::JsonValue) -> Self {
        Simulation {
            num_channels: num_channels,
            total_time: total_time,
            dt: dt,
            channels: vec![Channel::new(); num_channels as usize],
            stimulus: generate_stimulus(stim_intervals),
            emis_dists: EmisDist::from(emis_params),
            cum_emis_hist: vec![],
            occ_rates: generate_occ_rates(&c_r_p),
            cond_trans_matrix: generate_cond_trans_matrix(&c_r_p),
            c_r_p: c_r_p,
        }
    }

    fn run(&mut self) {
        for channel in &mut self.channels {
            let mut ts: f32 = 0.0;
            let mut stim_id: usize = 0;
            update_rates(&mut self.c_r_p, &self.stimulus, stim_id);
            //update_occ_rates(&self.c_r_p, &mut self.occ_rates);
            //self.update_cond_trans_matrix();
            while ts < (self.total_time as f32) {
                if ts > self.stimulus[stim_id].0 {
                    stim_id += 1;
                    update_rates(&mut self.c_r_p, &self.stimulus, stim_id); // only update rates once we've asynchonously
                                                // crossed a piec-wise voltage step boundary
                    //update_occ_rates(&self.c_r_p, &mut self.occ_rates);
                    //self.update_cond_trans_matrix();
                }
                channel.compute_occupancy_time(&self.occ_rates);
                channel.make_transition(&self.cond_trans_matrix);
                ts += channel.occ_time;
                println!("{}, {:?}", ts, self.occ_rates);
            }
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // initial values for delayed K+ rectifier conductance
    let params_json = json::parse(fs::read_to_string(PARAM_FILE_NAME)?.as_str()).unwrap();
    let mut model_sim = Simulation::from(params_json["num_channels"].as_u32().unwrap(),
                                         params_json["total_time"].as_u32().unwrap(),
                                         params_json["dt"].as_f32().unwrap(),
                                         ChannelRateParams::from(
                                            params_json["init_open_rate"].as_f32().unwrap(),
                                            params_json["open_exp_const"].as_f32().unwrap(),
                                            params_json["init_close_rate"].as_f32().unwrap(),
                                            params_json["close_exp_const"].as_f32().unwrap(),
                                        ),
                                         &params_json["stimulus"],
                                         &params_json["emissions"]);

    model_sim.run();

    // it's plotting time
    //let root = BitMapBackend::new(OUT_IMG_FILE_NAME, (1024, 600)).into_drawing_area();
    //root.fill(&WHITE)?;
    //let root = root.titled("State and Emission Sequences", ("sans-serif", 40))?;

    //let mut chart = ChartBuilder::on(&root)
    //    .margin(10)
    //    .set_label_area_size(LabelAreaPosition::Left, 45)
    //    .set_label_area_size(LabelAreaPosition::Bottom, 40)
    //    .build_cartesian_2d(-50f32..((model_sim.total_time as f32) / model_sim.dt + 50f32),
    //                                            -4.0f32..14.0f32)?;

    //chart
    //    .configure_mesh()
    //    .disable_x_mesh()
    //    .disable_y_mesh()
    //    .x_labels(6)
    //    .x_desc("time (microsecs)")
    //    .y_desc("Current (pA)")
    //    .draw()?;

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

    //root.present().expect("Unable to write result to file.");
    //println!("Result have been saved to {}", OUT_IMG_FILE_NAME);

    Ok(())
}

