use crate::util;
use crate::k_hmm;
/*
 * struct representing the simulation/ updating of the dynamic var n
 * in the HH scheme
 */
pub struct KHHSim {
    pub total_time: u32,
    pub dt: f32,
    pub k_g_max: f32,
    pub k_probs: Vec<f32>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_hist: Vec<f32>,
    pub c_p: k_hmm::KChannelParams
}

impl KHHSim {
    pub fn from(total_time: u32,
            dt: f32,
            k_g_max: f32,
            c_p: k_hmm::KChannelParams,
            stim_intervals: &json::JsonValue) -> Self {
        let mut k_probs_init = vec![0.0; (total_time as f32 / dt) as usize];
        k_probs_init[0] = 0.5; // initial prob at ts == 0
        KHHSim {
            total_time: total_time,
            dt: dt,
            k_g_max: k_g_max,
            k_probs: k_probs_init,
            stimulus: util::gen_stimulus(stim_intervals),
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

    pub fn run(&mut self) {
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

