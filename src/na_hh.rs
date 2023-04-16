use crate::util;
use crate::na_hmm;

pub struct NaHHSim {
    pub total_time: u32,
    pub dt: f32,
    pub na_g_max: f32,
    pub na_act_probs: Vec<f32>,
    pub na_inact_probs: Vec<f32>,
    pub stimulus: Vec<(f32, f32)>,
    pub emis_hist: Vec<f32>,
    pub c_p: na_hmm::NaChannelParams
}

impl NaHHSim {
    pub fn from(total_time: u32,
            dt: f32,
            na_g_max: f32,
            c_p: na_hmm::NaChannelParams,
            stim_intervals: &json::JsonValue) -> Self {
        let mut na_act_probs_init = vec![0.0; (total_time as f32 / dt) as usize];
        let mut na_inact_probs_init = vec![0.0; (total_time as f32 / dt) as usize];
        na_act_probs_init[0] = 0.5;
        na_inact_probs_init[0] = 0.5;
        NaHHSim {
            total_time: total_time,
            dt: dt,
            na_g_max: na_g_max,
            na_act_probs: na_act_probs_init,
            na_inact_probs: na_inact_probs_init,
            stimulus: util::gen_stimulus(stim_intervals),
            emis_hist: vec![0.0; (total_time as f32 / dt) as usize],
            c_p: c_p,
        }
    }

    fn update_rates(&mut self, voltage: f32) {
        self.c_p.m_open_prob = (self.c_p.init_m_open_prob * (voltage - self.c_p.m_v_offset_open))
                             / (1.0 - (-self.c_p.m_pre_v_fact_open * (voltage - self.c_p.m_v_offset_open)).exp());
        self.c_p.m_close_prob = self.c_p.init_m_close_prob * (self.c_p.close_m_exp_const * (voltage - self.c_p.m_v_offset_close)).exp();
        self.c_p.h_open_prob = self.c_p.init_h_open_prob * (self.c_p.open_h_exp_const * (voltage - self.c_p.h_v_offset_open)).exp();
        self.c_p.h_close_prob = self.c_p.init_h_close_prob / (1.0 + (self.c_p.h_pre_v_fact_close * (voltage - self.c_p.h_v_offset_close)).exp());
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
            self.na_act_probs[ts] = self.na_act_probs[ts-1]
                        + (1.0 - self.na_act_probs[ts-1]) * self.c_p.m_open_prob
                        - self.na_act_probs[ts-1] * self.c_p.m_close_prob;

            self.na_inact_probs[ts] = self.na_inact_probs[ts-1]
                        + (1.0 - self.na_inact_probs[ts-1]) * self.c_p.h_open_prob
                        - self.na_inact_probs[ts-1] * self.c_p.h_close_prob;


            self.emis_hist[ts] = self.na_g_max
                               * self.na_act_probs[ts].powf(3.0)
                               * self.na_inact_probs[ts]
                               * (self.stimulus[stim_id-1].1 - self.c_p.e_rev);
        }

    }
}
