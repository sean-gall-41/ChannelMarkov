pub mod k_hmm;
pub mod k_hh;
pub mod na_hmm;
pub mod util;

use k_hmm::{KHMMSim, KChannelParams};
use k_hh::KHHSim;
use na_hmm::{NaHMMSim, NaChannelParams};
use std::fs;
use std::io::Error;
use plotters::prelude::*;

const PARAM_FILE_NAME: &'static str = "params.json";
//const OUT_IMG_FILE_NAME: &'static str = "example_trace.png";
//const OUT_IMG_FILE_NAME: &'static str = "example_cum_emis_trace.png";
const OUT_IMG_FILE_NAME: &'static str = "example_nahmm_emis_trace.png";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // initial values for delayed K+ rectifier conductance
    let params_json = json::parse(fs::read_to_string(PARAM_FILE_NAME)?.as_str()).unwrap();
    //let k_model_params = &params_json["potassium"];
    //let mut model_hmm_sim = KHMMSim::from(k_model_params["num_channels"].as_u32().unwrap(),
    //                                      k_model_params["total_time"].as_u32().unwrap(),
    //                                      k_model_params["dt"].as_f32().unwrap(),
    //                                      KChannelParams::from(
    //                                         k_model_params["e_reverse"].as_f32().unwrap(),
    //                                         k_model_params["init_open_rate"].as_f32().unwrap()
    //                                         * k_model_params["dt"].as_f32().unwrap(),
    //                                         k_model_params["open_exp_const"].as_f32().unwrap(),
    //                                         k_model_params["init_close_rate"].as_f32().unwrap()
    //                                         * k_model_params["dt"].as_f32().unwrap(),
    //                                         k_model_params["close_exp_const"].as_f32().unwrap(),
    //                                      ),
    //                                      &k_model_params["stimulus"],
    //                                      &k_model_params["emissions"]);

    //model_hmm_sim.run();

    //let mut model_hhk_sim = KHHSim::from(k_model_params["total_time"].as_u32().unwrap(),
    //                                      k_model_params["dt"].as_f32().unwrap(),
    //                                      k_model_params["emissions"][4]["mu"].as_f32().unwrap(),
    //                                      KChannelParams::from(
    //                                         k_model_params["e_reverse"].as_f32().unwrap(),
    //                                         k_model_params["init_open_rate"].as_f32().unwrap()
    //                                         * k_model_params["dt"].as_f32().unwrap(),
    //                                         k_model_params["open_exp_const"].as_f32().unwrap(),
    //                                         k_model_params["init_close_rate"].as_f32().unwrap()
    //                                         * k_model_params["dt"].as_f32().unwrap(),
    //                                         k_model_params["close_exp_const"].as_f32().unwrap(),
    //                                      ),
    //                                      &k_model_params["stimulus"]);

    //model_hhk_sim.run();

    let na_model_params = &params_json["sodium"];
    let mut model_nahmm_sim = NaHMMSim::from(na_model_params["num_channels"].as_u32().unwrap(),
                                           na_model_params["total_time"].as_u32().unwrap(),
                                           na_model_params["dt"].as_f32().unwrap(),
                                           NaChannelParams::from(
                                              na_model_params["k_1"].as_f32().unwrap()
                                              * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["k_2"].as_f32().unwrap()
                                              * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["k_3"].as_f32().unwrap()
                                              * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["m_pre_v_fact_open"].as_f32().unwrap()
                                              * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["m_pre_v_fact_open"].as_f32().unwrap(),
                                              na_model_params["m_v_offset_open"].as_f32().unwrap(),
                                              na_model_params["init_m_close_rate"].as_f32().unwrap()
                                              * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["close_m_exp_const"].as_f32().unwrap(),
                                              na_model_params["m_v_offset_close"].as_f32().unwrap(),
                                              na_model_params["init_h_open_rate"].as_f32().unwrap()
                                              * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["h_v_offset_open"].as_f32().unwrap(),
                                              na_model_params["open_h_exp_const"].as_f32().unwrap(),
                                              1.0 * na_model_params["dt"].as_f32().unwrap(),
                                              na_model_params["h_pre_v_fact_close"].as_f32().unwrap(),
                                              na_model_params["h_v_offset_close"].as_f32().unwrap(),
                                              na_model_params["e_reverse"].as_f32().unwrap()
                                           ),
                                           &na_model_params["stimulus"],
                                           &na_model_params["emissions"]);
    model_nahmm_sim.run();
    // it's plotting time
    let root = BitMapBackend::new(OUT_IMG_FILE_NAME, (1024, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.titled(format!("HMM vs HH emissions num_channels={}", model_nahmm_sim.num_channels).as_str(),
                           ("sans-serif", 40))?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 45)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(-50f32..((model_nahmm_sim.total_time as f32) / model_nahmm_sim.dt + 50f32),
                                                (-30.0 * model_nahmm_sim.num_channels as f32)..(15.0 * model_nahmm_sim.num_channels as f32))?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_labels(6)
        .x_desc("time (microsecs)")
        .y_desc("Current (pA)")
        .draw()?;

    chart.draw_series(LineSeries::new(
        model_nahmm_sim.cum_emis.iter().enumerate().map(|(ts, e)| {
            (ts as f32, *e)
        }),
        &BLUE,
    ))?;
    //chart.draw_series(LineSeries::new(
    //    model_hhk_sim.emis_hist.iter().enumerate().map(|(ts, e)| {
    //        (ts as f32, model_hmm_sim.num_channels as f32 * *e)
    //    }),
    //    &BLACK,
    //))?;

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

