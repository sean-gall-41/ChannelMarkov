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
pub fn gen_stimulus(stim_intervals: &json::JsonValue) -> Vec<(f32, f32)> {
    assert!(stim_intervals.is_array());
    let stim: Vec<(f32, f32)> = stim_intervals.members().map(|p| {
        (p["step_t"].as_f32().unwrap(), p["step_v"].as_f32().unwrap())
    }).collect();
    stim
}

pub fn gen_emis_params(emis_params: &json::JsonValue) -> Vec<(f32, f32)> {
    assert!(emis_params.is_array());
    let emis: Vec<(f32, f32)> = emis_params.members().map(|p| {
        (p["mu"].as_f32().unwrap(), p["sigma"].as_f32().unwrap())
    }).collect();
    emis
}

