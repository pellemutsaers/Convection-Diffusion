// DEFINE YOUR OWN SOURCE FUNCTION
// CANVAS IS (0, 1) x (0, 1)

pub fn fill_sourcefunction(f_lx_vec: &mut Vec<f64>, m: usize, h_x: f64) {
    for i in 0..m - 1 {
        for j in 0..m - 1 {
            let x = (i as f64) * h_x;
            let y = (j as f64) * h_x;
    
            let center_x = 0.5;
            let center_y = 0.5;
            let radius = 0.05;
    
            let distance_squared = (x - center_x).powi(2) + (y - center_y).powi(2);
    
            if distance_squared <= radius * radius {
                f_lx_vec[i * (m - 1) + j] = 1.0;
            } else {
                f_lx_vec[i * (m - 1) + j] = 0.0;
            }
        }
    }
}