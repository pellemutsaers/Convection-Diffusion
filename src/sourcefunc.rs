// DEFINE YOUR OWN SOURCE FUNCTION
// CANVAS IS (0, 1) x (0, 1)

pub fn fill_sourcefunction(f_lx_vec: &mut Vec<f64>, m: usize, h_x: f64) {
    for i in 0..m - 1 {
        for j in 0..m - 1 {
            let x = (i as f64) * h_x;
            let y = (j as f64) * h_x;

            let mut source_value = 0.0;

            // // Dot 1 (center)
            // let center_x1 = 0.5;
            // let center_y1 = 0.5;
            // let radius1 = 0.05;
            // let distance_squared1 = (x - center_x1).powi(2) + (y - center_y1).powi(2);
            // if distance_squared1 <= radius1 * radius1 {
            //     source_value = 1.0;
            // }

            // Dot 2 (offset)
            let center_x2 = 0.2;
            let center_y2 = 0.8;
            let radius2 = 0.05;
            let distance_squared2 = (x - center_x2).powi(2) + (y - center_y2).powi(2);
            if distance_squared2 <= radius2 * radius2 {
                source_value = 1.0;
            }

            // Dot 3 (other offset)
            let center_x3 = 0.8;
            let center_y3 = 0.2;
            let radius3 = 0.05;
            let distance_squared3 = (x - center_x3).powi(2) + (y - center_y3).powi(2);
            if distance_squared3 <= radius3 * radius3 {
                source_value = 1.0;
            }

            f_lx_vec[i * (m - 1) + j] = source_value;
        }
    }
}

pub fn create_rotating_vector_field(m: usize) -> (Vec<f64>, Vec<f64>) {
    let mut vx = vec![0.0; (m - 1) * (m - 1)];
    let mut vy = vec![0.0; (m - 1) * (m - 1)];

    let center_x = (m - 1) as f64 / 2.0;
    let center_y = (m - 1) as f64 / 2.0;

    for i in 0..(m - 1) {
        for j in 0..(m - 1) {
            let index = i * (m - 1) + j;
            let x = j as f64;
            let y = i as f64;

            vx[index] = -(y - center_y) * 0.001; 
            vy[index] = (x - center_x) * 0.001;
        }
    }

    (vx, vy)
}