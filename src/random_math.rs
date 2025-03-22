use crate::matrix;

// Kronecker function, automatically adds result to "result"
pub fn kron(a: &matrix::Matrix, b: &matrix::Matrix, result: &mut matrix::Matrix) { 
    for (&(i, j), &a_val) in &a.data {
        for (&(k, l), &b_val) in &b.data {
            let idx = (i * b.rows + k, j * b.cols + l);
            let entry = result.data.entry(idx).or_insert(0.0);
            *entry += a_val * b_val;
        }
    }
}

pub fn forward_euler_step(
    ia_matrix: &matrix::Matrix,
    sol_vec: &Vec<f64>,
    src_rust: &Vec<f64>,
    vx: &Vec<f64>,
    vy: &Vec<f64>,
    m: usize,  // Grid resolution
) -> Vec<f64> {
    let mut new_sol: Vec<f64> = ia_matrix
        .mul_vector(sol_vec)
        .iter()
        .zip(src_rust.iter())
        .map(|(&val1, &val2)| val1 + val2)
        .collect();

    let h_x = 1.0 / m as f64;

    // Compute upwind finite differences for advection
    for i in 1..(m - 2) { // Exclude boundaries
        for j in 1..(m - 2) {
            let idx = i * (m - 1) + j;

            let u_ij = sol_vec[idx];

            // Upwind scheme for first-order convection term
            let du_dx = if vx[idx] > 0.0 {
                (u_ij - sol_vec[idx - 1]) / h_x  // Backward difference
            } else {
                (sol_vec[idx + 1] - u_ij) / h_x  // Forward difference
            };

            let du_dy = if vy[idx] > 0.0 {
                (u_ij - sol_vec[idx - (m - 1)]) / h_x
            } else {
                (sol_vec[idx + (m - 1)] - u_ij) / h_x
            };

            // Apply advection term
            new_sol[idx] -= vx[idx] * du_dx * h_x;
            new_sol[idx] -= vy[idx] * du_dy * h_x;
        }
    }

    new_sol
}


pub fn fill_dirichlet_advection(
    a_xx: &mut matrix::Matrix,
    a_yy: &mut matrix::Matrix,
    i_x: &mut matrix::Matrix,
    i_y: &mut matrix::Matrix,
    vx: &Vec<f64>,
    vy: &Vec<f64>,
    d: f64,
    hx: f64,
    hy: f64,
    m: usize
) {
    for i in 0..m - 2 {
        let peclet_x = vx[i] * hx / (2.0 * d);
        let peclet_y = vy[i] * hy / (2.0 * d);

        // Second-derivative matrix in X-direction (Diffusion + Advection)
        a_xx.insert(i, i, 2.0);
        a_xx.insert(i + 1, i, -1.0 - peclet_x);
        a_xx.insert(i, i + 1, -1.0 + peclet_x);
        i_x.insert(i, i, 1.0);

        // Second-derivative matrix in Y-direction (Diffusion + Advection)
        a_yy.insert(i, i, 2.0);
        a_yy.insert(i + 1, i, -1.0 - peclet_y);
        a_yy.insert(i, i + 1, -1.0 + peclet_y);
        i_y.insert(i, i, 1.0);
    }

    a_xx.insert(m - 2, m - 2, 2.0);
    i_x.insert(m - 2, m - 2, 1.0);

    a_yy.insert(m - 2, m - 2, 2.0);
    i_y.insert(m - 2, m - 2, 1.0);
}


pub fn fill_neumann_advection(
    a_xx: &mut matrix::Matrix,
    a_yy: &mut matrix::Matrix,
    i_x: &mut matrix::Matrix,
    i_y: &mut matrix::Matrix,
    vx: &Vec<f64>,
    vy: &Vec<f64>,
    d: f64,
    hx: f64,
    hy: f64,
    m: usize
) {
    for i in 0..m - 2 {
        let peclet_x = vx[i] * hx / (2.0 * d);
        let peclet_y = vy[i] * hy / (2.0 * d);

        // Second-derivative matrix in X-direction
        a_xx.insert(i, i, 2.0);
        a_xx.insert(i + 1, i, -1.0 - peclet_x);
        a_xx.insert(i, i + 1, -1.0 + peclet_x);
        i_x.insert(i, i, 1.0);

        // Second-derivative matrix in Y-direction
        a_yy.insert(i, i, 2.0);
        a_yy.insert(i + 1, i, -1.0 - peclet_y);
        a_yy.insert(i, i + 1, -1.0 + peclet_y);
        i_y.insert(i, i, 1.0);
    }

    // Adjust boundary conditions for Neumann
    a_xx.insert(0, 0, 1.0 - vx[0] * hx / (2.0 * d));
    a_xx.insert(m - 2, m - 2, 1.0 + vx[m - 2] * hx / (2.0 * d));
    i_x.insert(0, 0, 1.0);
    i_x.insert(m - 2, m - 2, 1.0);

    a_yy.insert(0, 0, 1.0 - vy[0] * hy / (2.0 * d));
    a_yy.insert(m - 2, m - 2, 1.0 + vy[m - 2] * hy / (2.0 * d));
    i_y.insert(0, 0, 1.0);
    i_y.insert(m - 2, m - 2, 1.0);
}
