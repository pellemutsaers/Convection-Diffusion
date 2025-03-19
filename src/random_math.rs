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

pub fn forward_euler_step(ia_matrix: &matrix::Matrix, sol_vec: &Vec<f64>, src_rust: &Vec<f64>) -> Vec<f64> {
    ia_matrix
        .mul_vector(sol_vec)
        .iter()
        .zip(src_rust.iter())
        .map(|(&val1, &val2)| val1 + val2)
        .collect()
}

pub fn fill_dirichlet(a_xx: &mut matrix::Matrix, i_x: &mut matrix::Matrix, peclet: f64, m: usize) {
    for i in 0..m - 2 {
        a_xx.insert(i, i, 2.0);
        a_xx.insert(i + 1, i, -1.0 - peclet);
        a_xx.insert(i, i + 1, -1.0 + peclet);
        i_x.insert(i, i, 1.0);
    }
    a_xx.insert(m - 2, m - 2, 2.0);
    i_x.insert(m - 2, m - 2, 1.0);
}

pub fn fill_neumann(a_xx: &mut matrix::Matrix, i_x: &mut matrix::Matrix, peclet: f64, m: usize) {
    for i in 0..m - 2 {
        a_xx.insert(i, i, 2.0);
        a_xx.insert(i + 1, i, -1.0 - peclet);
        a_xx.insert(i, i + 1, -1.0 + peclet);
        i_x.insert(i, i, 1.0);
    }
    a_xx.insert(0, 0, 1.0 - peclet);
    a_xx.insert(m - 2, m - 2, 1.0 + peclet);
    i_x.insert(0, 0, 1.0);
    i_x.insert(m - 2, m - 2, 1.0);
}