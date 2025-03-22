mod input;
mod matrix;
mod plotting;
mod random_math;
mod gif;
mod sourcefunc;

use std::time::SystemTime;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use random_math as math;

#[derive(Debug)]
enum BoundaryConditions {
    Dirichlet,
    Neumann,
}

impl BoundaryConditions {
    fn from_input(input: &str) -> Result<BoundaryConditions, String> {
        match input.to_lowercase().as_str() {
            "dirichlet" => Ok(BoundaryConditions::Dirichlet),
            "neumann" => Ok(BoundaryConditions::Neumann),
            _ => Err("Invalid boundary condition".to_string()),
        }
    }
}

struct ConvectionDiffusion {
    h_t: f64,
    matrix_a: matrix::Matrix,
    matrix_i: matrix::Matrix,
    f_lx_vec: Vec<f64>,
}

impl ConvectionDiffusion {
    fn new(d: f64, vx: Vec<f64>, vy: Vec<f64>, m: usize, h_t: f64, boundary_conditions: BoundaryConditions) -> Self {
        let h_x = 1.0 / m as f64;
        let size = (m - 1) * (m - 1);
        let mut f_lx_vec = vec![0.0; size];
        sourcefunc::fill_sourcefunction(&mut f_lx_vec, m, h_x);

        // Identity matrix
        let mut matrix_i = matrix::Matrix::new(size, size);
        for i in 0..size {
            matrix_i.insert(i, i, 1.0);
        }

        // Construct matrix A using vector fields
        let matrix_a = Self::construct_a(m, h_x, h_x, &vx, &vy, d, &boundary_conditions);

        Self {
            h_t,
            matrix_a,
            matrix_i,
            f_lx_vec,
        }
    }

    fn construct_a(
        m: usize,
        hx: f64,
        hy: f64,
        vx: &Vec<f64>,
        vy: &Vec<f64>,
        d: f64,
        boundary_conditions: &BoundaryConditions,
    ) -> matrix::Matrix {
        let mut a_xx = matrix::Matrix::new(m - 1, m - 1);
        let mut a_yy = matrix::Matrix::new(m - 1, m - 1);
        let mut i_x = matrix::Matrix::new(m - 1, m - 1);
        let mut i_y = matrix::Matrix::new(m - 1, m - 1);

        match boundary_conditions {
            BoundaryConditions::Dirichlet => {
                math::fill_dirichlet_advection(&mut a_xx, &mut a_yy, &mut i_x, &mut i_y, vx, vy, d, hx, hy, m);
            }
            BoundaryConditions::Neumann => {
                math::fill_neumann_advection(&mut a_xx, &mut a_yy, &mut i_x, &mut i_y, vx, vy, d, hx, hy, m);
            }
        }

        let mut a = matrix::Matrix::new((m - 1) * (m - 1), (m - 1) * (m - 1));

        // Diffusion + advection matrix (Kronecker products)
        let scaled_a_xx = a_xx.clone() * (d / (hx * hx));
        let scaled_a_yy = a_yy.clone() * (d / (hy * hy));
        math::kron(&scaled_a_xx, &i_y, &mut a);
        math::kron(&i_x, &scaled_a_yy, &mut a);

        a
    }

    fn solve(&self, t_end: f64, m: usize, vx: &Vec<f64>, vy: &Vec<f64>) -> Vec<f64> {
        let total_gif_duration = 5;
        let fps = 30;
        let desired_frame_count = total_gif_duration * fps;
        let n = (t_end / self.h_t).round() as usize; // Number of timesteps to take

        let mut frames = Vec::with_capacity(desired_frame_count as usize);
        let mut sol_vec = vec![0.0; (m-1) * (m-1)];

        let ia_matrix = &self.matrix_i - &(self.matrix_a.clone() * self.h_t);
        let src_rust: Vec<f64> = self.f_lx_vec.iter().map(|&val| val * self.h_t).collect();

        let pb = ProgressBar::new(n as u64);
        pb.set_style(ProgressStyle::default_bar().template("{elapsed_precise} [{bar:40}] {percent}%"));

        for step in 0..n {

            sol_vec = math::forward_euler_step(&ia_matrix, &sol_vec, &src_rust, vx, vy, m);

            if frames.len() < desired_frame_count as usize {
                let expected_step = (frames.len() as f64 / desired_frame_count as f64 * n as f64).round() as usize;
                if step == expected_step {
                    frames.push(sol_vec.clone());
                }
            }

            pb.set_position(step as u64);
        }
        pb.finish_with_message("Calculation complete!");
        gif::create_gif(frames, (m - 1) as usize, "Convection-diffusion.gif", total_gif_duration as u16, fps, 200);
        sol_vec
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let res = input::read_res();
    let bounds = input::read_bounds()?;
    let start = SystemTime::now();

    const D: f64 = 0.10;
    // Example: Uniform flow in x, no flow in y initially
    let (vx, vy) = sourcefunc::create_rotating_vector_field(res); // Create rotating vector field


    let h_x = 1.0 / res as f64;
    let h_y = h_x; // Assume uniform grid spacing

    let max_vx = vx.iter().cloned().fold(0.0, f64::max); // Maximum velocity in x
    let max_vy = vy.iter().cloned().fold(0.0, f64::max); // Maximum velocity in y

    let h_t = (1.0 / ((4.0 * D) / (h_x * h_x) + (4.0 * D) / (h_y * h_y) + max_vx / h_x + max_vy / h_y)) / 10.0; // Timestep, divide by more if unstable

    println!("Stable timestep: {}", h_t);

    let boundary_conditions = match BoundaryConditions::from_input(&bounds) {
        Ok(bc) => bc,
        Err(e) => {
            println!("{}", e);
            return Ok(());
        }
    };

    let solver = ConvectionDiffusion::new(D, vx.clone(), vy.clone(), res, h_t, boundary_conditions);
    let t_end = 1.0;
    let solution = solver.solve(t_end, res, &vx, &vy);

    plotting::plot_solution(&solution, res)?;

    println!("Calculation time: {} seconds.", start.elapsed().unwrap().as_secs_f64());

    Ok(())
}