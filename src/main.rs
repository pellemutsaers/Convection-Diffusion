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
    fn new(d: f64, v: f64, m: usize, h_t: f64, boundary_conditions: BoundaryConditions) -> Self {
        let h_x = 1.0 / m as f64;
        let peclet = v * h_x / (2.0 * d);
        let size = (m - 1) * (m - 1);
        let mut f_lx_vec = vec![0.0; size];
        sourcefunc::fill_sourcefunction(&mut f_lx_vec, m, h_x); // Define sourcefunction in sourcefunc.rs

        // Identity matrix
        let mut matrix_i = matrix::Matrix::new(size, size);
        for i in 0..size {
            matrix_i.insert(i, i, 1.0);
        }

        // Coefficient matrix
        let matrix_a = Self::construct_a(m, h_x, peclet, d, &boundary_conditions);

        Self {
            h_t,
            matrix_a,
            matrix_i,
            f_lx_vec,
        }
    }

    fn construct_a(m: usize, hx: f64, peclet: f64, d: f64, boundary_conditions: &BoundaryConditions) -> matrix::Matrix {
        let mut a_xx = matrix::Matrix::new(m - 1, m - 1);
        let mut i_x = matrix::Matrix::new(m - 1, m - 1);

        match boundary_conditions {
            BoundaryConditions::Dirichlet => { math::fill_dirichlet(&mut a_xx, &mut i_x, peclet, m); }
            BoundaryConditions::Neumann => { math::fill_neumann(&mut a_xx, &mut i_x, peclet, m); }
        }

        let mut a = matrix::Matrix::new((m - 1) * (m - 1), (m - 1) * (m - 1));

        // Magic math logic
        let scaled_a_xx = a_xx.clone() * (d / (hx * hx));
        math::kron(&scaled_a_xx, &i_x, &mut a);
        math::kron(&i_x, &scaled_a_xx, &mut a);

        a
    }

    // Currently only using Forward Euler: 
    // u_next = u + h_t*(-A*u + f), or rewritten:
    // u_next = (I − h_t*A)u + h_t*f
    fn solve(&self, t_end: f64, m: usize) -> Vec<f64> {
        let total_gif_duration= 5;
        let fps = 30;
        let desired_frame_count = total_gif_duration * fps;
        let n = (t_end / self.h_t).round() as usize; // Number of timesteps to take
        
        let mut frames = Vec::with_capacity(desired_frame_count as usize); // Used for creating simulation gif.
        let mut sol_vec = vec![0.0; (self.matrix_a.rows) as usize]; // The simulation data, sometimes copied into frames.

        let ia_matrix = &self.matrix_i - &(self.matrix_a.clone() * self.h_t); // Precalculating I − h_t*A
        let src_rust: Vec<f64> = self.f_lx_vec.iter().map(|&val| val * self.h_t).collect(); // Source function scaled by h_t
        
        let pb = ProgressBar::new(n as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{elapsed_precise} [{bar:40}] {percent}%")
                .progress_chars("=>-"),
        );
        
        // Time integrating using Forward Euler
        for step in 0..n {
            sol_vec = math::forward_euler_step(&ia_matrix, &sol_vec, &src_rust);
        
            // Decides if this frame should be added to the GIF
            if frames.len() < desired_frame_count as usize {
                let expected_step = (frames.len() as f64 / desired_frame_count as f64 * n as f64).round() as usize;
                if step == expected_step {
                    frames.push(sol_vec.clone());
                }
            }
            pb.set_position(step as u64);
        }
        pb.finish_with_message("Calculation complete!");
        gif::create_gif(frames, m - 1, "Convection-diffusion.gif", total_gif_duration as u16, fps);
        
        sol_vec // return the final solution at time t = t_end.
    }
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let res = input::read_res();
    // let res: usize = 100;
    let bounds = input::read_bounds()?;
    // let boundary_conditions = BoundaryConditions::Neumann;
    let start = SystemTime::now();

    const D: f64 = 0.10;
    const V: f64 = 0.80;
    let h_t = 1.0 / (4.0 * res as f64 * res as f64 * D + res as f64 * V);

    println!("Stable timestep: {}", h_t);

    let boundary_conditions = match BoundaryConditions::from_input(&bounds) {
        Ok(bc) => bc,
        Err(e) => {
            println!("{}", e);
            return Ok(());
        }
    };

    let solver = ConvectionDiffusion::new(D, V, res, h_t, boundary_conditions); // Initialize Convection-diffusion problem object
    let t_end = 1.0;
    let solution = solver.solve(t_end, res); // Solves the convection-diffusion-problem at t=t_end from t=0

    plotting::plot_solution(&solution, res)?;

    println!(
        "Calculation time: {} seconds.",
        start.elapsed().unwrap().as_secs_f64()
    );

    Ok(())
}
