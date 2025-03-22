use plotters::prelude::*;
use colorgrad::preset::magma;
use colorgrad::Gradient;

pub fn plot_solution(solution: &[f64], res: usize) -> Result<(), Box<dyn std::error::Error>> {
    // Create drawing area
    let root = BitMapBackend::new("Convection-diffusion.png", (1000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Heatmap", ("sans-serif", 30))
        // Add these lines to set margins to zero
        .margin_top(0)
        .margin_bottom(0)
        .margin_left(0)
        .margin_right(0)
        .build_cartesian_2d(0..res - 1, 0..res - 1)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    // Normalize the color mapping to the range [0, 1]
    let min_val = solution.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_val = solution.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let epsilon = 1e-16;
    let value_range = (max_val - min_val).max(epsilon);

    let grad = magma(); // Or choose another: inferno, magma, plasma, cividis, turbo
    // let grad = colorgrad::diverging_linear_bwr_55_98_c37_n256();  // For diverging colormaps

    for i in 0..res - 1 {
        for j in 0..res - 1 {
            let idx = i * (res - 1) + j;
            let value = solution[idx];

            // Normalize value to [0, 1] range
            let normalized_value = (value - min_val) / value_range;

            // Get the color from the gradient (CORRECTED)
            let color = grad.at(normalized_value as f32);
            let color = RGBColor(
                (color.r * 255.0) as u8,
                (color.g * 255.0) as u8,
                (color.b * 255.0) as u8,
            );

            chart.draw_series(std::iter::once(Rectangle::new(
                [(i, j), (i + 1, j + 1)],
                color.filled(),
            )))?;
        }
    }
    Ok(())
}