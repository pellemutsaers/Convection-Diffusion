use colorgrad::preset::magma;
use colorgrad::Gradient;
use gif::{Encoder, Frame, Repeat};
use std::fs::File;

/// Normalize the vector to fit in the 0-255 range
fn normalize(data: &[f64]) -> Vec<u8> {
    let min_val = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_val = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    data.iter()
        .map(|&v| ((v - min_val) / (max_val - min_val) * 255.0) as u8)
        .collect()
}

/// Generate a magma palette using the colorgrad crate
fn magma_palette() -> Vec<u8> {
    let grad = magma(); // Get the magma gradient

    let mut palette: Vec<u8> = Vec::new();

    for i in 0..256 {
        let rgba = grad.at(i as f32 / 255.0 as f32).to_rgba8(); // Get the color at each step
        palette.extend_from_slice(&rgba[0..3]); // Add RGB values to the palette
    }

    palette
}

/// Convert a solution vector into an indexed magma GIF frame (flipped) with interpolation
fn vector_to_frame(data: &[f64], size: usize, output_size: usize) -> Frame<'static> {
    let normalized = normalize(data);

    let mut interpolated_data = vec![0u8; output_size * output_size];

    for y_out in 0..output_size {
        for x_out in 0..output_size {
            let x = x_out as f64 * (size - 1) as f64 / (output_size - 1) as f64;
            let y = y_out as f64 * (size - 1) as f64 / (output_size - 1) as f64;

            let x_floor = x.floor() as usize;
            let y_floor = y.floor() as usize;
            let x_ceil = (x_floor + 1).min(size - 1);
            let y_ceil = (y_floor + 1).min(size - 1);

            let x_weight = x - x_floor as f64;
            let y_weight = y - y_floor as f64;

            let a = normalized[(size - 1 - y_floor) * size + x_floor] as f64;
            let b = normalized[(size - 1 - y_floor) * size + x_ceil] as f64;
            let c = normalized[(size - 1 - y_ceil) * size + x_floor] as f64;
            let d = normalized[(size - 1 - y_ceil) * size + x_ceil] as f64;

            let interpolated_value = (a * (1.0 - x_weight) * (1.0 - y_weight)
                + b * x_weight * (1.0 - y_weight)
                + c * (1.0 - x_weight) * y_weight
                + d * x_weight * y_weight) as u8;

            interpolated_data[y_out * output_size + x_out] = interpolated_value;
        }
    }

    let mut frame = Frame::from_indexed_pixels(
        output_size as u16,
        output_size as u16,
        &interpolated_data,
        None,
    );

    frame.delay = 5;
    frame
}

/// Create an animated GIF from solution frames with desired duration, FPS, and output size
pub fn create_gif(frames: Vec<Vec<f64>>, grid_size: usize, gif_name: &str, gif_duration_seconds: u16, fps: u8, output_size: usize) {
    let file = File::create(gif_name).unwrap();
    let magma_palette = magma_palette();

    let mut encoder = Encoder::new(file, output_size as u16, output_size as u16, &magma_palette[..]).unwrap();
    encoder.set_repeat(Repeat::Infinite).unwrap();

    let total_frames = frames.len();
    let desired_frame_count = gif_duration_seconds as usize * fps as usize;

    let mut frame_indices = Vec::new();
    if desired_frame_count >= total_frames {
        frame_indices = (0..total_frames).collect();
    } else {
        for i in 0..desired_frame_count {
            let index = (i as f64 / desired_frame_count as f64 * (total_frames - 1) as f64).round() as usize;
            frame_indices.push(index);
        }
    }

    let frame_delay = (100 / fps) as u16; // Delay in hundredths of a second

    for &index in frame_indices.iter() {
        let mut frame = vector_to_frame(&frames[index], grid_size, output_size);
        frame.delay = frame_delay;
        encoder.write_frame(&frame).unwrap();
    }
}