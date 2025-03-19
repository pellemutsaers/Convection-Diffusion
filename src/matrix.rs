use std::collections::HashMap;

#[derive(Clone)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: HashMap<(usize, usize), f64>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: HashMap::new(),
        }
    }

    pub fn mul_vector(&self, vec: &Vec<f64>) -> Vec<f64> {
        let mut result = vec![0.0; vec.len()];
        for (&(i, j), &val) in &self.data {
            result[i] += val * vec[j];
        }
        result
    }

    pub fn insert(&mut self, i: usize, j: usize, value: f64) {
        self.data.insert((i, j), value);
    }
}

impl std::ops::Mul<f64> for Matrix {
    type Output = Matrix;

    fn mul(mut self, scalar: f64) -> Matrix {
        for value in self.data.values_mut() {
            *value *= scalar;
        }
        self
    }
}

impl std::ops::Sub<&Matrix> for &Matrix {
    type Output = Matrix;

    fn sub(self, other: &Matrix) -> Matrix {
        let mut result = self.clone();

        for (&(i, j), &val) in &other.data {
            let entry = result.data.entry((i, j)).or_insert(0.0);
            *entry -= val;
        }
        result
    }
}