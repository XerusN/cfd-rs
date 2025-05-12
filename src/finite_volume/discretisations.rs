use std::ops::{Add, Div, Mul, Sub};

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::DMatrix;
use nalgebra_sparse::{csr::CsrRow, CsrMatrix};
use nalgebra::Vector2;

use super::base::ScalarVariable;

#[derive(Clone, Debug, PartialEq)]
pub struct Row<'a, 'b> {
    pub row: &'a CsrRow<'b, f64>,
    pub values: Vec<f64>,
}

impl<'a, 'b> Row<'a, 'b> {
    fn new(row: &'a CsrRow<'b, f64>, values: Vec<f64>) -> Row<'a, 'b> {
        Row {
            row,
            values,
        }
    }
    
    /// Creates same row structure but with zero values
    fn from_csr_row(csr_row: &'a CsrRow<'b, f64>) -> Self {
        let t = csr_row.col_indices();
        let a = t.len();
        Row::new(csr_row, vec![0.; a])
    }
}

impl<'a, 'b> Add for &Row<'a, 'b> {
    type Output = Row<'a, 'b>;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.row, rhs.row, "Trying to add Row not refering to the same row in the matrix ({self:?} and {rhs:?})");
        let mut values = self.values.clone();
        for (i, value) in values.iter_mut().enumerate() {
            *value += rhs.values[i];
        }
        Row::new(self.row, values)
    }
}

impl<'a, 'b> Sub for &Row<'a, 'b> {
    type Output = Row<'a, 'b>;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.row, rhs.row, "Trying to substract Row not refering to the same row in the matrix ({self:?} and {rhs:?})");
        let mut values = self.values.clone();
        for (i, value) in values.iter_mut().enumerate() {
            *value -= rhs.values[i];
        }
        Row::new(self.row, values)
    }
}

impl<'a, 'b> Add<f64> for &Row<'a, 'b> {
    type Output = Row<'a, 'b>;
    fn add(self, rhs: f64) -> Self::Output {
        let mut values = self.values.clone();
        for value in values.iter_mut() {
            *value += rhs;
        }
        Row::new(self.row, values)
    }
}

impl<'a, 'b> Sub<f64> for &Row<'a, 'b> {
    type Output = Row<'a, 'b>;
    fn sub(self, rhs: f64) -> Self::Output {
        let mut values = self.values.clone();
        for value in values.iter_mut() {
            *value -= rhs;
        }
        Row::new(self.row, values)
    }
}

impl<'a, 'b> Mul<f64> for &Row<'a, 'b> {
    type Output = Row<'a, 'b>;
    fn mul(self, rhs: f64) -> Self::Output {
        let mut values = self.values.clone();
        for value in values.iter_mut() {
            *value *= rhs;
        }
        Row::new(self.row, values)
    }
}

impl<'a, 'b> Div<f64> for &Row<'a, 'b> {
    type Output = Row<'a, 'b>;
    fn div(self, rhs: f64) -> Self::Output {
        let mut values = self.values.clone();
        for value in values.iter_mut() {
            *value /= rhs;
        }
        Row::new(self.row, values)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum SpaceDiscretizationConfig {
    CentralDifference,
    Upwind,
    SecondOrderUpwind,
    Fromm,
    Quick,
}

impl SpaceDiscretizationConfig {
    pub fn discretize<'a, 'b>(&self) -> fn(&Computational2DMesh, usize, Row<'a, 'b>, &ScalarVariable) -> Row<'a, 'b> {
        match *self {
            Self::CentralDifference => central_difference,
            _ => panic!("{self:?} not implemented"),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum TimeDiscretizationConfig {
    FirstOrderEuler,
}

pub fn central_difference<'a, 'b>(mesh: &Computational2DMesh, cell_id: usize, row: Row<'a, 'b>, value: &ScalarVariable) -> Row<'a, 'b> {
    todo!();
}
