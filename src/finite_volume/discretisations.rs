use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::DMatrix;
use nalgebra_sparse::{csr::CsrRow, CsrMatrix};
use nalgebra::Vector2;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct CellFaceArrays {
    pub cell: Option<Vec<f64>>,
    pub face: Option<Vec<f64>>
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct VectorCellFaceArrays {
    pub cell: Option<Vector2<Vec<f64>>>,
    pub face: Option<Vector2<Vec<f64>>>
}


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

#[derive(Clone, Debug, PartialEq)]
pub enum SpaceDiscretizationConfig {
    CentralDifference,
    Upwind,
    SecondOrderUpwind,
    Fromm,
    Quick,
}

impl SpaceDiscretizationConfig {
    fn discretize<'a, 'b>(&self) -> fn(&Computational2DMesh, usize, Row<'a, 'b>, &CellFaceArrays, Option<&VectorCellFaceArrays>) -> Row<'a, 'b> {
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

pub fn central_difference<'a, 'b>(mesh: &Computational2DMesh, cell_id: usize, row: Row<'a, 'b>, value: &CellFaceArrays, grad: Option<&VectorCellFaceArrays>) -> Row<'a, 'b> {
    todo!();
}