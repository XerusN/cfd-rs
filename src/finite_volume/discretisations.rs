use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::DMatrix;
use nalgebra_sparse::csr::CsrRow;

#[derive(Clone, Debug, PartialEq)]
pub enum SpaceDiscretizationConfig {
    CentralDifference,
    Upwind,
    SecondOrderUpwind,
    FROMM,
    Quick,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TimeDiscretizationConfig {
    FirstOrderEuler,
}

pub fn central_difference<'a>(mesh: &Computational2DMesh, cell_id: usize, row: CsrRow<'a, f64>) -> CsrRow<'a, f64> {
    todo!();
}