use cfd_rs_utils::mesh::computational_mesh::*;
use nalgebra::Vector2;

#[derive(Clone, Debug, PartialEq)]
pub struct CellFaceField {
    cell: Vec<f64>,
    face: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PhysicalQuantities {
    Speed(Vector2<Variables>),
    Pressure(Variables),
    KinematicViscosity(Variables),
    Density(Variables),
    Temperature(Variables),
    Other(Variables, String),
}

#[derive(Default, Clone, Debug, PartialEq)]
pub enum Variables {
    #[default]
    None,
    Constant(f64),
    CellField(Vec<f64>),
    FaceField(Vec<f64>),
    CellFaceField(CellFaceField),
}


pub trait Case {
    fn mesh(&self) -> Computational2DMesh;
    
    /// Gives the variables of the case (fields and constants)
    fn physical_quantities(&self) -> Vec<&PhysicalQuantities>;
}