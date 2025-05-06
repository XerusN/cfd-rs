use cfd_rs_utils::mesh::computational_mesh::*;
use nalgebra::Vector2;
use crate::finite_volume::{case::Case, config::CaseConfig};

#[derive(Clone, Debug, PartialEq)]
pub struct SimpleCase {
    pub name: String,
    pub step: usize,
    pub time: f64,
    pub time_step: f64,
    
    pub density: f64,
    pub kinematic_viscosity: f64,
    
    pub mesh: Computational2DMesh,
    
    pub speed: Vector2<Vec<f64>>,
    pub pressure: Vec<f64>,
    
    pub next_step: fn(&mut SimpleCase),
    
    /// Iterative solvers matrices and vectors
    pub a: nalgebra_sparse::CsrMatrix<f64>,
    pub b: nalgebra::DVector<f64>,
    pub x: nalgebra::DVector<f64>,
}

impl Case for SimpleCase {
    fn name(&self) -> &str {
        &self.name
    }
    
    fn step(&self) -> usize {
        self.step
    }
    
    fn time(&self) -> f64 {
        self.time
    }
    
    fn time_step(&self) -> f64 {
        self.time_step
    }
    
    fn export(&self, directory: &str) -> std::io::Result<()> {
        todo!()
    }
    
    fn import_from_file(file_name: &str) -> std::io::Result<()> {
        todo!()
    }
    
    fn next_step(&mut self) {
        (self.next_step)(self)
    }
}

pub fn setup(config: CaseConfig) -> SimpleCase {
    todo!()
}