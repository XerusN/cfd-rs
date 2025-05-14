use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::Vector2;

use super::gradients::GradientConfig;

#[derive(Debug, PartialEq, Clone)]
pub enum Field {
    CellScalarField(CellScalarField),
    CellVectorField(Vector2<CellScalarField>),
}

/// The grads will only be allocated if necessary
#[derive(Debug, PartialEq, Clone)]
pub struct CellScalarField {
    values: Vec<f64>,
    grads_cell: Vec<Vector2<f64>>,
    grads_face: Vec<Vector2<f64>>,
    gradients_up_to_date: bool,
}

impl CellScalarField {
    
    pub fn new(cells_num: usize, faces_num: usize, grads_cell_required: bool, grads_face_required: bool) -> Self {
        let values = vec![0.; cells_num];
        
        let grads_cell;
        if grads_cell_required {
            grads_cell = vec![Vector2::new(0., 0.); cells_num];
        } else {
            grads_cell = vec![];
        }
        
        let grads_face;
        if grads_face_required {
            grads_face = vec![Vector2::new(0., 0.); faces_num];
        } else {
            grads_face = vec![];
        }
        
        CellScalarField { values, grads_cell, grads_face, gradients_up_to_date: true }
    }
    
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    
    /// Borrowing mutably the values will set the gradients as not up to date no matter what you do with the values.
    pub fn values_mut(&mut self) -> &mut [f64] {
        self.gradients_up_to_date = false;
        &mut self.values
    }
    
    pub fn grads_cell(&self) -> &[Vector2<f64>] {
        &self.grads_cell
    }
    
    /// Will not set the gradients as up to date
    pub fn grads_cell_mut(&mut self) -> &mut [Vector2<f64>] {
        &mut self.grads_cell
    }
    
    pub fn grads_face(&self) -> &[Vector2<f64>] {
        &self.grads_face
    }
    
    /// Will not set the gradients as up to date
    pub fn grads_face_mut(&mut self) -> &mut [Vector2<f64>] {
        &mut self.grads_face
    }
    
    pub fn gradients_up_to_date(&self) -> bool {
        self.gradients_up_to_date
    }
    
    /// Updates the gradients to match the values
    pub fn update_grads(&mut self, mesh: Computational2DMesh, config: GradientConfig) {
        todo!();
        self.gradients_up_to_date = true;
    }
    
    /// Will set the gradients as updated with no check
    pub unsafe fn gradients_updated(&mut self) {
        self.gradients_up_to_date = true
    }
}