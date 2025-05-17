use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::{DVector, Vector2};

use super::{case::GradRequirements, gradients::{update_grads, GradientConfig}};

/// For now only support of scalar fields
#[derive(Debug, PartialEq, Clone)]
pub enum Field {
    Scalar(CellScalarField),
    //Vector2(Vector2<CellScalarField>),
}

/// The grads will only be allocated if necessary
#[derive(Debug, PartialEq, Clone)]
pub struct CellScalarField {
    values: DVector<f64>,
    grads_cell: DVector<Vector2<f64>>,
    grads_face: DVector<Vector2<f64>>,
    gradients_up_to_date: bool,
}

impl CellScalarField {
    pub fn new(
        cells_num: usize,
        faces_num: usize,
        grads_required: &GradRequirements,
    ) -> Self {
        let values = DVector::zeros(cells_num);

        let grads_cell;
        if grads_required.cell() {
            grads_cell = DVector::zeros(cells_num);
        } else {
            grads_cell = DVector::zeros(0);
        }

        let grads_face;
        if grads_required.face() {
            grads_face = DVector::zeros(faces_num);
        } else {
            grads_face = DVector::zeros(0);
        }

        CellScalarField {
            values,
            grads_cell,
            grads_face,
            gradients_up_to_date: true,
        }
    }

    pub fn values(&self) -> &DVector<f64> {
        &self.values
    }

    /// Borrowing mutably the values will set the gradients as not up to date no matter what you do with the values.
    pub fn values_mut(&mut self) -> &mut DVector<f64> {
        self.gradients_up_to_date = false;
        &mut self.values
    }

    pub fn grads_cell(&self) -> &DVector<Vector2<f64>> {
        &self.grads_cell
    }

    /// Will not set the gradients as up to date
    pub fn grads_cell_mut(&mut self) -> &mut DVector<Vector2<f64>> {
        &mut self.grads_cell
    }

    pub fn grads_face(&self) -> &DVector<Vector2<f64>> {
        &self.grads_face
    }

    /// Will not set the gradients as up to date
    pub fn grads_face_mut(&mut self) -> &mut DVector<Vector2<f64>> {
        &mut self.grads_face
    }

    pub fn gradients_up_to_date(&self) -> bool {
        self.gradients_up_to_date
    }

    /// Will set the gradients as updated with no check
    pub unsafe fn gradients_updated(&mut self) {
        self.gradients_up_to_date = true
    }
}

impl Field {
    /// Updates the gradients to match the values
    pub fn update_grads(&mut self, mesh: &Computational2DMesh, config: &GradientConfig) {
        update_grads(self, mesh, config);
        unsafe {
            match self {
                Field::Scalar(field) => field.gradients_updated(),
                // Field::Vector2(fields) => {
                //     fields.x.gradients_updated();
                //     fields.y.gradients_updated();
                // }
            }
        }
    }
}
