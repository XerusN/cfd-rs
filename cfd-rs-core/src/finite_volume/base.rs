use cfd_rs_utils::mesh::{
    assembled_mesh::{Mesh, MeshCore},
    computational_mesh::Computational2DMesh,
};
use nalgebra::{DVector, Point2, Vector2};

use crate::finite_volume::config::InitFunc;

use super::{
    boundary::BoundaryCondition,
    case::{CaseEquations, GradRequirements},
    config::{self, CaseConfig},
    equation::{Component, Variable},
    gradients::{update_grads, GradientConfig},
};

/// For now only support of scalar fields
#[derive(Debug, PartialEq, Clone)]
pub enum Field {
    Scalar(CellScalarField),
    Vector2(Vector2<CellScalarField>),
}

impl Field {}

/// The grads will only be allocated if necessary
#[derive(Debug, PartialEq, Clone)]
pub struct CellScalarField {
    values: DVector<f64>,
    face_values: DVector<f64>,
    grads_cell: DVector<Vector2<f64>>,
    grads_face: DVector<Vector2<f64>>,
    gradients_up_to_date: bool,
}

impl CellScalarField {
    pub fn new(
        mesh: &Computational2DMesh,
        grads_required: &GradRequirements,
        variable: &Variable,
        component: Component,
        config: &CaseConfig,
    ) -> Self {
        let init = match *config.initial_fields.get(variable).expect(&format![
            "No initialization defined for field {:?}",
            variable
        ]) {
            InitFunc::Scalar(func) => {
                if let Component::X = component {
                    func
                } else {
                    panic!("Trying to initialize scalar with vec function")
                }
            }
            InitFunc::Vector2(func_x, func_y) => match component {
                Component::X => func_x,
                Component::Y => func_y,
            },
        };
        let mut values = DVector::zeros(mesh.num_cells());
        for (i, cell) in mesh.cells().iter().enumerate() {
            values[i] = init(cell.centroid());
        }
        let face_values = DVector::zeros(mesh.num_faces());

        let grads_cell;
        if grads_required.cell() | grads_required.face() {
            grads_cell = DVector::zeros(mesh.num_cells());
        } else {
            grads_cell = DVector::zeros(0);
        }

        let grads_face;
        if grads_required.face() {
            grads_face = DVector::zeros(mesh.num_faces());
        } else {
            grads_face = DVector::zeros(0);
        }

        // println!(
        //     "cells: {} | faces: {} | grad_cells: {} | grad_faces: {}",
        //     values.len(),
        //     face_values.len(),
        //     grads_cell.len(),
        //     grads_face.len()
        // );

        CellScalarField {
            values,
            face_values,
            grads_cell,
            grads_face,
            gradients_up_to_date: false,
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

    pub fn face_values(&self) -> &DVector<f64> {
        &self.face_values
    }

    pub fn face_values_mut(&mut self) -> &mut DVector<f64> {
        &mut self.face_values
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

    pub fn get_deconstructed_field_mut(
        &mut self,
    ) -> (
        &mut DVector<f64>,
        &mut DVector<f64>,
        &mut DVector<Vector2<f64>>,
        &mut DVector<Vector2<f64>>,
    ) {
        (
            &mut self.values,
            &mut self.face_values,
            &mut self.grads_cell,
            &mut self.grads_face,
        )
    }
}

impl Field {
    /// Updates the gradients to match the values
    pub fn update_grads<M: MeshCore>(
        &mut self,
        grad_requirements: &GradRequirements,
        mesh: &Mesh<M>,
        config: &GradientConfig,
        bc: &Vec<BoundaryCondition>,
    ) {
        update_grads(self, grad_requirements, mesh, config, bc);
        unsafe {
            match self {
                Field::Scalar(field) => field.gradients_updated(),
                Field::Vector2(fields) => {
                    fields.x.gradients_updated();
                    fields.y.gradients_updated();
                }
            }
        }
    }
}
