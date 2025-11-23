use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};
use nalgebra::{DVector, Vector2};

use crate::finite_volume::{
    config::InitFunc,
    equation::variables::{ControlVolumeType, Variable},
};

use super::{
    boundary::BoundaryCondition,
    case::GradRequirements,
    config::CaseConfig,
    equation::Component,
    gradients::{update_grads, GradientConfig},
};

/// For now only support of scalar fields
#[derive(Debug, PartialEq, Clone)]
pub enum Field {
    Scalar(ScalarField),
    Vector2(Vector2<ScalarField>),
}

impl Field {}

/// The grads will only be allocated if necessary
#[derive(Debug, PartialEq, Clone)]
pub struct ScalarField {
    values: DVector<f64>,
    faces_values: DVector<f64>,
    grads_centers: DVector<Vector2<f64>>,
    grads_faces: DVector<Vector2<f64>>,
    gradients_up_to_date: bool,
    cvt: ControlVolumeType,
}

impl ScalarField {
    pub fn new<M: MeshCore>(
        mesh: &Mesh<M>,
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

        let cvt = variable.cvt();
        let n_values = match cvt {
            ControlVolumeType::Cells => mesh.cells.n,
            ControlVolumeType::Nodes => mesh.nodes.n,
        };
        let mut values = DVector::zeros(n_values);
        match cvt {
            ControlVolumeType::Cells => {
                for (i, centers) in mesh.cells.centers().iter().enumerate() {
                    values[i] = init(centers);
                }
            }
            ControlVolumeType::Nodes => {
                for (i, centers) in mesh.nodes.centers().iter().enumerate() {
                    values[i] = init(centers);
                }
            }
        }

        let faces_values = DVector::zeros(mesh.pairs.n);

        let grads_centers;
        if grads_required.cell() | grads_required.face() {
            grads_centers = DVector::zeros(n_values);
        } else {
            grads_centers = DVector::zeros(0);
        }

        let grads_faces;
        if grads_required.face() {
            grads_faces = DVector::zeros(mesh.pairs.n);
        } else {
            grads_faces = DVector::zeros(0);
        }

        // println!(
        //     "cells: {} | faces: {} | grad_cells: {} | grad_faces: {}",
        //     values.len(),
        //     face_values.len(),
        //     grads_cell.len(),
        //     grads_face.len()
        // );

        ScalarField {
            values,
            faces_values,
            grads_centers,
            grads_faces,
            gradients_up_to_date: false,
            cvt: cvt.clone(),
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

    pub fn faces_values(&self) -> &DVector<f64> {
        &self.faces_values
    }

    pub fn faces_values_mut(&mut self) -> &mut DVector<f64> {
        &mut self.faces_values
    }

    pub fn grads_centers(&self) -> &DVector<Vector2<f64>> {
        &self.grads_centers
    }

    /// Will not set the gradients as up to date
    pub fn grads_centers_mut(&mut self) -> &mut DVector<Vector2<f64>> {
        &mut self.grads_centers
    }

    pub fn grads_faces(&self) -> &DVector<Vector2<f64>> {
        &self.grads_faces
    }

    /// Will not set the gradients as up to date
    pub fn grads_faces_mut(&mut self) -> &mut DVector<Vector2<f64>> {
        &mut self.grads_faces
    }

    pub fn cvt(&self) -> &ControlVolumeType {
        &self.cvt
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
        &ControlVolumeType,
    ) {
        (
            &mut self.values,
            &mut self.faces_values,
            &mut self.grads_centers,
            &mut self.grads_faces,
            &self.cvt,
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
