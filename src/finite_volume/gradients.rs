use std::iter::Skip;

use super::{
    base::{CellScalarField, Field}, boundary::BoundaryCondition, case::{Case, GradRequirements}, equation::*, interpolations::GradientInterpConfig
};
use cfd_rs_utils::mesh::{computational_mesh::Computational2DMesh, indices::FaceIndex};
use nalgebra::Vector2;

#[derive(Clone, Debug, PartialEq)]
pub enum GradientScheme {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GradientConfig {
    pub scheme: GradientScheme,
    pub interp: GradientInterpConfig,
}

pub fn update_grads(field: &mut Field, grad_requirements: &GradRequirements, mesh: &Computational2DMesh, config: &GradientConfig, bc: &Vec<BoundaryCondition>) {
    match field {
        Field::Scalar(field) => update_grad_scalar(field, grad_requirements, mesh, config, bc),
        // Field::Vector2(field) => {
        //     update_grad_scalar(&mut field.x, mesh, config);
        //     update_grad_scalar(&mut field.y, mesh, config);
        // }
    }
}

fn update_grad_scalar(
    field: &mut CellScalarField,
    grad_requirements: &GradRequirements,
    mesh: &Computational2DMesh,
    config: &GradientConfig,
    bc: &Vec<BoundaryCondition>,
) {
    if !field.gradients_up_to_date() {
        if grad_requirements.cell() | grad_requirements.face() {
            match config.scheme {
                GradientScheme::GreenGaussCompact => green_gauss_compact(field, mesh),
                _ => unimplemented!("GradientScheme not implemented for {:?}", config.scheme),
            }
        }
        
        if grad_requirements.face() {
            match config.interp {
                _ => unimplemented!("GradientInterp not implemented for {:?}", config.interp),
            }
        }
    }
}

fn green_gauss_compact(field: &mut CellScalarField, mesh: &Computational2DMesh) {
    
    let (values, face_values, grads, face_grads) = field.get_deconstructed_field_mut();
    
    for (i, value) in face_values.iter_mut().enumerate() {
        let (id_1, id_2, g_c) = match mesh.geometric_weighting_factor(FaceIndex(i)) {
            None => continue,
            Some(value) => value,
        };
        *value = values[id_1.0]*g_c + values[id_2.0]*(1. - g_c);
    }
    for (i, grad) in grads.iter_mut().enumerate() {
        
    }
    
    
    
    
    todo!()
}