use super::{
    base::{CellScalarField, Field}, boundary::BoundaryCondition, case::{GradRequirements}, interpolations::GradientInterpConfig
};
use cfd_rs_utils::mesh::{computational_mesh::{Computational2DMesh, Patch}, indices::{CellIndex, FaceIndex}};
use log::info;
use nalgebra::Vector2;

/// Arbitrary value, has to be checked
const GREEN_GAUSS_COMPACT_ITER: usize = 3;

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
                GradientScheme::GreenGaussCompact => green_gauss_compact(field, mesh, bc),
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

fn green_gauss_compact(field: &mut CellScalarField, mesh: &Computational2DMesh, bc: &Vec<BoundaryCondition>) {
    
    let (values, face_values, grads, face_grads) = field.get_deconstructed_field_mut();
    
    for (i, value) in face_values.iter_mut().enumerate() {
        let (patch_1, patch_2, g_c) = mesh.geometric_weighting_factor(FaceIndex(i));
        
        let id_1 = match *patch_1 {
            Patch::Cell(id) => id,
            Patch::Boundary(id) => {
                let id_2 = match *patch_2 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(_) => panic!("Face with two boundaries as neighbors")
                };
                
                match bc[id.0] {
                    BoundaryCondition::Dirichlet(bc_value) => *value = bc_value,
                    BoundaryCondition::Neumann(bc_value) => todo!(),
                }
                
                continue;
            },
        };
        let id_2 = match *patch_2 {
            Patch::Cell(id) => id,
            Patch::Boundary(id) => {
                match bc[id.0] {
                    BoundaryCondition::Dirichlet(bc_value) => *value = bc_value,
                    BoundaryCondition::Neumann(bc_value) => todo!(),
                }
                continue;
            },
        };
        *value = values[id_1.0]*g_c + values[id_2.0]*(1. - g_c);
    }
    
    for (cell, grad) in grads.iter_mut().enumerate() {
        let (faces_id, normals) = mesh.surface_vectors_from_cell_with_faces_id(CellIndex(cell));
        
        *grad = Vector2::zeros();
        for (i, face_id) in faces_id.iter().enumerate() {
            *grad += face_values[face_id.0]*normals[i];
        }
        *grad /= mesh.cells()[cell].volume();
    }
    
    for j in 0..GREEN_GAUSS_COMPACT_ITER {    
        for (i, value) in face_values.iter_mut().enumerate() {
            let (patch_1, patch_2, g_c) = mesh.geometric_weighting_factor(FaceIndex(i));
            
            let id_1 = match *patch_1 {
                Patch::Cell(id) => id,
                Patch::Boundary(id) => {
                    let id_2 = match *patch_2 {
                        Patch::Cell(id) => id,
                        Patch::Boundary(_) => panic!("Face with two boundaries as neighbors")
                    };
                    
                    match bc[id.0] {
                        BoundaryCondition::Dirichlet(bc_value) => *value = bc_value,
                        BoundaryCondition::Neumann(bc_value) => todo!(),
                    }
                    
                    continue;
                },
            };
            let id_2 = match *patch_2 {
                Patch::Cell(id) => id,
                Patch::Boundary(id) => {
                    match bc[id.0] {
                        BoundaryCondition::Dirichlet(bc_value) => *value = bc_value,
                        BoundaryCondition::Neumann(bc_value) => todo!(),
                    }
                    continue;
                },
            };
            *value += g_c*grads[id_1.0].dot(&(mesh.middle_point_from_face(FaceIndex(i)) - mesh.cells()[id_1.0].centroid())) + (1. - g_c)*grads[id_2.0].dot(&(mesh.middle_point_from_face(FaceIndex(i)) - mesh.cells()[id_2.0].centroid()));
        }
        
        for (cell, grad) in grads.iter_mut().enumerate() {
            let (faces_id, normals) = mesh.surface_vectors_from_cell_with_faces_id(CellIndex(cell));
            
            *grad = Vector2::zeros();
            for (i, face_id) in faces_id.iter().enumerate() {
                *grad += face_values[face_id.0]*normals[i];
            }
            *grad /= mesh.cells()[cell].volume();
        }
    }
    
    println!("Gradiennts updated");
}