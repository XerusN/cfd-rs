use super::{
    base::{CellScalarField, Field},
    boundary::BoundaryCondition,
    case::GradRequirements,
    equation::Component,
    interpolations::GradientInterpConfig,
};
use cfd_rs_utils::mesh::{
    computational_mesh::{Computational2DMesh, Patch},
    indices::{CellIndex, FaceIndex},
};

/// Recommanded value in book: 2
const GREEN_GAUSS_COMPACT_ITER: usize = 2;

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

pub fn update_grads(
    field: &mut Field,
    grad_requirements: &GradRequirements,
    mesh: &Computational2DMesh,
    config: &GradientConfig,
    bc: &Vec<BoundaryCondition>,
) {
    match field {
        Field::Scalar(field) => {
            update_grad_scalar(field, grad_requirements, mesh, config, bc, &Component::X)
        }
        Field::Vector2(field) => {
            update_grad_scalar(
                &mut field.x,
                grad_requirements,
                mesh,
                config,
                bc,
                &Component::X,
            );
            update_grad_scalar(
                &mut field.y,
                grad_requirements,
                mesh,
                config,
                bc,
                &Component::X,
            );
        }
    }
}

fn update_grad_scalar(
    field: &mut CellScalarField,
    grad_requirements: &GradRequirements,
    mesh: &Computational2DMesh,
    config: &GradientConfig,
    bc: &Vec<BoundaryCondition>,
    component: &Component,
) {
    if !field.gradients_up_to_date() {
        println!("Updating");
        if grad_requirements.cell() | grad_requirements.face() {
            match config.scheme {
                GradientScheme::GreenGaussCompact => {
                    green_gauss_compact(field, mesh, bc, component)
                }
                _ => todo!("GradientScheme not implemented for {:?}", config.scheme),
            }
        }

        if grad_requirements.face() {
            match config.interp {
                GradientInterpConfig::AveragedCorrected => {
                    averaged_corrected_interp(field, mesh, bc, component)
                }
                _ => todo!("GradientInterp not implemented for {:?}", config.interp),
            }
        }
    }
}

fn green_gauss_compact(
    field: &mut CellScalarField,
    mesh: &Computational2DMesh,
    bc: &Vec<BoundaryCondition>,
    component: &Component,
) {
    let (values, face_values, grads, _) = field.get_deconstructed_field_mut();

    for (i, value) in face_values.iter_mut().enumerate() {
        let (patch_1, patch_2, _) = mesh.geometric_weighting_factor(FaceIndex(i));
        //println!("{g_c}");
        let id_1 = match *patch_1 {
            Patch::Cell(id) => id,
            Patch::Boundary(id) => {
                let id_2 = match *patch_2 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(_) => panic!("Face with two boundaries as neighbors"),
                };

                match &bc[id.0] {
                    BoundaryCondition::Dirichlet(bc_value) => {
                        *value = bc_value.get_value(component);
                        // println!("{:?}", bc_value.get_value(component));
                    }
                    // Check Neumann implementation
                    BoundaryCondition::Neumann(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        *value = values[id_2.0]
                            + bc_value
                                * (mesh.middle_point_from_face(FaceIndex(i))
                                    - mesh.cells()[id_2.0].centroid())
                                .dot(
                                    &mesh.faces()[i]
                                        .normal_from_cell(id_2)
                                        .expect("Mesh not coherent"),
                                );
                        //println!("{:?} {:?}", value, values[id_2.0]);
                    }
                }

                continue;
            }
        };

        let id_2 = match *patch_2 {
            Patch::Cell(id) => id,
            Patch::Boundary(id) => {
                match &bc[id.0] {
                    BoundaryCondition::Dirichlet(bc_value) => {
                        *value = bc_value.get_value(component)
                    }
                    BoundaryCondition::Neumann(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        *value = values[id_1.0]
                            + bc_value
                                * (mesh.middle_point_from_face(FaceIndex(i))
                                    - mesh.cells()[id_1.0].centroid())
                                .dot(
                                    &mesh.faces()[i]
                                        .normal_from_cell(id_1)
                                        .expect("Mesh not coherent"),
                                )
                    }
                }
                continue;
            }
        };
        *value = (values[id_1.0] + values[id_2.0]) * 0.5;
    }

    for (cell, grad) in grads.iter_mut().enumerate() {
        let (faces_id, normals) = mesh.normals_from_cell_with_faces_id(CellIndex(cell));

        grad.x = 0.;
        grad.y = 0.;
        for (i, face_id) in faces_id.iter().enumerate() {
            *grad += face_values[face_id.0] * normals[i] * mesh.faces()[face_id.0].area();
        }
        *grad /= mesh.cells()[cell].volume();
    }

    let mut i = 0;
    for _ in 0..GREEN_GAUSS_COMPACT_ITER {
        let mut norm = 0.;
        for (i, value) in face_values.iter_mut().enumerate() {
            let old = value.clone();
            let (patch_1, patch_2, _) = mesh.geometric_weighting_factor(FaceIndex(i));

            let id_1 = match *patch_1 {
                Patch::Cell(id) => id,
                Patch::Boundary(id) => {
                    let id_2 = match *patch_2 {
                        Patch::Cell(id) => id,
                        Patch::Boundary(_) => panic!("Face with two boundaries as neighbors"),
                    };

                    match &bc[id.0] {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            *value = bc_value.get_value(component)
                        }
                        BoundaryCondition::Neumann(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            *value = values[id_2.0]
                                + bc_value
                                    * (mesh.middle_point_from_face(FaceIndex(i))
                                        - mesh.cells()[id_2.0].centroid())
                                    .dot(
                                        &mesh.faces()[i]
                                            .normal_from_cell(id_2)
                                            .expect("Mesh not coherent"),
                                    )
                        }
                    }
                    continue;
                }
            };
            let id_2 = match *patch_2 {
                Patch::Cell(id) => id,
                Patch::Boundary(id) => {
                    match &bc[id.0] {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            *value = bc_value.get_value(component)
                        }
                        BoundaryCondition::Neumann(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            *value = values[id_1.0]
                                + bc_value
                                    * (mesh.middle_point_from_face(FaceIndex(i))
                                        - mesh.cells()[id_1.0].centroid())
                                    .dot(
                                        &mesh.faces()[i]
                                            .normal_from_cell(id_1)
                                            .expect("Mesh not coherent"),
                                    )
                        }
                    }
                    continue;
                }
            };
            *value = (values[id_1.0] + values[id_2.0]) * 0.5;
            *value += 0.5
                * (grads[id_1.0] + grads[id_2.0]).dot(
                    &(mesh.middle_point_from_face(FaceIndex(i))
                        - mesh.cells()[id_1.0]
                            .centroid()
                            .lerp(mesh.cells()[id_2.0].centroid(), 0.5)),
                );
            norm += (*value - old).abs();
        }

        for (cell, grad) in grads.iter_mut().enumerate() {
            let (faces_id, normals) = mesh.normals_from_cell_with_faces_id(CellIndex(cell));

            grad.x = 0.;
            grad.y = 0.;
            for (i, face_id) in faces_id.iter().enumerate() {
                *grad += face_values[face_id.0] * normals[i] * mesh.faces()[face_id.0].area();
            }

            *grad /= mesh.cells()[cell].volume();
        }

        i += 1;

        norm /= mesh.num_faces() as f64;
        // println!("norm = {:.9e}", norm);

        if norm < 1e-3 {
            break;
        }
    }

    // println!("Gradients updated {i}");
}

fn averaged_corrected_interp(
    field: &mut CellScalarField,
    mesh: &Computational2DMesh,
    bc: &Vec<BoundaryCondition>,
    component: &Component,
) {
    let (values, _, grads, face_grads) = field.get_deconstructed_field_mut();

    for (i, face_grad) in face_grads.iter_mut().enumerate() {
        let (patch_1, patch_2, g_c) = mesh.geometric_weighting_factor(FaceIndex(i));

        let id_1 = match *patch_1 {
            Patch::Cell(id) => id,
            Patch::Boundary(id) => {
                let id_2 = match *patch_2 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(_) => panic!("Face with two boundaries as neighbors"),
                };

                // Check!!
                match &bc[id.0] {
                    BoundaryCondition::Dirichlet(bc_value) => {
                        let d = mesh.middle_point_from_face(FaceIndex(i))
                            - mesh.cells()[id_2.0].centroid();
                        let normal = mesh.faces()[i]
                            .normal_from_cell(id_2)
                            .expect("Mesh not coherent");

                        let bc_value = bc_value.get_value(component);
                        *face_grad = (bc_value - values[id_2.0]) / (d.dot(&normal)) * normal;
                    }
                    BoundaryCondition::Neumann(bc_value) => {
                        let normal = mesh.faces()[i]
                            .normal_from_cell(id_2)
                            .expect("Mesh not coherent");
                        let bc_value = bc_value.get_value(component);
                        *face_grad =
                            grads[id_2.0] + (bc_value - grads[id_2.0].dot(&normal)) * normal;
                    }
                }

                continue;
            }
        };
        let id_2 = match *patch_2 {
            Patch::Cell(id) => id,
            Patch::Boundary(id) => {
                match &bc[id.0] {
                    BoundaryCondition::Dirichlet(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        let d = mesh.middle_point_from_face(FaceIndex(i))
                            - mesh.cells()[id_1.0].centroid();
                        // Assumption on the way geometric weighting factor behaves (same patch order as face)
                        let normal = mesh.faces()[i]
                            .normal_from_cell(id_1)
                            .expect("Mesh not coherent");

                        *face_grad = (bc_value - values[id_1.0]) / (d.dot(&(-normal))) * &normal;
                    }
                    BoundaryCondition::Neumann(bc_value) => {
                        let normal = mesh.faces()[i]
                            .normal_from_cell(id_1)
                            .expect("Mesh not coherent");
                        let bc_value = bc_value.get_value(component);
                        
                        // SUSPICIOUS
                        *face_grad =
                            grads[id_1.0] + (bc_value - grads[id_1.0].dot(&normal)) * normal;
                    }
                }
                continue;
            }
        };
        let d_cf = mesh.cells()[id_2.0].centroid() - mesh.cells()[id_1.0].centroid();
        let d_cf_norm = d_cf.norm();
        let mean_grad = g_c * grads[id_1.0] + (1. - g_c) * grads[id_2.0];
        // To check
        *face_grad = mean_grad
            + ((values[id_2.0] - values[id_1.0]) - mean_grad.dot(&d_cf)) * d_cf / d_cf_norm.powi(2);
    }

    // println!("Gradients interpolated");
}
