use crate::finite_volume::equation::variables::ControlVolumeType;

use super::{
    boundary::BoundaryCondition,
    case::GradRequirements,
    equation::Component,
    fields::{Field, ScalarField},
};
use cfd_rs_utils::mesh::{
    assembled_mesh::{Mesh, MeshCore, Patch},
    indices::{CellIndex, FaceIndex},
};
use nalgebra::Vector2;

/// Recommanded value in book: 2
const GREEN_GAUSS_COMPACT_ITER: usize = 2;

#[derive(Clone, Debug, PartialEq)]
pub enum GradientScheme {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub enum GradientInterpConfig {
    Averaged,
    #[default]
    AveragedCorrected,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GradientConfig {
    pub scheme: GradientScheme,
    pub interp: GradientInterpConfig,
}

pub fn update_grads<M: MeshCore>(
    field: &mut Field,
    grad_requirements: &GradRequirements,
    mesh: &Mesh<M>,
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
                &Component::Y,
            );
        }
    }
}

fn update_grad_scalar<M: MeshCore>(
    field: &mut ScalarField,
    grad_requirements: &GradRequirements,
    mesh: &Mesh<M>,
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

fn green_gauss_compact<M: MeshCore>(
    field: &mut ScalarField,
    mesh: &Mesh<M>,
    boundary_conditions: &Vec<BoundaryCondition>,
    component: &Component,
) {
    let (values, face_values, grads, _, cvt) = field.get_deconstructed_field_mut();
    
    let pairs = &mesh.pairs;
    let bnd = &mesh.boundaries;
    let cv_centers = match cvt {
        ControlVolumeType::Cells => mesh.cells.centers(),
        ControlVolumeType::Nodes => mesh.nodes.centers(),
    };
    
    for grad in grads.iter_mut() {
        grad.x = 0.;
        grad.y = 0.;
    }
    
    let volumes = match cvt {
        ControlVolumeType::Cells => mesh.cells.volumes(),
        ControlVolumeType::Nodes => mesh.nodes.volumes(),
    };
    
    for i_iter in 0..GREEN_GAUSS_COMPACT_ITER {
        for pair in 0..pairs.n {
            if pairs.on_bnd()[pair] {
                continue;
            }
            
            let i_cv = match cvt {
                ControlVolumeType::Cells => {
                    let patches = &pairs.neighboring_cells()[pair];
                    [
                        {
                            if let Patch::Cell(i_cell) = patches[0] { i_cell } else {panic!()}
                        },
                        {
                            if let Patch::Cell(i_cell) = patches[1] { i_cell } else {panic!()}
                        },
                    ]
                },
                ControlVolumeType::Nodes => pairs.nodes()[pair],
            };
            
            
            
            face_values[pair] = (values[i_cv[0]] + values[i_cv[1]]) * 0.5;
            face_values[pair] += 0.5
                * (grads[i_cv[0]] + grads[i_cv[1]]).dot(
                    &(pairs.centers()[pair]
                        - cv_centers[i_cv[0]]
                            .lerp(&cv_centers[i_cv[1]], 0.5)),
            );
            
            if i_iter != 0 {
                
            }
        }
        
        for (i_bnd, bc) in boundary_conditions.iter().enumerate() {
            match cvt {
                ControlVolumeType::Nodes => {
                    match bc {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for &pair in &bnd.faces()[i_bnd] {
                                face_values[pair] = bc_value;
                            }
                        },
                        BoundaryCondition::Neumann(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for &pair in &bnd.faces()[i_bnd] {
                                let nodes = pairs.nodes()[pair];
                                let normal = pairs.nodes_normals()[pair];
                                let normal_correction_angle = normal.angle(&pairs.cells_normals()[pair]).abs();
                                face_values[pair] = (values[nodes[0]] + values[nodes[1]]) * 0.5 * (normal_correction_angle.sin() + normal_correction_angle.cos() * bc_value * pairs.lengths()[pair]);
                            }
                        },
                    }
                },
                ControlVolumeType::Cells => {
                    match bc {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for &pair in &bnd.faces()[i_bnd] {
                                face_values[pair] = bc_value;
                            }
                        },
                        BoundaryCondition::Neumann(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for i in 0..bnd.faces()[i_bnd].len() {
                                let i_face = bnd.faces()[i_bnd][i];
                                let i_cell = bnd.faces()[i_bnd][i];
                                match pairs.neighboring_cells()[i_face][0] {
                                    Patch::Cell(_) => {
                                        // Check sign in .dot()
                                        face_values[i_face] = values[i_cell] + bc_value * (pairs.centers()[i_face] - mesh.cells.centers()[i_cell]).dot(&(- pairs.cells_normals()[i_face]))
                                    }
                                    Patch::Boundary(_) => {
                                        face_values[i_face] = values[i_cell] + bc_value * (pairs.centers()[i_face] - mesh.cells.centers()[i_cell]).dot(&pairs.cells_normals()[i_face])
                                    },
                                }
                            }
                        },
                    }
                }
            }
        }
        
        for grad in grads {
            grad.x = 0.;
            grad.y = 0.;
        }
        
        for pair in 0..pairs.n {
            let area = match cvt {
                ControlVolumeType::Cells => pairs.cells_areas()[pair],
                ControlVolumeType::Nodes => pairs.nodes_areas()[pair],
            };
            let normal = match cvt {
                ControlVolumeType::Cells => pairs.cells_normals()[pair],
                ControlVolumeType::Nodes => pairs.nodes_normals()[pair],
            };
            
            
            match cvt {
                ControlVolumeType::Cells => {
                    let patches = pairs.neighboring_cells()[pair];
                    let flux = area*normal*face_values[pair];
                    match patches[0] {
                        Patch::Boundary(_) => (),
                        Patch::Cell(cell) => {
                            grads[cell] = flux/volumes[cell];
                        },
                    };
                    match patches[1] {
                        Patch::Boundary(_) => (),
                        Patch::Cell(cell) => {
                            grads[cell] = - flux/volumes[cell];
                        },
                    };
                },
                ControlVolumeType::Nodes => {
                    let nodes = pairs.nodes()[pair];
                    let flux = area*normal*face_values[pair];
                    grads[nodes[0]] = flux/volumes[nodes[0]];   //ToCheck sign
                    grads[nodes[1]] = - flux/volumes[nodes[1]];
                },
            }
        }
    }
}

fn averaged_corrected_interp<M: MeshCore>(
    field: &mut ScalarField,
    mesh: &Mesh<M>,
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
                        let tangent = Vector2::new(-normal.y, normal.x);
                        let bc_value = bc_value.get_value(component);
                        *face_grad = grads[id_2.0].dot(&tangent) * tangent + bc_value * normal;
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
                        let tangent = Vector2::new(-normal.y, normal.x);
                        let bc_value = bc_value.get_value(component);

                        // SUSPICIOUS
                        *face_grad = grads[id_1.0].dot(&tangent) * tangent + bc_value * normal;
                    }
                }
                continue;
            }
        };
        let d_cf = mesh.cells()[id_2.0].centroid() - mesh.cells()[id_1.0].centroid();
        let d_cf_norm_2 = d_cf.norm_squared();
        let mean_grad = g_c * grads[id_1.0] + (1. - g_c) * grads[id_2.0];
        // To check
        *face_grad = mean_grad
            + ((values[id_2.0] - values[id_1.0]) - mean_grad.dot(&d_cf)) * d_cf / d_cf_norm_2;
    }

    // println!("Gradients interpolated");
}
