use crate::finite_volume::equation::variables::ControlVolumeType;

use super::{
    boundary::BoundaryCondition,
    case::GradRequirements,
    equation::Component,
    fields::{Field, ScalarField},
};
use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore, Patch};

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
                            if let Patch::Cell(i_cell) = patches[0] {
                                i_cell
                            } else {
                                panic!()
                            }
                        },
                        {
                            if let Patch::Cell(i_cell) = patches[1] {
                                i_cell
                            } else {
                                panic!()
                            }
                        },
                    ]
                }
                ControlVolumeType::Nodes => pairs.nodes()[pair],
            };

            face_values[pair] = (values[i_cv[0]] + values[i_cv[1]]) * 0.5;
            face_values[pair] += 0.5
                * (grads[i_cv[0]] + grads[i_cv[1]]).dot(
                    &(pairs.centers()[pair] - cv_centers[i_cv[0]].lerp(&cv_centers[i_cv[1]], 0.5)),
                );
        }

        for (i_bnd, bc) in boundary_conditions.iter().enumerate() {
            match cvt {
                ControlVolumeType::Nodes => match bc {
                    BoundaryCondition::Dirichlet(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        for &pair in &bnd.faces()[i_bnd] {
                            face_values[pair] = bc_value;
                            values[pairs.nodes()[pair][0]] = bc_value;
                            values[pairs.nodes()[pair][1]] = bc_value;
                        }
                    }
                    BoundaryCondition::Neumann(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        for &pair in &bnd.faces()[i_bnd] {
                            let nodes = pairs.nodes()[pair];
                            let normal = pairs.nodes_normals()[pair];
                            let normal_correction_angle =
                                normal.angle(&pairs.cells_normals()[pair]).abs();
                            face_values[pair] = (values[nodes[0]] + values[nodes[1]])
                                * 0.5
                                * (normal_correction_angle.sin()
                                    + normal_correction_angle.cos()
                                        * bc_value
                                        * pairs.lengths()[pair]);
                        }
                    }
                },
                ControlVolumeType::Cells => {
                    match bc {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for &pair in &bnd.faces()[i_bnd] {
                                face_values[pair] = bc_value;
                            }
                        }
                        BoundaryCondition::Neumann(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for i in 0..bnd.faces()[i_bnd].len() {
                                let i_face = bnd.faces()[i_bnd][i];
                                let i_cell = bnd.faces()[i_bnd][i];
                                match pairs.neighboring_cells()[i_face][0] {
                                    Patch::Cell(_) => {
                                        // Check sign in .dot()
                                        face_values[i_face] = values[i_cell]
                                            + bc_value
                                                * (pairs.centers()[i_face]
                                                    - mesh.cells.centers()[i_cell])
                                                    .dot(&(-pairs.cells_normals()[i_face]))
                                    }
                                    Patch::Boundary(_) => {
                                        face_values[i_face] = values[i_cell]
                                            + bc_value
                                                * (pairs.centers()[i_face]
                                                    - mesh.cells.centers()[i_cell])
                                                    .dot(&pairs.cells_normals()[i_face])
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for grad in grads.iter_mut() {
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
                    let patches = &pairs.neighboring_cells()[pair];
                    let flux = area * normal * face_values[pair];
                    match patches[0] {
                        Patch::Boundary(_) => (),
                        Patch::Cell(cell) => {
                            grads[cell] = flux / volumes[cell];
                        }
                    };
                    match patches[1] {
                        Patch::Boundary(_) => (),
                        Patch::Cell(cell) => {
                            grads[cell] = -flux / volumes[cell];
                        }
                    };
                }
                ControlVolumeType::Nodes => {
                    let nodes = pairs.nodes()[pair];
                    let flux = area * normal * face_values[pair];
                    grads[nodes[0]] = flux / volumes[nodes[0]]; //ToCheck sign
                    grads[nodes[1]] = -flux / volumes[nodes[1]];
                }
            }
        }

        // ToDo add convergence check
    }
}

fn averaged_corrected_interp<M: MeshCore>(
    field: &mut ScalarField,
    mesh: &Mesh<M>,
    boundary_conditions: &Vec<BoundaryCondition>,
    component: &Component,
) {
    let (values, _, grads, face_grads, cvt) = field.get_deconstructed_field_mut();

    let pairs = &mesh.pairs;
    let cv_centers = match cvt {
        ControlVolumeType::Cells => mesh.cells.centers(),
        ControlVolumeType::Nodes => mesh.nodes.centers(),
    };
    let bnd = &mesh.boundaries;

    for pair in 0..pairs.n {
        if pairs.on_bnd()[pair] {
            continue;
        }

        let i_cv = match cvt {
            ControlVolumeType::Cells => {
                let patches = &pairs.neighboring_cells()[pair];
                [
                    {
                        if let Patch::Cell(i_cell) = patches[0] {
                            i_cell
                        } else {
                            panic!()
                        }
                    },
                    {
                        if let Patch::Cell(i_cell) = patches[1] {
                            i_cell
                        } else {
                            panic!()
                        }
                    },
                ]
            }
            ControlVolumeType::Nodes => pairs.nodes()[pair],
        };

        let d_cf = cv_centers[i_cv[1]] - cv_centers[i_cv[0]];
        let d_cf_norm = d_cf.norm();
        let e_cf = d_cf / d_cf_norm;
        let geometric_weighting = (pairs.centers()[pair] - cv_centers[i_cv[1]]).norm() / d_cf_norm;

        face_grads[pair] =
            grads[i_cv[0]] * geometric_weighting + grads[i_cv[1]] * (1. - geometric_weighting);
        face_grads[pair] = face_grads[pair]
            + e_cf
                * ((values[i_cv[1]] - values[i_cv[0]]) / d_cf_norm - face_grads[pair].dot(&e_cf));
    }

    for (i_bnd, bc) in boundary_conditions.iter().enumerate() {
        match cvt {
            ControlVolumeType::Nodes => {
                match bc {
                    // ToCheck : not using the Dirichlet bc
                    BoundaryCondition::Dirichlet(_) => {
                        for &pair in &bnd.faces()[i_bnd] {
                            let i_cv = pairs.nodes()[pair];
                            let d_cf = cv_centers[i_cv[1]] - cv_centers[i_cv[0]];
                            let d_cf_norm = d_cf.norm();
                            let e_cf = d_cf / d_cf_norm;
                            let geometric_weighting =
                                (pairs.centers()[pair] - cv_centers[i_cv[1]]).norm() / d_cf_norm;

                            face_grads[pair] = grads[i_cv[0]] * geometric_weighting
                                + grads[i_cv[1]] * (1. - geometric_weighting);
                            face_grads[pair] = face_grads[pair]
                                + e_cf
                                    * ((values[i_cv[1]] - values[i_cv[0]]) / d_cf_norm
                                        - face_grads[pair].dot(&e_cf));
                        }
                    }
                    BoundaryCondition::Neumann(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        for &pair in &bnd.faces()[i_bnd] {
                            let i_cv = pairs.nodes()[pair];
                            let d_cf = cv_centers[i_cv[1]] - cv_centers[i_cv[0]];
                            let d_cf_norm = d_cf.norm();
                            let e_cf = d_cf / d_cf_norm;
                            let geometric_weighting =
                                (pairs.centers()[pair] - cv_centers[i_cv[1]]).norm() / d_cf_norm;

                            face_grads[pair] = grads[i_cv[0]] * geometric_weighting
                                + grads[i_cv[1]] * (1. - geometric_weighting);
                            face_grads[pair] = face_grads[pair]
                                + e_cf
                                    * ((values[i_cv[1]] - values[i_cv[0]]) / d_cf_norm
                                        - face_grads[pair].dot(&e_cf));
                            let cell_normal = match pairs.neighboring_cells()[pair][0] {
                                Patch::Cell(_) => pairs.cells_normals()[pair],
                                Patch::Boundary(_) => -pairs.cells_normals()[pair],
                            };
                            face_grads[pair] = -bc_value * cell_normal
                                + face_grads[pair].dot(&pairs.nodes_normals()[pair])
                                    * pairs.nodes_normals()[pair];
                        }
                    }
                }
            }
            ControlVolumeType::Cells => match bc {
                BoundaryCondition::Dirichlet(bc_value) => {
                    let bc_value = bc_value.get_value(component);
                    for i in 0..bnd.faces()[i_bnd].len() {
                        let i_face = bnd.faces()[i_bnd][i];
                        let i_cell = bnd.cells()[i_bnd][i];
                        let d = pairs.centers()[i_face] - cv_centers[i_cell];
                        let normal = match pairs.neighboring_cells()[i_face][0] {
                            Patch::Cell(_) => pairs.cells_normals()[i_face],
                            Patch::Boundary(_) => -pairs.cells_normals()[i_face],
                        };
                        face_grads[i_face] =
                            (bc_value - values[i_cell]) / (d.dot(&(-normal))) * normal;
                    }
                }
                BoundaryCondition::Neumann(bc_value) => {
                    let bc_value = bc_value.get_value(component);
                    for i in 0..bnd.faces()[i_bnd].len() {
                        let i_face = bnd.faces()[i_bnd][i];
                        let i_cell = bnd.cells()[i_bnd][i];
                        let normal = match pairs.neighboring_cells()[i_face][0] {
                            Patch::Cell(_) => pairs.cells_normals()[i_face],
                            Patch::Boundary(_) => -pairs.cells_normals()[i_face],
                        };
                        face_grads[i_face] = -bc_value * normal
                            + grads[i_cell].dot(&pairs.nodes_normals()[i_face])
                                * pairs.nodes_normals()[i_face];
                    }
                }
            },
        }
    }
    // println!("Gradients interpolated");
}
