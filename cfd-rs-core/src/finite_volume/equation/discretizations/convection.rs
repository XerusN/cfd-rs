use std::cell::RefCell;

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore, Patch};
use nalgebra::Vector2;

use crate::finite_volume::{
    boundary::BoundaryCondition,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{
        variables::ControlVolumeType, Component, EquationSolver, IntegrationCategory, Variable,
    },
    fields::Field,
};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum ConvectionScheme {
    UpwindSecondOrder,
}

impl ConvectionScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::UpwindSecondOrder => GradRequirements::new(true, true),
        }
    }

    pub fn discretize<M: MeshCore>(
        &self,
        var: &Variable,
        component: &Component,
        speed: &Variable,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
        equation_cvt: &ControlVolumeType,
    ) {
        let bc = config
            .bc
            .map
            .get(var)
            .expect("Missing boundary condition for field");
        let speed = find_var_in_fields(speed, fields);

        let var = find_var_in_fields(var, fields);

        match *self {
            Self::UpwindSecondOrder => upwind_second_order(
                var,
                component,
                speed,
                solver,
                mesh,
                bc,
                integration,
                coeff,
                equation_cvt,
            ),
        }
    }
}

/// Maybe wrong sign on rhs update
fn upwind_second_order<M: MeshCore>(
    field: &RefCell<Field>,
    component: &Component,
    speed: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Mesh<M>,
    boundary_conditions: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
    equation_cvt: &ControlVolumeType,
) {
    let speed = speed.borrow();
    let speed = match *speed {
        Field::Scalar(_) => panic!("Speed has to be a vector"),
        Field::Vector2(ref values) => values,
    };

    let (_, rhs) = solver.solver_borrow_mut();

    let field = field.borrow();
    let field = match *field {
        Field::Scalar(ref value) => value,
        Field::Vector2(ref value) => match *component {
            Component::X => &value.x,
            Component::Y => &value.y,
        },
    };
    let field_cvt = field.cvt();

    let pairs = &mesh.pairs;
    let (areas, normals) = match equation_cvt {
        ControlVolumeType::Cells => (pairs.cells_areas(), pairs.cells_normals()),
        ControlVolumeType::Nodes => (pairs.nodes_areas(), pairs.nodes_normals()),
    };
    let field_normals = match equation_cvt {
        ControlVolumeType::Cells => pairs.cells_normals(),
        ControlVolumeType::Nodes => pairs.nodes_normals(),
    };

    let field_cv_centers = match field_cvt {
        ControlVolumeType::Cells => mesh.cells.centers(),
        ControlVolumeType::Nodes => mesh.nodes.centers(),
    };

    let bnd = &mesh.boundaries;

    match integration {
        IntegrationCategory::Explicit => {
            for pair in 0..pairs.n {
                if pairs.on_bnd()[pair] {
                    continue;
                }

                let face_speed =
                    Vector2::new(speed.x.faces_values()[pair], speed.y.faces_values()[pair]);
                let flow_rate = areas[pair] * normals[pair].dot(&face_speed);

                let upwind_cv;
                // Compute Phi_f for the field control volume
                if field_normals[pair].dot(&face_speed) < 0. {
                    upwind_cv = match field_cvt {
                        ControlVolumeType::Cells => match pairs.neighboring_cells()[pair][1] {
                            Patch::Boundary(_) => panic!("Undected boundary"),
                            Patch::Cell(i) => i,
                        },
                        ControlVolumeType::Nodes => pairs.nodes()[pair][1],
                    }
                } else if field_normals[pair].dot(&face_speed) > 0. {
                    upwind_cv = match field_cvt {
                        ControlVolumeType::Cells => match pairs.neighboring_cells()[pair][0] {
                            Patch::Boundary(_) => panic!("Undetected boundary"),
                            Patch::Cell(i) => i,
                        },
                        ControlVolumeType::Nodes => pairs.nodes()[pair][0],
                    }
                } else {
                    continue;
                }

                let d_uf = pairs.centers()[pair] - field_cv_centers[upwind_cv];
                let face_value = field.values()[upwind_cv]
                    + (2. * field.grads_centers()[upwind_cv] - field.grads_faces()[pair])
                        .dot(&d_uf);
                let flux = flow_rate * face_value * coeff;

                // Sum over faces on equation control volume
                match equation_cvt {
                    ControlVolumeType::Cells => {
                        let cells = &pairs.neighboring_cells()[pair];
                        match cells[0] {
                            Patch::Cell(i_cell) => rhs[i_cell] -= flux,
                            Patch::Boundary(_) => panic!("Undetected boundary"),
                        }
                        match cells[1] {
                            Patch::Cell(i_cell) => rhs[i_cell] += flux,
                            Patch::Boundary(_) => panic!("Undetected boundary"),
                        }
                    }
                    ControlVolumeType::Nodes => {
                        let nodes = pairs.nodes()[pair];
                        rhs[nodes[0]] -= flux;
                        rhs[nodes[1]] += flux;
                    }
                }
            }

            for (i_bnd, bc) in boundary_conditions.iter().enumerate() {
                match field_cvt {
                    ControlVolumeType::Nodes => {
                        match bc {
                            BoundaryCondition::Dirichlet(_) => (), //ToCheck
                            BoundaryCondition::Neumann(_) => {
                                for &pair in &bnd.faces()[i_bnd] {
                                    let face_speed = Vector2::new(
                                        speed.x.faces_values()[pair],
                                        speed.y.faces_values()[pair],
                                    );
                                    let flow_rate = areas[pair] * normals[pair].dot(&face_speed);

                                    let upwind_cv;
                                    // Compute Phi_f for the field control volume
                                    if field_normals[pair].dot(&face_speed) < 0. {
                                        upwind_cv = pairs.nodes()[pair][1];
                                    } else if field_normals[pair].dot(&face_speed) > 0. {
                                        upwind_cv = pairs.nodes()[pair][1];
                                    } else {
                                        continue;
                                    }

                                    let d_uf = pairs.centers()[pair] - field_cv_centers[upwind_cv];
                                    let face_value = field.values()[upwind_cv]
                                        + (2. * field.grads_centers()[upwind_cv]
                                            - field.grads_faces()[pair])
                                            .dot(&d_uf);
                                    let flux = flow_rate * face_value * coeff;

                                    // Sum over faces on equation control volume
                                    match equation_cvt {
                                        ControlVolumeType::Cells => {
                                            let cells = &pairs.neighboring_cells()[pair];
                                            match cells[0] {
                                                Patch::Cell(i_cell) => rhs[i_cell] -= flux,
                                                Patch::Boundary(_) => (),
                                            }
                                            match cells[1] {
                                                Patch::Cell(i_cell) => rhs[i_cell] += flux,
                                                Patch::Boundary(_) => (),
                                            }
                                        }
                                        ControlVolumeType::Nodes => {
                                            let nodes = pairs.nodes()[pair];
                                            rhs[nodes[0]] -= flux;
                                            rhs[nodes[1]] += flux;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ControlVolumeType::Cells => match bc {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            for i in 0..bnd.faces()[i_bnd].len() {
                                let i_face = bnd.faces()[i_bnd][i];

                                let face_speed = Vector2::new(
                                    speed.x.faces_values()[i_face],
                                    speed.y.faces_values()[i_face],
                                );
                                let flow_rate = areas[i_face] * normals[i_face].dot(&face_speed);
                                let flux = flow_rate * bc_value * coeff;
                                match equation_cvt {
                                    ControlVolumeType::Cells => {
                                        let cells = &pairs.neighboring_cells()[i_face];
                                        match cells[0] {
                                            Patch::Cell(i_cell) => rhs[i_cell] -= flux,
                                            Patch::Boundary(_) => (),
                                        }
                                        match cells[1] {
                                            Patch::Cell(i_cell) => rhs[i_cell] += flux,
                                            Patch::Boundary(_) => (),
                                        }
                                    }
                                    ControlVolumeType::Nodes => {
                                        let nodes = pairs.nodes()[i_face];
                                        rhs[nodes[0]] -= flux;
                                        rhs[nodes[1]] += flux;
                                    }
                                }
                            }
                        }
                        BoundaryCondition::Neumann(_) => {
                            for i in 0..bnd.faces()[i_bnd].len() {
                                let i_face = bnd.faces()[i_bnd][i];
                                let i_cell = bnd.cells()[i_bnd][i];
                                let sign;
                                if let Patch::Cell(_) = pairs.neighboring_cells()[i_face][0] {
                                    sign = 1.
                                } else {
                                    sign = -1.;
                                }

                                let face_speed = Vector2::new(
                                    speed.x.faces_values()[i_face],
                                    speed.y.faces_values()[i_face],
                                );
                                let flow_rate = areas[i_face] * normals[i_face].dot(&face_speed);

                                let flux;
                                if field_normals[i_face].dot(&face_speed) * sign < 0. {
                                    let d_cf = pairs.centers()[i_face] - field_cv_centers[i_cell];
                                    let face_value = field.values()[i_cell]
                                        + field.grads_faces()[i_face].dot(&d_cf);
                                    flux = flow_rate * face_value * coeff;
                                } else if field_normals[i_face].dot(&face_speed) * sign > 0. {
                                    let d_uf = pairs.centers()[i_face] - field_cv_centers[i_cell];
                                    let face_value = field.values()[i_cell]
                                        + (2. * field.grads_centers()[i_cell]
                                            - field.grads_faces()[i_face])
                                            .dot(&d_uf);
                                    flux = flow_rate * face_value * coeff;
                                } else {
                                    continue;
                                }

                                match equation_cvt {
                                    ControlVolumeType::Cells => {
                                        let cells = &pairs.neighboring_cells()[i_face];
                                        match cells[0] {
                                            Patch::Cell(i_cell) => rhs[i_cell] -= flux,
                                            Patch::Boundary(_) => (),
                                        }
                                        match cells[1] {
                                            Patch::Cell(i_cell) => rhs[i_cell] += flux,
                                            Patch::Boundary(_) => (),
                                        }
                                    }
                                    ControlVolumeType::Nodes => {
                                        let nodes = pairs.nodes()[i_face];
                                        rhs[nodes[0]] -= flux;
                                        rhs[nodes[1]] += flux;
                                    }
                                }
                            }
                        }
                    },
                }
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}
