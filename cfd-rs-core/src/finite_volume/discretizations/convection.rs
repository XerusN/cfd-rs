use std::cell::{RefCell, RefMut};

use cfd_rs_utils::mesh::{
    computational_mesh::{Computational2DMesh, Patch},
    indices::{CellIndex, FaceIndex},
};
use nalgebra::{Scalar, Vector2};

use crate::finite_volume::{
    base::Field,
    boundary::BoundaryCondition,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{Component, EquationSolver, IntegrationCategory, Variable},
};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum ConvectionScheme {
    UpwindSecondOrder,
    /// Not recommanded, very unstable
    CentralDifference,
}

impl ConvectionScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::UpwindSecondOrder => GradRequirements::new(true, true),
            Self::CentralDifference => GradRequirements::new(false, true),
        }
    }

    pub fn discretize(
        &self,
        var: &Variable,
        component: &Component,
        speed: &Variable,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
    ) {
        let bc = config
            .bc
            .map
            .get(var)
            .expect("Missing boundary condition for field");
        let speed = find_var_in_fields(speed, fields);

        let var = find_var_in_fields(var, fields);

        match *self {
            Self::UpwindSecondOrder => {
                upwind_second_order(var, component, speed, solver, mesh, bc, integration, coeff)
            }
            Self::CentralDifference => {
                central_difference(var, component, speed, solver, mesh, bc, integration, coeff)
            }
        }
    }
}

/// Maybe wrong sign on rhs update
fn upwind_second_order(
    field: &RefCell<Field>,
    component: &Component,
    speed: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    boundary_condition: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
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

    match integration {
        IntegrationCategory::Explicit => {
            for (face_id, face) in mesh.faces().iter().enumerate() {
                let face_speed = Vector2::new(
                    speed.x.face_values()[face_id],
                    speed.y.face_values()[face_id],
                );
                let flow_rate = face.area() * face.normal().dot(&face_speed);
                let id_1 = match face.patches().0 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(id) => {
                        let id_2 = match face.patches().1 {
                            Patch::Cell(id) => id,
                            Patch::Boundary(_) => panic!("Face with two boundaries as neighbors"),
                        };

                        let d_cf_2 = mesh.middle_point_from_face(FaceIndex(face_id))
                            - mesh.cells()[id_2.0].centroid();

                        if flow_rate > 0. {
                            match &boundary_condition[id.0] {
                                BoundaryCondition::Dirichlet(bc_value) => {
                                    let face_field = bc_value.get_value(component)
                                        + field.grads_face()[face_id].dot(&d_cf_2);
                                    rhs[id_2.0] += coeff * face_field * flow_rate;
                                }
                                // Check Neumann implementation
                                BoundaryCondition::Neumann(bc_value) => {
                                    let face_field = field.face_values()[face_id]
                                        + bc_value.get_value(component)
                                            * face.normal().dot(&d_cf_2);
                                    rhs[id_2.0] += coeff * face_field * flow_rate;
                                }
                            }
                        } else {
                            let face_field = field.values()[id_2.0]
                                + (2. * field.grads_cell()[id_2.0] - field.grads_face()[face_id])
                                    .dot(&d_cf_2);
                            rhs[id_2.0] += coeff * face_field * flow_rate;
                        }
                        continue;
                    }
                };

                let id_2 = match face.patches().1 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(id) => {
                        let d_cf_1 = mesh.middle_point_from_face(FaceIndex(face_id))
                            - mesh.cells()[id_1.0].centroid();

                        if flow_rate < 0. {
                            match &boundary_condition[id.0] {
                                BoundaryCondition::Dirichlet(bc_value) => {
                                    let face_field = bc_value.get_value(component)
                                        + field.grads_face()[face_id].dot(&d_cf_1);
                                    rhs[id_1.0] -= coeff * face_field * flow_rate;
                                }
                                // Check Neumann implementation
                                BoundaryCondition::Neumann(bc_value) => {
                                    let face_field = field.face_values()[face_id]
                                        + bc_value.get_value(component)
                                            * face.normal().dot(&d_cf_1);
                                    rhs[id_1.0] -= coeff * face_field * flow_rate;
                                }
                            }
                        } else {
                            let face_field = field.values()[id_1.0]
                                + (2. * field.grads_cell()[id_1.0] - field.grads_face()[face_id])
                                    .dot(&d_cf_1);
                            rhs[id_1.0] -= coeff * face_field * flow_rate;
                        }
                        continue;
                    }
                };
                let face_field;
                if flow_rate < 0. {
                    let d_cf_2 = mesh.middle_point_from_face(FaceIndex(face_id))
                        - mesh.cells()[id_2.0].centroid();
                    face_field = field.values()[id_2.0]
                        + (2. * field.grads_cell()[id_2.0] - field.grads_face()[face_id])
                            .dot(&d_cf_2);
                } else {
                    let d_cf_1 = mesh.middle_point_from_face(FaceIndex(face_id))
                        - mesh.cells()[id_1.0].centroid();
                    face_field = field.values()[id_1.0]
                        + (2. * field.grads_cell()[id_1.0] - field.grads_face()[face_id])
                            .dot(&d_cf_1);
                }
                rhs[id_1.0] -= coeff * face_field * flow_rate;
                rhs[id_2.0] += coeff * face_field * flow_rate;
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}

fn _upwind_second_order(
    field: &RefCell<Field>,
    component: &Component,
    speed: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    boundary_condition: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
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

    match integration {
        IntegrationCategory::Explicit => {
            for (cell_id, cell) in mesh.cells().iter().enumerate() {
                let mut f = 0.;
                let (faces_id, normals) =
                    mesh.normal_vectors_from_cell_with_faces_id(CellIndex(cell_id));
                for i in 0..faces_id.len() {
                    let face_speed = Vector2::new(
                        speed.x.face_values()[faces_id[i].0],
                        speed.y.face_values()[faces_id[i].0],
                    );
                    let flow_rate =
                        mesh.faces()[faces_id[i].0].area() * normals[i].dot(&face_speed);
                    //println!("{:?} {:?}", flow_rate, face_speed);
                    let d_cf = mesh.middle_point_from_face(faces_id[i]) - cell.centroid();
                    let face_field = field.values()[cell_id]
                        + (2. * field.grads_cell()[cell_id] - field.grads_face()[faces_id[i].0])
                            .dot(&d_cf);
                    rhs[cell_id] -= coeff * flow_rate * face_field;
                }
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}

fn central_difference(
    field: &RefCell<Field>,
    component: &Component,
    speed: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    boundary_condition: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
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

    match integration {
        IntegrationCategory::Explicit => {
            println!("ok");
            for (cell_id, cell) in mesh.cells().iter().enumerate() {
                let (faces_id, normals) =
                    mesh.normal_vectors_from_cell_with_faces_id(CellIndex(cell_id));
                for i in 0..faces_id.len() {
                    let face_speed = Vector2::new(
                        speed.x.face_values()[faces_id[i].0],
                        speed.y.face_values()[faces_id[i].0],
                    );
                    let flow_rate =
                        mesh.faces()[faces_id[i].0].area() * normals[i].dot(&face_speed);
                    //println!("{:?} {:?}", flow_rate, face_speed);
                    let d_cf = mesh.middle_point_from_face(faces_id[i]) - cell.centroid();
                    let face_field =
                        field.values()[cell_id] + field.grads_face()[faces_id[i].0].dot(&d_cf);
                    rhs[cell_id] -= coeff * flow_rate * face_field;
                }
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}
