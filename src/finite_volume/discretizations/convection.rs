use std::cell::{RefCell, RefMut};

use cfd_rs_utils::mesh::{computational_mesh::{Computational2DMesh, Patch}, indices::CellIndex};
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
            },
            Self::CentralDifference => {
                central_difference(var, component, speed, solver, mesh, bc, integration, coeff)
            }
        }
    }
}

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
                let id_1 = match face.patches().0 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(id) => {
                        let id_2 = match face.patches().1 {
                            Patch::Cell(id) => id,
                            Patch::Boundary(_) => panic!("Face with two boundaries as neighbors"),
                        };

                        match &boundary_condition[id.0] {
                            BoundaryCondition::Dirichlet(bc_value) => {
                                rhs[id_2.0] += 0.;
                            }
                            // Check Neumann implementation
                            BoundaryCondition::Neumann(bc_value) => {
                                rhs[id_2.0] += 0.;
                            }
                        }
                        continue;
                    }
                };

                let id_2 = match face.patches().1 {
                    Patch::Cell(id) => id,
                    Patch::Boundary(id) => {
                        match &boundary_condition[id.0] {
                            BoundaryCondition::Dirichlet(bc_value) => {
                                rhs[id_1.0] += 0.;
                            }
                            BoundaryCondition::Neumann(bc_value) => {
                                rhs[id_1.0] += 0.;
                            }
                        }
                        continue;
                    }
                };
                rhs[id_1.0] += 0.;
                rhs[id_2.0] += 0.;
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
            println!("ok");
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
                    let face_field = field.values()[cell_id]
                        + field.grads_face()[faces_id[i].0]
                            .dot(&d_cf);
                    rhs[cell_id] -= coeff * flow_rate * face_field;
                }
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}
