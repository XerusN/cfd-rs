use std::{cell::RefCell, ops::Deref};

use cfd_rs_utils::mesh::{
    computational_mesh::{Computational2DMesh, Patch},
    indices::CellIndex,
};
use nalgebra_sparse::SparseEntryMut;

use crate::finite_volume::{
    base::Field,
    boundary::BoundaryCondition,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{
        variables::ControlVolume, Component, EquationSolver, IntegrationCategory, Variable,
    },
};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum LaplacianScheme {
    OrthogonalCorrection,
    MinimalCorrection,
}

impl LaplacianScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::OrthogonalCorrection => GradRequirements::new(false, true),
            Self::MinimalCorrection => GradRequirements::new(false, true),
        }
    }

    pub fn discretize(
        &self,
        var: &Variable,
        component: &Component,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
        equation_cv: &ControlVolume,
    ) {
        let field = find_var_in_fields(var, fields);
        let bc = config
            .bc
            .map
            .get(var)
            .expect("Missing boundary condition for field");
        match *self {
            Self::OrthogonalCorrection => orthogonal_correction(
                component,
                solver,
                field,
                mesh,
                bc,
                integration,
                coeff,
                equation_cv,
            ),
            Self::MinimalCorrection => minimal_correction(
                component,
                solver,
                field,
                mesh,
                bc,
                integration,
                coeff,
                equation_cv,
            ),
        }
    }
}

fn orthogonal_correction(
    component: &Component,
    solver: &mut EquationSolver,
    field: &RefCell<Field>,
    mesh: &Computational2DMesh,
    boundary_condition: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
    equation_cv: &ControlVolume,
) {
    let (matrix, rhs) = solver.solver_borrow_mut();
    let field = field.borrow();
    let field = match field.deref() {
        Field::Scalar(value) => value,
        Field::Vector2(value) => match *component {
            Component::X => &value.x,
            Component::Y => &value.y,
        },
    };

    match *integration {
        IntegrationCategory::Implicit => {
            for cell in 0..rhs.len() {
                let neighbors_and_faces = mesh.neighboring_patches_and_faces(CellIndex(cell));
                let mut row = matrix
                    .get_row_mut(cell)
                    .expect("Bad Initialization of matrix");
                let mut f_c = 0.;
                for (neighbor, face, face_id) in neighbors_and_faces {
                    match *neighbor {
                        Patch::Cell(id) => {
                            let normal = face
                                .normal_from_cell(CellIndex(cell))
                                .expect("Incoherence in face and cell connection");
                            let d_cf =
                                mesh.cells()[id.0].centroid() - mesh.cells()[cell].centroid();
                            let e_f = face.area() * d_cf.normalize();
                            let f_f = -face.area() / d_cf.magnitude();
                            f_c -= f_f;
                            match row
                                .get_entry_mut(id.0)
                                .expect("Bad Initialization of matrix")
                            {
                                SparseEntryMut::NonZero(value) => *value += f_f * coeff,
                                SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                            }
                            let t_f = face.area() * normal - e_f;
                            println!("{:?} {:?}", e_f, t_f);
                            rhs[cell] += coeff * field.grads_face()[face_id.0].dot(&t_f);
                        }
                        Patch::Boundary(id) => match &boundary_condition[id.0] {
                            BoundaryCondition::Dirichlet(bc_value) => {
                                let d_cb = face.middle_point(mesh.vertices())
                                    - mesh.cells()[cell].centroid();
                                let e_b = face.area() * d_cb.normalize();
                                let f_b = face.area() / d_cb.magnitude();
                                f_c += f_b;
                                let t_b = face.area()
                                    * face
                                        .normal_from_cell(CellIndex(cell))
                                        .expect("Incoherence in face and cell connection")
                                    - e_b;
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] += coeff
                                    * (f_b * bc_value + field.grads_face()[face_id.0].dot(&t_b));
                            }
                            BoundaryCondition::Neumann(bc_value) => {
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] -= coeff * bc_value * face.area();
                            }
                        },
                    }
                }

                match row
                    .get_entry_mut(cell)
                    .expect("Bad Initialization of matrix")
                {
                    SparseEntryMut::NonZero(value) => *value += f_c * coeff,
                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                }

                // println!("{:?}", row)
            }
            // panic!("Test")
        }
        IntegrationCategory::Explicit => {
            for cell in 0..rhs.len() {
                let neighbors_and_faces = mesh.neighboring_patches_and_faces(CellIndex(cell));
                let mut f_c = 0.;
                for (neighbor, face, face_id) in neighbors_and_faces {
                    match *neighbor {
                        Patch::Cell(id) => {
                            let d_cf =
                                mesh.cells()[id.0].centroid() - mesh.cells()[cell].centroid();
                            let e_f = face.area() * d_cf.normalize();
                            let f_f = -e_f.magnitude() / d_cf.magnitude();
                            f_c -= f_f;
                            rhs[cell] -= f_f * coeff * field.values()[id.0];
                            let t_f = face.area()
                                * face
                                    .normal_from_cell(CellIndex(cell))
                                    .expect("Incoherence in face and cell connection")
                                - e_f;
                            rhs[cell] += coeff * field.grads_face()[face_id.0].dot(&t_f);
                        }
                        Patch::Boundary(id) => match &boundary_condition[id.0] {
                            BoundaryCondition::Dirichlet(bc_value) => {
                                let d_cb = face.middle_point(mesh.vertices())
                                    - mesh.cells()[cell].centroid();
                                let e_b = face.area() * d_cb.normalize();
                                let f_b = e_b.magnitude() / d_cb.magnitude();
                                f_c += f_b;
                                let t_b = face.area()
                                    * face
                                        .normal_from_cell(CellIndex(cell))
                                        .expect("Incoherence in face and cell connection")
                                    - e_b;
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] += coeff
                                    * (f_b * bc_value + field.grads_face()[face_id.0].dot(&t_b));
                            }
                            BoundaryCondition::Neumann(bc_value) => {
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] -= coeff * bc_value * face.area();
                            }
                        },
                    }
                }
                rhs[cell] -= f_c * coeff * field.values()[cell];
            }
        }
    }
}

// Not prod ready
fn minimal_correction(
    component: &Component,
    solver: &mut EquationSolver,
    field: &RefCell<Field>,
    mesh: &Computational2DMesh,
    boundary_condition: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
    equation_cv: &ControlVolume,
) {
    let (matrix, rhs) = solver.solver_borrow_mut();
    let field = field.borrow();
    let field = match field.deref() {
        Field::Scalar(value) => value,
        Field::Vector2(value) => match *component {
            Component::X => &value.x,
            Component::Y => &value.y,
        },
    };

    match *integration {
        IntegrationCategory::Implicit => {
            for cell in 0..rhs.len() {
                let neighbors_and_faces = mesh.neighboring_patches_and_faces(CellIndex(cell));
                let mut row = matrix
                    .get_row_mut(cell)
                    .expect("Bad Initialization of matrix");
                let mut f_c = 0.;
                for (neighbor, face, face_id) in neighbors_and_faces {
                    match *neighbor {
                        Patch::Cell(id) => {
                            let normal = face
                                .normal_from_cell(CellIndex(cell))
                                .expect("Incoherence in face and cell connection");
                            let d_cf =
                                mesh.cells()[id.0].centroid() - mesh.cells()[cell].centroid();
                            let mut e_f = d_cf.normalize();
                            e_f *= face.area() * (e_f.angle(&normal)).cos();
                            let f_f = e_f.magnitude() / d_cf.magnitude();
                            f_c -= f_f;
                            match row
                                .get_entry_mut(id.0)
                                .expect("Bad Initialization of matrix")
                            {
                                SparseEntryMut::NonZero(value) => *value += f_f * coeff,
                                SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                            }
                            let t_f = face.area() * normal - e_f;
                            println!("{:?} {:?}", e_f, t_f);
                            rhs[cell] -= coeff * field.grads_face()[face_id.0].dot(&t_f);
                        }
                        Patch::Boundary(id) => match &boundary_condition[id.0] {
                            BoundaryCondition::Dirichlet(bc_value) => {
                                let d_cb = face.middle_point(mesh.vertices())
                                    - mesh.cells()[cell].centroid();
                                let normal = face
                                    .normal_from_cell(CellIndex(cell))
                                    .expect("Incoherence in face and cell connection");
                                let mut e_b = d_cb.normalize();
                                e_b *= face.area() * (e_b.angle(&normal)).cos();
                                let f_b = e_b.magnitude() / d_cb.magnitude();
                                f_c += f_b;
                                let t_b = face.area() * normal - e_b;
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] += coeff
                                    * (f_b * bc_value + field.grads_face()[face_id.0].dot(&t_b));
                            }
                            BoundaryCondition::Neumann(bc_value) => {
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] -= coeff * bc_value * face.area();
                            }
                        },
                    }
                }

                match row
                    .get_entry_mut(cell)
                    .expect("Bad Initialization of matrix")
                {
                    SparseEntryMut::NonZero(value) => *value += f_c * coeff,
                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                }

                // println!("{:?}", row)
            }
            // panic!("Test")
        }
        // Not yet minimal correction
        IntegrationCategory::Explicit => {
            for cell in 0..rhs.len() {
                let neighbors_and_faces = mesh.neighboring_patches_and_faces(CellIndex(cell));
                let mut f_c = 0.;
                for (neighbor, face, face_id) in neighbors_and_faces {
                    match *neighbor {
                        Patch::Cell(id) => {
                            let d_cf =
                                mesh.cells()[id.0].centroid() - mesh.cells()[cell].centroid();
                            let e_f = face.area() * d_cf.normalize();
                            let f_f = -e_f.magnitude() / d_cf.magnitude();
                            f_c -= f_f;
                            rhs[cell] -= f_f * coeff * field.values()[id.0];
                            let t_f = face.area()
                                * face
                                    .normal_from_cell(CellIndex(cell))
                                    .expect("Incoherence in face and cell connection")
                                - e_f;
                            rhs[cell] += coeff * field.grads_face()[face_id.0].dot(&t_f);
                        }
                        Patch::Boundary(id) => match &boundary_condition[id.0] {
                            BoundaryCondition::Dirichlet(bc_value) => {
                                let d_cb = face.middle_point(mesh.vertices())
                                    - mesh.cells()[cell].centroid();
                                let e_b = face.area() * d_cb.normalize();
                                let f_b = e_b.magnitude() / d_cb.magnitude();
                                f_c += f_b;
                                let t_b = face.area()
                                    * face
                                        .normal_from_cell(CellIndex(cell))
                                        .expect("Incoherence in face and cell connection")
                                    - e_b;
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] += coeff
                                    * (f_b * bc_value + field.grads_face()[face_id.0].dot(&t_b));
                            }
                            BoundaryCondition::Neumann(bc_value) => {
                                let bc_value = bc_value.get_value(component);
                                rhs[cell] -= coeff * bc_value * face.area();
                            }
                        },
                    }
                }
                rhs[cell] -= f_c * coeff * field.values()[cell];
            }
        }
    }
}
