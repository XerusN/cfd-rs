use std::{cell::{RefCell, RefMut}, ops::Deref};

use cfd_rs_utils::mesh::{computational_mesh::{Computational2DMesh, Patch}, indices::CellIndex};

use super::{
    base::Field,
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    case::{GradRequirements, VariableFields},
    config::{CaseConfig, Schemes},
    equation::{Component, Dimension, Equation, EquationSolver, IntegrationCategory, Variable},
};
use std::ops::DerefMut;

pub mod convection;
pub mod divergence;
pub mod laplacian;
pub mod time_schemes;

#[derive(Clone, Debug, PartialEq)]
pub enum DifferentialOperator {
    Laplacian(Variable, IntegrationCategory),
    Convection {
        var: Variable,
        speed: Variable,
        integration: IntegrationCategory,
    },
    Divergence(Variable, IntegrationCategory),
    TimeDerivative(Variable),
    RhieChowGrad{velocity: Variable, pressure: Variable, pressure_coeff: f64},
}

impl DifferentialOperator {
    pub fn required_grads(&self, schemes: &Schemes) -> (&Variable, GradRequirements) {
        match self {
            Self::Laplacian(var, _) => (var, schemes.laplacian.required_grads()),
            Self::Convection { var, .. } => (var, schemes.convection.required_grads()),
            Self::Divergence(var, _) => (var, schemes.divergence.required_grads()),
            Self::TimeDerivative(var) => (var, schemes.transient.required_grads()),
            Self::RhieChowGrad{pressure, ..} => (pressure, GradRequirements::new(true, true)),
        }
    }
    
    // Strange for convection and Rhie-Chow
    pub fn variable(&self) -> &Variable {
        match self {
            Self::Laplacian(var, _) => &var,
            Self::Convection { var, .. } => &var,
            Self::Divergence(var, _) => &var,
            Self::TimeDerivative(var) => &var,
            Self::RhieChowGrad{pressure, ..} => &pressure,
        }
    }

    pub fn integration(&self) -> Option<IntegrationCategory> {
        match self {
            Self::Laplacian(_, int) => Some(int.clone()),
            Self::Convection { integration, .. } => Some(integration.clone()),
            Self::Divergence(_, int) => Some(int.clone()),
            Self::TimeDerivative(_) => None,
            Self::RhieChowGrad{..} => Some(IntegrationCategory::Explicit),
        }
    }

    pub fn discretize(
        &self,
        component: &Component,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        time_step: f64,
        coeff: f64,
    ) {
        match self {
            Self::Laplacian(var, integration) => {
                config.schemes.laplacian.discretize(
                    &var,
                    component,
                    solver,
                    fields,
                    mesh,
                    config,
                    &integration,
                    coeff,
                );
            }
            Self::Convection {
                var,
                integration,
                speed,
            } => config.schemes.convection.discretize(
                &var,
                component,
                &speed,
                solver,
                fields,
                mesh,
                config,
                &integration,
                coeff,
            ),
            Self::Divergence(var, integration) => {
                config.schemes.divergence.discretize(
                    &var,
                    solver,
                    fields,
                    mesh,
                    config,
                    &integration,
                    coeff,
                );
            }
            Self::TimeDerivative(var) => config
                .schemes
                .transient
                .discretize(&var, component, solver, mesh, fields, time_step, coeff),
            Self::RhieChowGrad{velocity, pressure, pressure_coeff} => rhie_chow(&velocity, &pressure, *pressure_coeff, component, solver, fields, mesh, config, coeff),
        }
    }
}

pub fn find_var_in_fields<'a>(var: &'a Variable, fields: &'a VariableFields) -> &'a RefCell<Field> {
    &fields
        .map
        .get(var)
        .expect(&format!("Missing variable {var:?} in fields",))
        .0
}

pub fn rhie_chow(
    pressure: &Variable,
    velocity: &Variable,
    pressure_coeff: f64,
    component: &Component,
    solver: &mut EquationSolver,
    fields: &VariableFields,
    mesh: &Computational2DMesh,
    config: &CaseConfig,
    coeff: f64,
) {
    // let bc_velocity = config
    //     .bc
    //     .map
    //     .get(velocity)
    //     .expect("Missing boundary condition for field");
    
    // let (matrix, rhs) = solver.solver_borrow_mut();
    
    // let field_velocity = find_var_in_fields(velocity, fields);
    // let field_velocity = field_velocity.borrow();
    // let field_velocity = match field_velocity.deref() {
    //     Field::Scalar(value) => &value,
    //     Field::Vector2(value) => match *component {
    //         Component::X => &value.x,
    //         Component::Y => &value.y,
    //     },
    // };
    // let field_pressure = find_var_in_fields(pressure, fields);
    // let field_pressure = field_pressure.borrow();
    // let field_pressure = match field_pressure.deref() {
    //     Field::Scalar(value) => &value,
    //     Field::Vector2(value) => match *component {
    //         Component::X => &value.x,
    //         Component::Y => &value.y,
    //     },
    // };
    
    // for (cell, rhs_cell) in rhs.iter_mut().enumerate() {
        
    //     let neighbors_and_faces = mesh.neighboring_patches_and_faces(CellIndex(cell));
    //     let mut sum_u = 0.;
    //     for (neighbor, face, face_id) in neighbors_and_faces {
    //         match *neighbor {
    //             Patch::Cell(id) => {
    //                 let d_cf =
    //                     mesh.cells()[id.0].centroid() - mesh.cells()[cell].centroid();
    //                 let e_f = face.area() * d_cf.normalize();
    //                 let f_f = -e_f.magnitude() / d_cf.magnitude();
    //                 f_c -= f_f;
    //                 rhs[cell] -= f_f * coeff * field.values()[id.0];
    //                 let t_f = face.area()
    //                     * face
    //                         .normal_from_cell(CellIndex(cell))
    //                         .expect("Incoherence in face and cell connection")
    //                     - e_f;
    //                 rhs[cell] += coeff * field.grads_face()[face_id.0].dot(&t_f);
    //             }
    //             Patch::Boundary(id) => match &bc[id.0] {
    //                 BoundaryCondition::Dirichlet(bc_value) => {
    //                     let d_cb = face.middle_point(mesh.vertices())
    //                         - mesh.cells()[cell].centroid();
    //                     let e_b = face.area() * d_cb.normalize();
    //                     let f_b = e_b.magnitude() / d_cb.magnitude();
    //                     f_c += f_b;
    //                     let t_b = face.area()
    //                         * face
    //                             .normal_from_cell(CellIndex(cell))
    //                             .expect("Incoherence in face and cell connection")
    //                         - e_b;
    //                     let bc_value = bc_value.get_value(component);
    //                     rhs[cell] += coeff
    //                         * (f_b * bc_value + field.grads_face()[face_id.0].dot(&t_b));
    //                 }
    //                 BoundaryCondition::Neumann(bc_value) => {
    //                     let bc_value = bc_value.get_value(component);
    //                     rhs[cell] -= coeff * bc_value * face.area();
    //                 }
    //             },
    //         }
    //     }
        
        
        
        
        
    //     *rhs_cell += coeff;
    // }
    
    
}