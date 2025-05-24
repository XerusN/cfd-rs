use std::cell::{RefCell, RefMut};

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use crate::finite_volume::{
    base::Field,
    boundary::BoundaryCondition,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{Component, EquationSolver, Variable},
};

use super::find_var_in_fields;

/// Only explicit time schemes are usable for now
#[derive(Clone, Debug, PartialEq)]
pub enum TimeIntegration {
    ForwardEuler,
}

impl TimeIntegration {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::ForwardEuler => GradRequirements::new(false, false),
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
        coeff: f64,
    ) {
        let bc = config
            .bc
            .map
            .get(var)
            .expect("Missing boundary condition for field");
        let var = find_var_in_fields(var, fields);

        match *self {
            Self::ForwardEuler => forward_euler(var, solver, mesh, bc, coeff),
        }
    }
}

fn forward_euler(
    var: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    boundary_condition: &[BoundaryCondition],
    coeff: f64,
) {
    todo!()
}
