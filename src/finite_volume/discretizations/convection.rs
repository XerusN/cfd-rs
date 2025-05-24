use std::cell::{RefCell, RefMut};

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

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
}

impl ConvectionScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::UpwindSecondOrder => GradRequirements::new(false, false),
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
        }
    }
}

fn upwind_second_order(
    var: &RefCell<Field>,
    component: &Component,
    speed: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    boundary_condition: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
) {
    todo!()
}
