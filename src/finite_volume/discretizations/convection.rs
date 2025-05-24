use std::cell::RefMut;

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use crate::finite_volume::{
    base::Field, boundary::FieldsBoundaryConditions, case::{GradRequirements, VariableFields}, config::CaseConfig, equation::{Component, EquationSolver, IntegrationCategory, System, Variable}
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
        speed: &(Variable, Variable),
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
    ) {
        let speed = (
            find_var_in_fields(var, system, fields),
            find_var_in_fields(var, system, fields),
        );

        let var = find_var_in_fields(var, system, fields);

        match *self {
            Self::UpwindSecondOrder => {
                upwind_second_order(var, speed, system, integration, fields, coeff)
            }
        }
    }
}

fn upwind_second_order(
    var: &RefMut<Field>,
    speed: (&RefMut<Field>, &RefMut<Field>),
    system: &System,
    integration: &IntegrationCategory,
    fields: &Vec<RefMut<Field>>,
    coeff: f64,
) {
    todo!()
}
