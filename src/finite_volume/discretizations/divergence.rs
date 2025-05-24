use std::cell::RefMut;

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use crate::finite_volume::{
    base::Field,
    boundary::FieldsBoundaryConditions,
    case::GradRequirements,
    equation::{IntegrationCategory, System, Variable},
};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum DivergenceScheme {
    RhieAndChow,
}

impl DivergenceScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::RhieAndChow => GradRequirements::new(false, false),
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
    ) {
        let var = find_var_in_fields(var, system, fields);

        match *self {
            Self::RhieAndChow => rhie_and_chow(var, system, integration, coeff),
        }
    }
}

fn rhie_and_chow(
    var: &RefMut<Field>,
    system: &mut System,
    integration: &IntegrationCategory,
    coeff: f64,
) {
    todo!()
}
