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
    Centered,
}

impl DivergenceScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::Centered => GradRequirements::new(false, false),
        }
    }

    pub fn discretize(
        &self,
        var: &Variable,
        system: &mut System,
        mesh: &Computational2DMesh,
        bc: &FieldsBoundaryConditions,
        integration: &IntegrationCategory,
        fields: &Vec<RefMut<Field>>,
        coeff: f64,
    ) {
        let var = find_var_in_fields(var, system, fields);

        match *self {
            Self::Centered => centered(var, system, integration, coeff),
        }
    }
}

fn centered(
    var: &RefMut<Field>,
    system: &mut System,
    integration: &IntegrationCategory,
    coeff: f64,
) {
    todo!()
}
