use std::cell::RefMut;

use crate::finite_volume::{
    base::Field,
    case::GradRequirements,
    equation::{System, Variable},
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
        system: &mut System,
        fields: &Vec<RefMut<Field>>,
        coeff: f64,
    ) {
        let var = find_var_in_fields(var, system, fields);

        match *self {
            Self::ForwardEuler => forward_euler(var, system, coeff),
        }
    }
}

fn forward_euler(var: &RefMut<Field>, system: &mut System, coeff: f64) {}
