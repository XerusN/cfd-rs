use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{IntegrationCategory, System, Variable}};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum LaplacianScheme {
    Centered,
}

impl LaplacianScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::Centered => GradRequirements::new(false, false),
        }
    }
    
    pub fn discretize(&self, var: &Variable, system: &mut System, integration: &IntegrationCategory, fields: &Vec<RefMut<Field>>, coeff: f64) {
        
        let var = find_var_in_fields(var, system, fields);
        
        match *self {
            Self::Centered => centered(var, system, integration, coeff),
        }
        
    }
}

fn centered(var: &RefMut<Field>, system: &mut System, integration: &IntegrationCategory, coeff: f64) {
    
}