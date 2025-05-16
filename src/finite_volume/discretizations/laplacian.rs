use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{IntegrationCategory, System, Variable}};

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
    
    pub fn discretize(&self, var: &Variable, system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>, coeff: f64) {
        
        match *self {
            Self::Centered => centered(var, system, integration, fields, coeff),
        }
        
    }
}

fn centered(var: &Variable, system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>, coeff: f64) {
    
}