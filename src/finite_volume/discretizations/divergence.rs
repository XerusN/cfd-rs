use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{IntegrationCategory, System}};

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
    
    pub fn discretize(&self, system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>) {
        
        match *self {
            Self::Centered => centered(system, integration, fields),
        }
        
    }
}

fn centered(system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>) {
    
}
