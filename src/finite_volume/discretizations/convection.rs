use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{IntegrationCategory, System}};

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
    
    pub fn discretize(&self, system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>) {
        
        match *self {
            Self::UpwindSecondOrder => upwind_second_order(system, integration, fields),
        }
        
    }
}

fn upwind_second_order(system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>) {
    
}
