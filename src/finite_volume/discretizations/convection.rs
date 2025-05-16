use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{IntegrationCategory, System, Variable}};

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
    
    pub fn discretize(&self, var: &Variable, speed: &(Variable, Variable), system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>, coeff: f64) {
        
        match *self {
            Self::UpwindSecondOrder => upwind_second_order(var, system, integration, fields, coeff),
        }
        
    }
}

fn upwind_second_order(var: &Variable, system: &mut System, integration: &IntegrationCategory, fields: &RefMut<Field>, coeff: f64) {
    
}
