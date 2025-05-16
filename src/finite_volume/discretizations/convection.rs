use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{IntegrationCategory, System, Variable}};

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
    
    pub fn discretize(&self, var: &Variable, speed: &(Variable, Variable), system: &mut System, integration: &IntegrationCategory, fields: &Vec<RefMut<Field>>, coeff: f64) {
        
        let speed = (find_var_in_fields(var, system, fields),  find_var_in_fields(var, system, fields));
        
        let var = find_var_in_fields(var, system, fields);
        
        
        match *self {
            Self::UpwindSecondOrder => upwind_second_order(var, speed, system, integration, fields, coeff),
        }
        
    }
}

fn upwind_second_order(var: &RefMut<Field>, speed: (&RefMut<Field>, &RefMut<Field>), system: &mut System, integration: &IntegrationCategory, fields: &Vec<RefMut<Field>>, coeff: f64) {
    
}
