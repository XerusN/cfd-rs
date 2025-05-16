use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::{System, Variable}};

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

    pub fn discretize(&self, var: &Variable, system: &mut System, fields: &Vec<RefMut<Field>>, coeff: f64) {
        
        match *self {
            Self::ForwardEuler => forward_euler(var, system, fields, coeff),
        }
    }
}

fn forward_euler(var: &Variable, system: &mut System, fields: &Vec<RefMut<Field>>, coeff: f64) {
    
}