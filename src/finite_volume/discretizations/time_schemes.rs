use std::cell::RefMut;

use crate::finite_volume::{base::Field, case::GradRequirements, equation::System};

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

    pub fn discretize(&self, system: &mut System, fields: &RefMut<Field>) {
        
        match *self {
            Self::ForwardEuler => forward_euler(system, fields),
        }
    }
}

fn forward_euler(system: &mut System, fields: &RefMut<Field>) {
    
}