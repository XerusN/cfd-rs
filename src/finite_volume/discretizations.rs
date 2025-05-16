use std::cell::RefMut;

use super::{base::Field, case::GradRequirements, config::{CaseConfig, Schemes}, equation::{IntegrationCategory, System, Variable}};

pub mod convection;
pub mod divergence;
pub mod laplacian;
pub mod time_schemes;

#[derive(Clone, Debug, PartialEq)]
pub enum DifferentialOperator {
    Laplacian(Variable, IntegrationCategory),
    Convection {
        var: Variable,
        speed: (Variable, Variable),
        integration: IntegrationCategory,
    },
    Divergence(Variable, IntegrationCategory),
    TimeDerivative(Variable),
}

impl DifferentialOperator {
    pub fn required_grads(&self, schemes: &Schemes) -> (&Variable, GradRequirements) {
        match self {
            Self::Laplacian(var, _) => (var, schemes.laplacian.required_grads()),
            Self::Convection { var, .. } => (var, schemes.convection.required_grads()),
            Self::Divergence(var, _) => (var, schemes.divergence.required_grads()),
            Self::TimeDerivative(var) => (var, schemes.transient.required_grads()),
        }
    }

    pub fn variable(&self) -> &Variable {
        match self {
            Self::Laplacian(var, _) => &var,
            Self::Convection { var, .. } => &var,
            Self::Divergence(var, _) => &var,
            Self::TimeDerivative(var) => &var,
        }
    }
    
    pub fn integration(&self) -> Option<&IntegrationCategory> {
        match self {
            Self::Laplacian(_, int) => Some(&int),
            Self::Convection { integration, .. } => Some(&integration),
            Self::Divergence(_, int) => Some(&int),
            Self::TimeDerivative(_) => None,
        }
    }
    
    pub fn discretize(&self, system: &mut System, fields: &RefMut<Field>, schemes: &Schemes) {
        match self {
            Self::Laplacian(var, integration) => schemes.laplacian.discretize(&var, system, &integration, fields),
            Self::Convection {var, integration, speed } => schemes.convection.discretize(&var, &speed, system, &integration, fields),
            Self::Divergence(var, integration) => schemes.divergence.discretize(&var, system, &integration, fields),
            Self::TimeDerivative(var) => schemes.transient.discretize(&var, system, fields),
        }
    }
}