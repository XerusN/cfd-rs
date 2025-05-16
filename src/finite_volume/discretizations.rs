use std::cell::RefMut;

use super::{
    base::Field,
    case::GradRequirements,
    config::Schemes,
    equation::{IntegrationCategory, System, Variable},
};
use std::ops::DerefMut;

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

    pub fn discretize(
        &self,
        system: &mut System,
        fields: &Vec<RefMut<Field>>,
        schemes: &Schemes,
        coeff: f64,
    ) {
        match self {
            Self::Laplacian(var, integration) => {
                schemes
                    .laplacian
                    .discretize(&var, system, &integration, fields, coeff);
            }
            Self::Convection {
                var,
                integration,
                speed,
            } => schemes
                .convection
                .discretize(&var, &speed, system, &integration, fields, coeff),
            Self::Divergence(var, integration) => {
                schemes
                    .divergence
                    .discretize(&var, system, &integration, fields, coeff)
            }
            Self::TimeDerivative(var) => schemes.transient.discretize(&var, system, fields, coeff),
        }
    }
}

fn find_var_in_fields<'a, 'b>(
    var: &Variable,
    system: &System,
    fields: &'a Vec<RefMut<'b, Field>>,
) -> &'a RefMut<'b, Field> {
    let index = system
        .fields_required()
        .iter()
        .position(|var_i| var == var_i)
        .expect(&format!(
            "Missing variable {var:?} in fields for equation {:?}",
            system.equation()
        ));
    &fields[index]
}
