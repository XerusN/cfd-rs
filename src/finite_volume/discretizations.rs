use std::cell::RefMut;

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use super::{
    base::Field,
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    case::{GradRequirements, VariableFields},
    config::{CaseConfig, Schemes},
    equation::{Component, Dimension, Equation, EquationSolver, IntegrationCategory, Variable},
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
        speed: Variable,
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
        component: &Component,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        coeff: f64,
    ) {
        match self {
            Self::Laplacian(var, integration) => {
                config.schemes
                    .laplacian
                    .discretize(&var, component, solver, fields, mesh, config, &integration, coeff);
            }
            Self::Convection {
                var,
                integration,
                speed,
            } => config.schemes.convection.discretize(
                &var,
                &speed,
                solver,
                mesh,
                config,
                &integration,
                fields,
                coeff,
            ),
            Self::Divergence(var, integration) => {
                config.schemes
                    .divergence
                    .discretize(&var, component, solver, fields, mesh, config, &integration, coeff);
            }
            Self::TimeDerivative(var) => config.schemes
                .transient
                .discretize(&var, component, solver, fields, mesh, config, coeff),
        }
    }
}

fn find_var_in_fields<'a, 'b>(
    var: &Variable,
    fields: &VariableFields,
) -> &'a RefMut<'b, Field> {
    todo!();
    let index = equation
        .fields_required()
        .iter()
        .position(|var_i| var == var_i)
        .expect(&format!(
            "Missing variable {var:?} in fields for equation {:?}",
            equation
        ));
    &fields[index]
}
