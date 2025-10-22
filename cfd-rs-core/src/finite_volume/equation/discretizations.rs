use std::{
    cell::{RefCell, RefMut},
    ops::Deref,
};

use cfd_rs_utils::mesh::{
    computational_mesh::{Computational2DMesh, Patch},
    indices::CellIndex,
};

use crate::finite_volume::equation::variables::ControlVolume;

use super::{
    super::{
        base::Field,
        boundary::{BoundaryCondition, FieldsBoundaryConditions},
        case::{GradRequirements, VariableFields},
        config::{CaseConfig, Schemes},
    },
    Component, Dimension, Equation, EquationSolver, IntegrationCategory, Variable,
};

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

    // Strange for convection and Rhie-Chow
    pub fn variable(&self) -> &Variable {
        match self {
            Self::Laplacian(var, _) => &var,
            Self::Convection { var, .. } => &var,
            Self::Divergence(var, _) => &var,
            Self::TimeDerivative(var) => &var,
        }
    }

    pub fn integration(&self) -> Option<IntegrationCategory> {
        match self {
            Self::Laplacian(_, int) => Some(int.clone()),
            Self::Convection { integration, .. } => Some(integration.clone()),
            Self::Divergence(_, int) => Some(int.clone()),
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
        time_step: f64,
        coeff: f64,
        equation_cv: &ControlVolume,
    ) {
        match self {
            Self::Laplacian(var, integration) => {
                config.schemes.laplacian.discretize(
                    &var,
                    component,
                    solver,
                    fields,
                    mesh,
                    config,
                    &integration,
                    coeff,
                    equation_cv,
                );
            }
            Self::Convection {
                var,
                integration,
                speed,
            } => config.schemes.convection.discretize(
                &var,
                component,
                &speed,
                solver,
                fields,
                mesh,
                config,
                &integration,
                coeff,
                equation_cv,
            ),
            Self::Divergence(var, integration) => {
                config.schemes.divergence.discretize(
                    &var,
                    solver,
                    fields,
                    mesh,
                    config,
                    &integration,
                    coeff,
                    equation_cv,
                );
            }
            Self::TimeDerivative(var) => config.schemes.transient.discretize(
                &var,
                component,
                solver,
                mesh,
                fields,
                time_step,
                coeff,
                equation_cv,
            ),
        }
    }
}

pub fn find_var_in_fields<'a>(var: &'a Variable, fields: &'a VariableFields) -> &'a RefCell<Field> {
    &fields
        .map
        .get(var)
        .expect(&format!("Missing variable {var:?} in fields",))
        .0
}
