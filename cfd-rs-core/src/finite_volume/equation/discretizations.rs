use std::cell::RefCell;

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};

use crate::finite_volume::{boundary::FieldsBoundaryConditions, equation::variables::ControlVolumeType};

use super::{
    super::{
        solvers::{GradRequirements, VariableFields},
        config::{CaseConfig, Schemes},
        fields::Field,
    },
    Component, EquationSolver, IntegrationCategory, Variable,
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

    pub fn discretize<M: MeshCore>(
        &self,
        component: &Component,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        schemes: &Schemes,
        boundary_conditions: &FieldsBoundaryConditions,
        time_step: f64,
        coeff: f64,
        equation_cvt: &ControlVolumeType,
    ) {
        match self {
            Self::Laplacian(var, integration) => {
                schemes.laplacian.discretize(
                    &var,
                    component,
                    solver,
                    fields,
                    mesh,
                    boundary_conditions,
                    &integration,
                    coeff,
                    equation_cvt,
                );
            }
            Self::Convection {
                var,
                integration,
                speed,
            } => schemes.convection.discretize(
                &var,
                component,
                &speed,
                solver,
                fields,
                mesh,
                    boundary_conditions,
                &integration,
                coeff,
                equation_cvt,
            ),
            Self::Divergence(var, integration) => {
                schemes.divergence.discretize(
                    &var,
                    solver,
                    fields,
                    mesh,
                    boundary_conditions,
                    &integration,
                    coeff,
                    equation_cvt,
                );
            }
            Self::TimeDerivative(var) => schemes.transient.discretize(
                &var,
                component,
                solver,
                mesh,
                fields,
                time_step,
                coeff,
                equation_cvt,
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
