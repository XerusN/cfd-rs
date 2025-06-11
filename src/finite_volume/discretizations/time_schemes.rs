use std::{
    cell::{RefCell, RefMut},
    ops::Deref,
};

use cfd_rs_utils::mesh::{
    computational_mesh::{Computational2DMesh, Patch},
    indices::CellIndex,
};
use nalgebra_sparse::SparseEntryMut;

use crate::finite_volume::{
    base::Field,
    boundary::BoundaryCondition,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{Component, EquationSolver, Variable}, mesh::mesh,
};

use super::find_var_in_fields;

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

    pub fn discretize(
        &self,
        var: &Variable,
        component: &Component,
        solver: &mut EquationSolver,
        mesh: &Computational2DMesh,
        fields: &VariableFields,
        time_step: f64,
        coeff: f64,
    ) {
        let field = find_var_in_fields(var, fields);

        match *self {
            Self::ForwardEuler => forward_euler(component, solver, mesh, field, time_step, coeff),
        }
    }
}

fn forward_euler(
    component: &Component,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    field: &RefCell<Field>,
    time_step: f64,
    coeff: f64,
) {
    let (matrix, rhs) = solver.solver_borrow_mut();
    let field = field.borrow();
    let field = match field.deref() {
        Field::Scalar(value) => &value,
        Field::Vector2(value) => match *component {
            Component::X => &value.x,
            Component::Y => &value.y,
        },
    };

    let time_step_inv = 1. / time_step;

    for cell in 0..rhs.len() {
        let mut row = matrix
            .get_row_mut(cell)
            .expect("Bad Initialization of matrix");

        let f_c = time_step_inv * mesh.cells()[cell].volume();

        rhs[cell] += f_c * coeff * field.values()[cell];

        match row
            .get_entry_mut(cell)
            .expect("Bad Initialization of matrix")
        {
            SparseEntryMut::NonZero(value) => *value += f_c * coeff,
            SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
        }
    }
}
