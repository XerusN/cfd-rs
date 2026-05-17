use std::{cell::RefCell, ops::Deref};

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};
use nalgebra_sparse::SparseEntryMut;

use crate::finite_volume::{
    solvers::{GradRequirements, VariableFields},
    equation::{variables::ControlVolumeType, Component, EquationSolver, Variable},
    fields::Field,
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

    pub fn discretize<M: MeshCore>(
        &self,
        var: &Variable,
        component: &Component,
        solver: &mut EquationSolver,
        mesh: &Mesh<M>,
        fields: &VariableFields,
        time_step: f64,
        coeff: f64,
        equation_cvt: &ControlVolumeType,
    ) {
        let field = find_var_in_fields(var, fields);

        match *self {
            Self::ForwardEuler => forward_euler(
                component,
                solver,
                mesh,
                field,
                time_step,
                coeff,
                equation_cvt,
            ),
        }
    }
}

fn forward_euler<M: MeshCore>(
    component: &Component,
    solver: &mut EquationSolver,
    mesh: &Mesh<M>,
    field: &RefCell<Field>,
    time_step: f64,
    coeff: f64,
    equation_cvt: &ControlVolumeType,
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

    assert_eq!(
        field.cvt(),
        equation_cvt,
        "Implicit term, field ({:?}) and equation ({:?}) ControlVolumeTypes should be the same",
        field.cvt(),
        equation_cvt
    );

    let time_step_inv = 1. / time_step;
    let volume = equation_cvt.volumes(mesh);

    for i in 0..rhs.len() {
        let mut row = matrix.get_row_mut(i).expect("Bad Initialization of matrix");

        let f_c = time_step_inv * volume[i];

        rhs[i] += f_c * coeff * field.values()[i];

        match row.get_entry_mut(i).expect("Bad Initialization of matrix") {
            SparseEntryMut::NonZero(value) => *value += f_c * coeff,
            SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
        }
    }
}
