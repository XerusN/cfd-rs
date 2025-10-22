use std::cell::RefCell;

use cfd_rs_utils::mesh::{computational_mesh::Computational2DMesh, indices::CellIndex};
use nalgebra::Vector2;

use crate::finite_volume::{
    fields::Field,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{variables::ControlVolume, EquationSolver, IntegrationCategory, Variable},
};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum DivergenceScheme {
    Basic,
}

impl DivergenceScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            // Indirectly(face values needed)
            Self::Basic => GradRequirements::new(true, false),
        }
    }

    pub fn discretize(
        &self,
        var: &Variable,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
        equation_cv: &ControlVolume,
    ) {
        let field = find_var_in_fields(var, fields);

        match *self {
            Self::Basic => basic(field, solver, mesh, integration, coeff, equation_cv),
        }
    }
}

fn basic(
    field: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Computational2DMesh,
    integration: &IntegrationCategory,
    coeff: f64,
    equation_cv: &ControlVolume,
) {
    let field = field.borrow();
    let field = match *field {
        Field::Scalar(_) => panic!("No implemtation of divergence for a scalar field"),
        Field::Vector2(ref values) => values,
    };

    let (_, rhs) = solver.solver_borrow_mut();

    match integration {
        IntegrationCategory::Explicit => {
            for (cell_id, _) in mesh.cells().iter().enumerate() {
                let (faces_id, normals) =
                    mesh.normal_vectors_from_cell_with_faces_id(CellIndex(cell_id));
                for i in 0..faces_id.len() {
                    let face_values = Vector2::new(
                        field.x.face_values()[faces_id[i].0],
                        field.y.face_values()[faces_id[i].0],
                    );
                    let flow_rate =
                        mesh.faces()[faces_id[i].0].area() * normals[i].dot(&face_values);
                    rhs[cell_id] -= coeff * flow_rate;
                }
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}
