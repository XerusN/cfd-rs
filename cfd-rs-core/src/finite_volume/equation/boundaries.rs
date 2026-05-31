use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};
use nalgebra_sparse::SparseEntryMut;

use crate::finite_volume::{
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    equation::{
        variables::{ControlVolumeType, Variable},
        Component, EquationSolver,
    },
    solvers::{find_var_in_fields, VariableFields},
};

pub fn enforce_strong_bcs<M: MeshCore>(
    unknown: &Variable,
    component: &Component,
    solver: &mut EquationSolver,
    mesh: &Mesh<M>,
    fields: &VariableFields,
) {
    let field = find_var_in_fields(unknown, &fields);
    let field = field.borrow();
    let boundary_conditions = field.boundary_condition();

    let cvt = unknown.cvt();
    let (matrix, rhs) = solver.solver_borrow_mut();
    let bnd = &mesh.boundaries;

    match cvt {
        ControlVolumeType::Nodes => {
            for (i_bnd, bc) in boundary_conditions.iter().enumerate() {
                match bc {
                    BoundaryCondition::Dirichlet(bc_value) => {
                        let bc_value = bc_value.get_value(component);
                        for cv in &bnd.nodes()[i_bnd] {
                            let mut row = matrix
                                .get_row_mut(*cv)
                                .expect("Bad Initialization of matrix");
                            for value in row.values_mut() {
                                *value = 0.;
                            }
                            match row
                                .get_entry_mut(*cv)
                                .expect("Bad Initialization of matrix")
                            {
                                SparseEntryMut::NonZero(value) => *value = 1.,
                                SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                            }

                            rhs[*cv] = bc_value;
                        }
                    }
                    BoundaryCondition::Neumann(_) => (),
                }
            }
        }
        ControlVolumeType::Cells => (),
    }
}
