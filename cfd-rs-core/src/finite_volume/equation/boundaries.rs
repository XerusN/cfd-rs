use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};
use nalgebra_sparse::SparseEntryMut;

use crate::finite_volume::{
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    config::CaseConfig,
    equation::{
        Component, EquationSolver, variables::{ControlVolumeType, Variable}
    },
};

pub fn enforce_strong_bcs<M: MeshCore>(
    unknown: &Variable,
    component: &Component,
    solver: &mut EquationSolver,
    mesh: &Mesh<M>,
    boundary_conditions: &FieldsBoundaryConditions,
) {
    let boundary_conditions = boundary_conditions
        .map
        .get(unknown)
        .expect("Missing boundary condition for field");

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
