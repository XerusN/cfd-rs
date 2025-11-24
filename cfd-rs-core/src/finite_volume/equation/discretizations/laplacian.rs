use std::{cell::RefCell, ops::Deref};

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore, Patch};
use nalgebra_sparse::SparseEntryMut;

use crate::finite_volume::{
    boundary::BoundaryCondition,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{
        variables::ControlVolumeType, Component, EquationSolver, IntegrationCategory, Variable,
    },
    fields::Field,
};

use super::find_var_in_fields;

#[derive(Clone, Debug, PartialEq)]
pub enum LaplacianScheme {
    OrthogonalCorrection,
}

impl LaplacianScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::OrthogonalCorrection => GradRequirements::new(false, true),
        }
    }

    pub fn discretize<M: MeshCore>(
        &self,
        var: &Variable,
        component: &Component,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
        equation_cvt: &ControlVolumeType,
    ) {
        let field = find_var_in_fields(var, fields);
        let bc = config
            .bc
            .map
            .get(var)
            .expect("Missing boundary condition for field");
        match *self {
            Self::OrthogonalCorrection => orthogonal_correction(
                component,
                solver,
                field,
                mesh,
                bc,
                integration,
                coeff,
                equation_cvt,
            ),
        }
    }
}

fn orthogonal_correction<M: MeshCore>(
    component: &Component,
    solver: &mut EquationSolver,
    field: &RefCell<Field>,
    mesh: &Mesh<M>,
    boundary_conditions: &Vec<BoundaryCondition>,
    integration: &IntegrationCategory,
    coeff: f64,
    equation_cvt: &ControlVolumeType,
) {
    let (matrix, rhs) = solver.solver_borrow_mut();
    let field = field.borrow();
    let field = match field.deref() {
        Field::Scalar(value) => value,
        Field::Vector2(value) => match *component {
            Component::X => &value.x,
            Component::Y => &value.y,
        },
    };

    let field_cvt = field.cvt();
    let pairs = &mesh.pairs;
    let (areas, normals) = match equation_cvt {
        ControlVolumeType::Cells => (pairs.cells_areas(), pairs.cells_normals()),
        ControlVolumeType::Nodes => (pairs.nodes_areas(), pairs.nodes_normals()),
    };
    let centers = match field_cvt {
        ControlVolumeType::Cells => mesh.cells.centers(),
        ControlVolumeType::Nodes => mesh.nodes.centers(),
    };
    
    let bnd = &mesh.boundaries;

    match *integration {
        IntegrationCategory::Implicit => {
            assert_eq!(
                field_cvt,
                equation_cvt,
                "Implicit term, field ({:?}) and equation ({:?}) ControlVolumeTypes should be the same",
                field_cvt,
                equation_cvt
            );

            for pair in 0..pairs.n {
                if pairs.on_bnd()[pair] {
                    continue;
                }

                let s_f = areas[pair] * normals[pair];

                let cvs = match field_cvt {
                    ControlVolumeType::Nodes => pairs.nodes()[pair],
                    ControlVolumeType::Cells => {
                        let cells = pairs.neighboring_cells()[pair]
                            .iter()
                            .map(|cell| match cell {
                                Patch::Boundary(_) => panic!("Undetected boundary"),
                                Patch::Cell(value) => *value,
                            })
                            .collect::<Vec<usize>>();
                        // ToOptimize
                        [cells[0], cells[1]]
                    }
                };

                let e_f = areas[pair] * (centers[cvs[1]] - centers[cvs[0]]).normalize();
                let t_f = s_f - e_f;
                let d_cf = (centers[cvs[1]] - centers[cvs[0]]).norm();

                let flux_f = coeff * e_f.norm() / d_cf;

                let mut row = matrix
                    .get_row_mut(cvs[0])
                    .expect("Bad Initialization of matrix");
                match row
                    .get_entry_mut(cvs[1])
                    .expect("Bad Initialization of matrix")
                {
                    SparseEntryMut::NonZero(value) => *value += flux_f,
                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                }
                match row
                    .get_entry_mut(cvs[0])
                    .expect("Bad Initialization of matrix")
                {
                    SparseEntryMut::NonZero(value) => *value -= flux_f,
                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                }

                let mut row = matrix
                    .get_row_mut(cvs[1])
                    .expect("Bad Initialization of matrix");
                match row
                    .get_entry_mut(cvs[0])
                    .expect("Bad Initialization of matrix")
                {
                    SparseEntryMut::NonZero(value) => *value += flux_f,
                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                }
                match row
                    .get_entry_mut(cvs[1])
                    .expect("Bad Initialization of matrix")
                {
                    SparseEntryMut::NonZero(value) => *value -= flux_f,
                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                }

                rhs[cvs[0]] -= field.grads_faces()[pair].dot(&t_f);
                rhs[cvs[1]] += field.grads_faces()[pair].dot(&t_f);
            }

            for (i_bnd, bc) in boundary_conditions.iter().enumerate() {
                match field_cvt {
                    ControlVolumeType::Nodes => {
                        match bc {
                            BoundaryCondition::Dirichlet(_) => (), //ToCheck
                            BoundaryCondition::Neumann(_) => {
                                for &pair in &bnd.faces()[i_bnd] {
                                    let s_f = areas[pair] * normals[pair];

                                    let cvs = pairs.nodes()[pair];

                                    let e_f = areas[pair] * (centers[cvs[1]] - centers[cvs[0]]).normalize();
                                    let t_f = s_f - e_f;
                                    let d_cf = (centers[cvs[1]] - centers[cvs[0]]).norm();

                                    let flux_f = coeff * e_f.norm() / d_cf;

                                    let mut row = matrix
                                        .get_row_mut(cvs[0])
                                        .expect("Bad Initialization of matrix");
                                    match row
                                        .get_entry_mut(cvs[1])
                                        .expect("Bad Initialization of matrix")
                                    {
                                        SparseEntryMut::NonZero(value) => *value += flux_f,
                                        SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                                    }
                                    match row
                                        .get_entry_mut(cvs[0])
                                        .expect("Bad Initialization of matrix")
                                    {
                                        SparseEntryMut::NonZero(value) => *value -= flux_f,
                                        SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                                    }

                                    let mut row = matrix
                                        .get_row_mut(cvs[1])
                                        .expect("Bad Initialization of matrix");
                                    match row
                                        .get_entry_mut(cvs[0])
                                        .expect("Bad Initialization of matrix")
                                    {
                                        SparseEntryMut::NonZero(value) => *value += flux_f,
                                        SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                                    }
                                    match row
                                        .get_entry_mut(cvs[1])
                                        .expect("Bad Initialization of matrix")
                                    {
                                        SparseEntryMut::NonZero(value) => *value -= flux_f,
                                        SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                                    }

                                    rhs[cvs[0]] -= field.grads_faces()[pair].dot(&t_f);
                                    rhs[cvs[1]] += field.grads_faces()[pair].dot(&t_f);
                                }
                            }
                        }
                    }
                    ControlVolumeType::Cells => match bc {
                        BoundaryCondition::Dirichlet(bc_value) => {
                            let bc_value = bc_value.get_value(component);
                            
                            for i in 0..bnd.faces().len() {
                                let i_face = bnd.faces()[i_bnd][i];
                                let i_cell = bnd.cells()[i_bnd][i];
                                let s_f = areas[i_face] * normals[i_face];
                                
                                let sign;
                                if let Patch::Cell(_) = pairs.neighboring_cells()[i_face][0] {
                                    sign = 1.;
                                } else {
                                    sign = -1.;
                                }
                                let e_f = areas[i_face] * (pairs.centers()[i_face] - centers[i_cell]).normalize()*sign;
                                let t_f = s_f - e_f;
                                let d_cf = (pairs.centers()[i_face] - centers[i_cell]).norm();

                                let flux_f = coeff * e_f.norm() / d_cf;
                                
                                let mut row = matrix
                                    .get_row_mut(i_cell)
                                    .expect("Bad Initialization of matrix");
                                match row
                                    .get_entry_mut(i_cell)
                                    .expect("Bad Initialization of matrix")
                                {
                                    SparseEntryMut::NonZero(value) => *value += flux_f,
                                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                                }
                                
                                rhs[i_cell] -= sign*(field.grads_faces()[i_face].dot(&t_f) + flux_f*bc_value);
                            }
                        }
                        BoundaryCondition::Neumann(_) => {
                            for i in 0..bnd.faces().len() {
                                let i_face = bnd.faces()[i_bnd][i];
                                let i_cell = bnd.cells()[i_bnd][i];
                                let s_f = areas[i_face] * normals[i_face];
                                
                                let sign;
                                if let Patch::Cell(_) = pairs.neighboring_cells()[i_face][0] {
                                    sign = 1.;
                                } else {
                                    sign = -1.;
                                }
                                let e_f = areas[i_face] * (pairs.centers()[i_face] - centers[i_cell]).normalize()*sign;
                                let t_f = s_f - e_f;
                                let d_cf = (pairs.centers()[i_face] - centers[i_cell]).norm();

                                let flux_f = coeff * e_f.norm() / d_cf;
                                
                                let mut row = matrix
                                    .get_row_mut(i_cell)
                                    .expect("Bad Initialization of matrix");
                                match row
                                    .get_entry_mut(i_cell)
                                    .expect("Bad Initialization of matrix")
                                {
                                    SparseEntryMut::NonZero(value) => *value += flux_f,
                                    SparseEntryMut::Zero => panic!("Bad Initialization of matrix"),
                                }
                                
                                rhs[i_cell] -= sign*(field.grads_faces()[i_face].dot(&t_f) + flux_f*field.faces_values()[i_face]);
                            }
                        }
                    },
                }
            }
        }
        IntegrationCategory::Explicit => {
            assert_eq!(
                field_cvt,
                equation_cvt,
                "Implicit term, field ({:?}) and equation ({:?}) ControlVolumeTypes should be the same",
                field_cvt,
                equation_cvt
            );

            for pair in 0..pairs.n {
                if pairs.on_bnd()[pair] {
                    continue;
                }

                let s_f = areas[pair] * normals[pair];

                let cvs = match equation_cvt {
                    ControlVolumeType::Nodes => pairs.nodes()[pair],
                    ControlVolumeType::Cells => {
                        let cells = pairs.neighboring_cells()[pair]
                            .iter()
                            .map(|cell| match cell {
                                Patch::Boundary(_) => panic!("Undetected boundary"),
                                Patch::Cell(value) => *value,
                            })
                            .collect::<Vec<usize>>();
                        // ToOptimize
                        [cells[0], cells[1]]
                    }
                };

                rhs[cvs[0]] -= field.grads_faces()[pair].dot(&s_f);
                rhs[cvs[1]] += field.grads_faces()[pair].dot(&s_f);
            }

            for (i_bnd, _) in boundary_conditions.iter().enumerate() {
                match equation_cvt {
                    ControlVolumeType::Nodes => {
                        for &pair in &bnd.faces()[i_bnd] {
                            let s_f = areas[pair] * normals[pair];

                            let cvs = pairs.nodes()[pair];

                            rhs[cvs[0]] -= field.grads_faces()[pair].dot(&s_f);
                            rhs[cvs[1]] += field.grads_faces()[pair].dot(&s_f);
                        }
                    },
                    ControlVolumeType::Cells => {
                        for &pair in &bnd.faces()[i_bnd] {
                            
                            let s_f = areas[pair] * normals[pair];

                            let cvs = &pairs.neighboring_cells()[pair];
                            
                            if let Patch::Cell(cell) = cvs[0] {
                                rhs[cell] -= field.grads_faces()[pair].dot(&s_f);
                            }
                            if let Patch::Cell(cell) = cvs[1] {
                                rhs[cell] += field.grads_faces()[pair].dot(&s_f);
                            }
                        }
                    },
                }
            }
        }
    }
}