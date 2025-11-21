use std::cell::RefCell;

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore, Patch};

use crate::finite_volume::{
    fields::Field,
    case::{GradRequirements, VariableFields},
    config::CaseConfig,
    equation::{variables::ControlVolumeType, EquationSolver, IntegrationCategory, Variable},
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

    pub fn discretize<M: MeshCore>(
        &self,
        var: &Variable,
        solver: &mut EquationSolver,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        config: &CaseConfig,
        integration: &IntegrationCategory,
        coeff: f64,
        equation_cvt: &ControlVolumeType,
    ) {
        let field = find_var_in_fields(var, fields);

        match *self {
            Self::Basic => basic(field, solver, mesh, integration, coeff, equation_cvt),
        }
    }
}

fn basic<M: MeshCore>(
    field: &RefCell<Field>,
    solver: &mut EquationSolver,
    mesh: &Mesh<M>,
    integration: &IntegrationCategory,
    coeff: f64,
    equation_cvt: &ControlVolumeType,
) {
    let field = field.borrow();
    let field = match *field {
        Field::Scalar(_) => panic!("No implemtation of divergence for a scalar field"),
        Field::Vector2(ref values) => values,
    };
    let field_cvt = field.x.cvt();

    let (_, rhs) = solver.solver_borrow_mut();
    
    let pairs = &mesh.pairs;
    let (areas, normals) = match equation_cvt {
        ControlVolumeType::Cells => (pairs.cells_areas(), pairs.cells_normals()),
        ControlVolumeType::Nodes => (pairs.nodes_areas(), pairs.nodes_normals()),
    };
    
    match integration {
        IntegrationCategory::Explicit => {
            for pair in 0..pairs.n {
                let flow_rate = coeff*areas[pair]*(normals[pair].x*field.x.faces_values()[pair] + normals[pair].y*field.y.faces_values()[pair]);
                match equation_cvt {
                    ControlVolumeType::Cells => {
                        let cells = &pairs.neighboring_cells()[pair];
                        match cells[0] {
                            Patch::Cell(i_cell) => rhs[i_cell] -= flow_rate,
                            Patch::Boundary(_) => (),
                        }
                        match cells[1] {
                            Patch::Cell(i_cell) => rhs[i_cell] += flow_rate,
                            Patch::Boundary(_) => (),
                        }
                    },
                    ControlVolumeType::Nodes => {
                        let nodes = pairs.nodes()[pair];
                        rhs[nodes[0]] -= flow_rate;
                        rhs[nodes[1]] += flow_rate;
                    },
                }
            }
        }
        IntegrationCategory::Implicit => todo!(),
    }
}
