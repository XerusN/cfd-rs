use hashbrown::HashMap;
use log::warn;
use nalgebra_sparse_linalg::iteratives::{Amg, IterativeSolver, gauss_seidel};
use std::{ops::Deref, vec};

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};

use crate::finite_volume::{equation::operations::Op, linalg::easy_jacobi};

use super::{
    case::{Case, GradRequirements, VariableFields},
    config::{CaseConfig, Schemes},
    error::CfdError,
    fields::Field,
};
use discretizations::{find_var_in_fields, DifferentialOperator};

use variables::{ControlVolumeType, Dimension, Variable};

pub mod discretizations;
pub mod macros;
pub mod operations;
pub mod variables;

#[derive(Clone, Debug, PartialEq)]
pub enum IntegrationCategory {
    Implicit,
    Explicit,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum Component {
    X,
    Y,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FieldOperator {
    DifferentialOperator(DifferentialOperator),
    Field(Variable, IntegrationCategory),
    Gradient(Variable),
}

impl FieldOperator {
    pub fn variable(&self) -> &Variable {
        match self {
            FieldOperator::DifferentialOperator(diff_op) => diff_op.variable(),
            FieldOperator::Gradient(var) => &var,
            FieldOperator::Field(var, _) => &var,
        }
    }

    pub fn required_grads(&self, schemes: &Schemes) -> (&Variable, GradRequirements) {
        match self {
            FieldOperator::DifferentialOperator(diff_op) => diff_op.required_grads(schemes),
            FieldOperator::Gradient(var) => (&var, GradRequirements::new(true, false)),
            FieldOperator::Field(var, _) => (&var, GradRequirements::new(false, false)),
        }
    }
}

/// Implicit formulations are allowed on one variable only
#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    lhs: Op,
    rhs: Op,
    unknown: Variable,
    variables_requirements: HashMap<Variable, GradRequirements>,
}

impl Equation {
    pub fn lhs(&self) -> &Op {
        &self.lhs
    }

    pub fn rhs(&self) -> &Op {
        &self.rhs
    }

    pub fn unknown(&self) -> &Variable {
        &self.unknown
    }
    
    pub fn cvt(&self) -> &ControlVolumeType {
        &self.unknown.cvt()
    }

    pub fn fields_required(&self) -> &HashMap<Variable, GradRequirements> {
        &self.variables_requirements
    }

    pub fn collect_field_operators(&self) -> Vec<FieldOperator> {
        let mut collector = vec![];
        self.lhs.collect_field_operators(&mut collector);
        self.rhs.collect_field_operators(&mut collector);
        collector
    }

    pub fn new(lhs: Op, rhs: Op, schemes: &Schemes) -> Result<Equation, CfdError> {
        let mut unknown_var = None;

        let mut variables_requirements = HashMap::new();

        let mut collector = vec![];
        (lhs.clone() - rhs.clone()).collect_field_operators(&mut collector);

        for f_op in &collector {
            match f_op {
                FieldOperator::DifferentialOperator(diff_op) => match diff_op {
                    DifferentialOperator::TimeDerivative(var) => match unknown_var {
                        None => unknown_var = Some(var),
                        Some(current) => {
                            if current != var {
                                return Err(CfdError::EquationInconsistentUnknown {
                                    lhs,
                                    rhs,
                                    var1: current.clone(),
                                    var2: var.clone(),
                                });
                            }
                        }
                    },

                    DifferentialOperator::Convection {
                        var, integration, ..
                    }
                    | DifferentialOperator::Divergence(var, integration)
                    | DifferentialOperator::Laplacian(var, integration) => {
                        if let IntegrationCategory::Implicit = integration {
                            match unknown_var {
                                None => unknown_var = Some(var),
                                Some(current) => {
                                    if current != var {
                                        return Err(CfdError::EquationInconsistentUnknown {
                                            lhs,
                                            rhs,
                                            var1: current.clone(),
                                            var2: var.clone(),
                                        });
                                    }
                                }
                            }
                        }
                    }
                },
                FieldOperator::Gradient(_) => (),
                FieldOperator::Field(var, integration) => {
                    if let IntegrationCategory::Implicit = integration {
                        match unknown_var {
                            None => unknown_var = Some(var),
                            Some(current) => {
                                if current != var {
                                    return Err(CfdError::EquationInconsistentUnknown {
                                        lhs,
                                        rhs,
                                        var1: current.clone(),
                                        var2: var.clone(),
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }

        for f_op in &collector {
            let var = f_op.variable();
            let (_, required_grad) = f_op.required_grads(schemes);
            variables_requirements
                .entry(var.clone())
                .and_modify(|current: &mut GradRequirements| {
                    current.update_requirements(&required_grad)
                })
                .or_insert(required_grad);
        }

        match unknown_var {
            None => Err(CfdError::EquationNoUnknown { lhs: lhs, rhs: rhs }),
            Some(var) => Ok(Equation {
                lhs,
                rhs,
                unknown: var.clone(),
                variables_requirements,
            }),
        }
    }

    pub fn solve<M: MeshCore, T: Case<M>>(case: &mut T, name: &str) -> usize {
        let time_step = case.time_step();
        let (solvers, variable_fields, equations, mesh, config) = case.equation_solver_borrow();

        let equation = equations
            .map
            .get(name)
            .expect(&format!("This equation is not defined: {name:?}"));

        solve(solvers.get_from_cvt_mut(equation.cvt()), equation, variable_fields, mesh, config, time_step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EquationSolver {
    matrix: CsrMatrix<f64>,
    rhs: DVector<f64>,
    cvt: ControlVolumeType,
}

impl EquationSolver {
    pub fn new<M: MeshCore>(mesh: &Mesh<M>, cvt: ControlVolumeType) -> Self {
        let n = match cvt {
            ControlVolumeType::Cells => mesh.cells.n,
            ControlVolumeType::Nodes => mesh.nodes.n,
        };
        let neighbors = match cvt {
            ControlVolumeType::Cells => mesh.cells.neighboring_cells(),
            ControlVolumeType::Nodes => mesh.nodes.neighboring_nodes(),
        };

        let mut matrix = CooMatrix::new(n, n);

        for i in 0..n {
            matrix.push(i, i, 0.);
            for neighbor in &neighbors[i] {
                matrix.push(i, *neighbor, 0.)
            }
        }

        let matrix = CsrMatrix::from(&matrix);

        let rhs = DVector::zeros(n);

        EquationSolver { matrix, rhs, cvt }
    }

    pub fn clear(&mut self) {
        for v in self.matrix.values_mut() {
            *v = 0.;
        }
        for v in self.rhs.iter_mut() {
            *v = 0.;
        }
    }

    pub fn matrix(&self) -> &CsrMatrix<f64> {
        &self.matrix
    }

    pub fn matrix_mut(&mut self) -> &mut CsrMatrix<f64> {
        &mut self.matrix
    }

    pub fn cvt(&self) -> &ControlVolumeType {
        &self.cvt
    }

    pub fn rhs(&self) -> &DVector<f64> {
        &self.rhs
    }

    pub fn rhs_mut(&mut self) -> &mut DVector<f64> {
        &mut self.rhs
    }

    pub fn solver_borrow_mut(&mut self) -> (&mut CsrMatrix<f64>, &mut DVector<f64>) {
        (&mut self.matrix, &mut self.rhs)
    }

    fn add_scalar<M: MeshCore>(&mut self, scalar: f64, mesh: &Mesh<M>) {
        let volumes = self.cvt.volumes(mesh);
        for (i, v) in self.rhs.iter_mut().enumerate() {
            *v -= scalar * volumes[i];
        }
    }

    fn add_gradient<M: MeshCore>(
        &mut self,
        var: &Variable,
        component: &Component,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        coeff: f64,
    ) {
        let volumes = self.cvt.volumes(mesh);

        let field = find_var_in_fields(var, fields);
        let field = field.borrow();
        let field = match *field {
            Field::Scalar(ref values) => values,
            Field::Vector2(_) => panic!("Can't add the gradient of a vector field"),
        };
        let grads = match self.cvt {
            ControlVolumeType::Cells => match var.cvt() {
                ControlVolumeType::Cells => field.grads_centers(),
                ControlVolumeType::Nodes => {
                    todo!()
                }
            },
            ControlVolumeType::Nodes => match var.cvt() {
                ControlVolumeType::Cells => {
                    todo!()
                }
                ControlVolumeType::Nodes => field.grads_centers(),
            },
        };

        match component {
            Component::X => {
                for (i, v) in self.rhs.iter_mut().enumerate() {
                    *v -= grads[i].x * coeff * volumes[i];
                }
            }
            Component::Y => {
                for (i, v) in self.rhs.iter_mut().enumerate() {
                    *v -= grads[i].y * coeff * volumes[i];
                }
            }
        }
    }

    fn add_field<M: MeshCore>(
        &mut self,
        var: &Variable,
        component: &Component,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        coeff: f64,
        integration: &IntegrationCategory,
    ) {
        let volumes = self.cvt.volumes(mesh);

        match integration {
            IntegrationCategory::Explicit => {
                let field = find_var_in_fields(var, fields);
                let field = field.borrow();
                let field = match field.deref() {
                    Field::Scalar(value) => value,
                    Field::Vector2(value) => match *component {
                        Component::X => &value.x,
                        Component::Y => &value.y,
                    },
                };
                let values = match self.cvt {
                    ControlVolumeType::Cells => match var.cvt() {
                        ControlVolumeType::Cells => field.values(),
                        ControlVolumeType::Nodes => {
                            todo!()
                        }
                    },
                    ControlVolumeType::Nodes => match var.cvt() {
                        ControlVolumeType::Cells => {
                            todo!()
                        }
                        ControlVolumeType::Nodes => field.faces_values(),
                    },
                };

                for (i, v) in self.rhs.iter_mut().enumerate() {
                    *v -= values[i] * coeff * volumes[i];
                }
            }
            IntegrationCategory::Implicit => {
                if self.cvt != *var.cvt() {
                    panic!();
                }
                for (i, j, value) in self.matrix_mut().triplet_iter_mut() {
                    if i == j {
                        *value += coeff * volumes[i];
                    }
                }
            }
        }
    }

    pub fn apply_op<M: MeshCore>(
        &mut self,
        op: &Op,
        component: &Component,
        fields: &VariableFields,
        mesh: &Mesh<M>,
        config: &CaseConfig,
        time_step: f64,
        coeff: f64,
    ) {
        match op {
            Op::Add(op) => {
                self.apply_op(
                    &op.as_ref().0,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                );
                self.apply_op(
                    &op.as_ref().1,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                );
            }
            Op::Sub(op) => {
                self.apply_op(
                    &op.as_ref().0,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                );
                self.apply_op(
                    &op.as_ref().1,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    -coeff,
                );
            }
            Op::MulScalar(scalar, op) => {
                self.apply_op(
                    op.as_ref(),
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff * scalar,
                );
            }
            Op::DivScalar(scalar, op) => {
                self.apply_op(
                    op.as_ref(),
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff / scalar,
                );
            }
            Op::MulVector(vector, op) => {
                let value = match component {
                    Component::X => vector.x,
                    Component::Y => vector.y,
                };

                self.apply_op(
                    op.as_ref(),
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff * value,
                );
            }
            Op::FieldOperator(f_op) => match f_op {
                FieldOperator::DifferentialOperator(diff_op) => diff_op.discretize(
                    component,
                    self,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                    &self.cvt().clone(),
                ),
                FieldOperator::Field(var, integration) => {
                    self.add_field(var, component, fields, mesh, coeff, integration)
                }
                FieldOperator::Gradient(var) => {
                    self.add_gradient(var, component, fields, mesh, coeff)
                }
            },
            Op::Scalar(scalar) => {
                self.add_scalar(*scalar * coeff, mesh);
            }
            Op::Vector2(vector) => match component {
                Component::X => self.add_scalar(vector.x * coeff, mesh),
                Component::Y => self.add_scalar(vector.y * coeff, mesh),
            },
        }
    }
}

fn solve<M: MeshCore>(
    solver: &mut EquationSolver,
    equation: &Equation,
    fields: &mut VariableFields,
    mesh: &Mesh<M>,
    config: &CaseConfig,
    time_step: f64,
) -> usize {
    solver.clear();

    let lhs = equation.lhs.clone();
    let rhs = equation.rhs.clone();

    let eq = lhs - rhs;

    for (var, field) in &fields.map {
        println!("Update grads: {:?}", var);
        field.0.borrow_mut().update_grads(
            &field.1,
            mesh,
            &config.schemes.gradients,
            config
                .bc
                .map
                .get(var)
                .expect("Boundary Condition missing for field"),
        );
    }

    match equation.unknown().dim() {
        Dimension::Scalar => {
            
            let mut result = false;

            let component = Component::X;
            solver.apply_op(&eq, &component, &fields, mesh, config, time_step, 1.);
            
            // for row in solver.matrix.row_iter() {
            //     println!("{:?}", row);
            // }

            let field_cell = &fields
                .map
                .get_mut(equation.unknown())
                .expect("Missing field for equation");
            
            let mut field = (field_cell.0.borrow_mut(), field_cell.1.clone());
            for i in 0..1 {
                let scalar_field = match &mut *field.0 {
                    Field::Scalar(ref mut scalar_field) => scalar_field,
                    _ => panic!("Unknown should be scalar"),
                };

                // To change
                warn!("Matrix cloned for amg");
                warn!("Hard-coded tol and max_iter for solve");
                // let result = gauss_seidel::solve_with_initial_guess(&solver.matrix, &solver.rhs, scalar_field.values_mut(), 10000, 1e-4);
                
                // let mut linalg_solver = Amg::with_smoothing(1e-4, 0.8, 100, 4, 4);
                // linalg_solver.init(
                //     solver.matrix(),
                //     solver.rhs(),
                //     Some(scalar_field.values_mut()),
                // );
                // result = linalg_solver.solve_iterations(solver.matrix(), solver.rhs(), 100);

                // *scalar_field.values_mut() = linalg_solver.x.clone();

                easy_jacobi(&solver.matrix, &solver.rhs, scalar_field.values_mut());

                if !result {
                    println!("Update Grads");
                    field.0.update_grads(
                        &field.1,
                        mesh,
                        &config.schemes.gradients,
                        config
                            .bc
                            .map
                            .get(equation.unknown())
                            .expect("Boundary Condition missing for field"),
                    );
                } else {
                    break;
                }
            }

            // let result = nalgebra_sparse_linalg::iteratives::jacobi::solve_with_initial_guess(solver.matrix(), solver.rhs(), field.values_mut(), 1000, 1e-6);

            // let result = iteratives::amg::solve_with_initial_guess(
            //     solver.matrix().clone(),
            //     &solver.rhs,
            //     field.values_mut(),
            //     1000,
            //     1e-12,
            //     0.8,
            // );

            // if !result {
            //     panic!("Did not converge when solving {:?}", equation)
            // }
        }
        Dimension::Vector2 => {
            let buffer_cell = fields
                .map
                .get_mut(equation.unknown())
                .expect("Missing field for equation")
                .0
                .clone();

            for component in [Component::X, Component::Y] {
                solver.apply_op(&eq, &component, &fields, mesh, config, time_step, 1.);

                // for row in solver.matrix.row_iter() {
                //     println!("{:?}", row);
                // }

                let field_cell = &fields
                    .map
                    .get_mut(equation.unknown())
                    .expect("Missing field for equation")
                    .0;

                let mut buffer = buffer_cell.borrow_mut();

                let buffer = match &mut *buffer {
                    Field::Vector2(ref mut scalar_field) => match component {
                        Component::X => &mut scalar_field.x,
                        Component::Y => &mut scalar_field.y,
                    },
                    _ => panic!("Unknown should be vector"),
                };

                // let mut linalg_solver = Amg::with_smoothing(1e-4, 0.8, 100, 4, 4);
                // linalg_solver.init(solver.matrix(), solver.rhs(), Some(buffer.values_mut()));
                // let result = linalg_solver.solve_iterations(solver.matrix(), solver.rhs(), 100);
                // *buffer.values_mut() = linalg_solver.x.clone();

                easy_jacobi(&solver.matrix, &solver.rhs, buffer.values_mut());

                // let result = iteratives::amg::solve_with_initial_guess(
                //     solver.matrix().clone(),
                //     &solver.rhs,
                //     field.values_mut(),
                //     1000,
                //     1e-12,
                //     0.8,
                // );

                // if !result {
                //     panic!("Did not converge when solving {:?}", equation)
                // }
            }

            fields
                .map
                .get_mut(equation.unknown())
                .expect("Missing field for equation")
                .0 = buffer_cell.clone();
        }
    }

    0
}
