use hashbrown::HashMap;
use std::cell::{Ref, RefMut};

use super::super::equation::{EquationSolver, Variable};
use crate::finite_volume::{
    base::Field,
    case::{Case, CaseEquations, VariableFields},
    config::{CaseConfig, GeometryConfig, Schemes},
    discretizations::DifferentialOperator,
    equation::{Dimension, Equation, FieldOperator, IntegrationCategory, Op},
    mesh::mesh,
};

use cfd_rs_utils::{control::OutputControl, mesh::computational_mesh::*};

#[derive(Clone, Debug, PartialEq)]
pub struct PoissonCase {
    name: String,
    step: usize,
    time: f64,
    time_step: f64,

    config: CaseConfig,

    mesh: Computational2DMesh,

    fields: VariableFields,
    equations: CaseEquations,
    solver: EquationSolver,
}

impl Case for PoissonCase {
    fn name(&self) -> &str {
        &self.name
    }

    fn step(&self) -> usize {
        self.step
    }

    fn time(&self) -> f64 {
        self.time
    }

    fn time_step(&self) -> f64 {
        self.time_step
    }

    fn import_from_file(file_name: &str) -> std::io::Result<()> {
        todo!()
    }

    fn config(&self) -> &CaseConfig {
        &self.config
    }

    fn fields_list(&self) -> Vec<&Variable> {
        self.fields.map.keys().collect()
    }

    #[inline]
    fn field(&self, var: &Variable) -> Option<Ref<Field>> {
        match self.fields.map.get(var) {
            None => None,
            Some((field, _)) => Some(field.borrow()),
        }
    }

    #[inline]
    fn field_mut(&mut self, var: &Variable) -> Option<RefMut<Field>> {
        match self.fields.map.get_mut(var) {
            None => None,
            Some((field, _)) => Some(field.borrow_mut()),
        }
    }

    fn equations_list(&self) -> Vec<&String> {
        self.equations.map.keys().collect()
    }

    fn equation(&self, name: &str) -> Option<&Equation> {
        self.equations.map.get(name)
    }

    fn equation_mut(&mut self, name: &str) -> Option<&mut Equation> {
        self.equations.map.get_mut(name)
    }

    fn solver(&self) -> &EquationSolver {
        &self.solver
    }

    fn solver_mut(&mut self) -> &mut EquationSolver {
        &mut self.solver
    }

    fn schemes(&self) -> &Schemes {
        &self.config.schemes
    }

    fn mesh(&self) -> &Computational2DMesh {
        &self.mesh
    }

    fn equation_solver_borrow(
        &mut self,
    ) -> (
        &mut EquationSolver,
        &mut VariableFields,
        &CaseEquations,
        &Computational2DMesh,
        &CaseConfig,
    ) {
        (
            &mut self.solver,
            &mut self.fields,
            &self.equations,
            &self.mesh,
            &self.config,
        )
    }

    fn next_step(&mut self) {
        // Check if gradients are correctly updated
        Equation::solve(self, "Poisson");

        self.time += self.time_step;
        self.step += 1;
    }

    fn new(config: CaseConfig) -> Self {
        let mesh = mesh(&config.geometry);

        let mut equations = CaseEquations::new();

        let t = Variable::new("T".to_string(), Dimension::Scalar);

        let lhs = Op::FieldOperator(FieldOperator::DifferentialOperator(
            DifferentialOperator::Laplacian(t, IntegrationCategory::Implicit),
        ));
        let rhs = Op::Scalar(0.);

        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Poisson".to_string(), eq).unwrap();

        let fields = VariableFields::new(&equations, &mesh);

        let solver = EquationSolver::new(&mesh);

        Self {
            name: "Poisson-2D".to_string(),

            time: 0.,
            time_step: 1.,
            step: 0,

            config,
            mesh,
            fields,
            equations: equations,
            solver,
        }
    }
}
