use hashbrown::HashMap;
use std::cell::{Ref, RefMut};

use super::super::equation::{EquationSolver, Variable};
use crate::{
    convection, divergence,
    finite_volume::{
        base::Field,
        case::{Case, CaseEquations, VariableFields},
        config::{CaseConfig, GeometryConfig, Schemes},
        discretizations::DifferentialOperator,
        equation::{Dimension, Equation, FieldOperator, IntegrationCategory, Op},
        mesh::mesh,
    },
    gradient, laplacian, time_derivative,
};

use cfd_rs_utils::{
    control::OutputControl,
    mesh::{
        assembled_mesh::{Mesh, MeshCore},
        computational_mesh::*,
    },
};

#[derive(Clone, Debug, PartialEq)]
pub struct SimpleCase<T: MeshCore> {
    name: String,
    step: usize,
    time: f64,
    time_step: f64,

    config: CaseConfig,

    mesh: Mesh<T>,

    fields: VariableFields,
    equations: CaseEquations,
    solver: EquationSolver,

    density: f64,
    kinematic_viscosity: f64,
}

impl<T: MeshCore> Case<T> for SimpleCase<T> {
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

    fn mesh(&self) -> &Mesh<T> {
        &self.mesh
    }

    fn equation_solver_borrow(
        &mut self,
    ) -> (
        &mut EquationSolver,
        &mut VariableFields,
        &CaseEquations,
        &Mesh<T>,
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
        println!("Prediction");
        Equation::solve(self, "Prediction");

        println!("Div");
        Equation::solve(self, "Div");

        self.step += 1;
        self.export().unwrap();

        println!("Poisson");
        Equation::solve(self, "Poisson");
        println!("Grad");
        Equation::solve(self, "Grad");

        self.step += 1;
        self.export().unwrap();

        println!("Correction");
        Equation::solve(self, "Correction");

        println!("Div");
        Equation::solve(self, "Div");

        self.time += self.time_step;
        self.step += 1;
    }

    fn new(config: CaseConfig) -> Self {
        // let mesh = quad_mesh(&config.geometry);
        let mesh = mesh(&config.geometry);

        let density = 1.;
        let kinematic_viscosity = 1.;

        let time_step = 0.00001;

        let mut equations = CaseEquations::new();

        let p = Variable::new("P".to_string(), Dimension::Scalar);
        let u = Variable::new("U".to_string(), Dimension::Vector2);
        let grad_p = Variable::new("grad(P)".to_string(), Dimension::Vector2);
        let div_u = Variable::new("div(U)".to_string(), Dimension::Scalar);

        let lhs = Op::FieldOperator(FieldOperator::Field(
            div_u.clone(),
            IntegrationCategory::Implicit,
        ));
        let rhs = divergence!(&u, IntegrationCategory::Explicit);
        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Div".to_string(), eq).unwrap();

        let lhs = Op::FieldOperator(FieldOperator::Field(
            grad_p.clone(),
            IntegrationCategory::Implicit,
        ));
        let rhs = gradient!(&p);
        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Grad".to_string(), eq).unwrap();

        // ------------------------------------

        let lhs = time_derivative!(&u) + convection!(&u, &u, IntegrationCategory::Explicit);
        let rhs = kinematic_viscosity * laplacian!(&u, IntegrationCategory::Explicit);
        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Prediction".to_string(), eq).unwrap();

        let lhs = laplacian!(&p, IntegrationCategory::Implicit);
        let rhs = density / time_step * divergence!(&u, IntegrationCategory::Explicit);
        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Poisson".to_string(), eq).unwrap();

        let lhs = time_derivative!(&u);
        let rhs = -1. / density * gradient!(&p);
        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Correction".to_string(), eq).unwrap();

        // ------------------------------------

        let fields = VariableFields::new(&equations, &mesh, &config);

        let solver = EquationSolver::new(&mesh);

        Self {
            name: "Simple-2D".to_string(),

            time: 0.,
            time_step,
            step: 0,

            config,
            mesh,
            fields,
            equations: equations,
            solver,

            density,
            kinematic_viscosity,
        }
    }
}
