use std::cell::{Ref, RefMut};

use super::super::equation::{EquationSolver};
use crate::{
    finite_volume::{
        case::{Case, CaseEquations, SolversSet, VariableFields}, config::{CaseConfig, Schemes}, equation::{discretizations::DifferentialOperator, Equation, FieldOperator, IntegrationCategory, operations::Op, variables::{ControlVolumeType, Dimension, Variable}}, fields::Field, mesh::mesh
    },
    laplacian,
};

use cfd_rs_utils::{mesh::{assembled_mesh::{Mesh, MeshCore}}};

#[derive(Clone, Debug, PartialEq)]
pub struct PoissonCase<M: MeshCore> {
    name: String,
    step: usize,
    time: f64,
    time_step: f64,

    config: CaseConfig,

    mesh: Mesh<M>,

    fields: VariableFields,
    equations: CaseEquations,
    solvers: SolversSet,
}

impl<M: MeshCore> Case<M> for PoissonCase<M> {
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

    fn solver(&self, cvt: &ControlVolumeType) -> &EquationSolver {
        &self.solvers.get_from_cvt(cvt)
    }

    fn solver_mut(&mut self, cvt: &ControlVolumeType) -> &mut EquationSolver {
        self.solvers.get_from_cvt_mut(cvt)
    }

    fn schemes(&self) -> &Schemes {
        &self.config.schemes
    }

    fn mesh(&self) -> &Mesh<M> {
        &self.mesh
    }

    fn equation_solver_borrow(
        &mut self,
    ) -> (
        &mut SolversSet,
        &mut VariableFields,
        &CaseEquations,
        &Mesh<M>,
        &CaseConfig,
    ) {
        (
            &mut self.solvers,
            &mut self.fields,
            &self.equations,
            &self.mesh,
            &self.config,
        )
    }

    fn next_step(&mut self) {
        Equation::solve(self, "Poisson");

        self.time += self.time_step;
        self.step += 1;
    }

    fn new(config: CaseConfig, mesh: Mesh<M>) -> Self {

        let mut equations = CaseEquations::new();

        let t = Variable::new("T".to_string(), Dimension::Scalar, ControlVolumeType::Cells);

        let lhs = laplacian!(&t, IntegrationCategory::Implicit);

        let rhs = Op::Scalar(0.);

        let eq = Equation::new(lhs, rhs, &config.schemes).expect("Equation not valid");
        equations.add_eq("Poisson".to_string(), eq).unwrap();

        let fields = VariableFields::new(&equations, &mesh, &config);

        let solvers = SolversSet::new(&mesh);

        Self {
            name: "Poisson-2D".to_string(),

            time: 0.,
            time_step: 1.,
            step: 0,

            config,
            mesh,
            fields,
            equations: equations,
            solvers,
        }
    }
}
