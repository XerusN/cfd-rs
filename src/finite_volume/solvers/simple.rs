use hashbrown::HashMap;
use std::cell::{Ref, RefMut};

use super::super::equation::{System, Variable};
use crate::finite_volume::{
    base::Field,
    case::{Case, CaseSystems, GradRequirements, VariableFields},
    config::{CaseConfig, Schemes},
    discretizations::DifferentialOperator,
    equation::{Dimension, Equation, IntegrationCategory, Op},
    mesh::mesh,
};

use cfd_rs_utils::mesh::computational_mesh::*;

#[derive(Clone, Debug, PartialEq)]
pub struct SimpleCase {
    name: String,
    step: usize,
    time: f64,
    time_step: f64,

    config: CaseConfig,

    mesh: Computational2DMesh,

    fields: VariableFields,
    systems: CaseSystems,
}

impl Case for SimpleCase {
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
        self.systems.map.keys().collect()
    }

    fn equation(&self, name: &str) -> Option<&System> {
        self.systems.map.get(name)
    }

    fn equation_mut(&mut self, name: &str) -> Option<&mut System> {
        self.systems.map.get_mut(name)
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
        &mut CaseSystems,
        &mut VariableFields,
        &Computational2DMesh,
        &CaseConfig,
    ) {
        (
            &mut self.systems,
            &mut self.fields,
            &self.mesh,
            &self.config,
        )
    }

    fn next_step(&mut self) {
        println!("Ok");
        // Check if gradients are correctly updated
        System::solve(self, "Poisson");

        self.time += self.time_step;
        self.step += 1;
    }

    fn new(config: CaseConfig) -> Self {
        let mesh = mesh(&config.geometry);

        let u_x = Variable::new("U_x".to_string(), Dimension::Scalar);
        let u_y = Variable::new("U_y".to_string(), Dimension::Scalar);
        let u_x_temp = Variable::new("U_x_temp".to_string(), Dimension::Scalar);
        let u_y_temp = Variable::new("U_y_temp".to_string(), Dimension::Scalar);
        let p = Variable::new("P".to_string(), Dimension::Scalar);

        let mut systems = HashMap::new();

        let mut variables_glob = HashMap::new();

        let lhs = Op::Discretize(DifferentialOperator::Laplacian(
            p,
            IntegrationCategory::Implicit,
        ));
        let rhs = Op::Scalar(0.);
        let eq = Equation::new(lhs, rhs).expect("Equation not valid");
        let (system, variables) = eq.into_system(&mesh, &config.schemes);
        systems.insert("Poisson".to_string(), system);
        for (variable, grad) in &variables {
            let old_value = variables_glob.try_insert(variable.clone(), grad.clone());
            match old_value {
                Ok(_) => (),
                Err(mut old_value) => old_value.value.update_requirements(&grad),
            }
        }

        let fields = VariableFields::new(variables_glob, &mesh);

        SimpleCase {
            name: "Simple-2D".to_string(),

            time: 0.,
            time_step: 1.,
            step: 0,

            config,
            mesh,
            fields,
            systems: CaseSystems { map: systems },
        }
    }
}
