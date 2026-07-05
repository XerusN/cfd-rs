use std::cell::{Ref, RefMut};

use super::super::equation::EquationSolver;
use crate::{
    finite_volume::{
        boundary::BoundaryCondition,
        config::{GeometryConfig, GradientConfig, InitFunctionsTrait, Schemes, SchemesConfig},
        equation::{
            discretizations::{
                convection::ConvectionScheme, divergence::DivergenceScheme,
                laplacian::LaplacianScheme, time_schemes::TimeIntegration, DifferentialOperator,
            },
            operations::Op,
            variables::{ControlVolumeType, Dimension, Variable},
            Equation, FieldOperator, IntegrationCategory,
        },
        fields::Field,
        gradients::{GradientInterp, GradientMethods, GradientScheme},
        mesh::mesh,
        solvers::{
            EquationSolversSet, EquationsEnum, EquationsEnumConfigTrait, EquationsSet, Parameters,
            SolverCore, VariableFields, VariablesEnum, VariablesEnumConfigTrait, VariablesHashMaps,
        },
    },
    gradient, laplacian,
    post_processing::OutputConfig,
};

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};
use hashbrown::HashMap;
use log::{debug, info, trace, warn};

const CVT: ControlVolumeType = ControlVolumeType::Nodes;
const MAIN_SCHEMES: Schemes = Schemes {
    transient: TimeIntegration::ForwardEuler,
    convection: ConvectionScheme::UpwindSecondOrder,
    laplacian: LaplacianScheme::OrthogonalCorrection,
    divergence: DivergenceScheme::Basic,
};
const MAIN_GRADIENT_CONFIG: GradientMethods = GradientMethods {
    scheme: GradientScheme::GreenGaussCompact,
    interp: GradientInterp::AveragedCorrected,
};

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum MainVariables {
    T,
    Grad,
    Laplacian,
}

crate::generate_solver_variables! {
    MainVariables,
    T => {
        name: "T",
        dim: Dimension::Scalar,
        cvt: CVT,
        default_gradient_config: MAIN_GRADIENT_CONFIG,
    },
    Grad => {
        name: "Grad_T",
        dim: Dimension::Vector2,
        cvt: CVT,
        default_gradient_config: MAIN_GRADIENT_CONFIG,
    },
    Laplacian => {
        name: "Laplacian_T",
        dim: Dimension::Scalar,
        cvt: CVT,
        default_gradient_config: MAIN_GRADIENT_CONFIG,
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum MainEquations {
    Poisson,
    Laplacian,
    Gradient,
}

crate::generate_solver_equations! {
    MainEquations,
    Poisson => {
        name: "Poisson",
        lhs: laplacian!(MainVariables::T.var(), IntegrationCategory::Implicit),
        rhs: Op::Scalar(1.),
        default_schemes: MAIN_SCHEMES,
    },
    Gradient => {
        name: "Gradient",
        lhs: Op::FieldOperator(FieldOperator::Field(
            MainVariables::Grad.var(),
            IntegrationCategory::Implicit,
        )),
        rhs: gradient!(MainVariables::T.var()),
        default_schemes: MAIN_SCHEMES,
    },
    Laplacian => {
        name: "Laplacian",
        lhs: Op::FieldOperator(FieldOperator::Field(
            MainVariables::Laplacian.var(),
            IntegrationCategory::Implicit,
        )),
        rhs: laplacian!(MainVariables::T.var(), IntegrationCategory::Explicit),
        default_schemes: MAIN_SCHEMES,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    pub run_name: String,
    pub geometry: GeometryConfig,
    pub output: Vec<OutputConfig>,
}

pub trait PoissonUserFunctions {}

#[derive(Clone, Debug, PartialEq)]
pub struct PoissonSolver<M: MeshCore, U: PoissonUserFunctions> {
    pub config: Config,

    pub core: SolverCore<M>,

    pub user_functions: U,
}

impl<M: MeshCore, U: PoissonUserFunctions> PoissonSolver<M, U> {
    pub fn next_step(&mut self) {
        info!("Solving Poisson equation");
        Equation::solve(&mut self.core, MainEquations::Poisson.name());
        info!("Solving Laplacian equation");
        Equation::solve(&mut self.core, MainEquations::Laplacian.name());
        info!("Solving gradient equation");
        Equation::solve(&mut self.core, MainEquations::Gradient.name());

        warn!("Add User functions capability");

        self.core.increment();
    }

    pub fn new<
        I: InitFunctionsTrait,
        VC: VariablesEnumConfigTrait<MainVariables, I>,
        EC: EquationsEnumConfigTrait<MainEquations>,
    >(
        config: Config,
        mesh: Mesh<M>,
        main_eq_config: EC,
        main_var_config: VC,
        user_functions: U,
    ) -> Self {
        info!("In PoissonSolver init");
        
        let mut equations = EquationsSet::new();

        equations
            .add_eq_from_enum(&MainEquations::Poisson, &main_eq_config)
            .unwrap();
        equations
            .add_eq_from_enum(&MainEquations::Laplacian, &main_eq_config)
            .unwrap();
        equations
            .add_eq_from_enum(&MainEquations::Gradient, &main_eq_config)
            .unwrap();
        
        debug!("Equations built");
        // ---------

        let mut var_hashmaps = VariablesHashMaps::new();

        let var = MainVariables::T;
        var_hashmaps.add_var_from_enum(&var, &main_var_config);

        let var = MainVariables::Grad;
        var_hashmaps.add_var_from_enum(&var, &main_var_config);

        let var = MainVariables::Laplacian;
        var_hashmaps.add_var_from_enum(&var, &main_var_config);

        let fields = VariableFields::new(&equations, &mesh, var_hashmaps);
        
        debug!("Fields built");

        // ---------

        let solvers = EquationSolversSet::new(&mesh);

        let core_output = config.output.clone();

        let solver = Self {
            config,
            core: SolverCore {
                equation_solvers: solvers,
                variable_fields: fields,
                equations,
                mesh,
                parameters_set: Parameters::new(vec![]),
                name: "Poisson-2D".to_string(),
                step: 0,
                time: 0.,
                time_step: 0.,
                postproc: core_output,
                is_final_iter: false,
            },
            user_functions,
        };
        debug!("Core constructed");
        solver.core.try_output();
        solver
    }
}
