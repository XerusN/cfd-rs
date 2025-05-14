use super::case::Case;

/// Implementation of the creation of calculation graph for matrix creation (OpenFoam style)

#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Add(Box<(Op, Op)>),
    Sub(Box<(Op, Op)>),
    Discretize(DifferentialOperator),
    MulScalar(f64, Box<Op>),
    Equal(Box<(Op, Op)>),
}

#[derive(Clone, Debug, PartialEq)]
pub enum DifferentialOperator {
    Laplacian(Variable),
    Convection {var: Variable, speed: Variable},
    Divergence(Variable),
    TimeDerivative(Variable),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    pub op: Op,
    pub unknown: Variable,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub struct Variable {
    name: String,
}

impl Variable {
    pub fn new(name: String) -> Self {
        Variable{name}
    }
    
    pub fn name(&self) -> &str {
        &self.name
    }
}

impl Equation {
    
    pub fn new(op: Op, unknown: Variable) -> Equation {
        Equation { op, unknown }
    }
    
    pub fn construct_fn<T: Case>(eq: &Equation, case: &mut T) -> impl Fn(&mut T) -> () {
        
        
        
        |_| ()
    }
    
}