use super::case::Case;

/// Implementation of the creation of calculation graph for matrix creation (OpenFoam style)

#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Add(Box<(Op, Op)>),
    Sub(Box<(Op, Op)>),
    Discretize(Variable),
    MulScalar(f64, Box<Op>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    pub op: Op,
    pub unknown: Variable,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Variable {
    Pressure,
    Speed2D,
    Temperature,
}

impl Equation {
    
    pub fn new(op: Op, unknown: Variable) -> Equation {
        Equation { op, unknown }
    }
    
    pub fn construct_fn<T: Case>(eq: &Equation, case: &mut T) -> impl Fn(&mut T) -> () {
        
        
        
        |_| ()
    }
    
    
    
    
}