use hashbrown::HashMap;
use nalgebra::Vector2;

use super::equation::{variables::Variable, Component};

#[derive(PartialEq, Debug, Clone)]
pub enum BoundaryCondition {
    Neumann(BoundaryValue),
    Dirichlet(BoundaryValue),
    // Mixed
}

#[derive(PartialEq, Debug, Clone)]
pub enum BoundaryValue {
    Scalar(f64),
    Vector2(Vector2<f64>),
}

impl BoundaryValue {
    pub fn get_value(&self, component: &Component) -> f64 {
        match *self {
            BoundaryValue::Scalar(value) => {
                if let Component::X = component {
                    value
                } else {
                    panic!("Expecting a vectorized boundary condition")
                }
            }
            BoundaryValue::Vector2(value) => match component {
                Component::X => value.x,
                Component::Y => value.y,
            },
        }
    }
}

#[derive(PartialEq, Debug, Clone)]
pub struct FieldsBoundaryConditions {
    pub map: HashMap<Variable, Vec<BoundaryCondition>>,
}

impl FieldsBoundaryConditions {
    pub fn new(map: HashMap<Variable, Vec<BoundaryCondition>>) -> Self {
        FieldsBoundaryConditions { map }
    }
}
