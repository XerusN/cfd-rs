use std::ops::Index;

use cfd_rs_utils::mesh::indices::BoundaryPatchIndex;
use hashbrown::HashMap;

use super::equation::Variable;

#[derive(PartialEq, Debug, Clone)]
pub enum BoundaryCondition {
    Neumann(f64),
    Dirichlet(f64),
    // Mixed
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