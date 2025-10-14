use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct Boundaries {
    names: Vec<String>,
    faces: Vec<Vec<usize>>,
    nodes: Vec<Vec<usize>>,
}

impl Boundaries {
    pub fn new(names: Vec<String>, faces: Vec<Vec<usize>>, nodes: Vec<Vec<usize>>) -> Self {
        Self {
            names,
            faces,
            nodes,
        }
    }
}
