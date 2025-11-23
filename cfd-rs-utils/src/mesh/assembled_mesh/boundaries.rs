use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct Boundaries {
    pub n: usize,
    names: Vec<String>,
    faces: Vec<Vec<usize>>,
    nodes: Vec<Vec<usize>>,
    cells: Vec<Vec<usize>>,
}

impl Boundaries {
    pub fn new(
        names: Vec<String>,
        faces: Vec<Vec<usize>>,
        nodes: Vec<Vec<usize>>,
        cells: Vec<Vec<usize>>,
    ) -> Self {
        let n = names.len();
        assert_eq!(n, faces.len());
        assert_eq!(n, nodes.len());
        assert_eq!(n, cells.len());
        Self {
            n,
            names,
            faces,
            nodes,
            cells,
        }
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn faces(&self) -> &[Vec<usize>] {
        &self.faces
    }

    pub fn nodes(&self) -> &[Vec<usize>] {
        &self.nodes
    }

    pub fn cells(&self) -> &[Vec<usize>] {
        &self.cells
    }
}
