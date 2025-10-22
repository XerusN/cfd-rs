use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub struct Variable {
    name: String,
    dim: Dimension,
    cv: ControlVolume,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum Dimension {
    Scalar,
    Vector2,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum ControlVolume {
    Nodes,
    Cells,
}

impl ControlVolume {
    pub fn volumes<'a, M: MeshCore>(&self, mesh: &'a Mesh<M>) -> &'a [f64] {
        match self {
            ControlVolume::Cells => mesh.cells.volumes(),
            ControlVolume::Nodes => mesh.nodes.volumes(),
        }
    }
}

impl Variable {
    pub fn new(name: String, dim: Dimension, cv: ControlVolume) -> Self {
        Variable { name, dim, cv }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn dim(&self) -> &Dimension {
        &self.dim
    }

    pub fn cv(&self) -> &ControlVolume {
        &self.cv
    }
}
