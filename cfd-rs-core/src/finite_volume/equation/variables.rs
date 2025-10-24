use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub struct Variable {
    name: String,
    dim: Dimension,
    cvt: ControlVolumeType,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum Dimension {
    Scalar,
    Vector2,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum ControlVolumeType {
    Nodes,
    Cells,
}

impl ControlVolumeType {
    pub fn volumes<'a, M: MeshCore>(&self, mesh: &'a Mesh<M>) -> &'a [f64] {
        match self {
            ControlVolumeType::Cells => mesh.cells.volumes(),
            ControlVolumeType::Nodes => mesh.nodes.volumes(),
        }
    }
}

impl Variable {
    pub fn new(name: String, dim: Dimension, cvt: ControlVolumeType) -> Self {
        Variable { name, dim, cvt }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn dim(&self) -> &Dimension {
        &self.dim
    }

    pub fn cvt(&self) -> &ControlVolumeType {
        &self.cvt
    }
}
