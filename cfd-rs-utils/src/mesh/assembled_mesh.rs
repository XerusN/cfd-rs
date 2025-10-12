use serde::{Deserialize, Serialize};

pub use core::MeshCore;
pub use nodes::Nodes;
pub use cells::Cells;
pub use pairs::Pairs;

mod core;
mod nodes;
mod cells;
mod pairs;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(bound = "T: MeshCore")]
pub struct Mesh<T: MeshCore> {
    core: T,
    pub nodes: Nodes,
    pub cells: Cells,
    pub pairs: Pairs,
}

impl<T: MeshCore> Mesh<T> {
    
    
    
}

impl<T: MeshCore> From<T> for Mesh<T> {
    
    fn from(value: T) -> Self {
        
        
        
        
        
        
    }
    
}