use nalgebra::Vector2;
use serde::{de::DeserializeOwned, Serialize};
use std::fmt::Debug;

pub trait MeshCore: Debug + Clone + PartialEq + DeserializeOwned + Serialize {
    
}