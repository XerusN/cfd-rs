use cfd_rs_utils::mesh::{assembled_mesh::Mesh, computational_mesh::manual_meshes::quad_square};
use nalgebra::Vector2;

fn main() {
    let mesh = quad_square(&Vector2::new(1., 1.), &Vector2::new(20, 20));
    let assembled = Mesh::from(mesh);
    assembled.export("./out/test.vtu".to_owned()).unwrap();
}
