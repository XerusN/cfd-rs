use cfd_rs_utils::mesh::{
    assembled_mesh::Mesh,
    computational_mesh::{manual_meshes::quad_square, Computational2DMesh},
};

fn main() {
    // let mesh = quad_square(&Vector2::new(1., 1.), &Vector2::new(3, 3));
    // for face in mesh.faces() {
    //     println!("{:?}", face.vertices())
    // }
    // println!("{:?}", mesh.vertices());
    // println!("------------------------------");
    // let assembled = Mesh::from(mesh);
    // assembled.export_cell_centered("./out/test.vtu".to_owned()).unwrap();
    // assembled.export_node_centered("./out/test_node.vtu".to_owned()).unwrap();
    // println!("{:?}", assembled.nodes);
    // println!("------------------------------");
    // println!("{:?}", assembled.pairs);
    // println!("------------------------------");
    // println!("{:?}", assembled.cells);
    // println!("------------------------------");
    // println!("{:?}", assembled.boundaries);

    let mesh = Computational2DMesh::deserialize_file(
        "/home/xaviern/Documents/VSCode/CFD/cfd-rs-workspace/cfd-rs/target/exports/mesh.cfd",
    )
    .unwrap();
    let assembled = Mesh::from(mesh);
    assembled
        .export_cell_centered("./out/test.vtu".to_owned())
        .unwrap();
    assembled
        .export_node_centered("./out/test_node.vtu".to_owned())
        .unwrap();
}
