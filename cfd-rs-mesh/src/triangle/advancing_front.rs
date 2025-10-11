use cfd_rs_utils::control::OutputControl;
use cfd_rs_utils::*;
use errors::MeshError;
use mesh::indices::*;
use mesh::Modifiable2DMesh;
use std::fs;
use steps::*;

use crate::triangle::delaunay::delaunay_check_and_swap;

pub mod steps;
pub mod utils;

/// Advancing front algorithm.
/// The boundaries need to be a closed loop.
/// For now only single boundaries can be used safely (no hole and only one Cell).
///
/// # Safety
///
/// The function is at a really early stage of its implementation, so be careful with the input.
///
/// Only one cell is expected for the input.
pub fn advancing_front(
    mesh: &mut Modifiable2DMesh,
    element_size: f64,
    output_control: OutputControl,
) -> Result<(), MeshError> {
    let mut first_cell = None;

    let max_it = 200000;

    for (i, parent) in mesh.0.parents().iter().enumerate() {
        if parent == &mesh::Parent::Cell {
            first_cell = Some(ParentIndex(i))
        }
    }

    // With how the function are defined, this should remain the cell containing each edge from the front
    let first_cell = match first_cell {
        Some(value) => value,
        None => return Err(MeshError::WrongMeshInitialisation),
    };

    let mut front = vec![first_cell];

    refine_boundary(mesh, element_size);

    let mut i = 0;
    if let OutputControl::Iteration(_) = output_control {
        fs::remove_dir_all("./output").unwrap();
        fs::create_dir("output").unwrap();
        mesh.0
            .export_vtk(format!("./output/advancing_{}.vtk", i).as_str())
            .expect("");
    }
    if let OutputControl::Final = output_control {
        fs::remove_dir_all("./output").unwrap();
        fs::create_dir("output").unwrap();
    }

    loop {
        if front.is_empty() {
            break;
        }

        i += 1;
        if i > max_it {
            return Err(MeshError::MaxIterationReached { max_it });
        }
        println!(
            "{:?}) Front size : {:?}",
            i,
            front
                .iter()
                .map(|&parent| mesh.0.he_from_parent(parent).len())
                .collect::<Vec<usize>>()
        );

        let mut base_edge = select_base_edge(mesh, &front);
        let first_edge = base_edge;

        let mut result;

        loop {
            result = new_element(mesh, &mut front, base_edge, element_size);
            if let Err(MeshError::NoElementCreatable(_)) = result {
                base_edge = mesh.0.he_to_next_he()[base_edge];
                if base_edge == first_edge {
                    return Err(MeshError::NoElementCreatable(base_edge));
                }
            } else if let Err(err) = result {
                return Err(err);
            } else {
                break;
            }
        }

        if let Ok(mut stack) = result {
            println!("Delaunay");
            loop {
                if stack.is_empty() {
                    break;
                }

                delaunay_check_and_swap(mesh, &mut stack)?;

                println!("Stack size: {:?}", stack.len());
            }
        }

        println!();

        if let OutputControl::Iteration(step) = output_control {
            if i % step == 0 {
                mesh.0
                    .export_vtk(format!("./output/advancing_{}.vtk", i).as_str())
                    .expect("");
            }
        }

        let check = mesh.0.check_mesh();
        if let Err(error) = check {
            panic!("{:?}", error)
        }
    }

    if let OutputControl::Final = output_control {
        mesh.0
            .export_vtk(format!("./output/advancing_{}.vtk", i).as_str())
            .expect("");
    }

    println!("------------------------------------------------");
    println!("                Meshing complete                ");

    Ok(())
}
