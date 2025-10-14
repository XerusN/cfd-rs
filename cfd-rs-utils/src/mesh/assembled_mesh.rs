use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

pub use boundaries::Boundaries;
pub use cells::Cells;
pub use core::{MeshCore, Patch};
pub use nodes::Nodes;
pub use pairs::Pairs;
use std::{
    fs::File,
    io::{self, Write},
    path::PathBuf,
};

use crate::geometry::area;

mod boundaries;
mod cells;
mod core;
mod nodes;
mod pairs;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(bound = "T: MeshCore")]
pub struct Mesh<T: MeshCore> {
    core: T,
    pub nodes: Nodes,
    pub cells: Cells,
    pub pairs: Pairs,
    pub boundaries: Boundaries,
}

impl<T: MeshCore> Mesh<T> {}

impl<T: MeshCore> From<T> for Mesh<T> {
    fn from(core: T) -> Self {
        let n_nodes = core.n_nodes();
        let n_cells = core.n_cells();
        // Made for 2D!!
        let n_pairs = core.n_faces();
        let n_boundaries = core.n_boundaries();

        let mut nodes_centers = Vec::with_capacity(n_nodes);
        let mut nodes_volumes = Vec::with_capacity(n_nodes);
        let mut nodes_areas = Vec::with_capacity(n_nodes);
        let mut nodes_normals = Vec::with_capacity(n_nodes);
        let mut nodes_neighboring_nodes = Vec::with_capacity(n_nodes);
        let mut nodes_neighboring_cells = Vec::with_capacity(n_nodes);
        let mut nodes_neighboring_pairs = Vec::with_capacity(n_nodes);

        let mut cells_centers = Vec::with_capacity(n_cells);
        let mut cells_volumes = Vec::with_capacity(n_cells);
        let mut cells_areas = Vec::with_capacity(n_cells);
        let mut cells_normals = Vec::with_capacity(n_cells);
        let mut cells_neighboring_nodes = Vec::with_capacity(n_cells);
        let mut cells_neighboring_cells = Vec::with_capacity(n_cells);
        let mut cells_neighboring_pairs = Vec::with_capacity(n_cells);

        let mut pairs_centers = Vec::with_capacity(n_pairs);
        let mut pairs_lengths = Vec::with_capacity(n_pairs);
        let mut pairs_normals = Vec::with_capacity(n_pairs);
        let mut pairs_nodes = Vec::with_capacity(n_pairs);
        let mut pairs_neighboring_cells = Vec::with_capacity(n_pairs);

        let mut boundaries_names = Vec::with_capacity(n_boundaries);
        let mut boundaries_nodes = Vec::with_capacity(n_boundaries);
        let mut boundaries_faces = Vec::with_capacity(n_boundaries);

        // ---------------------------

        // Connectivity
        for i_bnd in 0..n_boundaries {
            boundaries_names.push(core.boundary(i_bnd));
            boundaries_faces.push(vec![]);
            boundaries_nodes.push(vec![]);
        }
        for i_node in 0..n_nodes {
            nodes_centers.push(core.node(i_node));
            nodes_neighboring_nodes.push(core.node_to_nodes(i_node));
            // need to check for consistency
            let mut i_pairs = core.node_to_faces(i_node);
            let i_cells = core.node_to_cells(i_node);
            nodes_neighboring_pairs.push(core.node_to_faces(i_node));
            nodes_neighboring_cells.push(core.node_to_cells(i_node));

            assert_eq!(
                i_pairs.len(),
                i_cells.len(),
                "Inconsistency in core connectivity"
            );

            let mut pairs = Vec::with_capacity(i_pairs.len());
            let mut cells = Vec::with_capacity(i_cells.len());

            let mut current_pair = i_pairs.swap_remove(0);
            pairs.push(current_pair);

            loop {
                let current_cell = {
                    if core.face_to_nodes(*pairs.last().unwrap())[0] == i_node {
                        core.face_to_neighbors(*pairs.last().unwrap())[1].clone()
                    } else {
                        core.face_to_neighbors(*pairs.last().unwrap())[0].clone()
                    }
                };
                cells.push(current_cell.clone());

                if i_pairs.len() == 0 {
                    break;
                }

                match current_cell {
                    Patch::Cell(i_cell) => 'pairs: {
                        let cell_faces = core.cell_to_faces(i_cell);
                        for (i, &pair) in i_pairs.iter().enumerate() {
                            if cell_faces.contains(&pair) {
                                current_pair = i_pairs.swap_remove(i);
                                break 'pairs;
                            }
                        }
                        panic!("Wrong connecvity in core mesh")
                    }
                    // Might mess up for very weird nodes connectivities (multiple unliked cells on boundaries for the same node)
                    // Should not be wrong for cfd usable meshes
                    Patch::Boundary(i_bnd) => 'pairs: {
                        boundaries_nodes[i_bnd].push(i_node);
                        for (i, &pair) in i_pairs.iter().enumerate() {
                            let face_neighbors = core.face_to_neighbors(pair);
                            for neighbor in face_neighbors {
                                if let Patch::Boundary(_) = neighbor {
                                    current_pair = i_pairs.swap_remove(i);
                                    break 'pairs;
                                }
                            }
                        }
                        panic!("Wrong connecvity in core mesh")
                    }
                }

                pairs.push(current_pair);
            }

            nodes_neighboring_pairs.push(pairs);
            nodes_neighboring_cells.push(cells);
        }
        for i_cell in 0..n_cells {
            cells_centers.push(core.node(i_cell));
            cells_neighboring_nodes.push(core.cell_to_nodes(i_cell));
            // Need to check for consitency face/neighbor
            cells_neighboring_pairs.push(core.cell_to_faces(i_cell));
            cells_neighboring_cells.push(core.cell_to_neighbors(i_cell));
        }
        for i_pair in 0..n_pairs {
            pairs_nodes.push(core.face_to_nodes(i_pair));

            let pair_to_neighbors = core.face_to_neighbors(i_pair);
            for patch in &pair_to_neighbors {
                match patch {
                    Patch::Boundary(i_bnd) => boundaries_faces[*i_bnd].push(i_pair),
                    Patch::Cell(_) => (),
                }
            }
            pairs_neighboring_cells.push(pair_to_neighbors);
        }

        // Geometry
        for i_pair in 0..n_pairs {
            let i_nodes = pairs_nodes[i_pair];
            let nodes = [core.node(i_nodes[0]), core.node(i_nodes[1])];
            pairs_centers.push(nodes[0].lerp(&nodes[1], 0.5));
            pairs_lengths.push((nodes[1] - nodes[0]).magnitude());
            pairs_normals.push((nodes[1] - nodes[0]).normalize());
        }
        for i_node in 0..n_nodes {
            let mut cv_nodes = vec![];
            for i in 0..nodes_neighboring_pairs[i_node].len() {
                cv_nodes.push(pairs_centers[nodes_neighboring_pairs[i_node][i]]);
                match nodes_neighboring_cells[i_node][i] {
                    Patch::Cell(i_cell) => cv_nodes.push(cells_centers[i_cell]),
                    Patch::Boundary(_) => cv_nodes.push(nodes_centers[i_node]),
                };
            }

            // For 2D only
            let mut areas = Vec::with_capacity(cv_nodes.len() / 2);
            let mut normals = Vec::with_capacity(cv_nodes.len() / 2);
            for i_node in 0..cv_nodes.len() / 2 {
                let vector;
                if i_node == 0 {
                    vector = cv_nodes[i_node] - *cv_nodes.last().unwrap();
                } else {
                    vector = cv_nodes[i_node * 2] - cv_nodes[i_node * 2 - 1];
                }
                let area_1 = vector.magnitude();
                let vector = vector.normalize();
                let normal_1 = Vector2::new(vector.y, -vector.x).normalize();
                let vector = cv_nodes[(i_node * 2 + 1) % cv_nodes.len()] - cv_nodes[i_node * 2];
                let area_2 = vector.magnitude();
                let vector = vector.normalize();
                let normal_2 = Vector2::new(vector.y, -vector.x).normalize();
                areas.push(area_1 + area_2);
                normals.push(normal_1.lerp(&normal_2, area_2 / (area_1 + area_2)));
            }
            nodes_areas.push(areas);
            nodes_normals.push(normals);
            nodes_volumes.push(area(&cv_nodes, &nodes_centers[i_node]));
        }
        for i_cell in 0..n_cells {
            let mut nodes = Vec::with_capacity(cells_neighboring_nodes[i_cell].len());
            for i_node in &cells_neighboring_nodes[i_cell] {
                nodes.push(nodes_centers[*i_node]);
            }
            // For 2D only
            let mut areas = Vec::with_capacity(nodes.len());
            let mut normals = Vec::with_capacity(nodes.len());
            for i_node in 0..nodes.len() {
                let vector;
                if i_node == 0 {
                    vector = nodes[i_node] - nodes[nodes.len() - 1];
                } else {
                    vector = nodes[i_node] - nodes[i_node - 1];
                }
                areas.push(vector.magnitude());
                normals.push(Vector2::new(vector.y, -vector.x).normalize());
            }

            cells_areas.push(areas);
            cells_normals.push(normals);
            cells_volumes.push(area(&nodes, &cells_centers[i_cell]));
        }

        // ---------------------------

        let nodes = Nodes::new(
            nodes_centers,
            nodes_volumes,
            nodes_areas,
            nodes_normals,
            nodes_neighboring_nodes,
            nodes_neighboring_cells,
            nodes_neighboring_pairs,
        );

        let cells = Cells::new(
            cells_centers,
            cells_volumes,
            cells_areas,
            cells_normals,
            cells_neighboring_nodes,
            cells_neighboring_cells,
            cells_neighboring_pairs,
        );

        let pairs = Pairs::new(
            pairs_centers,
            pairs_lengths,
            pairs_normals,
            pairs_nodes,
            pairs_neighboring_cells,
        );

        let boundaries = Boundaries::new(boundaries_names, boundaries_faces, boundaries_nodes);

        Mesh {
            core,
            nodes,
            cells,
            pairs,
            boundaries,
        }
    }
}

impl<T: MeshCore> Mesh<T> {
    pub fn export(&self, path: String) -> io::Result<()> {
        let path = PathBuf::from(path);

        let mut file = File::create(path)?;

        writeln!(
            file,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(file, "  <UnstructuredGrid>")?;
        writeln!(
            file,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            self.nodes.centers.len(),
            self.cells.centers.len()
        )?;
        writeln!(file, "      <Points>")?;
        writeln!(
            file,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for node in &self.nodes.centers {
            write!(file, "{} {} 0 ", node.x, node.y)?;
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(file, "      </Points>")?;

        writeln!(file, "      <Cells>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for nodes in &self.cells.neighboring_nodes {
            for i_node in nodes {
                write!(file, "{} ", i_node)?;
            }
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        let mut offset = 0;
        for faces in &self.cells.neighboring_pairs {
            offset += faces.len();
            write!(file, "{} ", offset)?;
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"types\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for faces in &self.cells.neighboring_pairs {
            if faces.len() == 3 {
                write!(file, "5 ")?;
            } else if faces.len() == 4 {
                write!(file, "7 ")?;
            }
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(file, "      </Cells>")?;

        writeln!(file, "      <CellData>")?;

        export_scalar(&mut file, &self.cells.volumes, "Cell volumes")?;

        writeln!(file, "      </CellData>")?;

        writeln!(file, "      <NodeData>")?;

        export_scalar(&mut file, &self.nodes.volumes, "Node volumes")?;

        writeln!(file, "      </NodeData>")?;

        writeln!(file, "    </Piece>")?;
        writeln!(file, "  </UnstructuredGrid>")?;
        writeln!(file, "</VTKFile>")?;

        Ok(())
    }
}

fn export_scalar(file: &mut File, data: &[f64], name: &str) -> io::Result<()> {
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
        name,
    )?;
    write!(file, "          ")?;
    for value in data {
        write!(file, "{} ", value)?;
    }
    writeln!(file)?;
    writeln!(file, "        </DataArray>")?;
    Ok(())
}
