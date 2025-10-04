# CFD-RS

This project aims at creating a finite-volume cfd solver in rust. The work is still in early development so there is no stable API.

## Design

The solver uses cell centered formulation and is mostly designed to handle unstructured and non-orthogonal meshes.
A strong emphasis is given on creating a tool that enables to easily define new equations to solve (openFoam like).

For now the goal is to achieve a lid-driven cavity computation using the SIMPLE algorithm (or similar).