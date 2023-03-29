# Notes

## Problem Statement

Limited options for stress line generation and finite element analysis natively within grasshooper, resulted in need for better stress lines.

## Technologies

SFEPY - For computing stress tensors across a given mesh in a generic way.
Docker[Runs containers, which are light weight abstractions of a full system OS (less heavy relative to a full VM)] - Encapsulated the (substantial) dependencies required for SFEPY.
TETGEN - Tetrahedralization of surface mesh
Numpy - Data array management
PyVista - Visualization, mesh wrangling
Pydantic - Parsing
F3D - Visualization of results

## Flow

1) Tetraheadralize input surface mesh into volume mesh. If surface mesh is not a closed manifold, error (couldn't get this to work).
- Tetrahedralization didn't work for open manifolds and there wasn't a simple way to do this
with the libraries we were already using. As a work-around we simply converted the open surface mesh to a closed manifold in grasshopper via extrusion.
2) Calculate displacement and stress vectors for newly created volume mesh.
3) Calculate cauchy stress and strain from stress vectors
4) Calculate principal component vectors for each node in the volume mesh. [Get the eigen vectors from the cauchy stress matrices after converting them from voight notation, done using online example using numpy]

## Inputs

* Vertices
* Faces
* Young Modulus
* Poisson Ratio
* Load Constraints
* Fixed Constraints

## Outputs

* Tetrahedralized Mesh
* Displacement Vectors
* Cauchy Stress
* Cauchy Strain


## Results

Unable to produce stress vector field superior to the one produced natively by karamba. FME is hard.

