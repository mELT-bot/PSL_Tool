import pyvista
import tetgen
import numpy

import logging
import logging
import pathlib
from sfepy.discrete.fem import Mesh as SFEPYMesh

logger= logging.getLogger()
debug,info,warn = logger.debug,logger.info,logger.warn


def sfepy_from_file(filename):
    reader=pyvista.get_reader(filename)
    polydata=reader.read()
    return pyvista_to_sfepy(polydata)



def data_to_pyvistamesh(points:numpy.ndarray,faces:numpy.ndarray):
    polydata=pyvista.PolyData(var_inp=points,faces=faces)
    return pyvista_to_sfepy(polydata)


def pyvista_to_sfepy(polydata:pyvista.PolyData):
    t=tetgen.TetGen(polydata)
    nodes,elements = t.tetrahedralize()
    # Always tetrahedron
    desc=["3_4"]
    element_vert_counts = numpy.unique([len(e) for e in elements])
    if len(element_vert_counts)>1:
        raise Exception('Heterogenous nodes not supported') 
    elif len(element_vert_counts)==0:
        raise Exception('No elements found')
    elif element_vert_counts[0]!=4:
        raise Exception('Tetrahedron must have 4 vertices')
    # Correct for vertex indices not starting from 0
    element_offset=elements.flatten().min()
    element_vertex_indices = [elements-element_offset]
    # 1 material per 3d element
    material_ids=[[0]*len(elements)]
    return SFEPYMesh.from_data('mesh',nodes,None,element_vertex_indices,material_ids,desc)

def sfepy_mesh_from_data(name:str,vertices:numpy.ndarray,tetrahedron_vertex_indices:numpy.ndarray):
    # Always tetrahedron
    desc=["3_4"]
    tetra_vert_counts = numpy.unique([len(tetra) for tetra in tetrahedron_vertex_indices])
    if len(tetra_vert_counts)>1:
        raise Exception('Heterogenous nodes not supported') 
    elif len(tetra_vert_counts)==0:
        raise Exception('No elements found')
    elif tetra_vert_counts[0]!=4:
        raise Exception('Tetrahedron must have 4 vertices')
    # Correct for vertex indices not starting from 0
    element_offset=tetrahedron_vertex_indices.flatten().min()
    tetrahedron_vertex_indices = [tetrahedron_vertex_indices-element_offset]
    # 1 material per 3d element
    material_ids=[[0]*len(tetrahedron_vertex_indices)]
    return SFEPYMesh.from_data(name,vertices,None,tetrahedron_vertex_indices,material_ids,desc)