
import pyvista
import numpy
import numpy.linalg
import tetgen
import logging
import vtk
from types import SimpleNamespace
from typing import List,Tuple,Union
from pydantic import BaseModel,validator,ValidationError,conlist,conint,root_validator,confloat,constr
from .stress import compute_stress

from .model import AnalysisRequest,AnalysisResponse

logger= logging.getLogger()
debug,info,warn = logger.debug,logger.info,logger.warn




def compute_principal_stresses(tensor:numpy.ndarray):
    """
    Compute principal stresses from a cauchy stress tensor.
    """
    tensor=tensor.flatten()
    if len(tensor)==6:
        # Assume Voigt notation ie: [s11 s22 s33 s23 s13 s12]
        tensor=numpy.array([
             [tensor[0],tensor[5],tensor[4]]
            ,[tensor[5],tensor[1],tensor[3]]
            ,[tensor[4],tensor[3],tensor[2]]
        ])
    tensor=tensor.reshape(3,3)
    values,vectors = numpy.linalg.eig(tensor)
    values=values.tolist()
    p1,p2,p3 = [vectors[:,values.index(v)] for v in sorted(values)]
    return p1,p2,p3
arc=pyvista.CircularArc([1,0,0],[-1,0,0],[0,0,-2])
extruded=arc.extrude([0,1,0],capping=False).extrude([0,0,.1],capping=True).clip_closed_surface(normal=[0,0,1],tolerance=.0001).fill_holes(1000)


def vector_angle(a:numpy.ndarray,b:numpy.ndarray):
    inp = numpy.inner(a,b)
    return inp/(numpy.linalg.norm(a)*numpy.linalg.norm(b))

def create_principal_stress_lines(
     sm_vertices:numpy.ndarray
    ,sm_faces:numpy.ndarray
    ,sm_face_stride:int
    ,vm_tetrahedrons:numpy.ndarray
    ,vm_vertices:numpy.ndarray
    ,vm_principal_stress_vectors:numpy.ndarray
):
    vm=pyvista.UnstructuredGrid(
        {vtk.VTK_TETRA:numpy.array(vm_tetrahedrons).reshape(-1,4)}
        ,numpy.array(vm_vertices).reshape(-1,3)
    )
    sm=pyvista.PolyData(
        var_inp=sm_vertices
        ,faces=sm_faces
        ,n_faces=int(len(sm_faces)/sm_face_stride)
    )
    sm.compute_normals(auto_orient_normals=True)
    starting_point_inx=numpy.where(sm.points[:,2] == max(sm.points[:,2]))[0][0]
    starting_point = sm.points[starting_point]
    starting_normal = sm.point_normals[starting_point_inx]

    vm_psv=vm_principal_stress_vectors.reshape((-1,3,3))


    increment=1.0

    cell_index=vm.find_closest_cell(starting_point)
    ps = vm_psv[cell_index]
    
    angles=numpy.array([
        vector_angle(ps[0],starting_normal)
        ,vector_angle(ps[1],starting_normal)
        ,vector_angle(ps[2],starting_normal)
    ])
    mask = numpy.where(angles != max(angles))[0][0]
    direction_vecs=ps[mask]

    for dvec in direction_vecs:
        unitv=dvec/numpy.linalg.norm(dvec)
        next_point=unitv*increment
        next_pt=sm.find_closest_point(next_point)










    
    


def tetrahedralize(polydata:pyvista.PolyData)->Tuple[numpy.ndarray,numpy.ndarray]:

    t=tetgen.TetGen(polydata)
    vertices,tetrahedrons = t.tetrahedralize()

    tetra_vert_counts = numpy.unique([len(tetra) for tetra in tetrahedrons])
    assert len(tetra_vert_counts)==1, 'Heterogenous nodes not supported'
    assert len(tetra_vert_counts)>0,'No elements found'
    assert tetra_vert_counts[0]==4, 'Tetrahedron must have 4 vertices'


    # Correct for vertex indices not starting from 0
    element_offset=tetrahedrons.flatten().min()
    tetrahedrons = (tetrahedrons-element_offset).flatten()
    assert len(tetrahedrons) % 4 == 0, 'Tetrahedralization produced incorrect number of indices'
    return vertices.flatten(),tetrahedrons


def analyze(request:AnalysisRequest):
    extrusion_vec=[0,0,.1]
    clip_plane_normal = extrusion_vec / numpy.linalg.norm(extrusion_vec)


    csm_vertex_list=numpy.array(request.vertices)
    csm_face_list = numpy.array(request.faces,dtype=numpy.int32)
    vm_face_stride=request.face_stride

    pvmesh=pyvista.PolyData(
        var_inp=csm_vertex_list
        ,faces=csm_face_list
        ,n_faces=int(len(csm_face_list)/vm_face_stride)
    )
    if not pvmesh.is_manifold:
        debug('Mesh is not manifold, attempting to make manifold')
        pvmesh=pvmesh.extrude(extrusion_vec,capping=True) \
            .clip_closed_surface(normal=clip_plane_normal,tolerance=.0001) \
            .fill_holes(1000) \
            .triangulate() \
            .compute_normals(auto_orient_normals=True)
        vm_face_stride=3

    vm_vertex_list,vm_tetrahedron_list = tetrahedralize(pvmesh)
    displacement,cauchy_stress,cauchy_strain=compute_stress(
        'target'
        ,vm_vertex_list
        ,vm_tetrahedron_list
        ,request
    )
    displacement=displacement
    cauchy_stress=cauchy_stress
    cauchy_strain=cauchy_strain
    principal_stresses = list(map(compute_principal_stresses,cauchy_stress.reshape((-1,6))))
    principal_stresses = numpy.array(principal_stresses).flatten().tolist()
    rc=AnalysisResponse(
        surface_mesh_vertices=pvmesh.verts.flatten().tolist()
        ,surface_mesh_faces=pvmesh.faces.flatten().tolist()
        ,surface_mesh_face_stride=vm_face_stride
        ,volume_mesh_vertices=vm_vertex_list.flatten().tolist()
        ,volume_mesh_tetrahedrons=vm_tetrahedron_list.flatten().tolist()
        ,volume_mesh_node_displacements=displacement.flatten().tolist()
        ,volume_mesh_node_cauchy_stress=cauchy_stress.flatten().tolist()
        ,volume_mesh_node_cauchy_strain=cauchy_strain.flatten().tolist()
        ,volume_mesh_node_principal_stress_vectors=principal_stresses
    )
    return rc
