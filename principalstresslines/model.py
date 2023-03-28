import numpy
import logging
import pyvista
import tetgen
import vtk
from scipy.spatial import KDTree
from typing import List,Tuple,Union,Iterable
from pydantic import BaseModel,validator,ValidationError,conlist,conint,root_validator,confloat,constr



logger= logging.getLogger("model")
debug,info,warn = logger.debug,logger.info,logger.warn



class ConstraintRegionBase:
    def create_region_func(self,domain:numpy.ndarray)->numpy.ndarray:
        raise NotImplementedError()


class SphereConstraintRegion(BaseModel):
    type:constr(regex="^sphere$")
    origin:conlist(float,min_items=3,max_items=3)
    radius:confloat(gt=0.0)

    def create_region_func(self,domain:numpy.ndarray):
        def get_coords(coors):
            distance=numpy.array(list(map(numpy.linalg.norm,coors-self.origin)))
            return numpy.where((distance < self.radius))[0]
        return get_coords

class BoxConstraintRegion(BaseModel):
    type:constr(regex="^box$")
    min:conlist(float,min_items=3,max_items=3)
    max:conlist(float,min_items=3,max_items=3)

    @root_validator(allow_reuse=True)
    def validate_min_max(cls,values):
        bmin=values.get('min')
        bmax = values.get('max')
        
        assert all([bmin[i] <= bmax[i] for i in range(3)]), "min <= max does not hold."
        return values

    def create_region_func(self,domain:numpy.ndarray):
        def get_coords(coors):
            x, y, z = coors[:, 0], coors[:, 1], coors[:, 2]
            
            return numpy.where((
                    (x >= self.min[0]) & (x <= self.max[0]) \
                  & (y >= self.min[1]) & (y <= self.max[1]) \
                  & (z >= self.min[2]) & (z <= self.max[2])
            ))[0]
        return get_coords

class VertexConstraintRegion(BaseModel):
    type:constr(regex="^vertex$")
    vertices:conlist(int,min_items=1)
    epsilon:float=.0001


    def create_region_func(self,domain:numpy.ndarray):
        tree_verts=domain.reshape((-1,3))[self.vertices]
        tree=KDTree(tree_verts)
        def get_coords(coors,domain=None):
            d=tree.query(coors.reshape((-1,3)),k=1,p=2)[0]
            return numpy.where((d < self.epsilon ))[0]
        return get_coords

class FixedConstraint(BaseModel):
    regions:conlist(Union[
        SphereConstraintRegion
        ,BoxConstraintRegion
        ,VertexConstraintRegion
    ],min_items=1)

    def get_regions(self)->Iterable[ConstraintRegionBase]:
        return self.regions

class LoadConstraint(FixedConstraint):
    regions:conlist(Union[
        SphereConstraintRegion
        ,BoxConstraintRegion
        ,VertexConstraintRegion
    ],min_items=1)
    load_vector:conlist(float,min_items=3,max_items=3)
    is_constant:bool

    def get_regions(self)->Iterable[ConstraintRegionBase]:
        return self.regions

class AnalysisRequest(BaseModel):
    vertices:conlist(float,min_items=9)
    face_stride:conint(ge=3,le=4)
    faces:conlist(float,min_items=3)
    young_modulus:conint(ge=1)
    poisson_ratio:confloat(ge=0.0,le=1.0)
    load_constraints:conlist(LoadConstraint,min_items=1)
    fixed_constraints:conlist(FixedConstraint,min_items=1)
    @validator('vertices',allow_reuse=True)
    def validate_vertices_modulus(cls,v):
        assert len(v) % 3 == 0, "Vertex list length is not a multiple of 3."
        return v
    @root_validator(allow_reuse=True)
    def validate_face_modulus(cls,values):
        stride=values.get('face_stride')
        faces = values.get('faces')
        assert len(faces) % stride == 0, "len(faces) is not a multiple of stride."
        return values


    def get_fixed_constraints(self)->Iterable[FixedConstraint]:
        return self.fixed_constraints

    def get_load_constraints(self)->Iterable[LoadConstraint]:
        return self.load_constraints

class AnalysisResponse(BaseModel):
    surface_mesh_vertices:List[float]
    surface_mesh_face_stride:int
    surface_mesh_faces:list[int]
    
    volume_mesh_vertices:List[float]
    volume_mesh_tetrahedrons:List[int]
    
    volume_mesh_node_displacements:List[float] # Per vertex
    volume_mesh_node_cauchy_strain:List[float] # Per Node
    volume_mesh_node_cauchy_stress:List[float] # Per Node
    volume_mesh_node_principal_stress_vectors:List[float] # Per Node



    @root_validator(allow_reuse=True)
    def validate_face_modulus(cls,values):

        sm_verts = values.get('surface_mesh_vertices')
        assert len(sm_verts) % 3 ==0,'surface mesh vertices must be multiple of 3.'

        vm_verts=values.get('volume_mesh_vertices')
        assert len(sm_verts) % 3 ==0,'volume mesh vertices must be multiple of 3.'
        vm_vert_cn = len(vm_verts)

        tetrahedrons=values.get('volume_mesh_tetrahedrons')
        assert len(tetrahedrons) % 4 == 0 ,'volume_mesh_tetrahedrons list must be multiple of 4.'
        node_cn=len(tetrahedrons)/4


        dsp=values.get('volume_mesh_node_displacements')
        assert len(dsp) == vm_vert_cn,"node displacement length must exactly equal vm vertex count"

        stress=values.get('volume_mesh_node_cauchy_stress')
        assert len(stress) == node_cn*6,"cauchy stress tensor length must be 6 times the node count."

        pstress=values.get('volume_mesh_node_principal_stress_vectors')
        assert len(pstress) == node_cn*9,"principal stress vector length must be 9 times the node count."

        return values

    def save(self,filename_wo_suffix):

        pv=pyvista.UnstructuredGrid(
            {vtk.VTK_TETRA:numpy.array(self.volume_mesh_tetrahedrons).reshape(-1,4)}
            ,numpy.array(self.volume_mesh_vertices).reshape(-1,3)
        )
        pv.cell_data["cauchy-stress"]=numpy.array(self.volume_mesh_node_cauchy_stress).reshape(-1,6)
        pv.cell_data["principal_stress_vector_1"]=numpy.array(self.volume_mesh_node_principal_stress_vectors).reshape(-1,3,3)[:,0]
        pv.cell_data["principal_stress_vector_2"]=numpy.array(self.volume_mesh_node_principal_stress_vectors).reshape(-1,3,3)[:,1]
        pv.cell_data["principal_stress_vector_3"]=numpy.array(self.volume_mesh_node_principal_stress_vectors).reshape(-1,3,3)[:,2]


        pv.point_data["displacement"]=numpy.array(self.volume_mesh_node_displacements).reshape(-1,3)

        pv.save(filename_wo_suffix+".vtk")


        pv=pyvista.UnstructuredGrid(
            {vtk.VTK_TETRA:numpy.array(self.volume_mesh_tetrahedrons).reshape(-1,4)}
            ,numpy.array(self.volume_mesh_vertices).reshape(-1,3) + numpy.array(self.volume_mesh_node_displacements).reshape(-1,3)
        )
        pv.save(filename_wo_suffix+".deformed.vtk")