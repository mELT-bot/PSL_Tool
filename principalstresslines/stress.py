import numpy
import logging
import pathlib
from scipy.spatial import KDTree
from typing import Callable, Tuple
from principalstresslines.model import AnalysisRequest
from principalstresslines.model import BoxConstraintRegion, ConstraintRegionBase, SphereConstraintRegion, VertexConstraintRegion

from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete import Problem,FieldVariable,Material,Integral,Equation, Equations
from sfepy.mechanics.matcoefs import stiffness_from_lame
from sfepy.terms import Term
from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.base.base import IndexedStruct
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete.fem.utils import refine_mesh
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.base.base import Struct
from sfepy.discrete.fem import Mesh
from sfepy.mesh.mesh_tools import expand2d


logger= logging.getLogger()
debug,info,warn = logger.debug,logger.info,logger.warn




def ComposedRegionFunc(funcs):
    def get_coords(coors,domain=None):
        r=numpy.array([])
        for func in funcs:
            r=numpy.union1d(r,func(coors))
        return r
    return get_coords


def _mesh_from_tetrahedron_data(name:str,vertices:numpy.ndarray,tetrahedron_vertex_indices:numpy.ndarray):
    # Always tetrahedron
    desc=["3_4"]
    assert len(tetrahedron_vertex_indices.flatten()) % 4 == 0,'tetrahedron index length not a multiple of 4!'

    tetrahedron_vertex_indices=tetrahedron_vertex_indices.reshape(-1,4)
    # Correct for vertex indices not starting from 0
    element_offset=tetrahedron_vertex_indices.flatten().min()
    tetrahedron_vertex_indices = tetrahedron_vertex_indices-element_offset
    # 1 material per 3d element
    material_ids=[[0]*len(tetrahedron_vertex_indices)]
    debug(f"vertices={len(vertices)},indices={len(tetrahedron_vertex_indices)},mats={len(material_ids)}")
    return Mesh.from_data(name,vertices.reshape(-1,3),None,[tetrahedron_vertex_indices],material_ids,desc)



def stress_strain(out, pb, state, extend=False):
    """
    Calculate and output strain and stress for given displacements.
    """
    ev = pb.evaluate
    strain = ev('ev_cauchy_strain.2.omega(u)', mode='el_avg')
    stress = ev('ev_cauchy_stress.2.omega(asphalt.D, u)', mode='el_avg',
                copy_materials=False)

    out['cauchy_strain'] = Struct(name='cauchy_strain', mode='cell',
                                  data=strain, dofs=None)
    out['cauchy_stress'] = Struct(name='cauchy_stress', mode='cell',
                                  data=stress, dofs=None)
    return out

def compute_stress(
    name:str, 
    tetrahedron_vertices:numpy.ndarray,
    tetrahedron_vertex_indices:numpy.ndarray,
    request:AnalysisRequest
)->Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
    """
    Compute displacement, cauchy stress, and cauchy strain.

    Returns:
        Displacement, cauchy stress, and cauchy strain
    """
    assert len(tetrahedron_vertices) >=0,"tetrahedron vertices not found!"
    assert len(tetrahedron_vertices.flatten()) % 3 == 0 ,"Tetrahedron vertex list length not multiple of 3!"
    assert len(tetrahedron_vertex_indices.flatten()) % 4 == 0,"Tetrahedron index list length not multiple of 4"

    mesh = _mesh_from_tetrahedron_data(name,tetrahedron_vertices,tetrahedron_vertex_indices)
    debug('Starting solve')
    debug(f"Mesh.cmesh={mesh.cmesh}")
    debug(f"Mesh.desc={mesh.descs},Mesh.name={mesh.name}")

    debug(f"Mesh.cmesh.cell_groups={mesh.cmesh.cell_groups}")

    domain = FEDomain('domain',mesh)

    omega = domain.create_region('omega','all')
    debug(f"omega={[len(e) for e in omega.entities]}")

    # Fields creation
    field = Field.from_args('displacement', numpy.float64, 'vector', omega,approx_order=2)
    # Variables
    u = FieldVariable('u', 'unknown', field)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    # Integrals
    integral = Integral('i', order=4)
    ## Thought this had to be 0 but seems to work w/ 4...
    integral0 = Integral('i', order=4)


    # Structural Material Baseline
    stiffness=stiffness_from_youngpoisson(3, request.young_modulus, request.poisson_ratio)
    asphalt = Material('asphalt', D=stiffness)
    t1 = Term.new('dw_lin_elastic(asphalt.D, v, u)',
                integral, omega, asphalt=asphalt, v=v, u=u)

    tbase=t1

    # Start region creation
    rverts=numpy.array(request.vertices)
    ## Handle fixed constraints
    ebcs=[]
    for ix,constraint in enumerate(request.get_fixed_constraints()):
        funcs=list([r.create_region_func(rverts) for r in constraint.get_regions()])
        r=domain.create_region(f"fixed{ix}",f"vertices by fixed{ix}_func",functions={
            f"fixed{ix}_func":ComposedRegionFunc(funcs)
        }) #allow_empty=false
        point_cn=len(field.get_dofs_in_region(r))
        debug(f"region 'fixed{ix}' has {point_cn} points.")
        ebcs.append(EssentialBC(f'r{ix}',r,{'u.all':0.0}))
    ## Handle load constraints
    for ix,constraint in enumerate(request.get_load_constraints()):
        funcs=list([r.create_region_func(rverts) for r in constraint.get_regions()])
        r=domain.create_region(f"rload{ix}",f"vertices by load{ix}_func",functions={
            f"load{ix}_func":ComposedRegionFunc(funcs)
        })
        loaded_point_cn=len(field.get_dofs_in_region(r))
        debug(f"region 'load{ix}' has {point_cn} points.")
        load_vec = constraint.load_vector 
        if hasattr(constraint,'is_constant') and not constraint.is_constant:
            load_vec = numpy.array(constraint.load_vector)/loaded_point_cn
        load = Material(f'load{ix}', values={'.val' : [load_vec]*loaded_point_cn})
        tbase -= Term.new(f'dw_point_load(load{ix}.val, v)', integral0, r, load0=load, v=v)
    # End region creation



    # Equation setup
    eq = Equation('balance', tbase)
    eqs = Equations([eq])

    # Solve
    ls = ScipyDirect({})
    nls_status = IndexedStruct()
    nls = Newton({}, lin_solver=ls, status=nls_status)

    pb = Problem('elasticity', equations=eqs)
    pb.set_bcs(ebcs=Conditions(ebcs))
    pb.set_solver(nls)
    status = IndexedStruct()
    variables = pb.solve(status=status,post_process_hook=stress_strain)



    debug(f'Nonlinear solver status: {nls_status}')
    debug(f'Stationary solver status: {status}')

    #if outfile:
    #    debug(f"Saving at path: {outfile}")
    #    pb.save_state(outfile, variables,post_process_hook=stress_strain)
    #    pb.save_regions_as_groups(str(pathlib.Path(outfile).with_suffix(''))+'_groups')

    extend=False
    o=variables.create_output(fill_value=None,extend=extend,linearization=pb.linearization)
    o= stress_strain(o,pb,variables,extend=extend)
    odc={k:v.to_dict() for k,v in o.items()}
    displacement=odc['u']['data']
    debug(str(odc))

    cauchy_stress = odc['cauchy_stress']['data']
    cauchy_stress=cauchy_stress.reshape((-1,6))

    cauchy_strain = odc['cauchy_strain']['data']
    cauchy_strain = cauchy_strain.reshape((-1,6))

    debug(f"displacement_vector_cn={len(displacement)}, array_lens={numpy.unique([len(d) for d in displacement])}")
    debug(f"cauchy_stress_tensor_cn={len(cauchy_stress)}, array_lens={numpy.unique([len(d) for d in cauchy_stress])}")
    debug(f"cauchy_strain_tensor_cn={len(cauchy_strain)}, array_lens={numpy.unique([len(d) for d in cauchy_strain])}")

    return displacement,cauchy_stress,cauchy_strain