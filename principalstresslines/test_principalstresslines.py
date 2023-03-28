from genericpath import isfile
import unittest
import os
import pathlib
import pyvista

from principalstresslines.convert import sfepy_from_file
from principalstresslines.analysis import analyze
from principalstresslines.model import AnalysisRequest, BoxConstraintRegion, FixedConstraint, LoadConstraint, SphereConstraintRegion, VertexConstraintRegion
from principalstresslines.materials import Concrete,Wood

class TestSuiteTestCase(unittest.TestCase):
    def test_testsuite(self):
        self.assertTrue(True)




class TestSolver(unittest.TestCase):

    def test_shell(self):
        rdr=pyvista.get_reader('/app/testdata/input/Mesh_2cm-thickness.obj')
        s=rdr.read()
        verts=[
            5.193932, -2.45694, -0.001323
            ,5.374502, -2.070572, 0.000048
            ,5.524318, -1.673636, -0.001607
            ,5.639372, -1.271812, -0.000016
            ,5.725725, -0.860975, -0.000839
            ,5.77686, -0.447481, -0.000011
            ,5.794504, -0.031619, 0.00117
            ,5.782077, 0.381038, 0.00005
            ,5.730146, 0.79063, 0.002421
            ,5.658635, 1.18964, 0.000076
            ,5.536599, 1.589711, 0.005684
            ,5.41686, 1.964905, 0.000273
            ,5.21348, 2.364896, 0.011406
            ,-0.599453, 5.693778, 0.012726
            ,-1.061087, 5.668629, 0.000299
            ,-1.450341, 5.578235, 0.006353
            ,-1.864834, 5.478006, 0.000068
            ,-2.249463, 5.334589, 0.00135
            ,-2.630085, 5.166564, 7.9891e-6
            ,-2.993922, 4.968738, -0.000614
            ,-3.340078, 4.739955, -0.000065
            ,-3.668172, 4.490176, -0.00108
            ,-3.646731, -4.498327, 0.00147
            ,-3.317407, -4.756066, 0.000088
            ,-2.963182, -4.980076, 0.001827
            ,-2.594924, -5.184251, 0.000082
            ,-2.204411, -5.351834, 0.001936
            ,-1.807452, -5.496095, 0.000131
            ,-1.38841, -5.600329, 0.002358
            ,-0.97197, -5.681598, 0.000173
            ,-0.527651, -5.717593, 0.00206
        ]
        s.scale(1000,inplace=True)
        rqst=AnalysisRequest(
             vertices=s.points.flatten().tolist()
             ,face_stride=int(len(s.faces)/s.n_faces)
             ,faces=s.faces.flatten().tolist()
             ,young_modulus=Concrete.young_modulus_mpa
             ,poisson_ratio=Concrete.poisson_ratio
             ,load_constraints=[
                LoadConstraint(
                    regions=[
                        BoxConstraintRegion(
                            type="box"
                            ,min=[-999999,-999999,-999999]
                            ,max=[999999,999999,999999]
                        )
                    ]
                    ,load_vector=[0,0,-50*46000]
                    ,is_constant=True
                )
             ]
             ,fixed_constraints=[
                FixedConstraint(
                    regions=[
                         BoxConstraintRegion(
                            type="box"
                            ,min=[-999999,-999999,-999999]
                            ,max=[999999,999999,1000]
                        )
                    ]
                )
             ]
        )
        rsp=analyze(rqst)
        self.assertIsNotNone(rsp)
        rsp.save('/app/testdata/output/shell')

    def test_sphere(self):
        s=pyvista.Sphere(radius=2.5,center=(0,0,2.5))

        rqst=AnalysisRequest(
             vertices=s.points.flatten().tolist()
             ,face_stride=int(len(s.faces)/s.n_faces)
             ,faces=s.faces.flatten().tolist()
             ,young_modulus=Wood.young_modulus_mpa
             ,poisson_ratio=Wood.poisson_ratio
             ,load_constraints=[
                LoadConstraint(
                    regions=[
                        BoxConstraintRegion(
                            type="box"
                            ,min=[-999,-999,4.8]
                            ,max=[999,999,5]
                        )
                    ]
                    ,load_vector=[0,0,-25000]
                    ,is_constant=False
                )
             ]
             ,fixed_constraints=[
                FixedConstraint(
                    regions=[
                         BoxConstraintRegion(
                            type="box"
                            ,min=[-999,-999,0]
                            ,max=[999,999,0.2]
                        )
                    ]
                )
             ]
        )
        rsp=analyze(rqst)
        self.assertIsNotNone(rsp)
        rsp.save('/app/testdata/output/sphere')
    '''
    def test_input(self):
        for file in pathlib.Path('/app/testdata/input').glob('*.*'):
            if file.is_file():
                reader=pyvista.get_reader(file)
                polydata:pyvista.PolyData=reader.read()
                verts=polydata.verts
                faces=polydata.faces
                polydata.n_faces
                compute_principal_stress_lines(
                    verts,faces.flatten(),len(faces[0])
                )
                

                mesh=sfepy_from_file(file)
                outpath=pathlib.Path('/app/testdata/output').joinpath(file.stem).with_suffix('.vtk')
                solve(mesh,outfile=str(outpath))
    '''