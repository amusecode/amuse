from amuse.units import units
from amuse.datamodel import new_cartesian_grid
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase

try:
    from omuse.ext import grid_remappers
except:
    from amuse.ext import grid_remappers

import numpy

from amuse.datamodel.grids import *
from amuse.datamodel.staggeredgrid import StaggeredGrid

class TestGridRemappers(TestCase):

    def setUp(self):
        if not grid_remappers.matplotlib_available:
            self.skip("matplotlib not available")

    def test1(self):
        source=new_cartesian_grid((10,20),1.)
        target=new_cartesian_grid((5,10),2.,offset=[0.,0.15])
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.interpolating_2D_remapper(source,target)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)

    def test2(self):
        source=new_cartesian_grid((10,20),1. | units.m)
        target=new_cartesian_grid((5,10),2. | units.m,offset=[0.,0.15] | units.m)
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.interpolating_2D_remapper(source,target)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)

class TestGridRemappers_bilinear(TestCase):

    def test1(self):
        source=new_cartesian_grid((10,20),1.)
        target=new_cartesian_grid((5,10),2.,offset=[0.,0.25])
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.bilinear_2D_remapper(source,target)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)

    def test2(self):
        source=new_cartesian_grid((10,20),1.)
        target=new_cartesian_grid((8,8),0.125,offset=[0.5,0.5])
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.bilinear_2D_remapper(source,target)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)


    def test3(self):
        source=new_cartesian_grid((10,20),1. | units.m)
        target=new_cartesian_grid((5,10),2. | units.m,offset=[0.,0.25] | units.m)
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.bilinear_2D_remapper(source,target)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)

    def test4(self):
        source=new_cartesian_grid((10,20),1.)
        target=new_cartesian_grid((10,20),1.5,offset=[-0.5,-0.5])
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.bilinear_2D_remapper(source,target, check_inside=False)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.xcopy,numpy.clip(target.x,0.5,9.5))
        self.assertEqual(target.ycopy,numpy.clip(target.y,0.5,19.5))

class TestGridRemappers_nearest(TestCase):

    def test1(self):
        source=new_cartesian_grid((10,20),1.)
        target=new_cartesian_grid((5,10),1.,offset=[0.,0.])
        
        source.xcopy=source.x
        source.ycopy=source.y
        
        remapper=grid_remappers.nearest_2D_remapper(source,target)
        remapper.forward_mapping(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)



class TestConservativeSphericalRemapper(TestCase):

    def setUp(self):
        try:
            from omuse.community.cdo.interface import CDORemapper
        except:
            self.skip("conservative spherical remapper requires omuse.community.cdo.interface")

        #this test creates a structured staggered grid and an unstructured staggered grid
        #and then uses the conservative_spherical_remapper to remap values between the grids

        #define nodal points and triangles of a small test grid
        #got this grid from http://matplotlib.org/examples/pylab_examples/triplot_demo.html
        xy = numpy.asarray([
            [-0.101, 0.872], [-0.080, 0.883], [-0.069, 0.888], [-0.054, 0.890],
            [-0.045, 0.897], [-0.057, 0.895], [-0.073, 0.900], [-0.087, 0.898],
            [-0.090, 0.904], [-0.069, 0.907], [-0.069, 0.921], [-0.080, 0.919],
            [-0.073, 0.928], [-0.052, 0.930], [-0.048, 0.942], [-0.062, 0.949],
            [-0.054, 0.958], [-0.069, 0.954], [-0.087, 0.952], [-0.087, 0.959],
            [-0.080, 0.966], [-0.085, 0.973], [-0.087, 0.965], [-0.097, 0.965],
            [-0.097, 0.975], [-0.092, 0.984], [-0.101, 0.980], [-0.108, 0.980],
            [-0.104, 0.987], [-0.102, 0.993], [-0.115, 1.001], [-0.099, 0.996],
            [-0.101, 1.007], [-0.090, 1.010], [-0.087, 1.021], [-0.069, 1.021],
            [-0.052, 1.022], [-0.052, 1.017], [-0.069, 1.010], [-0.064, 1.005],
            [-0.048, 1.005], [-0.031, 1.005], [-0.031, 0.996], [-0.040, 0.987],
            [-0.045, 0.980], [-0.052, 0.975], [-0.040, 0.973], [-0.026, 0.968],
            [-0.020, 0.954], [-0.006, 0.947], [ 0.003, 0.935], [ 0.006, 0.926],
            [ 0.005, 0.921], [ 0.022, 0.923], [ 0.033, 0.912], [ 0.029, 0.905],
            [ 0.017, 0.900], [ 0.012, 0.895], [ 0.027, 0.893], [ 0.019, 0.886],
            [ 0.001, 0.883], [-0.012, 0.884], [-0.029, 0.883], [-0.038, 0.879],
            [-0.057, 0.881], [-0.062, 0.876], [-0.078, 0.876], [-0.087, 0.872],
            [-0.030, 0.907], [-0.007, 0.905], [-0.057, 0.916], [-0.025, 0.933],
            [-0.077, 0.990], [-0.059, 0.993]])
        triangles = numpy.asarray([
            [67, 66,  1], [65,  2, 66], [ 1, 66,  2], [64,  2, 65], [63,  3, 64],
            [60, 59, 57], [ 2, 64,  3], [ 3, 63,  4], [ 0, 67,  1], [62,  4, 63],
            [57, 59, 56], [59, 58, 56], [61, 60, 69], [57, 69, 60], [ 4, 62, 68],
            [ 6,  5,  9], [61, 68, 62], [69, 68, 61], [ 9,  5, 70], [ 6,  8,  7],
            [ 4, 70,  5], [ 8,  6,  9], [56, 69, 57], [69, 56, 52], [70, 10,  9],
            [54, 53, 55], [56, 55, 53], [68, 70,  4], [52, 56, 53], [11, 10, 12],
            [69, 71, 68], [68, 13, 70], [10, 70, 13], [51, 50, 52], [13, 68, 71],
            [52, 71, 69], [12, 10, 13], [71, 52, 50], [71, 14, 13], [50, 49, 71],
            [49, 48, 71], [14, 16, 15], [14, 71, 48], [17, 19, 18], [17, 20, 19],
            [48, 16, 14], [48, 47, 16], [47, 46, 16], [16, 46, 45], [23, 22, 24],
            [21, 24, 22], [17, 16, 45], [20, 17, 45], [21, 25, 24], [27, 26, 28],
            [20, 72, 21], [25, 21, 72], [45, 72, 20], [25, 28, 26], [44, 73, 45],
            [72, 45, 73], [28, 25, 29], [29, 25, 31], [43, 73, 44], [73, 43, 40],
            [72, 73, 39], [72, 31, 25], [42, 40, 43], [31, 30, 29], [39, 73, 40],
            [42, 41, 40], [72, 33, 31], [32, 31, 33], [39, 38, 72], [33, 72, 38],
            [33, 38, 34], [37, 35, 38], [34, 38, 35], [35, 37, 36]])

        num_elems = len(triangles)
        elements = UnstructuredGrid(num_elems)
        elements.n1 = triangles[:,0] - 1
        elements.n2 = triangles[:,1] - 1
        elements.n3 = triangles[:,2] - 1

        lons = numpy.zeros(num_elems, dtype=numpy.double)
        lats = numpy.zeros(num_elems, dtype=numpy.double)
        for i in range(num_elems):
            for n in triangles[i]:
                lons[i] += xy[n,0]/3.0
                lats[i] += xy[n,1]/3.0
        elements.lon = (lons | units.rad)
        elements.lat = (lats | units.rad)

        nodes = UnstructuredGrid(len(xy))
        nodes.lon = (xy[:,0] | units.rad)
        nodes.lat = (xy[:,1] | units.rad)

        self.unstructured = StaggeredGrid(elements, nodes)

        #generate corners for a simple structured grid as source grid
        shape = [5,5]
        lon_range = xy[:,0].max() - xy[:,0].min() + 0.025
        lon_min = xy[:,0].min() -0.0125
        lat_range = xy[:,1].max() - xy[:,1].min() + 0.025
        lat_min = xy[:,1].min() -0.0125
        ind = numpy.indices( (shape[0]+1,shape[1]+1))
        lats = numpy.array( ind[1] , dtype=numpy.float)
        lats = lat_min + lats/shape[1] * lat_range
        lons = numpy.array( ind[0] , dtype=numpy.float)
        lons = lon_min + lons/shape[0] * lon_range 

        corners = numpy.array([lons,lats])
        elements = new_structured_grid(shape, corners, axes_names=['lon', 'lat'])
        nodes = StructuredGrid(*ind[0].shape)
        nodes.lat = (lats | units.rad)
        nodes.lon = (lons | units.rad)
        self.structured = StaggeredGrid(elements, nodes)

        corners += 0.01 #shift the grid
        nodes = StructuredGrid(*ind[0].shape)
        nodes.lat = (corners[1] | units.rad)
        nodes.lon = (corners[0] | units.rad)
        self.structured2 = StaggeredGrid(elements, nodes)

    def test1(self):
        target = self.unstructured
        source = self.structured

        #set some values on the grid
        constant_field = numpy.zeros(source.nodes.shape, dtype=numpy.double)
        constant_field += 1.0
        source.nodes.const = constant_field

        constant_field = numpy.zeros(source.elements.shape, dtype=numpy.double)
        constant_field += 1.0
        source.elements.const2 = constant_field

        #create remapper
        remapper = grid_remappers.conservative_spherical_remapper(source, target)

        #remap values        
        remapper.forward_mapping(["const", "const2"])

        self.assertTrue(numpy.all(target.nodes.const >= 0.0), msg="Expecting all remapped values values to be larger than zero")
        self.assertTrue(numpy.all(target.elements.const2 >= 0.0), msg="Expecting all remapped values values to be larger than zero")

    def test2(self):
        target = self.structured2
        source = self.structured

        #set some values on the grid
        constant_field = numpy.zeros(source.nodes.shape, dtype=numpy.double)
        constant_field += 1.0
        source.nodes.const = constant_field

        constant_field = numpy.zeros(source.elements.shape, dtype=numpy.double)
        constant_field += 1.0
        source.elements.const2 = constant_field

        #create remapper
        remapper = grid_remappers.conservative_spherical_remapper(source, target)

        #remap values        
        remapper.forward_mapping(["const", "const2"])

        self.assertTrue(numpy.all(target.nodes.const >= 0.0), msg="Expecting all remapped values values to be larger than zero")
        self.assertTrue(numpy.all(target.elements.const2 >= 0.0), msg="Expecting all remapped values values to be larger than zero")

    def test3(self):
        target = self.unstructured
        source = self.structured

        #set some values on the grid
        constant_field = numpy.zeros(source.nodes.shape, dtype=numpy.double)
        constant_field += 1.0
        source.nodes.const = constant_field

        constant_field = numpy.zeros(source.elements.shape, dtype=numpy.double)
        constant_field += 1.0
        source.elements.const2 = constant_field

        #use a remapping channel
        channel = source.new_remapping_channel_to(target, remapper=grid_remappers.conservative_spherical_remapper)

        #remap values        
        channel.copy_attributes(["const", "const2"])

        self.assertTrue(numpy.all(target.nodes.const >= 0.0), msg="Expecting all remapped values values to be larger than zero")
        self.assertTrue(numpy.all(target.elements.const2 >= 0.0), msg="Expecting all remapped values values to be larger than zero")




class TestGridRemappingChannel(TestCase):
    def setUp(self):
        if not grid_remappers.matplotlib_available:
            self.skip("matplotlib not available")

    def test1(self):
        source=new_cartesian_grid((10,20),1.)
        target=new_cartesian_grid((5,10),2.,offset=[0.,0.5])
        
        channel=source.new_remapping_channel_to(target,grid_remappers.interpolating_2D_remapper)
                
        source.xcopy=source.x
        source.ycopy=source.y
        
        channel.copy_attributes(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)

    def test2(self):
        source=new_cartesian_grid((10,20),1. | units.m)
        target=new_cartesian_grid((5,10),2. | units.m ,offset=[0.,0.5] | units.m)
        
        channel=source.new_remapping_channel_to(target,grid_remappers.interpolating_2D_remapper)
                
        source.xcopy=source.x
        source.ycopy=source.y
        
        channel.copy_attributes(["xcopy","ycopy"])
        self.assertEqual(target.x,target.xcopy)
        self.assertEqual(target.y,target.ycopy)
