from amuse.units import units
from amuse.datamodel import new_cartesian_grid
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase

from amuse.ext import grid_remappers 

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
