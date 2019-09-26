
try:
    import Blender
    import bpy
    import pylab as pl
    from Blender import Mesh
    BLENDER_AVAILABLE = True
except ImportError as ex:
    Mesh = None
    BLENDER_AVAILABLE = False

class Primitives(object):
    def __init__(self):
        pass

    @classmethod
    def sphere(cli, segments, rings, radius):
        me = Mesh.Primitives.UVsphere(segments, rings, radius)
        scn = bpy.data.scenes.active     # link object to current scene
        ob = scn.objects.new(me, 'sphere')
        return ob

    @classmethod
    def cube(cli, radius):
        me = Mesh.Primitives.Cube(radius)
        scn = bpy.data.scenes.active     # link object to current scene
        ob = scn.objects.new(me, 'cube')
        return ob

def Redraw():
    Blender.Redraw()
