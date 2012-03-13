from amuse.community.galactics.interface import GalactICs, GalactICsInterface
from amuse.datamodel import ParticlesWithUnitsConverted

def _new_galactics_model(halo_number_of_particles, unit_system_converter=None, do_scale=False, **keyword_arguments):
    instance = GalactICs(unit_system_converter=unit_system_converter)
    instance.parameters.halo_number_of_particles = halo_number_of_particles
    for (key, value) in keyword_arguments.iteritems():
        setattr(instance.parameters, key, value)
    
    instance.generate_particles()
    result = instance.particles.copy_to_memory()
    instance.stop()
    
    result.move_to_center()
    if do_scale:
        result.scale_to_standard()
    
    if not unit_system_converter is None:
        result = ParticlesWithUnitsConverted(result, unit_system_converter.as_converter_from_si_to_generic())
        result = result.copy_to_memory()
    return result

def _create_docstring():
    _tmp_instance = GalactICs()
    docstring = "\nGalactICs documentation:\n" + GalactICsInterface.__doc__ + \
        "\n" + _tmp_instance.parameters.__doc__ + "\n"
    _tmp_instance.stop()
    return docstring

class _NewModelMetaclass(type):
    def _get_doc(self):
        return _create_docstring()
    __doc__ = property(_get_doc)

class _NewModelWrapper(object):
    __metaclass__ = _NewModelMetaclass
    def _get_doc(self):
        return _create_docstring()
    __doc__ = property(_get_doc)
    def __call__(self, *args, **kwargs):
        return _new_galactics_model(*args, **kwargs)

new_galactics_model = _NewModelWrapper()
