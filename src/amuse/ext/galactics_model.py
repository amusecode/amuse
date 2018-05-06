from amuse.community.galactics.interface import GalactICs, GalactICsInterface
from amuse.community.galactics.gas_interface import GaslactICs, GaslactICsInterface
from amuse.datamodel import ParticlesWithUnitsConverted, Particles
from amuse.datamodel.particles import ParticlesSuperset

def _new_galactics_model(halo_number_of_particles, unit_system_converter=None, do_scale=False, verbose=False, **keyword_arguments):

    code=keyword_arguments.pop("code")
    instance = code(unit_converter=unit_system_converter, redirection="none" if verbose else "null")
    instance.parameters.halo_number_of_particles = halo_number_of_particles
    for (key, value) in keyword_arguments.iteritems():
        setattr(instance.parameters, key, value)
    
    if verbose:
        print "adopted galaxy model parameters:"
        print instance.parameters
    
    instance.generate_particles()
    result = instance.particles.copy()
    
    if hasattr(instance,"gas_particles") and len(instance.gas_particles)>0:
        resultgas=instance.gas_particles.copy()
    else:
        resultgas=Particles()
    
    instance.stop()

    if len(resultgas)>0:
        allpart=ParticlesSuperset([result, resultgas])
    else:
        allpart=result
    allpart.move_to_center()

    # do_scale *is* possible for case with unit converter
    # note that in case of scaling the output galaxy parameters may be very different from the model 
    # parameters input
    if do_scale:
        print "Warning: do_scale for a large galactics model may be very slow"
        if verbose:
            print "Warning: do_scale typically changes the galaxy scale parameters quite a lot from the input parameters"
        if len(resultgas)>0:
            # this is fixable
            raise Exception("scaling of galaxy models with gas currently not possible")        
        allpart.scale_to_standard(convert_nbody=unit_system_converter)

    if not unit_system_converter is None:
        result = ParticlesWithUnitsConverted(result, unit_system_converter.as_converter_from_si_to_generic())
        result = result.copy()
        if len(resultgas)>0:
            resultgas = ParticlesWithUnitsConverted(resultgas, unit_system_converter.as_converter_from_si_to_generic())
            resultgas = resultgas.copy()

    if len(resultgas)>0:
      # resultincludes the gas particles (but not under the same keys)
      return resultgas, result
    else:
      return result

def _create_docstring(code, codeInterface):
    _tmp_instance = code()
    docstring = "\nGalactICs documentation:\n" + codeInterface.__doc__ + \
        "\n" + _tmp_instance.parameters.__doc__ + "\n"
    _tmp_instance.stop()
    return docstring

def _newModelWrapperWrapper(code, codeInterface):
    class _NewModelMetaclass(type):
        def _get_doc(self):
            return _create_docstring(code, codeInterface)
        __doc__ = property(_get_doc)

    class _NewModelWrapper(object):
        __metaclass__ = _NewModelMetaclass
        def _get_doc(self):
            return _create_docstring(code, codeInterface)
        @property
        def __doc__(self):
            return self._get_doc() 
        def __call__(self, *args, **kwargs):
            return _new_galactics_model(*args, code=code,**kwargs)
    return _NewModelWrapper

new_galactics_model = _newModelWrapperWrapper(GalactICs, GalactICsInterface)()

# these return two particle sets (gas, other_particles)..
new_galactics_gas_model = _newModelWrapperWrapper(GaslactICs, GaslactICsInterface)()
new_gaslactics_model = new_galactics_gas_model
