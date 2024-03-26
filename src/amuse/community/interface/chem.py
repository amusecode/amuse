from amuse.support.interface import InCodeComponentImplementation
from amuse.units import nbody_system
from amuse.units import generic_unit_converter
from amuse.community.interface import common
from amuse.community import remote_function
from amuse.units import units
from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification

class ChemicalModelingInterface(object):
    @remote_function
    def commit_particles():
        #IN: particle set
        pass

    @remote_function
    def commit_parameters():
        pass

    @remote_function(can_handle_array=True)
    def new_particle(number_density='d'|units.cm**-3,temperature='d'|units.K,ionrate='d'|units.s**-1):
        """
        Add a new particle
        """
        returns (particle_index='i')

    @remote_function(can_handle_array=True)
    def delete_particle(index_of_particle='i'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_abundances(index_of_species='i', name_of_species='s'):
        """
        returns the abundance(s) of (an array of) species
        """
        returns (abundance = 'd')

    @remote_function(can_handle_array=True)
    def set_abundances(index_of_species='i', name_of_species='s', abundance='d'):
        """
        Set abundances for a given (array of) species
        """
        returns ()

    @remote_function
    def set_state(index_of_particle='i',number_density='d'|units.cm**-3,temperature='d'|units.K,ionrate='d'|units.s**-1):
        """
        Set number density, temperature, and CR ionization rate for particle with index
        """
        returns ()

    @remote_function
    def get_state(index_of_particle):
        """
        Get number density, temperature, and CR ionisation rate from particle with index
        """
        returns (number_density='d'|units.cm**-3,temperature='d'|units.K,ionrate='d'|units.s**-1)

    @remote_function(can_handle_array=True)
    def get_index_of_species(name_of_species='s'):
        """
        Get index of species with given name
        """
        returns (index_of_species='i')

    @remote_function(can_handle_array=True)
    def get_name_of_species(index_of_species='i'):
        """
        Get name of species with given index
        """
        returns (name_of_species='s')

    @remote_function
    def evolve_model(time_end='d'):
        returns ()

    @remote_function
    def get_time():
        returns (time='d')

    @remote_function
    def get_firstlast_abundances():
        """
        Get the indices of the first and last species of the network
        """
        returns (first='i', last='i')
        
class ChemicalModelingThermalFeedbackInterface(object):
    """
    Interface for chemical codes that include heating/cooling effects
    """
    @remote_function(can_handle_array=True)
    def set_cooling_methods(cooling_methods='s'):
        """
        turn on chemical cooling methods
        """
        returns ()
    
    @remote_function(can_handle_array=True)
    def set_heating_methods(heating_methods='s'):
        """
        Turn off chemical cooling methods
        """
        returns ()
    
class ChemicalModeling(common.CommonCode):
    def __init__(self,legacy_interface, **options):
        common.CommonCode.__init__(self, legacy_interface,**options)
        #Create dictionary of all species in the network
        #self.species = dict()
        #first, last = self.get_firstlast_abundance()
        #for i in range(first, last+1):
        #    self.species[self.get_name_of_species(i)] = i-1
    
    def define_methods(self, handler):
        common.CommonCode.define_methods(self, handler)
        """
        TODO: Add all the methods that should be available on the user-end of the ChemicalModeling class
        """
        handler.add_method(
            "new_particle",
            (units.cm**-3, units.K, units.s**-1),
            (handler.INDEX, handler.ERROR_CODE)
        )
        handler.add_method(
            "delete_particle",
            (handler.INDEX),
            (handler.ERROR_CODE)
        )
        handler.add_method(
            "get_abundances",
            (handler.INDEX, handler.INDEX),
            (handler.NO_UNIT, handler.ERROR_CODE)
        )
        handler.add_method(
            "set_abundances",
            (handler.INDEX, handler.INDEX, handler.NO_UNIT),
            (handler.ERROR_CODE)
        )
        handler.add_method(
            "get_state",
            (handler.INDEX),
            (units.cm**-3, units.K, units.s**-1, handler.ERROR_CODE)
        )
        handler.add_method(
            "set_state",
            (handler.INDEX, units.cm**-3, units.K, units.s**-1),
            (handler.ERROR_CODE)
        )
        handler.add_method(
            "get_firstlast_abundance",
            (),
            (handler.NO_UNIT, handler.NO_UNIT, handler.ERROR_CODE),
        )
        handler.add_method(
            "evolve_model",
            (units.s),
            (handler.ERROR_CODE)
        )
        handler.add_method(
            "get_time",
            (),
            (units.s, handler.ERROR_CODE)
        )

    def define_particle_sets(self, handler):
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_state')
        handler.add_getter('particles', 'get_state')
        handler.add_gridded_getter('particles', 'get_abundance','get_firstlast_abundance', names = ('abundances',))
        handler.add_gridded_setter('particles', 'set_abundance','get_firstlast_abundance', names = ('abundances',))

