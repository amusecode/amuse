from amuse.community.interface.chem import ChemicalModelingInterface
from amuse.community.interface.chem import ChemicalModeling
from amuse.community.interface.common import CommonCode, CommonCodeInterface
from amuse.community import *
from pathlib import Path

class UCLchemInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
    UCLCHEM: A Gas-Grain Chemical Code for astrochemical modelling

    .. [#] ADS:2017AJ....154...38H (Holdship, J. ; Viti, S, ; Jiménez-Serra, I.; Makrymallis, A. ; Priestley, F. , 2017, AJ)
    """
    def __init__(self, mode = 'cpu', **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker='uclchem_worker',
            **options
        )
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function
    
    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    def recommit_parameters():
        return self.commit_parameters()
    
    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['number_density','temperature','ionrate','radfield']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['number_density','temperature','ionrate','radfield']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('aid', dtype='i', direction=function.IN)
        function.addParameter('abundance', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('aid', dtype='i', direction=function.IN)
        function.addParameter('abundance', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_firstlast_abundance():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('first', dtype='i', direction=function.OUT)
        function.addParameter('last', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_name_of_species():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('name', dtype='s', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_index_of_species():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('name', dtype='s', direction=function.IN)
        function.addParameter('index', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.OUT)
        for x in ['number_density','temperature','ionrate', 'radfield']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function   
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64',direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def run_model():
        function = LegacyFunctionSpecification()
        #function.must_handle_array=True
        #function.addParameter('n_particles', dtype='i', direction=function.LENGTH)
        function.addParameter('dictionary', dtype='s', direction=function.IN)
        #function.addParameter('out_species', dtype='s', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_species():
        function = LegacyFunctionSpecification()
        function.addParameter('species', dtype='s',direction=function.IN)
        function.result_type='i'
        return function

    def _reform_inputs(self,param_dict, out_species):
        if param_dict is None:
            param_dict = {}
        else:
            # lower case (and conveniently copy so we don't edit) the user's dictionary
            # this is key to UCLCHEM's "case insensitivity"
            new_param_dict = {}
            for k, v in param_dict.items():
                assert k.lower() not in new_param_dict, f"Lower case key {k} is already in the dict, stopping"
                if isinstance(v, Path):
                    v = str(v)
                new_param_dict[k.lower()] = v
            param_dict = new_param_dict.copy()
            del new_param_dict
        if out_species is not None:
            n_out = len(out_species)
            param_dict["outspecies"] = n_out
            out_species = " ".join(out_species)
        else:
            out_species = ""
            n_out = 0
        return n_out, param_dict, out_species
    
class UCLchem(CommonCode):
    def __init__(self, convert_nbody=None, **options):
        legacy_interface = UCLchemInterface(**options)
        self.uclchem_time = 0.0|units.yr
        InCodeComponentImplementation.__init__(self,legacy_interface)

    def evolve_model(self,tend):
        assert tend > self.uclchem_time, 'end time must be larger than uclchem_time'
        dictionary, out_species= self._build_dict(tend=tend)
        self.set_species(out_species)
        output = self.run_model(dictionary)
        self.uclchem_time = tend
        return output


    def _build_dict(self, tend):
        dictionary_list = []
        outSpecies = self.out_species
        attributes = self.particles.get_attribute_names_defined_in_store()
        #print(self.particles[0])
        for i in range(len(self.particles.key)):
            dictionary = dict()
            if 'temperature' in attributes:
                dictionary['initialTemp'] = self.particles.temperature.value_in(units.K)[i]
            if 'number_density' in attributes:
                dictionary['initialDens'] = self.particles.number_density.value_in(units.cm**-3)[i]
            if 'ionrate' in attributes:
                dictionary['zeta'] = self.particles.ionrate.value_in(units.cr_ion)[i]
            if 'radfield' in attributes:
                dictionary['radfield'] = self.particles.radfield.value_in(units.habing)[i]
            dictionary['finalTime'] = tend.value_in(units.yr)-self.uclchem_time.value_in(units.yr)
            #dictionary['currentTime'] = self.uclchem_time
            _, dictionary, outSpecies_out = self._reform_inputs(dictionary, outSpecies)
            dictionary_list.append(str(dictionary))
        return dictionary_list, str(outSpecies_out)
    
    def define_parameters(self, handler):
        handler.add_interface_parameter(
            "out_species",
            "Array of molecules to use",
            default_value = ['H','H2']
        )
    
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method(
            "set_state",
            (
                handler.INDEX,
                units.cm**-3,
                units.K,
                units.s**-1,
                units.habing,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state",
            (
                handler.INDEX,
            ),
            (
                units.cm**-3,
                units.K,
                units.s**-1,
                units.habing,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_abundance",
            (
                handler.INDEX,
                handler.INDEX,
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_abundance",
            (
                handler.INDEX,
                handler.INDEX,
                handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_firstlast_abundance",
            (
            ),
            (
                handler.NO_UNIT,
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "new_particle",
            (
                units.cm**-3,
                units.K,
                units.s**-1,
                units.habing
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "delete_particle",
            (
                handler.INDEX,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "evolve_model",
            (
                units.yr
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_time",
            (
            ),
            (
                units.yr,
                handler.ERROR_CODE
            )
        )
        
    
    def define_particle_sets(self, handler):
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_state')
        handler.add_getter('particles', 'get_state')
        handler.add_gridded_getter('particles', 'get_abundance','get_firstlast_abundance', names = ('abundances',))
        handler.add_gridded_setter('particles', 'set_abundance','get_firstlast_abundance', names = ('abundances',))
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','PARAMETER_CHANGE_A','invoke_state_change2')
        handler.add_transition('EDIT','PARAMETER_CHANGE_B','invoke_state_change2')
        handler.add_transition('PARAMETER_CHANGE_A','RUN','recommit_parameters')
        handler.add_transition('PARAMETER_CHANGE_B','EDIT','recommit_parameters')
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_method('RUN', 'evolve_model')
        handler.add_method('RUN', 'get_state')
        handler.add_method('RUN', 'get_abundance')


