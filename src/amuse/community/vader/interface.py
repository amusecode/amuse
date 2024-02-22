from amuse.community import *
from amuse.community.interface import common


class VaderInterface(CodeInterface,
					 LiteratureReferencesMixIn,
					 common.CommonCodeInterface):

	"""
	VADER is a code simulating the evolution of viscous thin accretion
	disks. It is developed by Mark Krumholz and John Forbes [1]. 

    .. [#] ADS:2015A&C....11....1K (Krumholz, M. R. and Forbes, J. C., Astronomy and Computing, Vol. 11 (2015):
    .. [#] ... VADER: A Flexible, Robust, Open-Source Code for Simulating Viscous Thin Accretion Disks)
	"""
		
	include_headers = ['worker_code.h']

	def __init__(self, mode = 'none', **keyword_arguments):
		CodeInterface.__init__(self,
			name_of_the_worker="vader_worker_" + mode,
			**keyword_arguments)
		LiteratureReferencesMixIn.__init__(self)

	@legacy_function
	def initialize_code():
		function = LegacyFunctionSpecification()
		function.result_type = 'int32'
		return function

	@legacy_function
	def evolve_model():
		function = LegacyFunctionSpecification()
		function.addParameter('tlim', dtype='float64',
							  direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def initialize_keplerian_grid():
		function = LegacyFunctionSpecification()
		function.addParameter('n', dtype='int32',
							  direction=function.IN)
		function.addParameter('linear', dtype='bool',
							  direction=function.IN)
		function.addParameter('rmin', dtype='float64', 
							  direction=function.IN)
		function.addParameter('rmax', dtype='float64', 
							  direction=function.IN)
		function.addParameter('m', dtype='float64', 
							  direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def update_keplerian_grid():
		function = LegacyFunctionSpecification()
		function.addParameter('m', dtype='float64',
							  direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def initialize_flat_grid():
		function = LegacyFunctionSpecification()
		function.addParameter('n', dtype='int32', 
							  direction=function.IN)
		function.addParameter('linear', dtype='bool', 
							  direction=function.IN)
		function.addParameter('rmin', dtype='float64', 
							  direction=function.IN)
		function.addParameter('rmax', dtype='float64', 
							  direction=function.IN)
		function.addParameter('vphi', dtype='float64', 
							  direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def update_flat_grid():
		function = LegacyFunctionSpecification()
		function.addParameter('vphi', dtype='float64',
							  direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def initialize_tabulated_grid():
		function = LegacyFunctionSpecification()
		function.addParameter('n', dtype='int32', 
							  direction=function.IN)
		function.addParameter('linear', dtype='int32', 
							  direction=function.IN)
		function.addParameter('rmin', dtype='float64', 
							  direction=function.IN)
		function.addParameter('rmax', dtype='float64', 
							  direction=function.IN)
		function.addParameter('bspline_degree', dtype='int32', 
							  direction=function.IN)
		function.addParameter('bspline_breakpoints', 	
							  dtype='int32',
							  direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def update_tabulated_grid():
		function = LegacyFunctionSpecification()
		function.addParameter('bspline_degree', dtype='int32',
							  direction=function.IN)
		function.addParameter('bspline_breakpoints', dtype='int32',
							  direction=function.IN)
		function.result_type = 'int32'
		return function


	@legacy_function
	def get_tabulated_size():
		function = LegacyFunctionSpecification()
		function.addParameter('nTab', dtype='int32', direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_tabulated_size():
		function = LegacyFunctionSpecification()
		function.addParameter('nTab', dtype='int32', direction=function.IN)
		function.result_type = 'int32'
		return function

	def get_table_range(self):
		nTab = self.get_tabulated_size()['nTab']
		return (0, nTab-1)

	@legacy_function
	def get_tabulated_radius():
		function = LegacyFunctionSpecification()
		function.can_handle_array = True
		function.addParameter('i', dtype='int32', direction=function.IN)
		function.addParameter('rTab', dtype='float64', direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_tabulated_radius():
		function = LegacyFunctionSpecification()
		function.can_handle_array = True
		function.addParameter('i', dtype='int32', direction=function.IN)
		function.addParameter('rTab', dtype='float64', direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_tabulated_velocity():
		function = LegacyFunctionSpecification()
		function.can_handle_array = True
		function.addParameter('i', dtype='int32', direction=function.IN)
		function.addParameter('vTab', dtype='float64', direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_tabulated_velocity():
		function = LegacyFunctionSpecification()
		function.can_handle_array = True
		function.addParameter('i', dtype='int32', direction=function.IN)
		function.addParameter('vTab', dtype='float64', direction=function.IN)
		function.result_type = 'int32'
		return function

	#grid getters&setters
	@legacy_function
	def get_position_of_index():
		function = LegacyFunctionSpecification()
		function.can_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('r', dtype='float64',
							  direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_index_of_position():
		function = LegacyFunctionSpecification()
		function.addParameter('r', dtype='float64', 
							  direction=function.IN)
		function.addParameter('i', dtype='int32',
							  direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_grid_column_density():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('sigma', dtype='float64',
							  direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_grid_pressure():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('pressure', dtype='float64',
							  direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_grid_internal_energy():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('internal_energy', dtype='float64',
							  direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_grid_state():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
							   direction=function.IN)
		for x in ['sigma', 'pressure', 'internal_energy']:
			function.addParameter(x, dtype='float64',
								   direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_grid_user_output():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('n', dtype='int32',
							   direction=function.IN)
		function.addParameter('i', dtype='int32',
							   direction=function.IN)
		function.addParameter('user_output', dtype='float64',
							   direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function


	@legacy_function
	def set_grid_column_density():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('sigma', dtype='float64',
							  direction=function.IN)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_grid_pressure():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('pressure', dtype='float64',
							  direction=function.IN)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_grid_internal_energy():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32', 
							  direction=function.IN)
		function.addParameter('internal_energy', dtype='float64',
							  direction=function.IN)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_grid_state():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
							   direction=function.IN)
		for x in ['sigma', 'pressure', 'internal_energy']:
			function.addParameter(x, dtype='float64',
								   direction=function.IN)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_grid_user_output():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('n', dtype='int32',
							   direction=function.IN)
		function.addParameter('i', dtype='int32',
							   direction=function.IN)
		function.addParameter('user_output', dtype='float64',
							   direction=function.IN)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	#viscous parameters getters&setters
	@legacy_function
	def get_alpha_function():
		function = LegacyFunctionSpecification()
		function.addParameter('alpha_func', dtype='bool',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_alpha_function():
		function = LegacyFunctionSpecification()
		function.addParameter('alpha_func', dtype='bool',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_alpha():
		function = LegacyFunctionSpecification()
		function.addParameter('alpha', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_alpha():
		function = LegacyFunctionSpecification()
		function.addParameter('alpha', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_eos_function():
		function = LegacyFunctionSpecification()
		function.addParameter('eos_func', dtype='bool',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_eos_function():
		function = LegacyFunctionSpecification()
		function.addParameter('eos_func', dtype='bool',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_gamma():
		function = LegacyFunctionSpecification()
		function.addParameter('gamma', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_gamma():
		function = LegacyFunctionSpecification()
		function.addParameter('gamma', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_delta():
		function = LegacyFunctionSpecification()
		function.addParameter('delta', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_delta():
		function = LegacyFunctionSpecification()
		function.addParameter('delta', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	#source getters&setters
	@legacy_function
	def get_mass_source_function():
		function = LegacyFunctionSpecification()
		function.addParameter('mass_source_func', dtype='bool',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_mass_source_function():
		function = LegacyFunctionSpecification()
		function.addParameter('mass_source_func', dtype='bool',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_internal_energy_source_function():
		function = LegacyFunctionSpecification()
		function.addParameter('internal_energy_source_func', 
						   dtype='bool', direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_internal_energy_source_function():
		function = LegacyFunctionSpecification()
		function.addParameter('internal_energy_source_func', 
						   dtype='bool', direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_mass_source_value():
		function = LegacyFunctionSpecification()
		function.addParameter('mass_source_value', 
						dtype='float64', direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_mass_source_value():
		function = LegacyFunctionSpecification()
		function.addParameter('mass_source_value', 
						dtype='float64', direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_internal_energy_source_value():
		function = LegacyFunctionSpecification()
		function.addParameter('internal_energy_source_value', 
							   dtype='float64', 
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_internal_energy_source_value():
		function = LegacyFunctionSpecification()
		function.addParameter('internal_energy_source_value', 
							   dtype='float64', 
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	#boundary getters&setters
	@legacy_function
	def get_inner_pressure_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres', dtype='int32',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_pressure_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres', dtype='int32',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_pressure_boundary_mass_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_pressure_boundary_torque_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_pressure_boundary_torque():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_pressure_boundary_mass_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_pressure_boundary_torque_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_pressure_boundary_torque():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_pres_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_enthalpy_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_enth', dtype='int32',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_enthalpy_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_enth', dtype='int32',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_enthalpy_boundary_enthalpy():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_enth_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_enthalpy_boundary_enthalpy_gradient():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_enth_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_enthalpy_boundary_enthalpy():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_enth_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_enthalpy_boundary_enthalpy_gradient():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_enth_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_boundary_function():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_func', dtype='bool',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_inner_boundary_function():
		function = LegacyFunctionSpecification()
		function.addParameter('ibc_func', dtype='bool',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_pressure_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres', dtype='int32',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_pressure_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres', dtype='int32',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_pressure_boundary_mass_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_pressure_boundary_torque_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_pressure_boundary_torque():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_pressure_boundary_mass_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_pressure_boundary_torque_flux():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_pressure_boundary_torque():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_pres_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_enthalpy_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_enth', dtype='int32',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_enthalpy_boundary_type():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_enth', dtype='int32',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_enthalpy_boundary_enthalpy():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_enth_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_enthalpy_boundary_enthalpy_gradient():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_enth_val', dtype='float64',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_enthalpy_boundary_enthalpy():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_enth_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_enthalpy_boundary_enthalpy_gradient():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_enth_val', dtype='float64',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_boundary_function():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_func', dtype='bool',
							   direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_outer_boundary_function():
		function = LegacyFunctionSpecification()
		function.addParameter('obc_func', dtype='bool',
							   direction=function.IN)
		function.result_type = 'int32'
		return function

	#various
	@legacy_function
	def get_number_of_cells():
		function = LegacyFunctionSpecification()
		function.addParameter('i', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	def get_index_range_for_grid(self):
		i, error = self.get_number_of_cells()
		return (0, i-1)

	def get_index_range_for_user_outputs(self):
		i, error = self.get_number_of_cells()
		j, error = self.get_nUserOut()
		return (0, j-1, 0, i-1)

	@legacy_function
	def get_area_of_index():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('area', dtype='float64',
								direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
								function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_rotational_velocity_of_index():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('vphi', 
								dtype='float64',
								direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
								function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_effective_potential_of_index():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('psiEff', 
								dtype='float64',
								direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
								function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_gravitational_potential_of_index():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('psi_grav', 
								dtype='float64',
								direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
								function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_number_of_user_parameters():
		function = LegacyFunctionSpecification()
		function.addParameter('n', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_number_of_user_parameters():
		function = LegacyFunctionSpecification()
		function.addParameter('n', dtype='int32',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_parameter():
		function = LegacyFunctionSpecification()
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('param', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_parameter():
		function = LegacyFunctionSpecification()
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('param', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_time():
		function = LegacyFunctionSpecification()
		function.addParameter('time', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_boundary_mass_out():
		function = LegacyFunctionSpecification()
		function.addParameter('mSrcOut', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_boundary_mass_out():
		function = LegacyFunctionSpecification()
		function.addParameter('mSrcOut', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_inner_boundary_energy_out():
		function = LegacyFunctionSpecification()
		function.addParameter('eSrcOut', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_outer_boundary_energy_out():
		function = LegacyFunctionSpecification()
		function.addParameter('eSrcOut', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_mass_source_out():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('mSrcOut', dtype='float64',
								direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							 	function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_energy_source_out():
		function = LegacyFunctionSpecification()
		function.must_handle_array = True
		function.addParameter('i', dtype='int32',
								direction=function.IN)
		function.addParameter('eSrcOut', dtype='float64',
								direction=function.OUT)
		function.addParameter('number_of_points', 'int32',
							  function.LENGTH)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_dtStart():
		function = LegacyFunctionSpecification()
		function.addParameter('dtStart', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_dtStart():
		function = LegacyFunctionSpecification()
		function.addParameter('dtStart', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_dtMin():
		function = LegacyFunctionSpecification()
		function.addParameter('dtMin', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_dtMin():
		function = LegacyFunctionSpecification()
		function.addParameter('dtMin', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_dtTol():
		function = LegacyFunctionSpecification()
		function.addParameter('dtTol', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_dtTol():
		function = LegacyFunctionSpecification()
		function.addParameter('dtTol', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_errTol():
		function = LegacyFunctionSpecification()
		function.addParameter('errTol', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_errTol():
		function = LegacyFunctionSpecification()
		function.addParameter('errTol', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_maxDtIncrease():
		function = LegacyFunctionSpecification()
		function.addParameter('maxDtIncrease', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_maxDtIncrease():
		function = LegacyFunctionSpecification()
		function.addParameter('maxDtIncrease', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_maxIter():
		function = LegacyFunctionSpecification()
		function.addParameter('maxIter', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_maxIter():
		function = LegacyFunctionSpecification()
		function.addParameter('maxIter', dtype='int32',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_interpOrder():
		function = LegacyFunctionSpecification()
		function.addParameter('interpOrder', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_interpOrder():
		function = LegacyFunctionSpecification()
		function.addParameter('interpOrder', dtype='int32',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_maxStep():
		function = LegacyFunctionSpecification()
		function.addParameter('maxStep', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_maxStep():
		function = LegacyFunctionSpecification()
		function.addParameter('maxStep', dtype='int32',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_useBE():
		function = LegacyFunctionSpecification()
		function.addParameter('useBE', dtype='bool',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_useBE():
		function = LegacyFunctionSpecification()
		function.addParameter('useBE', dtype='bool',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_verbosity():
		function = LegacyFunctionSpecification()
		function.addParameter('verbosity', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_verbosity():
		function = LegacyFunctionSpecification()
		function.addParameter('verbosity', dtype='int32',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_PreTimestep():
		function = LegacyFunctionSpecification()
		function.addParameter('PreTimestep', dtype='bool',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_PreTimestep():
		function = LegacyFunctionSpecification()
		function.addParameter('PreTimestep', dtype='bool',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_PostTimestep():
		function = LegacyFunctionSpecification()
		function.addParameter('PostTimestep', dtype='bool',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_PostTimestep():
		function = LegacyFunctionSpecification()
		function.addParameter('PostTimestep', dtype='bool',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_nUserOut():
		function = LegacyFunctionSpecification()
		function.addParameter('nUserOut', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_nUserOut():
		function = LegacyFunctionSpecification()
		function.addParameter('nUserOut', dtype='int32',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	@legacy_function
	def get_begin_time():
		function = LegacyFunctionSpecification()
		function.addParameter('begin_time', dtype='float64',
								direction=function.OUT)
		function.result_type = 'int32'
		return function

	@legacy_function
	def set_begin_time():
		function = LegacyFunctionSpecification()
		function.addParameter('begin_time', dtype='float64',
								direction=function.IN)
		function.result_type = 'int32'
		return function

	'''
	@legacy_function
	def get_nFail():
		function = LegacyFunctionSpecification()
		function.addParameter('nFail', dtype='int32',
								direction=function.OUT)
		function.result_type = 'int32'
		return function
	'''

	@legacy_function
	def cleanup_code():
		function = LegacyFunctionSpecification()
		function.result_type = 'int32'
		return function


class Vader(common.CommonCode):

	def __init__(self, unit_converter = None, mode = 'none', **options):

		self.mode = mode
		self.unit_converter = unit_converter

		common.CommonCode.__init__(self,  VaderInterface(mode = mode, **options), 
			**options)


	def define_state(self, handler): 
		common.CommonCode.define_state(self, handler)
		handler.add_transition('END', 'INITIALIZED', 
							'initialize_code', False)

		handler.add_transition('INITIALIZED','EDIT','commit_parameters')

		handler.add_transition('RUN', 'CHANGE_PARAMETERS_RUN', 
							   'before_set_parameter', False)
		handler.add_transition('EDIT', 'CHANGE_PARAMETERS_EDIT', 
							   'before_set_parameter', False)

		handler.add_transition('CHANGE_PARAMETERS_RUN', 'RUN', 
							   'recommit_parameters')
		handler.add_transition('CHANGE_PARAMETERS_EDIT', 'EDIT',
							   'recommit_parameters')

		handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
		handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
		
		handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
		handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')

		handler.add_method('RUN', 'before_get_parameter')
		handler.add_method('EDIT', 'before_get_parameter')

		handler.add_transition('EDIT', 'RUN', 
							   'initialize_keplerian_grid', False)
		handler.add_transition('EDIT', 'RUN', 
							   'initialize_flat_grid', False)
		handler.add_transition('EDIT', 'RUN', 
							   'initialize_tabulated_grid', False)

		handler.add_method('RUN', 'update_keplerian_grid')
		handler.add_method('RUN', 'update_flat_grid')
		handler.add_method('RUN', 'update_tabulated_grid')


		grid_write_properties = ['state', 'column_density', 'pressure',
			'internal_energy', 'user_output']

		for prop in grid_write_properties:
			handler.add_method('RUN', 'set_grid_'+prop)
			handler.add_method('RUN', 'get_grid_'+prop)
			handler.add_method('EDIT', 'set_grid_'+prop)
			handler.add_method('EDIT', 'get_grid_'+prop)


		grid_read_properties = ['mass_source_out', 'energy_source_out', 
			'position_of_index', 'area_of_index', 'effective_potential_of_index',
			'gravitational_potential_of_index', 'index_range_for_grid',
			'rotational_velocity_of_index']

		for prop in grid_read_properties:
			handler.add_method('RUN', 'get_'+prop)
			handler.add_method('EDIT', 'get_'+prop)


		code_properties = ['inner_boundary_mass_out', 'outer_boundary_mass_out',
			'inner_boundary_energy_out', 'outer_boundary_energy_out']

		for prop in code_properties:
			handler.add_method('RUN', 'get_'+prop)
			handler.add_method('EDIT', 'get_'+prop)


		handler.add_method('RUN', 'evolve_model')

		handler.add_method('EDIT', 'get_time')
		handler.add_method('RUN', 'get_time')


	def commit_parameters(self):
		self.parameters.send_not_set_parameters_to_code()
		self.parameters.send_cached_parameters_to_code()
		self.overridden().commit_parameters()


	def define_methods(self, builder):

		length = units.cm
		time   = units.s
		mass   = units.g

		velocity 	 = length / time
		energy   	 = mass * velocity**2.
		pressure 	 = mass / length / time**2.
		force		 = mass * velocity / time
		col_density  = mass / length**2.
		Eint_density = energy / length**2.


		builder.add_method(
			"evolve_model",
			(time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"initialize_code",
			(),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"initialize_keplerian_grid",
			(builder.NO_UNIT, builder.NO_UNIT, 
			 length, length, mass,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"update_keplerian_grid",
			(mass,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"initialize_flat_grid",
			(builder.NO_UNIT, builder.NO_UNIT,
			 length, length, velocity,),
			(builder.ERROR_CODE)
		)

		builder.add_method(
			"update_flat_grid",
			(velocity,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"initialize_tabulated_grid",
			(builder.NO_UNIT, builder.NO_UNIT,
			 length, length,
			 builder.NO_UNIT, builder.NO_UNIT,),
			(builder.ERROR_CODE)
		)

		builder.add_method(
			"update_tabulated_grid",
			(builder.NO_UNIT, builder.NO_UNIT),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_tabulated_size",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_tabulated_size",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_tabulated_radius",
			(builder.INDEX,),
			(length, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_tabulated_radius",
			(builder.INDEX, length),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_tabulated_velocity",
			(builder.INDEX,),
			(velocity, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_tabulated_velocity",
			(builder.INDEX, velocity),
			(builder.ERROR_CODE,)
		)

		'''
		builder.add_method(
			"initialize_grid",
			(builder.NO_UNIT, builder.NO_UNIT,
			 length, length, velocity, velocity
			 #beta units
			 #beta units
			 energy/mass,
			 energy/mass,
			 #whatever the hell g_h is
			 ,),
			(builder.ERROR_CODE)
		)
		'''

		#grid getters&setters
		builder.add_method(
			"get_position_of_index",
			(builder.INDEX,),
			(length, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_index_of_position",
			(length,),
			(builder.INDEX, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_area_of_index",
			(builder.INDEX,),
			(length**2., builder.ERROR_CODE)
		)

		builder.add_method(
			"get_rotational_velocity_of_index",
			(builder.INDEX,),
			(length / time, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_effective_potential_of_index",
			(builder.INDEX,),
			(energy / mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_gravitational_potential_of_index",
			(builder.INDEX,),
			(energy / mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_grid_column_density",
			(builder.INDEX,),
			(col_density, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_grid_pressure",
			(builder.INDEX,),
			(pressure * length, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_grid_internal_energy",
			(builder.INDEX,),
			(Eint_density, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_grid_state",
			(builder.INDEX,),
			(col_density, pressure * length, Eint_density, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_grid_user_output",
			(builder.INDEX, builder.INDEX),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_grid_column_density",
			(builder.INDEX, col_density),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_grid_pressure",
			(builder.INDEX, pressure * length),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_grid_internal_energy",
			(builder.INDEX, Eint_density),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_grid_state",
			(builder.INDEX, col_density, pressure * length, Eint_density),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_grid_user_output",
			(builder.INDEX, builder.INDEX, builder.NO_UNIT),
			(builder.ERROR_CODE,)
		)

		#viscous parameters getters&setters
		builder.add_method(
			"get_alpha_function",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_alpha_function",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_alpha",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_alpha",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_eos_function",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_eos_function",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_gamma",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_gamma",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_delta",
			(),
			(velocity**2, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_delta",
			(velocity**2,),
			(builder.ERROR_CODE,)
		)

		#source getters&setters
		builder.add_method(
			"get_mass_source_function",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_mass_source_function",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_internal_energy_source_function",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_internal_energy_source_function",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_mass_source_value",
			(),
			(col_density / time, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_mass_source_value",
			(col_density / time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_internal_energy_source_value",
			(),
			(Eint_density / time, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_internal_energy_source_value",
			(Eint_density / time,),
			(builder.ERROR_CODE,)
		)

		#boundary condition getters&setters
		builder.add_method(
			"get_inner_pressure_boundary_type",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_inner_pressure_boundary_type",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_inner_pressure_boundary_mass_flux",
			(),
			(mass / time, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_inner_pressure_boundary_torque_flux",
			(),
			(force * length / time, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_inner_pressure_boundary_torque",
			(),
			(force * length, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_inner_pressure_boundary_mass_flux",
			(mass / time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_inner_pressure_boundary_torque_flux",
			(force * length / time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_inner_pressure_boundary_torque",
			(force * length,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_inner_enthalpy_boundary_type",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_inner_enthalpy_boundary_type",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_inner_enthalpy_boundary_enthalpy",
			(),
			(energy / mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_inner_enthalpy_boundary_enthalpy_gradient",
			(),
			(energy / length / mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_inner_enthalpy_boundary_enthalpy",
			(energy / mass,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_inner_enthalpy_boundary_enthalpy_gradient",
			(energy / length / mass,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_inner_boundary_function",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_inner_boundary_function",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_outer_pressure_boundary_type",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_outer_pressure_boundary_type",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_outer_pressure_boundary_mass_flux",
			(),
			(mass / time, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_outer_pressure_boundary_torque_flux",
			(),
			(force * length / time, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_outer_pressure_boundary_torque",
			(),
			(force * length, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_outer_pressure_boundary_mass_flux",
			(mass / time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_outer_pressure_boundary_torque_flux",
			(force * length / time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_outer_pressure_boundary_torque",
			(force * length,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_outer_enthalpy_boundary_type",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_outer_enthalpy_boundary_type",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_outer_enthalpy_boundary_enthalpy",
			(),
			(energy / mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_outer_enthalpy_boundary_enthalpy_gradient",
			(),
			(energy / length / mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_outer_enthalpy_boundary_enthalpy",
			(energy / mass,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_outer_enthalpy_boundary_enthalpy_gradient",
			(energy / length / mass,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_outer_boundary_function",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_outer_boundary_function",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		#other outputs
		builder.add_method(
			"get_inner_boundary_mass_out",
			(),
			(mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_outer_boundary_mass_out",
			(),
			(mass, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_inner_boundary_energy_out",
			(),
			(energy, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_outer_boundary_energy_out",
			(),
			(energy, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_mass_source_out",
			(builder.INDEX,),
			(col_density, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_energy_source_out",
			(builder.INDEX,),
			(Eint_density, builder.ERROR_CODE)
		)

		#control parameters
		builder.add_method(
			"get_dtStart",
			(),
			(time, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_dtStart",
			(time,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_dtMin",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_dtMin",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_dtTol",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_dtTol",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_errTol",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_errTol",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_maxDtIncrease",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_maxDtIncrease",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_maxIter",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_maxIter",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_interpOrder",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_interpOrder",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_maxStep",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_maxStep",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_useBE",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_useBE",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_verbosity",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_verbosity",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_PreTimestep",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_PreTimestep",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_PostTimestep",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_PostTimestep",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		#various
		builder.add_method(
			"get_time",
			(),
			(time, builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_number_of_user_parameters",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_number_of_user_parameters",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_max_index",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"get_nUserOut",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)

		builder.add_method(
			"set_nUserOut",
			(builder.NO_UNIT,),
			(builder.ERROR_CODE,)
		)

		builder.add_method(
			"get_begin_time",
			(),
			(time, builder.ERROR_CODE,)
		)

		builder.add_method(
			"set_begin_time",
			(time,),
			(builder.ERROR_CODE,)
		)

		'''
		builder.add_method(
			"get_nFail",
			(),
			(builder.NO_UNIT, builder.ERROR_CODE)
		)
		'''

	def define_parameters(self, builder):

		length = units.cm
		time   = units.s
		mass   = units.g

		velocity 	 = length / time
		energy   	 = mass * velocity**2.
		pressure 	 = mass / length / time**2.
		force		 = mass * velocity / time
		col_density  = mass / length**2.
		Eint_density = energy / length**2.

		builder.add_method_parameter(
			"get_alpha_function",
			"set_alpha_function",
			"alpha_function",
			"True if user-defined function is to be used for the viscosity parameter",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_alpha",
			"set_alpha",
			"alpha",
			"viscosity parameter, see Krumholz 2015",
			default_value = 1.,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_eos_function",
			"set_eos_function",
			"equation_of_state_function",
			"nonzero if user-defined function is to be used for the equation of state",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_gamma",
			"set_gamma",
			"gamma",
			"equation of state parameter, see Krumholz 2015",
			default_value = 1.000001,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_delta",
			"set_delta",
			"delta",
			"equation of state parameter, see Krumholz 2015",
			default_value = 0. | velocity**2,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_pressure_boundary_type",
			"set_inner_pressure_boundary_type",
			"inner_pressure_boundary_type",
			"type for pressure inner boundary condition \n1: fixed mass flux \n2: fixed torque flux \n3: fixed torque",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_pressure_boundary_mass_flux",
			"set_inner_pressure_boundary_mass_flux",
			"inner_pressure_boundary_mass_flux",
			"constant value for pressure inner boundary condition, ignored if inner boundary function is nonzero",
			default_value = 0. | mass / time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_pressure_boundary_torque_flux",
			"set_inner_pressure_boundary_torque_flux",
			"inner_pressure_boundary_torque_flux",
			"constant value for pressure inner boundary condition, ignored if inner boundary function is nonzero",
			default_value = 0. | force * length / time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_pressure_boundary_torque",
			"set_inner_pressure_boundary_torque",
			"inner_pressure_boundary_torque",
			"constant value for pressure inner boundary condition, ignored if inner boundary function is nonzero",
			default_value = 0. | force * length,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_enthalpy_boundary_type",
			"set_inner_enthalpy_boundary_type",
			"inner_enthalpy_boundary_type",
			"type for enthalpy inner boundary condition \n1: fixed enthalpy value \n2: fixed enthalpy gradient",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_enthalpy_boundary_enthalpy",
			"set_inner_enthalpy_boundary_enthalpy",
			"inner_enthalpy_boundary_enthalpy",
			"constant value for enthalpy inner boundary condition, ignored if inner boundary function is nonzero",
			default_value = 0. | energy / mass,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_enthalpy_boundary_enthalpy_gradient",
			"set_inner_enthalpy_boundary_enthalpy_gradient",
			"inner_enthalpy_boundary_enthalpy_gradient",
			"constant value for enthalpy inner boundary condition, ignored if inner boundary function is nonzero",
			default_value = 0. | energy / length / mass,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_inner_boundary_function",
			"set_inner_boundary_function",
			"inner_boundary_function",
			"nonzero if user-defined function is to be used for the inner boundary",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_pressure_boundary_type",
			"set_outer_pressure_boundary_type",
			"outer_pressure_boundary_type",
			"type for pressure outer boundary condition \n1: fixed mass flux \n2: fixed torque flux \n3: fixed torque",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_pressure_boundary_mass_flux",
			"set_outer_pressure_boundary_mass_flux",
			"outer_pressure_boundary_mass_flux",
			"constant value for pressure outer boundary condition, ignored if outer boundary function is nonzero",
			default_value = 0. | mass / time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_pressure_boundary_torque_flux",
			"set_outer_pressure_boundary_torque_flux",
			"outer_pressure_boundary_torque_flux",
			"constant value for pressure outer boundary condition, ignored if outer boundary function is nonzero",
			default_value = 0. | force * length / time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_pressure_boundary_torque",
			"set_outer_pressure_boundary_torque",
			"outer_pressure_boundary_torque",
			"constant value for pressure outer boundary condition, ignored if outer boundary function is nonzero",
			default_value = 0. | force * length,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_enthalpy_boundary_type",
			"set_outer_enthalpy_boundary_type",
			"outer_enthalpy_boundary_type",
			"type for enthalpy outer boundary condition \n1: fixed enthalpy value \n2: fixed enthalpy gradient",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_enthalpy_boundary_enthalpy",
			"set_outer_enthalpy_boundary_enthalpy",
			"outer_enthalpy_boundary_enthalpy",
			"constant value for enthalpy outer boundary condition, ignored if outer boundary function is nonzero",
			default_value = 0. | energy / mass,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_enthalpy_boundary_enthalpy_gradient",
			"set_outer_enthalpy_boundary_enthalpy_gradient",
			"outer_enthalpy_boundary_enthalpy_gradient",
			"constant value for enthalpy outer boundary condition, ignored if outer boundary function is nonzero",
			default_value = 0. | energy / length / mass,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_outer_boundary_function",
			"set_outer_boundary_function",
			"outer_boundary_function",
			"nonzero if user-defined function is to be used for the outer boundary",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_mass_source_function",
			"set_mass_source_function",
			"mass_source_function",
			"nonzero if user-defined function is to be used for mass source",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_mass_source_value",
			"set_mass_source_value",
			"mass_source_value",
			"constant mass source value, ignored if mass source function is nonzero",
			default_value = 0. | col_density / time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_internal_energy_source_function",
			"set_internal_energy_source_function",
			"internal_energy_source_function",
			"nonzero if user-defined function is to be used for internal energy source",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_internal_energy_source_value",
			"set_internal_energy_source_value",
			"internal_energy_source_value",
			"constant internal energy source value, ignored if energy source function is nonzero",
			default_value = 0. | Eint_density / time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_number_of_user_parameters",
			"set_number_of_user_parameters",
			"number_of_user_parameters",
			"the number of user-defined parameters for use in custom vader problems",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_dtStart",
			"set_dtStart",
			"initial_timestep",
			"the initial timestep of the simulation",
			default_value = 1. | time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_dtMin",
			"set_dtMin",
			"minimum_timestep",
			"the minimum timestep of the simulation",
			default_value = 1E-15,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_dtTol",
			"set_dtTol",
			"maximum_tolerated_change",
			"the allowed relative change in column density/pressure/internal energy",
			default_value = 0.1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_errTol",
			"set_errTol",
			"error_tolerance",
			"implicit solver tolerance",
			default_value = 1E-6,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_maxDtIncrease",
			"set_maxDtIncrease",
			"maximum_timestep_increase",
			"the maximum allowed relative time increase",
			default_value = 1.5,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_maxIter",
			"set_maxIter",
			"maximum_iterations",
			"the maximum number of implicit iterations",
			default_value = 40,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_interpOrder",
			"set_interpOrder",
			"interpolation_order",
			"the interpolation order; piecewise 1) constant 2) linear 3) parabolic",
			default_value = 2,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_maxStep",
			"set_maxStep",
			"maximum_steps",
			"the maximum number of steps. negative indicates unbound",
			default_value = -1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_useBE",
			"set_useBE",
			"use_backwards_euler",
			"uses backwards Euler integration if true, Crank Nicolson if false",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_verbosity",
			"set_verbosity",
			"verbosity",
			"the verbosity of the simulation. 0 is silent, 3 very loud. note that this is only shown on console if redirection='none' is added to the vader initialization",
			default_value = 0,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_PreTimestep",
			"set_PreTimestep",
			"pre_timestep_function",
			"if True, the user-defined PreTimestep function is used",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_PostTimestep",
			"set_PostTimestep",
			"post_timestep_function",
			"if True, the user-defined PostTimestep function is used",
			default_value = False,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_nUserOut",
			"set_nUserOut",
			"number_of_user_outputs",
			"number of user-defined grid properties",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_begin_time",
			"set_begin_time",
			"begin_time",
			"time to start integration at",
			default_value = 0.|time,
			must_set_before_get = True
		)

		builder.add_method_parameter(
			"get_tabulated_size",
			"set_tabulated_size",
			"table_size",
			"length of the interpolation table",
			default_value = 1,
			must_set_before_get = True
		)

		builder.add_array_parameter(
			"get_tabulated_radius",
			"set_tabulated_radius",
			"get_table_range",
			"radius_table",
			"radii of the velocity curve interpolation table"
		)

		builder.add_array_parameter(
			"get_tabulated_velocity",
			"set_tabulated_velocity",
			"get_table_range",
			"velocity_table",
			"velocities of the velocity curve interpolation table"
		)


	def define_properties(self, builder):
		builder.add_property('get_time', 
					public_name = 'model_time')
		#builder.add_property('get_nFail',
		#			public_name = 'failed_iterations')
		builder.add_property('get_inner_boundary_mass_out',
					public_name = 'inner_boundary_mass_out')
		builder.add_property('get_outer_boundary_mass_out',
					public_name = 'outer_boundary_mass_out')
		builder.add_property('get_inner_boundary_energy_out',
					public_name = 'inner_boundary_energy_out')
		builder.add_property('get_outer_boundary_energy_out',
					public_name = 'outer_boundary_energy_out')


	def define_particle_sets(self, builder):
		builder.define_grid('grid')
		builder.set_grid_range('grid',
							   'get_index_range_for_grid')
		builder.add_getter('grid', 'get_position_of_index', 
						   names=('r',))
		builder.add_getter('grid', 'get_area_of_index',
						   names=('area',))
		builder.add_getter('grid', 
							'get_rotational_velocity_of_index',
						   names=('rotational_velocity',))
		builder.add_getter('grid',
							'get_effective_potential_of_index',
						   names=('effective_potential',))
		builder.add_getter('grid',
						'get_gravitational_potential_of_index',
						   names=('gravitational_potential',))

		builder.add_getter('grid', 'get_grid_column_density', 
						   names=('column_density',))
		builder.add_setter('grid', 'set_grid_column_density', 
						   names=('column_density',))

		builder.add_getter('grid', 'get_grid_pressure', 
						   names=('pressure',))
		builder.add_setter('grid', 'set_grid_pressure', 
						   names=('pressure',))
		
		builder.add_getter('grid', 'get_grid_internal_energy', 
						   names=('internal_energy',))
		builder.add_setter('grid', 'set_grid_internal_energy', 
						   names=('internal_energy',))

		builder.add_getter('grid', 'get_mass_source_out',
						   names=('mass_source_difference',))
		builder.add_getter('grid', 'get_energy_source_out',
				names=('internal_energy_source_difference',))


		builder.define_grid('grid_user')
		builder.set_grid_range('grid_user',
							   'get_index_range_for_user_outputs')

		builder.add_getter('grid_user', 'get_grid_user_output', 
						   names=('value',))
		builder.add_setter('grid_user', 'set_grid_user_output', 
						   names=('value',))


	@property
	def mass(self):
		return (self.grid.area*self.grid.column_density).sum()
