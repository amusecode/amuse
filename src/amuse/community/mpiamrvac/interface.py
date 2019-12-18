from amuse.community import *
from amuse.community.interface.hydro import HydrodynamicsInterface
from amuse.support.options import OptionalAttributes, option
from amuse.units import generic_unit_system
from amuse.community.interface.common import CommonCode

import os

class MpiAmrVacInterface(CodeInterface, HydrodynamicsInterface, StoppingConditionInterface,
        CodeWithDataDirectories):
    
    use_modules = ['mpiamrvac_interface', 'StoppingConditions']
    
    MODE_NORMAL = 'normal'
    MODE_3D = '3d'
    MODE_3D_ACC = '3d-acc'
    MODE_2D   = '2d'
    MODE_2D_ACC   = '2d-acc'
    MODE_1D   = '1d'
    MODE_1D_ACC   = '1d-acc'
    
    def __init__(self, mode = MODE_NORMAL, **options):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **options)
        CodeWithDataDirectories.__init__(self)
        self._mode = mode
        
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL or mode == self.MODE_3D:
            return 'mpiamrvac_worker'
        elif mode == self.MODE_3D_ACC:
            return 'mpiamrvac_worker_acc'
        elif mode == self.MODE_2D:
            return 'mpiamrvac_worker_2d'
        elif mode == self.MODE_2D_ACC:
            return 'mpiamrvac_worker_2dacc'
        elif mode == self.MODE_1D:
            return 'mpiamrvac_worker_1d'
        elif mode == self.MODE_1D_ACC:
            return 'mpiamrvac_worker_1dacc'
        else:
            return 'mpiamrvac_worker'
    
    #
    # options
    #
    
    @option(type="string")
    def default_parameters_filename(self):
        """
        Default parameter file for amrvac, has empty lists for
        all parameters.
        """
        if self._mode == self.MODE_2D:
            return os.path.join(self.data_directory, 'amrvac_2d.par')
        elif self._mode == self.MODE_2D_ACC:
            return os.path.join(self.data_directory, 'amrvac_2d-acc.par')
        elif self._mode == self.MODE_1D:
            return os.path.join(self.data_directory, 'amrvac_1d.par')
        elif self._mode == self.MODE_1D_ACC:
            return os.path.join(self.data_directory, 'amrvac_1d-acc.par')
        elif self._mode == self.MODE_3D_ACC:
            return os.path.join(self.data_directory, 'amrvac-acc.par')
        else:
            return os.path.join(self.data_directory, 'amrvac.par')
    
    #
    # parameters
    #
    
    @legacy_function
    def set_typeentropy():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeentropy():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_typefull1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typefull1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_typepred1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typepred1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_gamma():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gamma():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_dt():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dt():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nbufferx1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nbufferx1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nbufferx2():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nbufferx2():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nbufferx3():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nbufferx3():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mxnest():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mxnest():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dixb():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dixb():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_levmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_levmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_levmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_levmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_skipfinestep():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_skipfinestep():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time_advance():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_advance():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_courantpar():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_courantpar():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtpar():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dtpar():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtdiffpar():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dtdiffpar():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_t():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_t():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dtmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dtmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_residmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_residmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_residmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_residmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_residual():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_residual():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tfixgrid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tfixgrid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tvdlfeps():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tvdlfeps():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mcbeta():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mcbeta():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_divbdiff():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_divbdiff():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_smallp():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_smallp():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_smallrho():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_smallrho():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dmaxvel():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dmaxvel():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tolernr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tolernr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_absaccnr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_absaccnr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_cfrac():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_cfrac():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_x1ptms():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_x1ptms():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_x2ptms():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_x2ptms():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_x3ptms():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_x3ptms():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ptmass():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ptmass():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ratebdflux():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ratebdflux():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_normt():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_normt():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time_bc():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_bc():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_it():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_it():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_itmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_itmax():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_itmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_itmin():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_slowsteps():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_slowsteps():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typepario():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typepario():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_itfixgrid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_itfixgrid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nwauxio():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nwauxio():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_istep():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_istep():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nstep():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nstep():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_errorestimate():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_errorestimate():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nxdiffusehllc():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nxdiffusehllc():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typespherical():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typespherical():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_maxitnr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_maxitnr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nflatgetaux():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nflatgetaux():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_level_io():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_level_io():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ncool():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ncool():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_cmulti():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_cmulti():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_snapshotini():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_snapshotini():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ixtest1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ixtest1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ixtest2():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ixtest2():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ixtest3():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_ixtest3():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_iwtest():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_iwtest():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_idimtest():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_idimtest():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_saveigrid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_saveigrid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typecourant():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typecourant():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeresid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeresid():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeadvance():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeadvance():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typelimited():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typelimited():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typesourcesplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typesourcesplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typelimiter():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typelimiter():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typegradlimiter():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typegradlimiter():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeprolonglimit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeprolonglimit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typetvd():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typetvd():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typetvdlf():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typetvdlf():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeaverage():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeaverage():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typedimsplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typedimsplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeaxial():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeaxial():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typepoly():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typepoly():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typedivbdiff():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typedivbdiff():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typedivbfix():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typedivbfix():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typediv():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typediv():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typegrad():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typegrad():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeglm():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeglm():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_coolcurve():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_coolcurve():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_coolmethod():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_coolmethod():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typeghostfill():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typeghostfill():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typegridfill():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typegridfill():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_filenameout():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_filenameout():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_filenameini():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_filenameini():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_filenamelog():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_filenamelog():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_fileheadout():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_fileheadout():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_wnames():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_wnames():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_primnames():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_primnames():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_typefilelog():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_typefilelog():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_convert_type():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_convert_type():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dxfiletype():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dxfiletype():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_teststr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_teststr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time_accurate():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_accurate():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_addmpibarrier():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_addmpibarrier():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tmaxexact():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tmaxexact():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_treset():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_treset():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_itreset():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_itreset():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_firstprocess():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_firstprocess():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_fixprocess():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_fixprocess():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_flathllc():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_flathllc():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_flatcd():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_flatcd():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_flatsh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_flatsh():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_flatppm():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_flatppm():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_sourcesplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_sourcesplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_sourceunsplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_sourceunsplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_useprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_useprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dimsplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dimsplit():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_restrictprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_restrictprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_prolongprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_prolongprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_coarsenprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_coarsenprimitive():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_useprimitiverel():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_useprimitiverel():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amrentropy():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_amrentropy():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_divbfix():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_divbfix():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_divbwave():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_divbwave():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_compactres():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_compactres():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bnormlf():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_bnormlf():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_strictnr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_strictnr():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_strictsmall():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_strictsmall():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_strictzero():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_strictzero():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_strictgetaux():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_strictgetaux():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_usecovariant():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_usecovariant():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nocartesian():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nocartesian():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tfix():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tfix():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_convert():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_convert():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_saveprim():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_saveprim():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_uselimiter():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_uselimiter():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    # parameters file
    
    
    @legacy_function
    def set_parameters_filename():
        """
        Update name of the parameters file
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', 
            dtype='string',
            direction=function.IN,
            description = "filename of the parameters file"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
        """
        return function
        
    @legacy_function
    def get_parameters_filename():
        """
        Retrieve name of the parameters file
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', 
            dtype='string',
            direction=function.OUT,
            description = "filename of the parameters file"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
        """
        return function
        
    @legacy_function
    def get_current_error():
        """When a function returns an error, this will retrieve
        a description (if possible)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('string', 
            dtype='string',
            direction=function.OUT,
            description = "description of the error"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
        """
        return function
    #
    #
    
    @legacy_function    
    def initialize_grid():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function
        
    @legacy_function  
    def refine_grid():
        function = LegacyFunctionSpecification()  
        function.addParameter('must_advance', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_mesh_size():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('nmeshx', dtype='i', direction=function.OUT)
        function.addParameter('nmeshy', dtype='i', direction=function.OUT)
        function.addParameter('nmeshz', dtype='i', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_position_of_index():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        
        function.result_type = 'i'
        return function
        
        
    @legacy_function    
    def get_index_of_position():
        """
        Retrieves the i,j and k index of the grid cell containing the
        given x, y and z position. The cell is looked up
        in the grid specified by index_of_grid.
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        
        for x in ['i','j','k']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function
     
        
    @legacy_function
    def get_level_of_grid():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('level', dtype='i', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_cell_size_of_grid():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('dx', dtype='d', direction=function.OUT)
        function.addParameter('dy', dtype='d', direction=function.OUT)
        function.addParameter('dz', dtype='d', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

        
    @legacy_function   
    def setup_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('nmeshx', dtype='i', direction=function.IN)
        function.addParameter('nmeshy', dtype='i', direction=function.IN)
        function.addParameter('nmeshz', dtype='i', direction=function.IN)
        function.addParameter('xlength', dtype='d', direction=function.IN)
        function.addParameter('ylength', dtype='d', direction=function.IN)
        function.addParameter('zlength', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    
    

    #
    #
    #
    
    
    @legacy_function    
    def set_boundary_in_code():
        function = LegacyFunctionSpecification()  
        for x in ["xbound1","xbound2","ybound1","ybound2","zbound1","zbound2"]:
            function.addParameter(x, dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
        
    
    def set_boundary(self, xbound1, xbound2, ybound1, ybound2, zbound1, zbound2):
        map_from_amuse_to_mpiamrvac= {
            "reflective": "symm", 
            "outflow":"limitinflow", 
            "periodic":"periodic",
            "interface": "special",
        }
        
        return self.set_boundary_in_code(
            map_from_amuse_to_mpiamrvac.setdefault(xbound1, xbound1),
            map_from_amuse_to_mpiamrvac.setdefault(xbound2, xbound2),
            map_from_amuse_to_mpiamrvac.setdefault(ybound1, ybound1),
            map_from_amuse_to_mpiamrvac.setdefault(ybound2, ybound2),
            map_from_amuse_to_mpiamrvac.setdefault(zbound1, zbound1),
            map_from_amuse_to_mpiamrvac.setdefault(zbound2, zbound2)
        )
        
    #
    #
    #
    

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    #
    #
    #
    
    
    @legacy_function    
    def get_acceleration_grid_position_of_index():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        #function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_acceleration_grid_acceleration():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        #function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['a1','a2', 'a3']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def set_acceleration_grid_acceleration():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        #function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['a1','a2', 'a3']:
            function.addParameter(x, dtype='d', direction=function.IN)
        
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        
        function.result_type = 'i'
        return function
        
    
    @legacy_function
    def get_acceleration_grid_size():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('nmeshx', dtype='i', direction=function.OUT)
        function.addParameter('nmeshy', dtype='i', direction=function.OUT)
        function.addParameter('nmeshz', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
        
    

    @legacy_function
    def get_grid_acceleration():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['ax','ay','az']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        
        return function
    

    @legacy_function
    def set_grid_acceleration():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['ax','ay','az']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_hydro_state_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, default = 0)
        for x in ['rho','rhovx','rhovy','rhovz','rhoe']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'i' 
        function.must_handle_array = True
        return function
    
    
    
    
class MpiAmrVac(CommonCode):

    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        self.stopping_conditions = StoppingConditions(self)
        
        CommonCode.__init__(self,  MpiAmrVacInterface(**options), **options)
        
        self.set_parameters_filename(self.default_parameters_filename)
        
    def define_converter(self, handler):
        if self.unit_converter is None:
            return
        
        handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())

    
    
    def get_index_range_inclusive(self, index_of_grid = 1):
        nx, ny, nz = self.get_mesh_size(index_of_grid)
        
        return (0, nx-1, 0, ny-1, 0, nz-1)

    def define_properties(self, handler):
        handler.add_property('get_time', public_name = "model_time")
        
    def define_methods(self, handler):
        handler.add_method(
            'evolve_model',
            (generic_unit_system.time,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'commit_parameters',
            (),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_position_of_index',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_acceleration_grid_position_of_index',
            (handler.INDEX, handler.INDEX, handler.INDEX),
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, handler.ERROR_CODE,)
        )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.momentum_density
        energy =  generic_unit_system.energy_density
        acceleration =  generic_unit_system.length / generic_unit_system.time ** 2
        
        handler.add_method(
            'get_acceleration_grid_size',
            (),
            (handler.NO_UNIT,handler.NO_UNIT,handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_acceleration_grid_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX),
            (acceleration, acceleration, acceleration, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'set_acceleration_grid_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX, acceleration, acceleration, acceleration,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grid_energy_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, energy, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'set_grid_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, density, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'set_grid_momentum_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, momentum, momentum, momentum, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        
        
        handler.add_method(
            'get_grid_energy_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            ( energy,
            handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_grid_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (density,
            handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_grid_momentum_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            ( momentum, momentum, momentum, 
            handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'refine_grid',
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_level_of_grid',
            (handler.INDEX),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_gamma",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_gamma",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_typeentropy",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_typeentropy",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_typefull1",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_typefull1",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_typepred1",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_typepred1",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_typeadvance",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_typeadvance",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_courantpar",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_courantpar",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_mxnest",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_mxnest",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_time',
            (),
            (generic_unit_system.time, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'setup_mesh',
            (handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,  generic_unit_system.length,  generic_unit_system.length,  generic_unit_system.length, ),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_boundary',
            (handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'set_grid_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX, acceleration, acceleration, acceleration, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_grid_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (acceleration, acceleration, acceleration, handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_hydro_state_at_point',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,
                generic_unit_system.speed, generic_unit_system.speed, generic_unit_system.speed),
            (generic_unit_system.density, generic_unit_system.momentum_density, generic_unit_system.momentum_density, 
                generic_unit_system.momentum_density, generic_unit_system.energy_density, handler.ERROR_CODE)
        )
        self.stopping_conditions.define_methods(handler)
        
    def define_parameters(self, handler):
        
        
        handler.add_method_parameter(
            "get_gamma", 
            "set_gamma",
            "gamma", 
            "ratio of specific heats used in equation of state", 
            default_value = 1.6666666666666667
        )
        
        handler.add_method_parameter(
            "get_typeentropy", 
            "set_typeentropy",
            "entropy_type", 
            "type of the entropy", 
            default_value = 'nul'
        )
        handler.add_method_parameter(
            "get_typefull1", 
            "set_typefull1",
            "spatial_discretization_method", 
            "the spatial discretization method used for the time integration per activated grid leve", 
            default_value = 'tvdmu'
        )
        handler.add_method_parameter(
            "get_typepred1", 
            "set_typepred1",
            "predictor_step_discretization_method", 
            "the precitor step discretization method (only used when integration procedure is twostep')", 
            default_value = 'tvdmu'
        )
        handler.add_method_parameter(
            "get_typeadvance", 
            "set_typeadvance",
            "time_integration_procedure", 
            "time integration procedure", 
            default_value = 'twostep'
        )
        
        handler.add_method_parameter(
            "get_courantpar", 
            "set_courantpar",
            "courant_number", 
            "CFL number", 
            default_value = 0.7
        )
        
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "nmeshx",
            "nx", 
            "number of cells in the x direction", 
            10,
        )
        
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "nmeshy",
            "ny", 
            "number of cells in the y direction", 
            10,
        )
        
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "nmeshz",
            "nz", 
            "number of cells in the z direction", 
            10,
        )
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "xlength",
            "length_x", 
            "length of model in the x direction", 
            10 | generic_unit_system.length,
        )
        handler.add_caching_parameter(
            "setup_mesh", 
            "ylength",
            "length_y", 
            "length of model in the x direction", 
            10 | generic_unit_system.length,
        )
        handler.add_caching_parameter(
            "setup_mesh", 
            "zlength",
            "length_z", 
            "length of model in the z direction", 
            10 | generic_unit_system.length,
        )
        
        handler.add_vector_parameter(
            "mesh_size",
            "number of cells in the x, y and z directions",
            ("nx", "ny", "nz")
        )
        
        handler.add_vector_parameter(
            "mesh_length",
            "length of the model in the x, y and z directions",
            ("length_x", "length_y", "length_z")
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "xbound1",
            "xbound1", 
            "boundary conditions on first (inner, left) X boundary",
            "reflective",
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "xbound2",
            "xbound2", 
            "boundary conditions on second (outer, right) X boundary",
            "reflective",
        )
        
        handler.add_caching_parameter(
            "set_boundary", 
            "ybound1",
            "ybound1", 
            "boundary conditions on first (inner, front) Y boundary",
            "reflective",
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "ybound2",
            "ybound2", 
            "boundary conditions on second (outer, back) Y boundary", 
            "reflective",
        )
        
        handler.add_caching_parameter(
            "set_boundary", 
            "zbound1",
            "zbound1", 
            "boundary conditions on first (inner, bottom) Z boundary",
            "reflective",
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "zbound2",
            "zbound2", 
            "boundary conditions on second (outer, top) Z boundary",
            "reflective",
        )
        
        
        
        handler.add_vector_parameter(
            "x_boundary_conditions",
            "boundary conditions for the X directorion",
            ("xbound1", "xbound2")
        )
        
        
        handler.add_vector_parameter(
            "y_boundary_conditions",
            "boundary conditions for the Y directorion",
            ("ybound1", "ybound2")
        )
        
        
        handler.add_vector_parameter(
            "z_boundary_conditions",
            "boundary conditions for the Z directorion",
            ("zbound1", "zbound2")
        )
        
        
        handler.add_method_parameter(
            "get_mxnest", 
            "set_mxnest",
            "maximum_number_of_grid_levels", 
            "the maximum number of grid levels that can be used during the simulation, including the base grid level", 
            default_value = 3
        )
        handler.add_method_parameter(
            "get_time_accurate", 
            "set_time_accurate",
            "time_accurate", 
            "if false will evolve to the given time, if true will take accurate steps using courant timesteps", 
            default_value = 3
        )
        handler.add_method_parameter(
            "get_dtpar", 
            "set_dtpar",
            "timestep", 
            "if greater than zero will fix the timestep to the given value", 
            default_value = 3
        )
        
        self.stopping_conditions.define_parameters(handler)
        
        

    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
    
    def get_acceleration_grid_index_range_inclusive(self):
        nx, ny, nz = self.get_acceleration_grid_size()
        return (1, nx, 1, ny, 1, nz)
        
    def define_particle_sets(self, handler):
        
        handler.define_grid('acceleration_grid')
        handler.set_grid_range('acceleration_grid', 'get_acceleration_grid_index_range_inclusive')
        handler.add_getter('acceleration_grid', 'get_acceleration_grid_position_of_index', names=('x','y','z'))
        handler.add_getter('acceleration_grid', 'get_acceleration_grid_acceleration', names=('ax','ay','az'))
        handler.add_setter('acceleration_grid', 'set_acceleration_grid_acceleration', names=('ax','ay','az'))
        
    def itergrids(self):
        n = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield self._create_new_grid(self.specify_grid, index_of_grid = x)
    
    def specify_grid(self, definition, index_of_grid = 1):
        definition.set_grid_range('get_index_range_inclusive')
        
        definition.add_getter('get_position_of_index', names=('x','y','z'))
        
        definition.add_setter('set_grid_density', names=('rho',))
        
        definition.add_setter('set_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_setter('set_grid_energy_density', names=('energy',))
        
        definition.add_getter('get_grid_density', names=('rho',))
        definition.add_getter('get_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_getter('get_grid_energy_density', names=('energy',))
        
        
        definition.add_getter('get_grid_acceleration', names=('ax','ay','az'))
        definition.add_setter('set_grid_acceleration', names=('ax','ay','az'))
        
        definition.define_extra_keywords({'index_of_grid':index_of_grid})
        
    
    
    def define_state(self, handler): 
        CommonCode.define_state(self, handler)       
        #handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        
        handler.add_transition('EDIT', 'RUN', 'initialize_grid')
        handler.add_method('RUN', 'evolve_model')
        handler.add_method('RUN', 'get_hydro_state_at_point')
        for state in ['EDIT', 'RUN']:
            for methodname in [
                    'get_grid_density',
                    'set_grid_density',
                    'set_grid_energy_density',
                    'get_grid_energy_density',
                    'get_grid_momentum_density',
                    'set_grid_momentum_density', 
                    'get_position_of_index',
                    'get_index_of_position',
                    'set_grid_scalar',
                    'get_grid_scalar',
                    'get_mesh_size',
                    'get_acceleration_grid_acceleration',
                    'set_acceleration_grid_acceleration',
                    'get_acceleration_grid_size',
                    'get_mesh_size',
                    'get_number_of_grids',
                    'get_level_of_grid',
                    'refine_grid'
                ]:
                handler.add_method(state, methodname)     
                
        self.stopping_conditions.define_state(handler)    
    

