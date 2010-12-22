from amuse.legacy import *
from amuse.legacy.interface.common import CommonCodeInterface

from amuse.support.options import OptionalAttributes, option

import os

class MpiAmrVacInterface(LegacyInterface, CommonCodeInterface):
    
    use_modules = ['mpiamrvac_interface']
    
    def __init__(self, **options):
        LegacyInterface.__init__(self, name_of_the_worker="mpiamrvac_worker", **options)
    
    #
    # options
    #
    
    @option(type="string")
    def default_parameters_filename(self):
        """
        Default parameter file for amrvac, has empty lists for
        all parameters.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'mpiamrvac', 'input', 'amrvac.par')
    
    #
    # parameters
    #
    
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
        
    #
    #
    
    @legacy_function    
    def initialize_grid():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
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
    def get_number_of_grids():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='i', direction=function.OUT)
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
        
    
    
    @legacy_function
    def get_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['rho',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    

    @legacy_function
    def get_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['en',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    

    @legacy_function
    def get_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['rhovx', 'rhovy', 'rhovz',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function


    @legacy_function
    def set_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho',]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['en',]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    

    @legacy_function
    def set_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx', 'rhovy', 'rhovz',]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    #
    #
    #
    
    
    @legacy_function    
    def set_boundary():
        function = LegacyFunctionSpecification()  
        for x in ["xbound1","xbound2","ybound1","ybound2","zbound1","zbound2"]:
            function.addParameter(x, dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
        
    #
    #
    #
    
    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    
    
class MpiAmrVac(CodeInterface):

    def __init__(self, **options):
        CodeInterface.__init__(self,  MpiAmrVacInterface(**options), **options)
    
