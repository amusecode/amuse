import os
import os.path
import pickle
import random
import numpy
import hashlib
from amuse.community import *
from amuse.community.interface.common import CommonCode, CommonCodeInterface

from amuse.support.options import option
from subprocess import Popen, PIPE

from amuse.rfi.core import PythonCodeInterface


class GalactICsImplementation(object):
    
    def __init__(self):
        self._output_directory = "./"
        self._particles_generated = False
        self._particle_data = numpy.array([])
        self._bin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "bin")
        
        
    def get_output_path(self, output_directory):
        output_directory.value = self._output_directory
        return 0
    
    def set_output_path(self, output_directory):
        self._output_directory = output_directory
        return 0
    
    def set_src_bin_path(self, src_bin_path):
        self._bin_path = src_bin_path
        return 0
    
    def initialize_code(self):
        self.set_default_parameter_values()
        return 0
    
    #parameter getters & setters
    for par in ["_generate_halo_flag", "_halo_outer_radius", "_scale_velocity", 
            "_scale_radius", "_truncation_delta_r", "_inner_cusp_slope", "_outer_slope", 
            "_generate_disk_flag", "_disk_mass", "_disk_scale_length", "_disk_outer_radius", 
            "_disk_scale_height_sech2", "_disk_truncation_dr", "_generate_bulge_flag", 
            "_Sersic_index_n", "_bulge_velocity", "_bulge_scale_radius", "_radial_grid_delta_r", 
            "_number_of_grid_intervals", "_order_of_multipole_expansion", 
            "_central_radial_vel_dispersion", "_scale_length_of_sigR2", 
            "_number_of_radial_steps_correction_fns", "_number_of_iterations", 
            "_halo_streaming_fraction", "_halo_number_of_particles", "_halo_random_seed", 
            "_halo_do_center_flag", "_bulge_streaming_fraction", "_bulge_number_of_particles", 
            "_bulge_random_seed", "_bulge_do_center_flag", "_disk_number_of_particles", 
            "_disk_random_seed", "_disk_do_center_flag"]:
        exec("def get"+par+"(self, value): value.value = self."+par+"; return 0")
        exec("def set"+par+"(self, value): self."+par+" = value; return 0")
    
    def set_default_parameter_values(self):
        # Halo parameters:
        # do you want a halo (y or n)
        self._generate_halo_flag = True
        self._halo_outer_radius  = 300.0
        self._scale_velocity     = 3.26331115
        self._scale_radius = 6.06699419
        self._truncation_delta_r = 100.0
        self._inner_cusp_slope = 1.0
        self._outer_slope = 2.3
        
        # Disk parameters:
        # do you want a disk (y or n)
        self._generate_disk_flag = True
        self._disk_mass = 25.0
        self._disk_scale_length = 5.8097949
        self._disk_outer_radius = 40.5
        self._disk_scale_height_sech2 = 0.5
        self._disk_truncation_dr = 1.5
        
        # Bulge parameters:
        # do you want a Sersic bulge (y or n)
        self._generate_bulge_flag = True
        self._Sersic_index_n = 0.937324703
        self._ppp_ = -1.0
        self._bulge_velocity = 3.21182013
        self._bulge_scale_radius = 1.50395405
        
        # do you want a blackhole (y or n) (NOT SUPPORTED IN THIS VERSION)
        #n
        self._radial_grid_delta_r = 0.01
        self._number_of_grid_intervals = 90000
        self._order_of_multipole_expansion = 10
        # order of multipole expansion - even number - should be l=10 for models
        # with disks - l=0 for purely spherical models without a disk
        
        # parameters for in.diskdf
        self._central_radial_vel_dispersion = 0.73 # central radial vel dispersion (in vzdisp)
        self._scale_length_of_sigR2 = 5.8097949 # scale length of sig_r^2
        self._number_of_radial_steps_correction_fns = 10 # number of intervals for correction functions (min. 6)
        self._number_of_iterations = 50
        #self._psfile = "psfile" # doesn't seem to be used...
        
        # parameters for in.halo
        self._halo_streaming_fraction = 0.50 # [0.0, 1.0]; 0.5 means no net rotation
        self._halo_number_of_particles = 200000
        self._halo_random_seed = -1
        self._halo_do_center_flag = True
        # parameters for in.bulge
        self._bulge_streaming_fraction = 0.80 # [0.0, 1.0]; 0.5 means no net rotation
        self._bulge_number_of_particles = 50000
        self._bulge_random_seed = -1
        self._bulge_do_center_flag = True
        # parameters for in.disk
        self._disk_number_of_particles = 100000
        self._disk_random_seed = -1
        self._disk_do_center_flag = True
    
    def cleanup_code(self):
        return 0
    
    def generate_in_dbh_string(self):
        if self._generate_halo_flag:
            in_dbh = "y\n"
            in_dbh += "{0:.15} ".format(self._halo_outer_radius)
            in_dbh += "{0:.15} ".format(self._scale_velocity)
            in_dbh += "{0:.15} ".format(self._scale_radius)
            in_dbh += "{0:.15} ".format(self._truncation_delta_r)
            in_dbh += "{0:.15} ".format(self._inner_cusp_slope)
            in_dbh += "{0:.15}\n".format(self._outer_slope)
        else:
            in_dbh = "n\n"
        
        if self._generate_disk_flag:
            in_dbh += "y\n"
            in_dbh += "{0:.15} ".format(self._disk_mass)
            in_dbh += "{0:.15} ".format(self._disk_scale_length)
            in_dbh += "{0:.15} ".format(self._disk_outer_radius)
            in_dbh += "{0:.15} ".format(self._disk_scale_height_sech2)
            in_dbh += "{0:.15}\n".format(self._disk_truncation_dr)
        else:
            in_dbh += "n\n"
        
        if self._generate_bulge_flag:
            in_dbh += "y\n"
            in_dbh += "{0:.15} ".format(self._Sersic_index_n)
            in_dbh += "{0:.15} ".format(self._ppp_)
            in_dbh += "{0:.15} ".format(self._bulge_velocity)
            in_dbh += "{0:.15}\n".format(self._bulge_scale_radius)
        else:
            in_dbh += "n\n"
        
        in_dbh += "n\n"
        in_dbh += "{0:.15} ".format(self._radial_grid_delta_r)
        in_dbh += "{0}\n".format(self._number_of_grid_intervals)
        in_dbh += "{0}\n".format(self._order_of_multipole_expansion)
        return in_dbh
    
    def generate_in_diskdf_string(self):
        in_diskdf = "{0:.15} ".format(self._central_radial_vel_dispersion) # central radial vel dispersion (in vzdisp)
        in_diskdf += "{0:.15}\n".format(self._scale_length_of_sigR2) # scale length of sig_r^2
        in_diskdf += "{0}\n".format(self._number_of_radial_steps_correction_fns) # number of intervals for correction functions (min. 6)
        in_diskdf += "{0}\n".format(self._number_of_iterations)
        #in_diskdf += "{0}\n".format(self._psfile) # doesn't seem to be used...
        return in_diskdf
    
    def generate_in_halo_string(self):
        in_halo = "{0:.15}\n".format(self._halo_streaming_fraction)
        in_halo += "{0}\n".format(self._halo_number_of_particles)
        in_halo += "{0}\n".format(self._halo_random_seed)
        in_halo += "{0}\n".format(1 if self._halo_do_center_flag else 0)
        return in_halo
    
    def generate_in_bulge_string(self):
        in_bulge = "{0:.15}\n".format(self._bulge_streaming_fraction)
        in_bulge += "{0}\n".format(self._bulge_number_of_particles)
        in_bulge += "{0}\n".format(self._bulge_random_seed)
        in_bulge += "{0}\n".format(1 if self._bulge_do_center_flag else 0)
        return in_bulge
    
    def generate_in_disk_string(self):
        in_disk = "{0}\n".format(self._disk_number_of_particles)
        in_disk += "{0}\n".format(self._disk_random_seed)
        in_disk += "{0}\n".format(1 if self._disk_do_center_flag else 0)
        return in_disk
    
    def _new_dbh_dir(self, data_directory,in_dbh,in_diskdf):
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)
        with open(os.path.join(data_directory, "in.gendenspsi"), "w") as f:
            f.write("2000 40\n")
        # for clarity, also store the used input parameters in this directory:
        with open(os.path.join(data_directory, "in.dbh"), "w") as f:
            f.write(in_dbh)
        with open(os.path.join(data_directory, "in.diskdf"), "w") as f:
            f.write(in_diskdf)
        # remove finished-step files
        for f in ['dbh.finished','getfreqs.finished','diskdf.finished']:
          try:
            os.remove(os.path.join( data_directory, f))    
          except:
            pass

    def _data_directory(self,in_dbh,in_diskdf):
        modelhash=hashlib.sha1((in_dbh+in_diskdf).encode()).hexdigest()
        return os.path.join(self._output_directory, "model_"+modelhash)
        
    def model_present(self,x):
        in_dbh = self.generate_in_dbh_string()
        in_diskdf = self.generate_in_diskdf_string()
        data_directory=self._data_directory(in_dbh,in_diskdf)
        x.value=self._directory_contains_valid_model(data_directory)      
        return 0

    def _directory_contains_valid_model(self,data_directory):
        if os.path.exists(os.path.join( data_directory)) and \
           os.path.exists(os.path.join( data_directory, 'dbh.dat')) and \
           os.path.exists(os.path.join( data_directory, 'dbh.finished')) and \
           os.path.exists(os.path.join( data_directory, 'getfreqs.finished')) and \
           (os.path.exists(os.path.join( data_directory, 'diskdf.finished')) or not self._generate_disk_flag):
             return True
        else:
             return False
    
    def _location_dbh_dat(self, in_dbh,in_diskdf):
        data_directory=self._data_directory(in_dbh,in_diskdf)
        
        if self._directory_contains_valid_model(data_directory):
            is_new = False
        else:
            is_new = True
            self._new_dbh_dir(data_directory,in_dbh,in_diskdf)
        return data_directory, is_new
    
    def commit_parameters(self):
        try:
            in_dbh = self.generate_in_dbh_string()
            in_diskdf = self.generate_in_diskdf_string()
            dbh_dir, is_new = self._location_dbh_dat(in_dbh, in_diskdf)
            self._cwd = dbh_dir
            print(dbh_dir)
            if not is_new:
                return 0
            print("Writing output to:", self._cwd)
            
            proc=Popen([os.path.join(self._bin_path, "dbh")], 
                cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
            stdout,stderr=proc.communicate(in_dbh.encode('UTF-8'))
            print("(stdout, stderr) =", stdout,stderr)
            if proc.returncode==0:
              open(os.path.join(dbh_dir,"dbh.finished"),'a').close()
            
            proc=Popen([os.path.join(self._bin_path, "getfreqs")], 
                cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
            stdout,stderr=proc.communicate()
            print("(stdout, stderr) =", stdout,stderr)
            if proc.returncode==0:
              open(os.path.join(dbh_dir,"getfreqs.finished"),'a').close()

            if self._generate_disk_flag:
              proc=Popen([os.path.join(self._bin_path, "diskdf")], 
                      cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
              stdout,stderr=proc.communicate(in_diskdf.encode('UTF-8'))
              print("(stdout, stderr) =", stdout,stderr)
              if proc.returncode==0:
                open(os.path.join(dbh_dir,"diskdf.finished"),'a').close()
                                
            return 0
        except Exception as ex:
            print("Exception occurred in commit_parameters:", ex)
            raise
    
    def recommit_parameters(self):
        return self.commit_parameters()
    
    def get_number_of_particles_updated(self, number_of_particles_updated):
        if self._particles_generated:
            number_of_particles_updated.value = self._number_of_particles_updated
            self._particles_generated = False
        else:
            number_of_particles_updated.value = 0
        return 0
    
    def generate_particles(self):
        try:
            if self._generate_disk_flag:
                in_disk  = self.generate_in_disk_string()
                process = Popen([os.path.join(self._bin_path, "gendisk")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out,err=process.communicate(in_disk.encode('UTF-8'))
                if process.returncode != 0:
                    print("error:", err)
                    return -2
                disk_data=numpy.frombuffer(out,dtype="float32")
            else:
                disk_data=numpy.array([])    

            if self._generate_bulge_flag:
                in_bulge = self.generate_in_bulge_string()
                process = Popen([os.path.join(self._bin_path, "genbulge")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out,err=process.communicate(in_bulge.encode('UTF-8'))
                if process.returncode != 0:
                    print("error:", err)
                    return -3 
                bulge_data=numpy.frombuffer(out,dtype="float32")
            else:
                bulge_data=numpy.array([])    
                            
            if self._generate_halo_flag:
                in_halo  = self.generate_in_halo_string()
                process = Popen([os.path.join(self._bin_path, "genhalo")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out, err = process.communicate(in_halo.encode('UTF-8'))
                if process.returncode != 0:
                    print("error:", err)
                    return -4 
                halo_data=numpy.frombuffer(out,dtype="float32")
            else:
                halo_data=numpy.array([])    
            
            self._number_of_particles_updated = (len(halo_data)+len(bulge_data)+len(disk_data))//7
            self._number_of_halo_particles=len(halo_data)//7
            self._number_of_bulge_particles=len(bulge_data)//7
            self._number_of_disk_particles=len(disk_data)//7
            data=numpy.concatenate((disk_data,bulge_data,halo_data))              
            self._particle_data = numpy.reshape(data,( self._number_of_particles_updated,7))
            self._particles_generated = True
            return 0
        except Exception as ex:
            print("Exception occurred in generate_particles:", ex)
            return -1
    
    def get_number_of_particles(self,nhalo,nbulge,ndisk):
        try:
            nhalo.value=self._number_of_halo_particles
            nbulge.value=self._number_of_bulge_particles
            ndisk.value=self._number_of_disk_particles
            return 0
        except:
            return -1
    
    def get_mass(self, index_of_the_particle, mass, length):
        try:
            mass.value = self._particle_data[index_of_the_particle, 0]
            return 0
        except:        
            return -1
    
    def get_position(self, index_of_the_particle, x, y, z, length):
        try:
            x.value = self._particle_data[index_of_the_particle, 1]
            y.value = self._particle_data[index_of_the_particle, 2]
            z.value = self._particle_data[index_of_the_particle, 3]
            return 0
        except:        
            return -1
    
    def get_velocity(self, index_of_the_particle, vx, vy, vz, length):
        try:
            vx.value = self._particle_data[index_of_the_particle, 4]
            vy.value = self._particle_data[index_of_the_particle, 5]
            vz.value = self._particle_data[index_of_the_particle, 6]
            return 0
        except:        
            return -1
    


class GalactICsInterface(PythonCodeInterface, CommonCodeInterface, LiteratureReferencesMixIn,
        CodeWithDataDirectories):
    """
    GalactICs allows to generate self-consistent disc-bulge-halo galaxy models. 
    The bulge and halo distribution functions (DFs) are functions of E and L_z 
    only. The halo's flattening and rotation can be specified. The disc DF is a 
    function of E and L_z and a third 'integral', E_z, the vertical energy, 
    which is approximately conserved in a warm disc with vertical extent. A 
    simulation of a sample model shows that in practice the models are very 
    close to equilibrium, making them ideal for experiments on instabilities in 
    galactic discs.
    
    Relevant references:
        .. [#] Kuijken K., Dubinski J., 1995, MNRAS, 277, 1341 (original version)
        .. [#] Widrow L.M., Dubinski J., 2005, ApJ, 631, 838 (2nd version)
        .. [#] Widrow L.M., Pym B., Dubinski J., 2008, ApJ, 679, 1239 (current version)
    """
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, GalactICsImplementation, **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    def _check_if_worker_is_up_to_date(self):
        if not os.path.exists(os.path.join(GalactICsImplementation()._bin_path, "dbh")):
            raise exceptions.CodeException(
                "The worker code of the '{0}' interface class is not up to date.\n"
                "Please do a 'make clean; make' in the root directory.".format(type(self).__name__))
    
    new_particle = None
    
    def delete_particle(self, index_of_the_particle):
        return 0
    
    @legacy_function
    def get_output_path():
        function = LegacyFunctionSpecification()
        function.addParameter('output_directory', dtype='string', direction=function.OUT,
            description = "The path to the output directory.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_output_path():
        function = LegacyFunctionSpecification()
        function.addParameter('output_directory', dtype='string', direction=function.IN,
            description = "The path to the output directory.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_src_bin_path():
        function = LegacyFunctionSpecification()
        function.addParameter('src_bin_path', dtype='string', direction=function.IN,
            description = "The path to the Galactics binaries.")
        function.result_type = 'int32'
        return function
    
    #parameter getters & setters
    # boolean parameters
    for par in ["_generate_halo_flag", "_generate_disk_flag", "_generate_bulge_flag", 
            "_halo_do_center_flag", "_bulge_do_center_flag", "_disk_do_center_flag"]:
        exec("@legacy_function\ndef get"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='bool', direction=function.OUT)\n"
            "  function.result_type = 'int32'\n  return function")
        exec("@legacy_function\ndef set"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='bool', direction=function.IN)\n"
            "  function.result_type = 'int32'\n  return function")
    
    # integer parameters
    for par in [ "_number_of_grid_intervals", "_order_of_multipole_expansion", 
            "_number_of_radial_steps_correction_fns", "_number_of_iterations", 
            "_halo_number_of_particles", "_halo_random_seed", 
            "_bulge_number_of_particles", "_bulge_random_seed", 
            "_disk_number_of_particles", "_disk_random_seed"]:
        exec("@legacy_function\ndef get"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='int32', direction=function.OUT)\n"
            "  function.result_type = 'int32'\n  return function")
        exec("@legacy_function\ndef set"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='int32', direction=function.IN)\n"
            "  function.result_type = 'int32'\n  return function")
    
    # float parameters
    for par in ["_halo_outer_radius", "_scale_velocity", "_scale_radius", 
            "_truncation_delta_r", "_inner_cusp_slope", "_outer_slope", 
            "_disk_mass", "_disk_scale_length", "_disk_outer_radius", 
            "_disk_scale_height_sech2", "_disk_truncation_dr", "_Sersic_index_n", 
            "_bulge_velocity", "_bulge_scale_radius", "_radial_grid_delta_r", 
            "_central_radial_vel_dispersion", "_scale_length_of_sigR2", 
            "_halo_streaming_fraction", "_bulge_streaming_fraction"]:
        exec("@legacy_function\ndef get"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='float64', direction=function.OUT)\n"
            "  function.result_type = 'int32'\n  return function")
        exec("@legacy_function\ndef set"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='float64', direction=function.IN)\n"
            "  function.result_type = 'int32'\n  return function")
    
    def invoke_state_change2(self):
        pass
    
    def invoke_state_change_updated(self):
        pass
    
    @legacy_function
    def generate_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_number_of_particles_updated():
        """
        Return the number of particles added during the last generate_particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles_updated', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_particles():
        """
        Return the number of halo/bulge/disk particles of the last generate_particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_halo_particles', dtype='int32', direction=function.OUT)
        function.addParameter('number_of_bulge_particles', dtype='int32', direction=function.OUT)
        function.addParameter('number_of_disk_particles', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current x component of the position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current y component of the position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current z component of the position vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x component of the velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y component of the velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z component of the velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def model_present():
        """
        Return whether a valid galaxy model is present.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('model_present', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        return function


    def get_code_src_directory(self):
        return os.path.join(os.path.dirname(os.path.normpath(os.path.abspath(__file__))),'src')


class GalactICs(CommonCode):
    
    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        InCodeComponentImplementation.__init__(self, GalactICsInterface(**options), **options)
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        self.parameters.set_defaults()
        self.parameters.output_directory = self.get_output_directory()
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_output_path", 
            "set_output_path",
            "output_directory", 
            "The path to the output directory", 
            default_value = "./"
        )
        
        # boolean parameters
        for par in ["generate_halo_flag", "generate_disk_flag", "generate_bulge_flag"]:
            handler.add_boolean_parameter(
                "get_"+par,
                "set_"+par,
                par,
                "Flag specifying whether to generate a "+par[9:-5],
                True
            )
        for par in ["halo_do_center_flag", "bulge_do_center_flag", "disk_do_center_flag"]:
            handler.add_boolean_parameter(
                "get_"+par,
                "set_"+par,
                par,
                "Flag specifying whether to center the "+par[:-15]+" at the origin",
                True
            )
        
        # integer parameters
        handler.add_method_parameter(
            "get_number_of_grid_intervals",
            "set_number_of_grid_intervals",
            "number_of_grid_intervals",
            "Number of gridpoints in the radial direction",
            default_value = 90000
        )
        handler.add_method_parameter(
            "get_order_of_multipole_expansion",
            "set_order_of_multipole_expansion",
            "order_of_multipole_expansion",
            "order of multipole expansion - even number - should be l=10 for models with disks - l=0 for purely spherical models without a disk",
            default_value = 10
        )
        handler.add_method_parameter(
            "get_number_of_radial_steps_correction_fns",
            "set_number_of_radial_steps_correction_fns",
            "number_of_radial_steps_correction_fns_disk_df",
            "The number of intervals for correction functions (min. 6); used in calculation of the DF of the disk",
            default_value = 10
        )
        handler.add_method_parameter(
            "get_number_of_iterations",
            "set_number_of_iterations",
            "number_of_iterations_disk_df",
            "The number of iterations in calculation of the DF of the disk",
            default_value = 50
        )
        handler.add_method_parameter(
            "get_halo_number_of_particles",
            "set_halo_number_of_particles",
            "halo_number_of_particles",
            "The number of halo particles to generate",
            default_value = 200000
        )
        handler.add_method_parameter(
            "get_bulge_number_of_particles",
            "set_bulge_number_of_particles",
            "bulge_number_of_particles",
            "The number of bulge particles to generate",
            default_value = 50000
        )
        handler.add_method_parameter(
            "get_disk_number_of_particles",
            "set_disk_number_of_particles",
            "disk_number_of_particles",
            "The number of disk particles to generate",
            default_value = 100000
        )
        handler.add_method_parameter(
            "get_halo_random_seed",
            "set_halo_random_seed",
            "halo_random_seed",
            "The seed to the random number generator used to generate the halo particles",
            default_value = -1
        )
        handler.add_method_parameter(
            "get_bulge_random_seed",
            "set_bulge_random_seed",
            "bulge_random_seed",
            "The seed to the random number generator used to generate the bulge particles",
            default_value = -1
        )
        handler.add_method_parameter(
            "get_disk_random_seed",
            "set_disk_random_seed",
            "disk_random_seed",
            "The seed to the random number generator used to generate the disk particles",
            default_value = -1
        )
        
        # float parameters
        handler.add_method_parameter(
            "get_halo_outer_radius",
            "set_halo_outer_radius",
            "halo_outer_radius",
            "The halo is smoothly truncated at this radius",
            default_value = 300.0 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_scale_velocity",
            "set_scale_velocity",
            "halo_scale_velocity",
            "The velocity scale of the halo",
            default_value = 3.26331115 | nbody_system.speed
        )
        handler.add_method_parameter(
            "get_scale_radius",
            "set_scale_radius",
            "halo_scale_radius",
            "The length scale of the halo",
            default_value = 6.06699419 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_truncation_delta_r",
            "set_truncation_delta_r",
            "halo_truncation_width",
            "The width of the smooth truncation at halo_outer_radius",
            default_value = 100.0 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_inner_cusp_slope",
            "set_inner_cusp_slope",
            "halo_inner_cusp_slope",
            "The slope of inner cusp of the halo density profile",
            default_value = 1.0
        )
        handler.add_method_parameter(
            "get_outer_slope",
            "set_outer_slope",
            "halo_outer_slope",
            "The outer slope of the halo density profile",
            default_value = 2.3
        )
        handler.add_method_parameter(
            "get_disk_mass",
            "set_disk_mass",
            "disk_mass",
            "The mass of the disk",
            default_value = 25.0 | nbody_system.mass
        )
        handler.add_method_parameter(
            "get_disk_scale_length",
            "set_disk_scale_length",
            "disk_scale_length",
            "The length scale of the disk",
            default_value = 5.8097949 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_disk_outer_radius",
            "set_disk_outer_radius",
            "disk_outer_radius",
            "The disk is smoothly truncated at this radius",
            default_value = 40.5 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_disk_scale_height_sech2",
            "set_disk_scale_height_sech2",
            "disk_scale_height_sech2",
            "The vertical scale length of the disk. The disk falls off as sech^2 in the z-direction.",
            default_value = 0.5 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_disk_truncation_dr",
            "set_disk_truncation_dr",
            "disk_truncation_width",
            "The width of the smooth truncation at disk_outer_radius",
            default_value = 1.5 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_Sersic_index_n",
            "set_Sersic_index_n",
            "Sersic_index",
            "The Sersic index of the bulge (1.0 for a classical bulge)",
            default_value = 0.937324703
        )
        handler.add_method_parameter(
            "get_bulge_velocity",
            "set_bulge_velocity",
            "bulge_scale_velocity",
            "The velocity scale of the bulge",
            default_value = 3.21182013 | nbody_system.speed
        )
        handler.add_method_parameter(
            "get_bulge_scale_radius",
            "set_bulge_scale_radius",
            "bulge_scale_radius",
            "The length scale of the bulge",
            default_value = 1.50395405 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_radial_grid_delta_r",
            "set_radial_grid_delta_r",
            "radial_grid_delta_r",
            "Spacing of the grid in the radial direction",
            default_value = 0.01 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_central_radial_vel_dispersion",
            "set_central_radial_vel_dispersion",
            "disk_central_radial_velocity_dispersion",
            "The velocity dispersion of the disk in the radial direction at the center (in units of vertical velocity dispersion)",
            default_value = 0.73
        )
        handler.add_method_parameter(
            "get_scale_length_of_sigR2",
            "set_scale_length_of_sigR2",
            "disk_scale_length_of_sigR2",
            "The length scale of the exponential decline of the velocity dispersion of the disk in the radial direction.",
            default_value = 5.8097949 | nbody_system.length
        )
        handler.add_method_parameter(
            "get_halo_streaming_fraction",
            "set_halo_streaming_fraction",
            "halo_streaming_fraction",
            "Control for rotating halo: distribution function is split in positive and negative angular momentum, and recombined with this parameter (F = aF+ + (1-a)F-); 0.5 means no rotation",
            default_value = 0.50
        )
        handler.add_method_parameter(
            "get_bulge_streaming_fraction",
            "set_bulge_streaming_fraction",
            "bulge_streaming_fraction",
            "Control for rotating bulge: distribution function is split in positive and negative angular momentum, and recombined with this parameter (F = aF+ + (1-a)F-); 0.5 means no rotation",
            default_value = 0.80
        )
    
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method("generate_particles", (), (handler.ERROR_CODE,))
        handler.add_method("get_number_of_particles_updated", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler.add_method("get_mass", (handler.INDEX,), 
            (nbody_system.mass, handler.ERROR_CODE)
        )
        handler.add_method("get_position", (handler.INDEX,), 
            (nbody_system.length, nbody_system.length, nbody_system.length, handler.ERROR_CODE)
        )
        handler.add_method("get_velocity", (handler.INDEX,), 
            (nbody_system.speed, nbody_system.speed, nbody_system.speed, handler.ERROR_CODE)
        )
        
        handler.add_method("get_output_path", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_output_path", (handler.NO_UNIT,), (handler.ERROR_CODE,))
        
        for par in ["_number_of_grid_intervals", "_order_of_multipole_expansion", 
                "_number_of_radial_steps_correction_fns", "_number_of_iterations", 
                "_halo_number_of_particles", "_halo_random_seed", 
                "_bulge_number_of_particles", "_bulge_random_seed", 
                "_disk_number_of_particles", "_disk_random_seed",
                "_inner_cusp_slope", "_outer_slope", "_Sersic_index_n", "_central_radial_vel_dispersion", 
                "_halo_streaming_fraction", "_bulge_streaming_fraction"]:
            handler.add_method("get"+par, (), (handler.NO_UNIT, handler.ERROR_CODE,))
            handler.add_method("set"+par, (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        for par in ["_halo_outer_radius", "_scale_radius", "_truncation_delta_r", 
                "_disk_scale_length", "_disk_outer_radius", "_disk_scale_height_sech2", 
                "_disk_truncation_dr", "_bulge_scale_radius", "_radial_grid_delta_r", 
                "_scale_length_of_sigR2"]:
            handler.add_method("get"+par, (), (nbody_system.length, handler.ERROR_CODE,))
            handler.add_method("set"+par, (nbody_system.length, ), (handler.ERROR_CODE,))
        
        for par in ["_scale_velocity", "_bulge_velocity"]:
            handler.add_method("get"+par, (), (nbody_system.speed, handler.ERROR_CODE,))
            handler.add_method("set"+par, (nbody_system.speed, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_disk_mass", (), (nbody_system.mass, handler.ERROR_CODE,))
        handler.add_method("set_disk_mass", (nbody_system.mass, ), (handler.ERROR_CODE,))
    
    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_getter('particles', 'get_mass')
        handler.add_getter('particles', 'get_position')
        handler.add_getter('particles', 'get_velocity')
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')

        handler.add_method('CHANGE_PARAMETERS_RUN', 'model_present')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'model_present')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','model_present')
        handler.add_method('INITIALIZED','model_present')

        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        handler.add_method('UPDATE','before_get_parameter')
        
        handler.add_transition('EDIT', 'UPDATE', 'generate_particles', False)
        handler.add_transition('UPDATE', 'RUN', 'update_particle_set')
        handler.add_transition('RUN', 'EDIT', 'clear_particle_set')
        handler.add_method('RUN', 'invoke_state_change_updated')
        handler.add_method('EDIT', 'get_number_of_particles_updated')
        handler.add_method('UPDATE', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_number_of_particles')
        handler.add_method('RUN', 'get_mass')
        handler.add_method('RUN', 'get_position')
        handler.add_method('RUN', 'get_velocity')
    
    def generate_particles(self):
        result = self.overridden().generate_particles()
        self.invoke_state_change_updated()
    
    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.datamodel.incode_storage module, as
        it uses a lot of internal methods and info!
        """
        number_of_updated_particles = self.get_number_of_particles_updated()
        if number_of_updated_particles:
            self.particles._private.attribute_storage._add_indices(
                list(range(number_of_updated_particles))
            )
    
    def clear_particle_set(self):
        if len(self.particles):
            self.particles.remove_particles(self.particles)
    
    @property
    def halo_particles(self):
        nhalo,nbulge,ndisk=self.get_number_of_particles()
        return self.particles[ndisk+nbulge:]

    @property
    def bulge_particles(self):
        nhalo,nbulge,ndisk=self.get_number_of_particles()
        return self.particles[ndisk:ndisk+nbulge]

    @property
    def disk_particles(self):
        nhalo,nbulge,ndisk=self.get_number_of_particles()
        return self.particles[:ndisk]


Galactics = GalactICs
