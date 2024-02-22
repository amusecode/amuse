import os
import sys
import os.path
import pickle
import random
import numpy
import hashlib
from subprocess import Popen, PIPE

from amuse.units.core import *
from amuse.community import *
from amuse.community.interface.common import CommonCode, CommonCodeInterface
from amuse.support.options import option
from amuse.rfi.core import PythonCodeInterface

# KD95 = 1995MNRAS.277.1341K

parameters={
            "halo_type_parameter": dict(dtype="int32", default=3 ,description="type of halo 0=no,1=KD95,2=fixed (file),3=spherical df"), 
            "halo_scale_radius": dict(dtype="float64", default=13.7 | nbody_system.length ,description="halo scale parameter"), 
            "halo_virial_radius": dict(dtype="float64", default=40. | nbody_system.length ,description="halo 'virial' radius, beyond this the density is tapered smoothly"), 
            "halo_density_parameter": dict(dtype="float64", default=0.01 | nbody_system.density ,description="halo density parameter"),
            "halo_einasto_nindex": dict(dtype="float64", default=5. ,description="einasto profile n-index"),
            "halo_einasto_mass": dict(dtype="float64", default=400. ,description="einasto profile total mass"),
            "halo_central_potential": dict(dtype="float64", default=-28. | nbody_system.potential ,description="KD95 central potential"), 
            "halo_v0": dict(dtype="float64", default=2.5 | nbody_system.speed ,description="KD95 v0 parameter"),
            "halo_q": dict(dtype="float64", default=1.  ,description="KD95 q parameter"),
            "halo_coreparam": dict(dtype="float64", default=0.0 ,description="KD95 core parameter"), 
            "halo_Ra": dict(dtype="float64", default=4.| nbody_system.length ,description="KD95 Ra parameter"),
            "disk_type_parameter": dict(dtype="int32", default=3 ,description="type of disk (0=no,1=KD95,2=gas disk,3=KD95 + gas)"), 
            "disk_mass": dict(dtype="float64", default=50. ,description="(approx) stellar disk mass"), 
            "disk_scale_length": dict(dtype="float64", default=3. | nbody_system.length,description="stellar disk exponential scale length"),
            "disk_outer_radius": dict(dtype="float64", default=13. | nbody_system.length,description="stellar disk outer radius"), 
            "disk_scale_height": dict(dtype="float64", default=0.6 | nbody_system.length,description="stellar disk scale height"), 
            "disk_truncation_dr": dict(dtype="float64", default=1. | nbody_system.length,description="stellar disk truncation width"), 
            "gas_disk_mass": dict(dtype="float64", default=5. ,description="gas disk mass"), 
            "gas_disk_scale_length": dict(dtype="float64", default=3. | nbody_system.length,description="gas disk scale length (1/r profile)"), 
            "gas_disk_outer_radius": dict(dtype="float64", default=26. | nbody_system.length,description="gas disk outer radius"), 
            "gas_disk_truncation_dr": dict(dtype="float64", default=1. | nbody_system.length,description="gas disk truncation width"), 
            "gas_disk_sound_speed": dict(dtype="float64", default=0.12 | nbody_system.speed ,description="gas disk sound speed"),
            "gas_disk_gamma": dict(dtype="float64", default=1. ,description="gas disk polytropic index"),
            "gas_disk_number_of_radial_bins": dict(dtype="int32", default=50 ,description="gas disk nr of radial bins"), 
            "gas_disk_max_z": dict(dtype="float64", default=5. | nbody_system.length,description="gas disk max z of grid"), 
            "bulge_type_parameter": dict(dtype="int32", default=3 ,description="type of bulge (0=no, 1=KD95 (untested), 3=spherical df)"), 
            "bulge_cutoff_potential": dict(dtype="float64", default= -20 | nbody_system.potential ,description="KD95 cutoff potential"),
            "bulge_velocity_dispersion": dict(dtype="float64", default= 1. | nbody_system.speed ,description="KD95 velocity dispersion"),
            "bulge_scale_radius": dict(dtype="float64", default=0.8 | nbody_system.length,description="bulge scale radius"),
            "bulge_cutoff_radius": dict(dtype="float64", default=7. | nbody_system.length,description="bulge cutoff radius"),
            "bulge_density_parameter": dict(dtype="float64", default=5.5 | nbody_system.density,description="bulge density parameter"), 
            "radial_grid_delta_r": dict(dtype="float64", default=0.025 | nbody_system.length,description="radial grid cell size"), 
            "number_of_grid_intervals": dict(dtype="int32", default=20000 ,description="number of radial grid points"),
            "order_of_multipole_expansion": dict(dtype="int32", default=12 ,description="order of multipole expansion"), 
            "disk_central_radial_vdisp_over_z_vdisp": dict(dtype="float64", default=1. ,description="disk central rdial velocity dispersion in units of z dispersion"),
            "disk_scale_length_of_radial_vdisp": dict(dtype="float64", default=3. | nbody_system.length ,description="scale length of radial vel dispersion"), 
            "number_of_radial_steps_correction_fns": dict(dtype="int32", default=200 ,description="number of radial steps correction functions"),
            "number_of_iterations": dict(dtype="int32", default=12 ,description="number of diskdf iterations"), 
            "halo_streaming_fraction": dict(dtype="float64", default=0.5 ,description="halo streaming factor (to impose rotation, 0.5 is no rotation)"),
            "halo_number_of_particles": dict(dtype="int32", default=10000 ,description="number of particles in halo"),
            "halo_random_seed": dict(dtype="int64", default=0 ,description="halo random seed"), 
            "halo_do_center_flag": dict(dtype="bool", default=False ,description="whether to independently center halo"), 
            "bulge_streaming_fraction": dict(dtype="float64", default=0.5 ,description="bulge streaming fraction (to impose rotation)"),
            "bulge_number_of_particles": dict(dtype="int32", default=1000 ,description="number of bulge particles"), 
            "bulge_random_seed": dict(dtype="int64", default=12345678 ,description="bulge random seed"),
            "bulge_do_center_flag": dict(dtype="bool", default=False ,description="whether to independently center bulge "),
            "disk_number_of_particles": dict(dtype="int32", default=1000 ,description="stellar disk number of particles"), 
            "disk_random_seed": dict(dtype="int64", default=98765432 ,description="stellar disk random seed"),
            "disk_do_center_flag": dict(dtype="bool", default=False ,description="whether to independently center the stellar disk"),
            "gas_disk_number_of_particles": dict(dtype="int32", default=1000 ,description="number of gas particles"),
            "gas_disk_velocity_dispersion": dict(dtype="float64", default=0. ,description="velocity dispersion of gas particles"),
            "gas_disk_random_seed": dict(dtype="int64", default=543212345 ,description="random seed for gas disk"),
            "reuse_cached_model" : dict(dtype="bool", default=True, description="if True reuse cached preexisting model with the same parameters"),
}

class GaslactICsImplementation(object):
    
    def __init__(self):
        self._output_directory = "./"
        self._particles_generated = False
        self._particle_data = numpy.array([])
        self._bin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gbin")
    
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
    for par in parameters:
        exec("def get_"+par+"(self, value): value.value = self._"+par+"; return 0")
        exec("def set_"+par+"(self, value): self._"+par+" = value; return 0")
    
    def set_default_parameter_values(self):
        for par in parameters:
            if hasattr(parameters[par]["default"],"unit"):
              exec("self._"+par+"=parameters[par]['default'].number")
            else:
              exec("self._"+par+"=parameters[par]['default']")
            
    def cleanup_code(self):
        return 0
    
    def generate_in_dbh_string(self):
        in_dbh=[]
        
        in_dbh.append( str(self._halo_type_parameter) )
        if self._halo_type_parameter==1:
          in_dbh.append( str(self._halo_central_potential)+" "+
                         str(self._halo_v0)+" "+
                         str(self._halo_q)+" "+
                         str(self._halo_coreparam)+" "+
                         str(self._halo_Ra) )        
        if self._halo_type_parameter==3:
          in_dbh.append( str(self._halo_scale_radius)+" "+str(self._halo_virial_radius)+" "+str(self._halo_density_parameter) )

        in_dbh.append(str(self._disk_type_parameter))
        if self._disk_type_parameter in [1,3]:
          in_dbh.append( str(self._disk_mass)+" "+
                         str(self._disk_scale_length)+" "+
                         str(self._disk_outer_radius)+" "+
                         str(self._disk_scale_height)+" "+
                         str(self._disk_truncation_dr) )
        if self._disk_type_parameter in [2,3]:
          in_dbh.append( str(self._gas_disk_mass)+" "+
                         str(self._gas_disk_scale_length)+" "+
                         str(self._gas_disk_outer_radius)+" "+
                         str(self._gas_disk_truncation_dr) )
          in_dbh.append( str(self._gas_disk_sound_speed)+" "+
                         str(self._gas_disk_gamma)+" "+
                         str(self._gas_disk_number_of_radial_bins)+" "+
                         str(self._gas_disk_max_z) )

        in_dbh.append(str(self._bulge_type_parameter))
        if self._bulge_type_parameter==1:
          in_dbh.append( str(self._bulge_density_parameter)+" "+
                         str(self._bulge_cutoff_potential)+" "+
                         str(self._bulge_velocity_dispersion) )
        if self._bulge_type_parameter==3:
          in_dbh.append( str(self._bulge_scale_radius)+" "+
                         str(self._bulge_cutoff_radius)+" "+
                         str(self._bulge_density_parameter) )

        in_dbh.append( str(self._radial_grid_delta_r)+" "+
                       str(self._number_of_grid_intervals) )
        in_dbh.append( str(self._order_of_multipole_expansion))
        in_dbh.append( "")
          
        in_dbh="\n".join(in_dbh)

#        print "in_dbh:\n", in_dbh
        return in_dbh
    
    def generate_in_diskdf_string(self):
        in_diskdf=[]
        
        in_diskdf.append( str(self._disk_central_radial_vdisp_over_z_vdisp)+" "+str(self._disk_scale_length_of_radial_vdisp))
        in_diskdf.append( str(self._number_of_radial_steps_correction_fns))
        in_diskdf.append( str(self._number_of_iterations))

        in_diskdf.append("diskdf.ps/ps")

        in_diskdf="\n".join(in_diskdf)
#        print "in_diskdf:\n", in_diskdf
        return in_diskdf
    
    def generate_in_halo_string(self):
        in_halo = "{0:.15}\n".format(self._halo_streaming_fraction)
        in_halo += "{0}\n".format(self._halo_number_of_particles)
        in_halo += "{0}\n".format(self._halo_random_seed)
        in_halo += "{0}\n".format(1 if self._halo_do_center_flag else 0)
        in_halo += "dbh.dat\n"
#        print "in_halo:\n", in_halo
        return in_halo
    
    def generate_in_bulge_string(self):
        in_bulge = "{0:.15}\n".format(self._bulge_streaming_fraction)
        in_bulge += "{0}\n".format(self._bulge_number_of_particles)
        in_bulge += "{0}\n".format(self._bulge_random_seed)
        in_bulge += "{0}\n".format(1 if self._bulge_do_center_flag else 0)
        in_bulge += "dbh.dat\n"
#        print "in_bulge:\n", in_bulge
        return in_bulge
    
    def generate_in_disk_string(self):
        in_disk = "{0}\n".format(self._disk_number_of_particles)
        in_disk += "{0}\n".format(self._disk_random_seed)
        in_disk += "{0}\n".format(1 if self._disk_do_center_flag else 0)
        in_disk += "dbh.dat\n"
#        print "in_disk:\n", in_disk
        return in_disk

    def generate_in_gas_string(self):
        in_gas = "{0}\n".format(self._gas_disk_number_of_particles)
        in_gas += "{0}\n0.\n".format(self._gas_disk_velocity_dispersion)
        in_gas += "{0}\n".format(self._gas_disk_random_seed)
#        print "in_gas:\n", in_gas
        return in_gas

    def _new_dbh_dir(self, data_directory,in_dbh,in_diskdf):
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)
#        with open(os.path.join(data_directory, "in.gendenspsi"), "w") as f:
#            f.write("2000 40\n")
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
        
    def model_present(self,model_present):
        in_dbh = self.generate_in_dbh_string()
        in_diskdf = self.generate_in_diskdf_string()
        data_directory=self._data_directory(in_dbh,in_diskdf)
        model_present.value=self._directory_contains_valid_model(data_directory)      
        return 0

    def _directory_contains_valid_model(self,data_directory):
        if os.path.exists(os.path.join( data_directory)) and \
           os.path.exists(os.path.join( data_directory, 'dbh.dat')) and \
           os.path.exists(os.path.join( data_directory, 'dbh.finished')) and \
           os.path.exists(os.path.join( data_directory, 'getfreqs.finished')) and \
           (os.path.exists(os.path.join( data_directory, 'diskdf.finished')) or self._disk_type_parameter in [0,2]):
             return True
        else:
             return False
    
    def _location_dbh_dat(self, in_dbh,in_diskdf):
        data_directory=self._data_directory(in_dbh,in_diskdf)
        
        if self._directory_contains_valid_model(data_directory) and self._reuse_cached_model:
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
            if not is_new:
                return 0
            print("Writing output to:", self._cwd)

            print("\n running dbh \n\n")
            
            proc=Popen([os.path.join(self._bin_path, "dbh")], 
                cwd = self._cwd, stdin = PIPE, stdout = sys.stdout, stderr = sys.stderr)
            proc.communicate(in_dbh.encode())
            if proc.returncode==0:
              open(os.path.join(dbh_dir,"dbh.finished"),'a').close()
            else:
              raise Exception("dbh fail")

            print("\n running getfreqs \n\n")
            
            proc=Popen([os.path.join(self._bin_path, "getfreqs")], 
                cwd = self._cwd, stdin = PIPE, stdout = sys.stdout, stderr = sys.stderr)
            proc.communicate()
            if proc.returncode==0:
              open(os.path.join(dbh_dir,"getfreqs.finished"),'a').close()
            else:
              raise Exception("getfreqs fail")

            print("\n running diskdf \n\n")

            if self._disk_type_parameter in [1,3]:
              proc=Popen([os.path.join(self._bin_path, "diskdf")], 
                      cwd = self._cwd, stdin = PIPE, stdout = sys.stdout, stderr = sys.stderr)
              proc.communicate(in_diskdf.encode())
              if proc.returncode==0:
                open(os.path.join(dbh_dir,"diskdf.finished"),'a').close()
              else:
                raise Exception("diskdf fail")
                                
            return 0
        except Exception as ex:
            print("Exception occurred in commit_parameters:", ex)
            return -1
    
    def recommit_parameters(self):
        return self.commit_parameters()
    
    def get_number_of_particles_updated(self, number_of_particles_updated,number_of_gas_particles_updated):
        if self._particles_generated:
            number_of_particles_updated.value = self._number_of_particles_updated
            number_of_gas_particles_updated.value = self._number_of_gas_particles_updated
            self._particles_generated = False
        else:
            number_of_particles_updated.value = 0
        return 0
    
    def generate_particles(self):
        try:
            if self._disk_type_parameter in [1,3] and self._disk_number_of_particles>0:
                in_disk  = self.generate_in_disk_string()
                process = Popen([os.path.join(self._bin_path, "gendisk")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out,err=process.communicate(in_disk.encode())
                if process.returncode != 0:
                    print("error:", err)
                    return -2
                print(err ," ****")
                disk_data=numpy.frombuffer(out,dtype="float32")
            else:
                disk_data=numpy.array([])    

            if self._disk_type_parameter in [2,3] and self._gas_disk_number_of_particles>0:
                in_gas  = self.generate_in_gas_string()
                process = Popen([os.path.join(self._bin_path, "gengas")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out,err=process.communicate(in_gas.encode())
                if process.returncode != 0:
                    print("error:", err)
                    return -2
                gas_data=numpy.frombuffer(out,dtype="float32")
            else:
                gas_data=numpy.array([])

            if self._bulge_type_parameter in [1,3] and self._bulge_number_of_particles>0:
                in_bulge = self.generate_in_bulge_string()
                process = Popen([os.path.join(self._bin_path, "genbulge")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out,err=process.communicate(in_bulge.encode())
                if process.returncode != 0:
                    print("error:", err)
                    return -3 
                bulge_data=numpy.frombuffer(out,dtype="float32")
            else:
                bulge_data=numpy.array([])    
                            
            if self._halo_type_parameter in [1,3] and self._halo_number_of_particles>0:
                in_halo  = self.generate_in_halo_string()
                process = Popen([os.path.join(self._bin_path, "genhalo")], 
                    cwd = self._cwd, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                out, err = process.communicate(in_halo.encode())
                if process.returncode != 0:
                    print("error:", err)
                    return -4 
                halo_data=numpy.frombuffer(out,dtype="float32")
            else:
                halo_data=numpy.array([])    
                        
            self._number_of_particles_updated = len(gas_data)//8+(len(halo_data)+len(bulge_data)+len(disk_data))//7
            self._number_of_gas_particles_updated = len(gas_data)//8
            self._number_of_halo_particles=len(halo_data)//7
            self._number_of_bulge_particles=len(bulge_data)//7
            self._number_of_disk_particles=len(disk_data)//7
            self._number_of_gas_particles=len(gas_data)//8
            gas_posdata=gas_data[8*numpy.arange(self._number_of_gas_particles*7)//7]
            
            gammafactor=1. if self._gas_disk_gamma==1 else 1/self._gas_disk_gamma/(self._gas_disk_gamma-1)
            self._gas_internal_energy=gammafactor*gas_data[8*numpy.arange(self._number_of_gas_particles)+7]**2
            data=numpy.concatenate((gas_posdata,disk_data,bulge_data,halo_data))  
            self._particle_data = numpy.reshape(data,( self._number_of_particles_updated,7))
            self._particles_generated = True
            return 0
        except Exception as ex:
            print("Exception occurred in generate_particles:", ex)
            return -1
    
    def get_number_of_particles(self,nhalo,nbulge,ndisk,ngas):
        try:
            nhalo.value=self._number_of_halo_particles
            nbulge.value=self._number_of_bulge_particles
            ndisk.value=self._number_of_disk_particles
            ngas.value=self._number_of_gas_particles
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

    def get_internal_energy(self, index_of_the_particle, u, length):
        try:
            u.value = self._gas_internal_energy[index_of_the_particle]
            return 0
        except:        
            return -1

class GaslactICsInterface(PythonCodeInterface, CommonCodeInterface, LiteratureReferencesMixIn,
        CodeWithDataDirectories):
    """
    GaslactICs is a variant of GalactICs which allows for the inclusion of a gas 
    disk. It generates self-consistent disc-bulge-halo-gas galaxy models. 
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
        .. [#] Pelupessy, F. I. et al., 2013, The Astrophysical Multipurpose Software Environment, 
               Astronomy and Astrophysics 557, 84 [2013A&A...557A..84P] (gas version)
    """
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, GaslactICsImplementation, **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    def _check_if_worker_is_up_to_date(self):
        if not os.path.exists(os.path.join(GaslactICsImplementation()._bin_path, "dbh")):
            raise exceptions.CodeException(
                "The worker code of the '{0}' interface class is not up to date.\n"
                "Please do a 'make clean; make' in the root directory.".format(type(self).__name__))
    
    @option(type="string", sections=('data',))
    def output_data_root_directory(self):
        """
        The root directory of the output data,
        read - write directory
        """
        return os.path.join(get_amuse_root_dir(), 'data')
        
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'gaslactics', 'output')
    
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
    
    #parameter getters & setters
    # boolean parameters
    for par in parameters:
        dtype=parameters[par]["dtype"]
        if hasattr(parameters[par]["default"],"unit"):
          unit=parameters[par]["default"].unit.reference_string()
        else:
          unit="None"
        exec("@legacy_function\ndef get_"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n"
            "  function.result_type = 'int32'\n  return function")
        exec("@legacy_function\ndef set_"+par+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('value', dtype='"+dtype+"', direction=function.IN, unit="+unit+")\n"
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
        function.addParameter('number_of_gas_particles_updated', dtype='int32', direction=function.OUT)
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
        function.addParameter('number_of_gas_particles', dtype='int32', direction=function.OUT)
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
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.OUT, 
          description = "internal energy of gas particle", unit=nbody_system.speed**2)
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

class GaslactICs(CommonCode):
    
    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        InCodeComponentImplementation.__init__(self, GaslactICsInterface(**options), **options)
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        self.parameters.set_defaults()
        self.ensure_data_directory_exists(self.get_output_directory())
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
        for par in parameters:
            if parameters[par]["dtype"]=="bool":
              handler.add_boolean_parameter(
                  "get_"+par,
                  "set_"+par,
                  par,
                  parameters[par]["description"],
                  default_value=parameters[par]["default"]
              )
            else:
              handler.add_method_parameter(
                  "get_"+par,
                  "set_"+par,
                  par,
                  parameters[par]["description"],
                  default_value=parameters[par]["default"]
              )
        
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method("generate_particles", (), (handler.ERROR_CODE,))
        handler.add_method("get_number_of_particles_updated", (), (handler.NO_UNIT,handler.NO_UNIT, handler.ERROR_CODE,))
        
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
        
        for par in parameters:
            if hasattr(parameters[par]["default"],"unit"):
              unit=parameters[par]["default"].unit
            else:
              unit=handler.NO_UNIT          
            handler.add_method("get_"+par, (), (unit, handler.ERROR_CODE,))
            handler.add_method("set_"+par, (unit, ), (handler.ERROR_CODE,))
    
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

        handler.define_set('gas_particles', 'index_of_the_particle')
        handler.set_new('gas_particles', 'new_particle')
        handler.set_delete('gas_particles', 'delete_particle')
        handler.add_getter('gas_particles', 'get_mass')
        handler.add_getter('gas_particles', 'get_position')
        handler.add_getter('gas_particles', 'get_velocity')
        handler.add_getter('gas_particles', 'get_internal_energy')
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','PARAMETER_CHANGE_A','invoke_state_change2')
        handler.add_transition('EDIT','PARAMETER_CHANGE_B','invoke_state_change2')
        handler.add_transition('PARAMETER_CHANGE_A','RUN','recommit_parameters')
        handler.add_transition('PARAMETER_CHANGE_B','EDIT','recommit_parameters')
        handler.add_transition('EDIT', 'UPDATE', 'generate_particles', False)
        handler.add_transition('UPDATE', 'RUN', 'update_particle_set')
        handler.add_transition('RUN', 'EDIT', 'clear_particle_set')
        handler.add_method('RUN', 'invoke_state_change_updated')
        handler.add_method('EDIT', 'get_number_of_particles_updated')
        handler.add_method('UPDATE', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_mass')
        handler.add_method('RUN', 'get_position')
        handler.add_method('RUN', 'get_velocity')

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
        handler.add_method('RUN', 'get_internal_energy')

    def commit_parameters(self):
        if not self.model_present():
          print("generating galaxy model, this may take a while...")
        self.overridden().commit_parameters()  

    def recommit_parameters(self):
        if not self.model_present():
          print("(re)generating galaxy model, this may take a while...")
        self.overridden().recommit_parameters()  

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
        number_of_updated_particles,number_of_updated_gas_particles = self.get_number_of_particles_updated()
        if number_of_updated_particles:
            self.particles._private.attribute_storage._add_indices(
                list(range(number_of_updated_gas_particles,number_of_updated_particles))
            ) # this should generate disjoint sets (gas_particles not in particles) 
        if number_of_updated_gas_particles:
            self.gas_particles._private.attribute_storage._add_indices(
                list(range(number_of_updated_gas_particles))
            )

    def clear_particle_set(self):
        if len(self.particles):
            self.particles.remove_particles(self.particles)
        if len(self.gas_particles):
            self.gas_particles.remove_particles(self.gas_particles)

    @property
    def halo_particles(self):
        nhalo,nbulge,ndisk,ngas=self.get_number_of_particles()
        return self.particles[ngas+ndisk+nbulge:]

    @property
    def bulge_particles(self):
        nhalo,nbulge,ndisk,ngas=self.get_number_of_particles()
        return self.particles[ngas+ndisk:ngas+ndisk+nbulge]

    @property
    def disk_particles(self):
        nhalo,nbulge,ndisk,ngas=self.get_number_of_particles()
        return self.particles[ngas:ngas+ndisk]
