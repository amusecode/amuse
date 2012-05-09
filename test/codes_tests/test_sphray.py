import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.sphray.interface import SPHRayInterface,SPHRay
from amuse.units import units
from amuse.datamodel import Particles,create_particle_set
from amuse.community import ensure_data_directory_exists

from amuse.io import read_set_from_file

#default_options = dict(redirection="none")
default_options={}

class TestSPHRayInterface(TestWithMPI):

    def is_fortan_version_up_to_date(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'compilers')
            if is_configured:
                is_configured = hasattr(config.compilers, 'gfortran_version')
        except ImportError:
            is_configured = False
    
        if is_configured:
            if not config.compilers.gfortran_version:
                return True
            
            try:
                parts = [int(x) for x in config.compilers.gfortran_version.split('.')]
            except:
                parts = []
                
            if len(parts) < 2:
                return True
                
            return parts[0] >= 4 and parts[1] > 1
        else:
            return True
            
    def setUp(self):
        super(TestWithMPI, self).setUp()
        self.check_fortran_version()
        
    def check_fortran_version(self):
        if not self.is_fortan_version_up_to_date():
            self.skip('cannot compile, fortran module names cannot be resolved correctly in this gfortran version')
            
            
    def test1(self):
        print "Test 1: initialization"
        
        instance = SPHRayInterface(**default_options)
        self.assertEqual(0, instance.initialize_code())
        ensure_data_directory_exists(instance.output_directory())
        self.assertEquals(0, instance.set_data_directory(instance.data_directory()))        
        self.assertEquals(0, instance.set_output_directory(instance.output_directory()))        
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        print "Test 2: add, commit_particles"
        instance = SPHRayInterface(**default_options)
        ensure_data_directory_exists(instance.output_directory())
        self.assertEqual(0, instance.set_data_directory(instance.data_directory()))
        self.assertEquals(0, instance.set_output_directory(instance.output_directory()))        
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory(), 'sphray_4K')
        mass, hsml, x, y, z, rho, xe, u = self.read_gas_file(input_file)
        number_of_gas_particles = len(x)
        indices, errors = instance.new_gas_particle(mass, hsml, x, y, z, rho, xe, u)
        self.assertEqual(errors, [0]*number_of_gas_particles)
        self.assertEqual(indices, range(1,number_of_gas_particles+1))
        
        input_file = os.path.join(instance.data_directory(), 'test1_sources_001.1')
        L, xs, ys, zs, spctype = self.read_src_file(input_file)
        number_of_src_particles = len(xs)
        s_indices, errors = instance.new_src_particle(L, xs, ys, zs, spctype)
        self.assertEqual(errors, [0]*number_of_src_particles)
        self.assertEqual(s_indices, range(number_of_gas_particles+1,number_of_src_particles+number_of_gas_particles+1))
       
        self.assertEqual(0, instance.commit_particles())
        mass2, hsml2, x2, y2, z2, rho2, xe2, u2 , error = instance.get_state_gas(indices)
        self.assertAlmostEqual((mass-mass2), numpy.zeros_like(x), 7)
        self.assertAlmostEqual((hsml-hsml2), numpy.zeros_like(x), 7)
        self.assertAlmostEqual((x-x2)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((y-y2)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((z-z2)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((rho-rho2), numpy.zeros_like(x), 7)
#        self.assertAlmostEqual((xe-xe2), numpy.zeros_like(x), 7)
        self.assertAlmostEqual((u-u2), numpy.zeros_like(x), 7)

        L2, xs2, ys2, zs2, spctype2 , error = instance.get_state_src(s_indices)
        self.assertAlmostEqual((L-L2), numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((xs-xs2)/13200., numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((ys-ys2)/13200., numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((zs-zs2)/13200., numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((spctype-spctype2), numpy.zeros_like(xs), 7)

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test3(self):
        print "Test 3: add, commit_particles, setters, remove"
        instance = SPHRayInterface(**default_options)
        ensure_data_directory_exists(instance.output_directory())
        self.assertEqual(0, instance.set_data_directory(instance.data_directory()))
        self.assertEquals(0, instance.set_output_directory(instance.output_directory()))        
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory(), 'sphray_4K')
        mass, hsml, x, y, z, rho, xe, u = self.read_gas_file(input_file)
        number_of_gas_particles = len(x)
        indices, errors = instance.new_gas_particle(mass, hsml, x, y, z, rho, xe, u)
        self.assertEqual(errors, [0]*number_of_gas_particles)
        self.assertEqual(indices, range(1,number_of_gas_particles+1))
        
        input_file = os.path.join(instance.data_directory(), 'test1_sources_001.1')
        L, xs, ys, zs, spctype = self.read_src_file(input_file)
        number_of_src_particles = len(xs)
        s_indices, errors = instance.new_src_particle(L, xs, ys, zs, spctype)
        self.assertEqual(errors, [0]*number_of_src_particles)
        self.assertEqual(s_indices, range(number_of_gas_particles+1,number_of_src_particles+number_of_gas_particles+1))
       
        self.assertEqual(0, instance.commit_particles())
        mass2, hsml2, x2, y2, z2, rho2, xe2, u2 , error = instance.get_state_gas(indices)
        self.assertAlmostEqual((mass-mass2), numpy.zeros_like(x), 7)
        self.assertAlmostEqual((hsml-hsml2), numpy.zeros_like(x), 7)
        self.assertAlmostEqual((x-x2)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((y-y2)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((z-z2)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((rho-rho2), numpy.zeros_like(x), 7)
#        self.assertAlmostEqual((xe-xe2), numpy.zeros_like(x), 7)
        self.assertAlmostEqual((u-u2), numpy.zeros_like(x), 7)

        L2, xs2, ys2, zs2, spctype2 , error = instance.get_state_src(s_indices)
        self.assertAlmostEqual((L-L2), numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((xs-xs2)/13200., numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((ys-ys2)/13200., numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((zs-zs2)/13200., numpy.zeros_like(xs), 7)
        self.assertAlmostEqual((spctype-spctype2), numpy.zeros_like(xs), 7)

        error=instance.set_state_gas(100,1.,2.,3.,4.,5.,6.,7.,8.)
        self.assertEqual(0, error)
        error=instance.set_state_gas(10000,1.,2.,3.,4.,5.,6.,7.,8.)
        self.assertEqual(-1, error)
 
        m,h,x,y,z,rho,xe,u,error=instance.get_state_gas(100)
        self.assertEqual(0, error)
        self.assertAlmostEqual(m, 1., 7)
        self.assertAlmostEqual(h, 2., 7)
        self.assertAlmostEqual(x, 3., 7)
        self.assertAlmostEqual(y, 4., 7)
        self.assertAlmostEqual(z, 5., 7)
        self.assertAlmostEqual(rho, 6., 7)
#        self.assertAlmostEqual(xe, 7., 7)
        self.assertAlmostEqual(u, 8., 7)

        error=instance.remove_gas_particle(100)
        self.assertEqual(0, error)
        error=instance.remove_gas_particle(100)
        self.assertEqual(-4, error)


        error=instance.set_state_src(4097,1.,2.,3.,4.,5.)
        self.assertEqual(0, error)
        error=instance.set_state_src(10000,1.,2.,3.,4.,5.)
        self.assertEqual(-1, error)
 
        L,x,y,z,s,error=instance.get_state_src(4097)
        self.assertEqual(0, error)
        self.assertAlmostEqual(L, 1., 7)
        self.assertAlmostEqual(x, 2., 7)
        self.assertAlmostEqual(y, 3., 7)
        self.assertAlmostEqual(z, 4., 7)
        self.assertAlmostEqual(s, 5., 7)

        error=instance.remove_src_particle(4097)
        self.assertEqual(0, error)
        error=instance.remove_src_particle(4097)
        self.assertEqual(-4, error)

        instance.recommit_particles()
        error=instance.remove_gas_particle(100)
        self.assertEqual(-1, error)        
        error=instance.remove_gas_particle(4097)
        self.assertEqual(-1, error)


        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test4(self):
        print "Test 2: set, get time"
        instance = SPHRayInterface(**default_options)
        ensure_data_directory_exists(instance.output_directory())
        self.assertEqual(0, instance.set_data_directory(instance.data_directory()))
        self.assertEquals(0, instance.set_output_directory(instance.output_directory()))        
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory(), 'sphray_4K')
        mass, hsml, x, y, z, rho, xe, u = self.read_gas_file(input_file)
        number_of_gas_particles = len(x)
        indices, errors = instance.new_gas_particle(mass, hsml, x, y, z, rho, xe, u)
        self.assertEqual(errors, [0]*number_of_gas_particles)
        self.assertEqual(indices, range(1,number_of_gas_particles+1))
        
        input_file = os.path.join(instance.data_directory(), 'test1_sources_001.1')
        L, xs, ys, zs, spctype = self.read_src_file(input_file)
        number_of_src_particles = len(xs)
        s_indices, errors = instance.new_src_particle(L, xs, ys, zs, spctype)
        self.assertEqual(errors, [0]*number_of_src_particles)
        self.assertEqual(s_indices, range(number_of_gas_particles+1,number_of_src_particles+number_of_gas_particles+1))
       
        self.assertEqual(0, instance.commit_particles())
        
        time,err= instance.get_time()
        self.assertEqual(0,err)
        self.assertEqual(0.,time)
        err= instance.set_time(123.)
        self.assertEqual(-1,err)   # because we don't know whether it is safe yet
        time,err= instance.get_time()
        self.assertEqual(0,err)
        self.assertEqual(123.,time)

    def test2(self):
        instance=SPHRayInterface()
        instance.initialize_code()
        
        for x,l in [('isothermal',0), ('H_caseA',1),('He_caseA',1)]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,l))
            err=getattr(instance, 'set_'+x)(1)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,1))
            err=getattr(instance, 'set_'+x)(0)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals((x,result),(x,0))

        for x,l in [('raynumber',30000),('iontempsolver',2),('boundary',0)]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,l))
            err=getattr(instance, 'set_'+x)(1)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,1))
            err=getattr(instance, 'set_'+x)(0)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals((x,result),(x,0))

        for x,l in [('boxsize',13.2),("defaultspectype",-1.)]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertAlmostEqual( result,l ,6)
            err=getattr(instance, 'set_'+x)(1.)
            result,err=getattr(instance, 'get_'+x)()
            self.assertAlmostEqual( result,1. ,6)
            err=getattr(instance, 'set_'+x)(0)
            result,err=getattr(instance, 'get_'+x)()
            self.assertAlmostEqual( result,0. ,6)

        result,err=instance.get_globalHefraction()
        self.assertEqual(err, 0)
        self.assertEqual(result,0.)
        err=instance.set_globalHefraction(0.1)
        self.assertEqual(err, -2)

    def read_gas_file(self,filename):
        p=read_set_from_file(filename,'amuse')
        mass=p.mass.number
        hsml=p.smoothing_length.number
        x=p.x.number
        y=p.y.number
        z=p.z.number
        rho=p.rho.number
        u=p.internal_energy.number
        xe=numpy.zeros_like(x)
        return mass, hsml, x, y, z, rho, xe, u
        
    def read_src_file(self,filename):
        f=open(filename)
        lines=f.readlines()
        f.close()
        L=[]
        x=[]
        y=[]
        z=[]
        spctype=[]
        for line in lines:
          l=line.split()
          if len(l) == 9:
            L.append(float(l[6]))
            x.append(float(l[0]))
            y.append(float(l[1]))
            z.append(float(l[2]))
            spctype.append(float(l[7]))
        return numpy.array(L), numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(spctype)

class TestSPHRay(TestWithMPI):

    def test0(self):
        print "test1: basic startup and flow"
        instance=SPHRay()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.box_size = 100 | units.parsec
        self.assertEquals(instance.parameters.box_size, 100 | units.parsec)
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')

        self.assertEquals(instance.parameters.box_size, 100 | units.parsec)
        instance.cleanup_code()
        instance.stop()

    def test1(self):
        print "test1: adding particles"

        instance=SPHRay()

        gasparts=self.read_gas_file(os.path.join(instance.data_directory(), 'sphray_4K'))
        srcparts=self.read_src_file(os.path.join(instance.data_directory(), 'test1_sources_001.1'))
                
        self.assertEqual(len(instance.gas_particles),0)
        self.assertEqual(len(instance.src_particles),0)
        instance.gas_particles.add_particles(gasparts)
        instance.src_particles.add_particles(srcparts)
        self.assertEqual(len(instance.gas_particles),len(gasparts))
        self.assertEqual(len(instance.src_particles),len(srcparts))

        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')

        gaspart2=instance.gas_particles.copy()

        self.assertAlmostRelativeEquals(gasparts.position,gaspart2.position,6)
        self.assertAlmostRelativeEquals(gasparts.rho,gaspart2.rho,6)
        self.assertAlmostRelativeEquals(gasparts.u,gaspart2.u,6)

        instance.cleanup_code()
        instance.stop()

    def test2(self):
        print "test2: test parameters"
        instance=SPHRay()

        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')

        for par,val in [("isothermal_flag", False),
            ("hydrogen_case_A_flag", True),("helium_case_A_flag", True)]:
            val1=getattr(instance.parameters,par)
            self.assertEqual(val,val1)            
            setattr(instance.parameters, par, False)
            val1=getattr(instance.parameters,par)
            self.assertEqual(val1,False)            
            setattr(instance.parameters, par, True)
            val1=getattr(instance.parameters,par)
            self.assertEqual(val1,True)            
            setattr(instance.parameters, par, False)
            val1=getattr(instance.parameters,par)
            self.assertEqual(val1,False)            

        for par,val in [("number_of_rays", 30000),
            ("ionization_temperature_solver", 2),("boundary_condition", 0)]:
            val1=getattr(instance.parameters,par)
            self.assertEqual(val,val1)            
            setattr(instance.parameters, par, 123)
            val1=getattr(instance.parameters,par)
            self.assertEqual(val1,123)            

        for par,val,tval in [ ("box_size",13.2 | units.kpc,  10. | units.kpc),
                              ("default_spectral_type", -1.,1.)]:
            val1=getattr(instance.parameters,par)
            self.assertAlmostRelativeEqual(val,val1,6)            
            setattr(instance.parameters, par, tval)
            val1=getattr(instance.parameters,par)
            self.assertAlmostRelativeEqual(val1,tval,6)            

        
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        
        
        

    def read_gas_file(self,filename):
        p=read_set_from_file(filename,'amuse')
        mass=p.mass.number
        hsml=p.smoothing_length.number
        x=p.x.number
        y=p.y.number
        z=p.z.number
        rho=p.rho.number
        u=p.internal_energy.number
        xe=numpy.zeros_like(x)
        return create_particle_set(mass=mass | (10**10*units.MSun), hsml=hsml | (units.kpc), 
            x=x | (units.kpc), y=y| (units.kpc), z=z| (units.kpc), rho=rho | ((10**10*units.MSun) /(units.kpc)**3),
            xe=xe, u=u| (10**5 * units.cm/units.s)**2)
        
    def read_src_file(self,filename):
        f=open(filename)
        lines=f.readlines()
        f.close()
        L=[]
        x=[]
        y=[]
        z=[]
        spctype=[]
        for line in lines:
          l=line.split()
          if len(l) == 9:
            L.append(float(l[6]))
            x.append(float(l[0]))
            y.append(float(l[1]))
            z.append(float(l[2]))
            spctype.append(float(l[7]))
        return create_particle_set( L=numpy.array(L) | (10**48 * units.s**-1), 
            x=numpy.array(x) | (units.kpc), y=numpy.array(y)| (units.kpc), 
                  z=numpy.array(z)| (units.kpc), SpcType=numpy.array(spctype))
