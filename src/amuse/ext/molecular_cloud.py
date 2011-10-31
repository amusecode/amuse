import numpy

from math import sqrt

from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.ext.evrard_test import uniform_unit_cube
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import generic_unit_converter
from amuse.units import units

from amuse import datamodel
def make_ifft_real(nf,vi):
    if vi.ndim==3:
# body of cube
        vi[1:nf,1:2*nf,1:2*nf]=numpy.conj(vi[2*nf-1:nf:-1,2*nf-1:0:-1,2*nf-1:0:-1])

# 3 lower + middle planes
        vi[0,1:nf,1:2*nf]=numpy.conj(vi[0,2*nf-1:nf:-1,2*nf-1:0:-1])
        vi[1:nf,0,1:2*nf]=numpy.conj(vi[2*nf-1:nf:-1,0,2*nf-1:0:-1])
        vi[1:nf,1:2*nf,0]=numpy.conj(vi[2*nf-1:nf:-1,2*nf-1:0:-1,0])
        vi[nf,1:nf,1:2*nf]=numpy.conj(vi[nf,2*nf-1:nf:-1,2*nf-1:0:-1])

# 7 lines
        vi[0,0,1:nf]=numpy.conj(vi[0,0,2*nf-1:nf:-1])
        vi[0,1:nf,0]=numpy.conj(vi[0,2*nf-1:nf:-1,0])
        vi[1:nf,0,0]=numpy.conj(vi[2*nf-1:nf:-1,0,0])
        vi[0,nf,1:nf]=numpy.conj(vi[0,nf,2*nf-1:nf:-1])
        vi[nf,0,1:nf]=numpy.conj(vi[nf,0,2*nf-1:nf:-1])
        vi[nf,nf,1:nf]=numpy.conj(vi[nf,nf,2*nf-1:nf:-1])
        vi[nf,1:nf,0]=numpy.conj(vi[nf,2*nf-1:nf:-1,0])

# 8 points
        vi[0,0,0]=2*numpy.real(vi[0,0,0])
        vi[nf,0,0]=2*numpy.real(vi[nf,0,0])
        vi[0,nf,0]=2*numpy.real(vi[0,nf,0])
        vi[nf,nf,0]=2*numpy.real(vi[nf,nf,0])

        vi[0,0,nf]=2*numpy.real(vi[0,0,nf])
        vi[nf,0,nf]=2*numpy.real(vi[nf,0,nf])
        vi[0,nf,nf]=2*numpy.real(vi[0,nf,nf])
        vi[nf,nf,nf]=2*numpy.real(vi[nf,nf,nf])
        return vi
    
    return -1  
  
def random_field(nf=32, power=-3., seed=None):
    if seed is not None:
        numpy.random.seed(seed)
    
    freq=numpy.mgrid[-nf:nf,-nf:nf,-nf:nf]   

    fi,fj,fk=freq

    fi=fi.flatten()
    fj=fj.flatten()
    fk=fk.flatten()
    
    norm=-numpy.log(numpy.random.uniform(0.,1.,len(fi)))*(fi**2+fj**2+fk**2+1.e-30)**(power/4)
    phase=numpy.random.uniform(0.,1.,len(fi))*2*numpy.pi
    vi=norm*numpy.exp(phase*1j)

    vi=vi.reshape(nf*2,nf*2,nf*2)

    vi[nf,nf,nf]=0.
    
    vi=make_ifft_real(nf,vi)

    vi=numpy.fft.ifftshift( vi)

    vi=numpy.fft.ifftn(vi)

    if vi.imag.max()>1.e-16:
        print "check random field"
    return vi

def make_div_free(nf,vx,vy,vz):

    vx=numpy.fft.fftn(vx)
    vx=vx.flatten()
    vy=numpy.fft.fftn(vy)
    vy=vy.flatten()
    vz=numpy.fft.fftn(vz)
    vz=vz.flatten()

    freq=numpy.mgrid[-nf:1.*nf,-nf:1.*nf,-nf:1.*nf] 
    fi,fj,fk=freq
    fi=numpy.fft.fftshift( fi)
    fj=numpy.fft.fftshift( fj)
    fk=numpy.fft.fftshift( fk)

    fi=fi.flatten()
    fj=fj.flatten()
    fk=fk.flatten()
    ff=fi*fi+fj*fj+fk*fk+1.e-30
    
    vdotf=(vx*fi+vy*fj+vz*fk)
    vx=vx-fi*vdotf/ff
    vy=vy-fj*vdotf/ff
    vz=vz-fk*vdotf/ff

    del fi,fj,fk,ff

    vx=vx.reshape(2*nf,2*nf,2*nf)
    vy=vy.reshape(2*nf,2*nf,2*nf)
    vz=vz.reshape(2*nf,2*nf,2*nf)

# zero out nyquist freq planes: strictly speaking this is too drastic....
# inside the nyquist planes only v// f x f_mirror needs to be enforced (methinks) 
    vx[nf,0:2*nf,0:2*nf]=0.
    vx[0:2*nf,nf,0:2*nf]=0.
    vx[0:2*nf,0:2*nf,nf]=0.
    vy[nf,0:2*nf,0:2*nf]=0.
    vy[0:2*nf,nf,0:2*nf]=0.
    vy[0:2*nf,0:2*nf,nf]=0.
    vz[nf,0:2*nf,0:2*nf]=0.
    vz[0:2*nf,nf,0:2*nf]=0.
    vz[0:2*nf,0:2*nf,nf]=0.

    vx=numpy.fft.ifftn(vx)
    vy=numpy.fft.ifftn(vy)
    vz=numpy.fft.ifftn(vz)

    if vx.imag.max()>1.e-16:
        print "check div-free field"
    if vy.imag.max()>1.e-16:
        print "check div-free field"
    if vz.imag.max()>1.e-16:
        print "check div-free field"

    return vx.real,vy.real,vz.real

def interpolate_trilinear(x,y,z,farray):

    if farray.ndim!=3:
        return -1
      
    nx,ny,nz=farray.shape    
    dx=2./nx
    dy=2./ny
    dz=2./nz

    fx,xint=numpy.modf((x+1)/dx)
    fy,yint=numpy.modf((y+1)/dy)
    fz,zint=numpy.modf((z+1)/dz)

    xint=xint.astype('i')
    yint=yint.astype('i')
    zint=zint.astype('i')

    xint1=numpy.mod(xint+1,nx)
    yint1=numpy.mod(yint+1,nx)
    zint1=numpy.mod(zint+1,nx)

    q111 = farray[xint, yint, zint]
    q211 = farray[xint1, yint, zint]
    q221 = farray[xint1, yint1, zint]
    q121 = farray[xint, yint1, zint]
    q112 = farray[xint, yint, zint1]
    q212 = farray[xint1, yint, zint1]
    q222 = farray[xint1, yint1, zint1]
    q122 = farray[xint, yint1, zint1]

    return (q222* fx*fy*fz +  
      q122* (1-fx)*fy*fz +  
      q212* fx*(1-fy)*fz +  
      q112* (1-fx)*(1-fy)*fz +  
      q221* fx*fy*(1-fz) +  
      q121* (1-fx)*fy*(1-fz) +  
      q211* fx*(1-fy)*(1-fz) +  
      q111* (1-fx)*(1-fy)*(1-fz))

class molecular_cloud(object):
    def __init__(self,nf=32,power=-3.,targetN=10000, ethep_ratio=0.01,
                   convert_nbody=None,ekep_ratio=1.,seed=None,base_grid=None):
        self.nf=nf
        self.power=power
        self.targetN=targetN
        self.convert_nbody=convert_nbody
        self.seed=seed
        self.base_grid=base_grid
        self.ethep_ratio=ethep_ratio
        self.ekep_ratio=ekep_ratio
      
    def new_model(self):
        if self.seed is not None:
            numpy.random.seed(self.seed)
        vx_field=random_field(self.nf,self.power)
        vy_field=random_field(self.nf,self.power)
        vz_field=random_field(self.nf,self.power)

        vx_field,vy_field,vz_field=make_div_free(self.nf,vx_field,vy_field,vz_field)

        base_sphere=uniform_unit_sphere(self.targetN,base_grid=self.base_grid)
        x,y,z=base_sphere.make_xyz()
        self.actualN=len(x)
        vx=interpolate_trilinear(x,y,z,vx_field)
        vy=interpolate_trilinear(x,y,z,vy_field)
        vz=interpolate_trilinear(x,y,z,vz_field)
        mass=numpy.ones_like(x)/self.actualN


        Ep=3./5
        self.internalE=Ep*self.ethep_ratio
        Ek=0.5*mass[0]*(vx**2+vy**2+vz**2).sum()
        vfac=sqrt(1/self.ekep_ratio*Ep/Ek)
        vx=vx*vfac
        vy=vy*vfac
        vz=vz*vfac
        Ek=0.5*mass[0]*(vx**2+vy**2+vz**2).sum()
        print self.internalE/Ep
        print Ek/Ep

        internal_energy=numpy.ones_like(x)*self.internalE
      
        return (mass,x,y,z,vx,vy,vz,internal_energy)

    @property
    def result(self):
        mass,x,y,z,vx,vy,vz,u = self.new_model()
        result = datamodel.Particles(self.actualN)
        result.mass = nbody_system.mass.new_quantity(mass)
        result.x = nbody_system.length.new_quantity(x)
        result.y = nbody_system.length.new_quantity(y)
        result.z = nbody_system.length.new_quantity(z)
        result.vx = nbody_system.speed.new_quantity(vx)
        result.vy = nbody_system.speed.new_quantity(vy)
        result.vz = nbody_system.speed.new_quantity(vz)
        result.u = (nbody_system.speed**2).new_quantity(u)

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()

        return result

class constant_density_div_free_power_law_v_ism_cube(object):
    def __init__(self,nf=32,power=-3.,targetN=10000, eketh_ratio=1.,
                   convert=None,seed=None,base_grid=None):
        self.nf=nf
        self.power=power
        self.targetN=targetN
        self.convert=convert
        self.seed=seed
        self.base_grid=base_grid
        self.eketh_ratio=eketh_ratio
      
    def new_model(self):
        if self.seed is not None:
            numpy.random.seed(self.seed)
        vx_field=random_field(self.nf,self.power)
        vy_field=random_field(self.nf,self.power)
        vz_field=random_field(self.nf,self.power)

        vx_field,vy_field,vz_field=make_div_free(self.nf,vx_field,vy_field,vz_field)

        base_cube=uniform_unit_cube(self.targetN,base_grid=self.base_grid)
        x,y,z=base_cube.make_xyz()
        self.actualN=len(x)
        vx=interpolate_trilinear(x,y,z,vx_field)
        vy=interpolate_trilinear(x,y,z,vy_field)
        vz=interpolate_trilinear(x,y,z,vz_field)
        mass=numpy.ones_like(x)/self.actualN

        Ek=0.5*mass[0]*(vx**2+vy**2+vz**2).sum()
        vfac=sqrt(1/self.eketh_ratio/Ek)
        vx=vx*vfac
        vy=vy*vfac
        vz=vz*vfac
        Ek=0.5*mass[0]*(vx**2+vy**2+vz**2).sum()
        print Ek

        internal_energy=numpy.ones_like(x)
      
        return (mass,x,y,z,vx,vy,vz,internal_energy)

    @property
    def result(self):
        mass,x,y,z,vx,vy,vz,u = self.new_model()
        result = datamodel.Particles(self.actualN)
        result.mass = nbody_system.mass.new_quantity(mass)
        result.x = nbody_system.length.new_quantity(x)
        result.y = nbody_system.length.new_quantity(y)
        result.z = nbody_system.length.new_quantity(z)
        result.vx = nbody_system.speed.new_quantity(vx)
        result.vy = nbody_system.speed.new_quantity(vy)
        result.vz = nbody_system.speed.new_quantity(vz)
        result.u = (nbody_system.speed**2).new_quantity(u)

        if not self.convert is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert.as_converter_from_si_to_generic())
            result = result.copy_to_memory()

        return result

def ism_cube(L=10| units.parsec,density=(1.14 | units.amu/units.cm**3), u=50 | (units. kms)**2):
    
    total_mass=density*L**3
    internalE=total_mass*u
    convert = generic_unit_converter.ConvertBetweenGenericAndSiUnits(total_mass, L, internalE)
    return constant_density_div_free_power_law_v_ism_cube(convert=convert)


if __name__=="__main__":
    cloud=ism_cube()
    parts=cloud.result
    print parts[0].u**0.5
    print len(parts)*parts[0].mass.in_(units.MSun)

    mu=1.4 | units.amu
    gamma1=1.6667-1
    print 'Temp:', (gamma1*min(parts.u)*mu/constants.kB).in_(units.K)

    total_mass=10000. | units.MSun
    radius=10. | units.parsec
    print 'dens:',(total_mass*3/4./3.1415/radius**3).in_(units.amu/units.cm**3) 
  
  
  
  
