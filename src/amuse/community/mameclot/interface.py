import numpy
import os

from amuse.community import LiteratureReferencesMixIn
from amuse.datamodel import Particles,ParticlesWithUnitsConverted

from subprocess import Popen, PIPE

from amuse.units import units,nbody_system

class mameclot(LiteratureReferencesMixIn):
    """  
    MAMECLOT: MAke ME a CLuster Or Two
  
    Construct initial conditions of 1 (or 2) cluster(s) for N-body simulations 
  
    A 2 cluster system with orbit in the x-y plane is made if the mass ratio q>0 
    Five equilibrium models are available, a cut-off radius can be provided
    The velocities can be isotropic or radially anisotropic (a la Osipkov-Merritt)
    Optional angular momentum by aligning the angular momentum vectors along z
    System is scaled to N-body units with G = M = -4E = 1 (Heggie & Mathieu 1986)
    In the case of a two-body system E = mu*e_orb + E_1 + E_2 (mu = reduced mass) 
    The 2-body orbit is computed from input q, e_hat and l_hat (and eta if q<1):
    
    e_hat = e_orb/0.5<sigma^2> 
    l_hat = l_orb/<r_vir*sigma>
    eta => defines the relative radii: r2/r1 = q^eta, for example:
    eta=[-0.33/0/0.33/0.5] for equal [t_rh/r_h/rho_h/Sigma_h]
    
    Relevant references:
        .. [#] Gieles, M., https://github.com/mgieles/mameclot

    :argument targetN: Number of particles
    :argument convert_nbody: When given will convert the resulting set to SI units
    :argument cluster_model: cluster model to generate ["Henon"]
    :argument mass_ratio: mass ratio, zero means only 1 cluster [0]
    :argument size_eta_parameter: eta, relative cluster sizes: r2/r1 = q^eta [0.333]
    :argument imf: IMF [single mass]
    :argument angular_momentum_signs: Angular momentum in z-direction, -, + or x ["xx"]
    :argument fraction_of_max_rot_energy: Fraction of maximum rotational energy   [1],
    :argument OsipkovMerritt_anisotropy_radius_r0: Osipkov-Merritt anisotropy radius in units of r_0 [999]
    :argument cutoff_radius_halfmass_radius: Cut-off radius in units of r_h [20]
    :argument physical_radius: Physical scale [1| units.parsec]
    :argument distance: Distance between 2 clusters in N-body units [20 | nbody_system.length]
    :argument orbital_energy: Dimensionless orbital energy of two-cluster system [0]
    :argument orbital_angular_momentum: Dimensionless orbital angular momentum two-cluster system [4]
    :argument seed: random seed 123456
    """
    
    available_cluster_models={"Dehnen" : 0, "Hernquist":1 ,"Jaffe":2,"Henon": 3, "Plummer" :4 }
    available_imf_model={"single mass" : 0, "Kroupa" : 1}
    angular_momentum_sign_key={"xx":0,"++":1,"+x":2,"+-":3,"x+":4,"x-":5,"--":-1,"-x":-2,"-+":-3}
    
    def __init__(self, targetN=10000, cluster_model="Henon", mass_ratio=0, size_eta_parameter=0.3333333,
        imf="single mass", angular_momentum_signs="xx", fraction_of_max_rot_energy=1,
        OsipkovMerritt_anisotropy_radius_r0=999,cutoff_radius_halfmass_radius=20,physical_radius=1| units.parsec,
        distance=20 | nbody_system.length,orbital_energy=0,orbital_angular_momentum=0,seed=123456,
        convert_to_physical=False):
        
        LiteratureReferencesMixIn.__init__(self)

        self._bin_path = os.path.dirname(os.path.abspath(__file__))
        self._particles=None
        self.convert_to_physical=convert_to_physical
        self._exec="mameclot_worker"
        
        self.targetN=targetN
        self.cluster_model=cluster_model
        self.mass_ratio=mass_ratio
        self.size_eta_parameter=size_eta_parameter
        self.imf=imf
        self.angular_momentum_signs=angular_momentum_signs
        self.fraction_of_max_rot_energy=fraction_of_max_rot_energy
        self.OsipkovMerritt_anisotropy_radius_r0=OsipkovMerritt_anisotropy_radius_r0
        self.cutoff_radius_halfmass_radius=cutoff_radius_halfmass_radius
        self.physical_radius=physical_radius
        self.distance=distance
        self.orbital_energy=orbital_energy
        self.orbital_angular_momentum=orbital_angular_momentum
        self.seed=seed        

    def arguments(self):
        arguments=[]
        arguments.extend([ "-N", str(self.targetN) ])
        arguments.extend([ "-m", str(self.available_cluster_models[self.cluster_model]) ])
        arguments.extend([ "-q", str(self.mass_ratio) ])
        arguments.extend([ "-e", str(self.size_eta_parameter) ])
        arguments.extend([ "-i", str(self.available_imf_model[self.imf]) ])        
        arguments.extend([ "-l", self.angular_momentum_signs ])
        arguments.extend([ "-f", str(self.fraction_of_max_rot_energy) ])
        arguments.extend([ "-a", str(self.OsipkovMerritt_anisotropy_radius_r0) ])
        arguments.extend([ "-c", str(self.cutoff_radius_halfmass_radius) ])
        arguments.extend([ "-r", str(self.physical_radius.value_in(units.parsec)) ])
        arguments.extend([ "-d", str(self.distance) ])
        arguments.extend([ "-E", str(self.orbital_energy) ])        
        arguments.extend([ "-L", str(self.orbital_angular_momentum) ])
        arguments.extend([ "-s", str(self.seed) ])                              
        return arguments
                
    def make_model(self):
    
        call=[self._exec]+self.arguments()
        
        print(call)
        
        mameclot=Popen(call, stdout=PIPE,stderr=PIPE,executable=os.path.join(self._bin_path,self._exec))
        
        (out,err)=mameclot.communicate()
 
        print(err)
        
        outsplit=out.decode().strip().split("\n")
        errsplit=err.decode().strip().split("\n")
        
        if self.mass_ratio==0:
          nline=errsplit[6].split()
          mline=errsplit[7].split()
          rline=errsplit[8].split()
          N1=int(nline[2])
          N2=0
          mscale=(float(mline[4])/float(mline[2])) | units.MSun
          rscale=(float(rline[4])/float(mline[2])) | units.parsec
        else:
          nline=errsplit[8].split()
          n2line=errsplit[22].split()
          mline=errsplit[9].split()
          rline=errsplit[10].split()
          N1=int(nline[2])
          N2=int(n2line[2])
          mscale=(float(mline[4])/float(mline[2])) | units.MSun
          rscale=(float(rline[4])/float(mline[2])) | units.parsec
        print(N1,N2)
        
        N=len( outsplit)
        
        parts=Particles(N)
        
        masses=numpy.zeros((N,))
        energy=numpy.zeros((N,))
        
        position=numpy.zeros((N,3))
        velocity=numpy.zeros((N,3))
        
        for i,line in enumerate(outsplit):
            l=line.split()
            masses[i]=float(l[0])
            position[i,0:3]=[float(l[1]),float(l[2]),float(l[3])]
            velocity[i,0:3]=[float(l[4]),float(l[5]),float(l[6])]
            energy[i]=float(l[7])  
            
        parts.mass=masses | nbody_system.mass
        parts.position=position | nbody_system.length
        parts.velocity=velocity | nbody_system.speed
        parts.specific_energy=energy| nbody_system.specific_energy
  
        parts.move_to_center()
  
        if self.convert_to_physical:
            print("mass scale:", mscale)
            print("length scale:", rscale)  

            convert_nbody=nbody_system.nbody_to_si(mscale,rscale)
            parts = ParticlesWithUnitsConverted(parts, convert_nbody.as_converter_from_si_to_generic())
            parts = parts.copy()
  
        self._all=parts
        self._cluster1=parts[:N1]
        self._cluster2=parts[N1:]


    @property
    def result(self):
        self.make_model()
        return self._all

    @property
    def result_split(self):
        self.make_model()
        return self._all,self._cluster1,self._cluster2

def new_mameclot_model(*args,**keyword_arguments):
    """

    MAMECLOT: MAke ME a CLuster Or Two
  
    Construct initial conditions of 1 (or 2) cluster(s) for N-body simulations 
  
    A 2 cluster system with orbit in the x-y plane is made if the mass ratio q>0 
    Five equilibrium models are available, a cut-off radius can be provided
    The velocities can be isotropic or radially anisotropic (a la Osipkov-Merritt)
    Optional angular momentum by aligning the angular momentum vectors along z
    System is scaled to N-body units with G = M = -4E = 1 (Heggie & Mathieu 1986)
    In the case of a two-body system E = mu*e_orb + E_1 + E_2 (mu = reduced mass) 
    The 2-body orbit is computed from input q, e_hat and l_hat (and eta if q<1):
    
    e_hat = e_orb/0.5<sigma^2> 
    l_hat = l_orb/<r_vir*sigma>
    eta => defines the relative radii: r2/r1 = q^eta, for example:
    eta=[-0.33/0/0.33/0.5] for equal [t_rh/r_h/rho_h/Sigma_h]

    :argument targetN: Number of particles
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument cluster_model: cluster model to generate ["Henon"]
    :argument mass_ratio: mass ratio, zero means only 1 cluster [0]
    :argument size_eta_parameter:eta: Relative cluster sizes: r2/r1 = q^eta [0.333]
    :argument imf: IMF ["single mass"]
    :argument angular_momentum_signs: Angular momentum in z-direction, -, + or x ["xx"]
    :argument fraction_of_max_rot_energy: Fraction of maximum rotational energy   [1],
    :argument OsipkovMerritt_anisotropy_radius_r0: Osipkov-Merritt anisotropy radius in units of r_0 [999]
    :argument cutoff_radius_halfmass_radius: Cut-off radius in units of r_h [20]
    :argument physical_radius: Physical scale [1| units.parsec]
    :argument distance: Distance between 2 clusters in N-body units [20 | nbody_system.length]
    :argument orbital_energy: Dimensionless orbital energy of two-cluster system [0]
    :argument orbital_angular_momentum: Dimensionless orbital angular momentum two-cluster system [4]
    :argument seed: random seed 123456
    """
    uc = mameclot(*args, **keyword_arguments)
    return uc.result

if __name__=="__main__":
    from matplotlib import pyplot
    
    clusters,cluster1,cluster2=mameclot(mass_ratio=0.25,convert_to_physical=False).result_split
    
    print(cluster1.total_mass())
    print(len(cluster2))
    
    pyplot.plot(cluster1.x.number,cluster1.y.number,"r.")
    pyplot.plot(cluster2.x.number,cluster2.y.number,"g.")
    pyplot.show()


Mameclot = mameclot
