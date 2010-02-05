import numpy
cimport numpy
cimport cython

cdef extern from "BHTC/nbody_particle.h":
    ctypedef struct c_vec "vec":
        double index "operator[]" (int x)
        pass
    c_vec *new_vec "new vec" (double x, double y, double z)
    void del_vec "delete" (c_vec *instance)
    
    ctypedef struct c_nbody_particle "nbody_particle":
        void set_pos(c_vec new_pos)   
        c_vec get_pos()   
        void set_vel(c_vec new_pos)   
        c_vec get_vel()    
        void set_acc_gravity(c_vec new_pos) 
        c_vec get_acc_gravity()      
        void set_phi_gravity(double set_phi_gravity) 
        void set_mass(double mass) 
        double get_mass()
        double get_phi_gravity()  
        int get_index()
        void set_index(int x)
        void dump()
    c_nbody_particle *new_nbody_particle "new nbody_particle" ()
    void del_nbody_particle "delete" (c_nbody_particle *instance)
    
    ctypedef struct c_nbody_system "nbody_system":
        int n
        void set_particle_pointer(c_nbody_particle * p)
        c_nbody_particle * get_particle_pointer()
        c_nbody_particle * new_particle_pointer(unsigned long n)
        void set_nsize(int n)
        void calculate_gravity()
        void integrate(double dt)
        void dump()
    c_nbody_system *new_nbody_system "new nbody_system" ()
    void del_nbody_system "delete" (c_nbody_system *instance)
    
      

DTYPE = numpy.float64
ctypedef numpy.float64_t DTYPE_t
  
        
cdef class NBodySystem:
    cdef c_nbody_system *thisptr
    def __cinit__(self):
        self.thisptr = new_nbody_system()
    def __dealloc__(self):
        cdef c_nbody_particle *particles=self.thisptr.get_particle_pointer()
        del_nbody_particle(particles)
        del_nbody_system(self.thisptr)
    
    @cython.boundscheck(False)
    def set_state(self, numpy.ndarray[DTYPE_t, ndim=2] position, numpy.ndarray[DTYPE_t, ndim=2] velocity, numpy.ndarray[DTYPE_t, ndim=2] mass):
        assert position.dtype == DTYPE and velocity.dtype == DTYPE and mass.dtype == DTYPE

        cdef unsigned int n = position.shape[0]
        cdef c_vec * vec = NULL
        cdef c_nbody_particle *particles=self.thisptr.new_particle_pointer(n)
        
        for i from 0 <= i < n:
            particles[i].set_index(i)
            
            vec = new_vec(position[i,0],position[i,1],position[i,2])
            particles[i].set_pos(vec[0])
            del_vec(vec)
            
            vec = new_vec(velocity[i,0],velocity[i,1],velocity[i,2])
            particles[i].set_vel(vec[0])
            del_vec(vec)
            
            particles[i].set_mass(mass[i,0])
            
        self.thisptr.set_particle_pointer(particles)
        self.thisptr.n = n
        self.thisptr.set_nsize(n)
    def get_state(self, ):
        cdef unsigned int n = self.thisptr.n
        cdef numpy.ndarray[DTYPE_t, ndim=2] position = numpy.zeros([n, 3], dtype=DTYPE)
        cdef numpy.ndarray[DTYPE_t, ndim=2] velocity = numpy.zeros([n, 3], dtype=DTYPE)
        cdef numpy.ndarray[DTYPE_t, ndim=2] mass = numpy.zeros([n, 1], dtype=DTYPE)
        cdef c_vec vec
        cdef c_nbody_particle *particles=self.thisptr.get_particle_pointer()
        for i from 0 <= i < n:
            particles[i].set_index(i)
            
            vec = particles[i].get_pos()
            position[i,0] = vec.index(0)
            position[i,1] = vec.index(1)
            position[i,2] = vec.index(2)
            
            
            vec = particles[i].get_vel()
            velocity[i,0] = vec.index(0)
            velocity[i,1] = vec.index(1)
            velocity[i,2] = vec.index(2)
            
            mass[i,0] = particles[i].get_mass()

        return position, velocity, mass
    property n:
        def __get__(self): return self.n
    def dump(self):
        self.thisptr.dump()
    def calculate_gravity(self):
        self.thisptr.calculate_gravity()
    def integrate(self, double dt):
        self.thisptr.integrate(dt)
         
       
            


   