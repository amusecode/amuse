import numpy as np

class PolarSupport(object):
    def __init__(self):
        pass

    def homogeneous_sphere_N(self, N):
        #6/pi is the volumetric ratio of a cube
        #and a fitting sphere
        size_is_not_N = True
        while (size_is_not_N):
            M = 2.0*(np.random.random([N*6/np.pi, 3]) - 0.5 * np.ones([N*6/np.pi, 3]))
            norms = (M[::,0]*M[::,0]+\
                     M[::,1]*M[::,1]+\
                     M[::,2]*M[::,2])**0.5
    
            selection = np.where(norms<1.0)
    
            if len(selection[0])>=N:
                size_is_not_N = False
      
        Msphere = M[selection,0:3][0]
        inv_norms_fit = 1.0/norms[selection]
  
        Shell = np.array(np.diag(inv_norms_fit)*np.matrix(Msphere))

        x = Shell[0:N,0]
        y = Shell[0:N,1]
        z = Shell[0:N,2]
        return x, y, z

    def phase_to_polar(self, x, y, z, vx, vy, vz):
        r = (x**2 + y**2 + z**2)**0.5
        vr = (vx*x + vy*y + vz*z) / (x**2 + y**2 + z**2)**0.5
        vt = ((vx**2 + vy**2 + vz**2) - vr**2)**0.5
        return r, vr, vt

    def position_to_polar(self, x, y, z):
        return (x**2 + y**2 + z**2)**0.5

    def phase_to_cartesian(self, ra, vr, vt):
        ex, ey, ez = self.homogeneous_sphere_N(len(ra))
        
        x = ex * ra
        y = ey * ra
        z = ez * ra
        
        n1 = 1.0/np.sqrt(x**2+y**2)
        n2 = 1.0/np.sqrt(x**2+z**2)
        t1x = -y*n1
        t1y = x*n1
        t1z = np.zeros(len(x))
        t2x = -z*n2
        t2y = np.zeros(len(x))
        t2z = x*n2
        
        theta = np.random.random(len(ra))
        s = np.cos(theta)
        t = np.sin(theta)
        vx = vt * (s*t1x + t*t2x)
        vy = vt * (s*t1y + t*t2y)
        vz = vt * (s*t1z + t*t2z)
        return x, y, z, vx, vy, vz, ex, ey, ez

    def position_to_cartesian(self, ra):
        return ra * self.homogeneous_sphere_N(len(ra))


