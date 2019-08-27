"""
Stefan Umbreit, python version

.. [#] Casertano, S., Hut, P., *The Astrophysical Journal*, **298**, 80-94 (1985)
"""

import numpy as N

def density_estimators(n_points, r, mass, dims=3):
  """
  Calculate the density estimators that are used for the calculation of the
  core quantities.

  n_points - number of points to average over
  r, mass  - radial positions and masses of the stars 
             (both arrays must be sorted radially)
  """
  # calculate the total mass of all stars
  m_tot= N.sum(mass)

  nave= N.argmin(N.abs(0.5*m_tot-N.cumsum(mass)))
  rhoj= N.zeros((nave,2), N.float64)

  for i in range(nave):
    jmin= max(i-n_points/2, 0)
    jmax= min(jmin + n_points, nave-1)
    mrho= 0.
    #this is equivalent to their J-1 factor for the case of equal masses,
    #and seems like a good generalization for unequal masses 
    mrho= N.sum(mass[jmin+1:jmax])

    if dims==3:
      Vrj = 4.0/3.0 * N.pi * (r[jmax]**3 - r[jmin]**3)
    elif dims==2:
      Vrj = N.pi*(r[jmax]**2- r[jmin]**2)

    rhoj[i,1]= mrho/Vrj
    rhoj[i,0]= r[i]

  return(rhoj)

def core_radius(rhoj):
  """
  Calculates the core size using estimators of the density.

  rhoj - Density estimators (use density_estimator function)
  """

  return N.sum(rhoj[:,0]*rhoj[:,1])/N.sum(rhoj[:,1])


