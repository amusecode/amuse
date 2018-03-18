import numpy
cimport numpy

cimport mpi4py.MPI
cdef extern from "mpi.h":
    pass
cdef extern from "amuse_mpi.h":
    int c_set_comm_world "set_comm_world" (mpi4py.MPI.MPI_Comm world)

def set_comm_world(mpi4py.MPI.Comm comm not None):
    return c_set_comm_world(comm.ob_mpi)

cdef extern:
  
  
  void c_get_mass_loss_wind "aaa_get_mass_loss_wind" (int, double, double, double, double, double *);
  
  
  void c_get_gyration_radius "aaa_get_gyration_radius" (int, double, double, double, double, double, double, double, double *);
  
  
  void c_initialize "aaa_initialize" (double, double, double, double, double, double, int, int, int, int, int, double, int, double, double, double, int *);
  
  
  void c_get_time_step "aaa_get_time_step" (int, double, double, double, double, double, double *);
  
  
  void c_evolve0 "aaa_evolve0" (int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
  
  
def get_mass_loss_wind(kw, lum, r, mt, mc, mlout):
  
  cdef double output_mlout
  c_get_mass_loss_wind(kw, lum, r, mt, mc, &output_mlout);
  mlout.value = output_mlout


def get_gyration_radius(kw, mass, mt, r, lum, epoch, tm, tphys, rg):
  
  cdef double output_rg
  c_get_gyration_radius(kw, mass, mt, r, lum, epoch, tm, tphys, &output_rg);
  rg.value = output_rg


def initialize(z_in, neta_in, bwind_in, hewind_in, sigma1_in, sigma2_in, ifflag_in, wdflag_in, bhflag_in, nsflag_in, piflag_in, mxns_in, idum_in, pts1_in, pts2_in, pts3_in, status):
  
  cdef int output_status
  c_initialize(z_in, neta_in, bwind_in, hewind_in, sigma1_in, sigma2_in, ifflag_in, wdflag_in, bhflag_in, nsflag_in, piflag_in, mxns_in, idum_in, pts1_in, pts2_in, pts3_in, &output_status);
  status.value = output_status


def get_time_step(kw, mass, age, mt, tm, epoch, dt):
  
  cdef double output_dt
  c_get_time_step(kw, mass, age, mt, tm, epoch, &output_dt);
  dt.value = output_dt


def evolve0(kw, mass, mt, r, lum, mc, rc, menv, renv, ospin, epoch, tm, tphys, tphysf):
  
  cdef int inout_kw = kw.value
  cdef double inout_mass = mass.value
  cdef double inout_mt = mt.value
  cdef double inout_r = r.value
  cdef double inout_lum = lum.value
  cdef double inout_mc = mc.value
  cdef double inout_rc = rc.value
  cdef double inout_menv = menv.value
  cdef double inout_renv = renv.value
  cdef double inout_ospin = ospin.value
  cdef double inout_epoch = epoch.value
  cdef double inout_tm = tm.value
  cdef double inout_tphys = tphys.value
  cdef double inout_tphysf = tphysf.value
  c_evolve0(&inout_kw, &inout_mass, &inout_mt, &inout_r, &inout_lum, &inout_mc, &inout_rc, &inout_menv, &inout_renv, &inout_ospin, &inout_epoch, &inout_tm, &inout_tphys, &inout_tphysf);
  kw.value = inout_kw
  mass.value = inout_mass
  mt.value = inout_mt
  r.value = inout_r
  lum.value = inout_lum
  mc.value = inout_mc
  rc.value = inout_rc
  menv.value = inout_menv
  renv.value = inout_renv
  ospin.value = inout_ospin
  epoch.value = inout_epoch
  tm.value = inout_tm
  tphys.value = inout_tphys
  tphysf.value = inout_tphysf

