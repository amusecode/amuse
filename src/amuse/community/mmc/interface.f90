
MODULE MMC

CONTAINS

FUNCTION nonstandard_init()
  INTEGER :: nonstandard_init
  INTEGER :: init_sequence
  INTEGER :: res
  ! read initial parameters
  PRINT*,'calling input'
  res = init_sequence()
  nonstandard_init = res
  PRINT*,'init done'
END FUNCTION

FUNCTION set_mmc_data_directory(data_directory) 
  INTEGER :: set_mmc_data_directory
  CHARACTER(len=200) :: data_directory
  CALL amuse_set_mmc_data_directory(data_directory)
  set_mmc_data_directory = 0
END FUNCTION

FUNCTION get_time(time_)
  INTEGER :: get_time
  COMMON /SYSTEM/ TIME
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: time_
  time_ = time
  get_time = 0
END FUNCTION

FUNCTION set_time(time_)
  INTEGER :: set_time
  COMMON /SYSTEM/ TIME
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: time_
  time = time_
  set_time = 0
END FUNCTION


FUNCTION get_kinetic_energy(Ek)
  IMPLICIT NONE
  INTEGER :: res
  INTEGER :: total_kinetic_energy
  INTEGER :: get_kinetic_energy
  DOUBLE PRECISION :: Ek

  res = total_kinetic_energy(Ek)
  get_kinetic_energy = 0
END FUNCTION


FUNCTION new_particle( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: new_particle
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  new_particle=0
END FUNCTION

FUNCTION delete_particle( index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: delete_particle
  INTEGER :: index_of_the_particle
  delete_particle=0
END FUNCTION

FUNCTION set_state( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: set_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  set_state=0
END FUNCTION

FUNCTION get_state( index_of_the_particle, mass_, r_, vr_, vt_,x1,x2,x3)
  INTEGER :: get_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass_, r_, vr_, vt_, x1,x2,x3
  
  COMMON /BODY/ BODY(1000),VR(1000),VT(1000), &
                U(1000),XESCP(1000),XESCT(1000),VRR(1000),R(1000)

  COMMON /POSVEL/ X(1000,3), XDOT(1000,3)
  REAL*8  BODY, VR, VT, R
  REAL*8  X,XDOT

  mass_ = BODY(index_of_the_particle)
  r_ = R(index_of_the_particle)
  vr_ = VR(index_of_the_particle)
  vt_ = VT(index_of_the_particle)
  x1 = X(index_of_the_particle,1)
  y1 = X(index_of_the_particle,2)
  z1 = X(index_of_the_particle,3)

  get_state=0
END FUNCTION

FUNCTION get_number_of_particles( VALUE)
  IMPLICIT NONE
  INTEGER :: get_number_of_particles
  INTEGER :: VALUE
  get_number_of_particles=0
END FUNCTION

FUNCTION internal__redirect_outputs( stdoutfile, stderrfile)
  IMPLICIT NONE
  INTEGER :: internal__redirect_outputs
  CHARACTER(LEN=*) :: stdoutfile, stderrfile
  internal__redirect_outputs=0
END FUNCTION

FUNCTION run( )
  IMPLICIT NONE
  INTEGER :: run
  run=0
END FUNCTION

FUNCTION commit_parameters( )
  IMPLICIT NONE
  INTEGER commit_parameters
  commit_parameters=0
END FUNCTION

END MODULE
FUNCTION set_irun(tmp_)
  COMMON /IPARAM/ irun
  INTEGER irun
  INTEGER tmp_
  INTEGER set_irun
  irun = tmp_
  set_irun = 0
END FUNCTION set_irun

FUNCTION set_nt(tmp_)
  COMMON /IPARAM/ nt
  INTEGER nt
  INTEGER tmp_
  INTEGER set_nt
  nt = tmp_
  set_nt = 0
END FUNCTION set_nt

FUNCTION set_istart(tmp_)
  COMMON /IPARAM/ istart
  INTEGER istart
  INTEGER tmp_
  INTEGER set_istart
  istart = tmp_
  set_istart = 0
END FUNCTION set_istart

FUNCTION set_ncor(tmp_)
  COMMON /IPARAM/ ncor
  INTEGER ncor
  INTEGER tmp_
  INTEGER set_ncor
  ncor = tmp_
  set_ncor = 0
END FUNCTION set_ncor

FUNCTION set_nmin(tmp_)
  COMMON /IPARAM/ nmin
  INTEGER nmin
  INTEGER tmp_
  INTEGER set_nmin
  nmin = tmp_
  set_nmin = 0
END FUNCTION set_nmin

FUNCTION set_nz_naught(tmp_)
  COMMON /IPARAM/ nz_naught
  INTEGER nz_naught
  INTEGER tmp_
  INTEGER set_nz_naught
  nz_naught = tmp_
  set_nz_naught = 0
END FUNCTION set_nz_naught

FUNCTION set_nzonc(tmp_)
  COMMON /IPARAM/ nzonc
  INTEGER nzonc
  INTEGER tmp_
  INTEGER set_nzonc
  nzonc = tmp_
  set_nzonc = 0
END FUNCTION set_nzonc

FUNCTION set_nminzo(tmp_)
  COMMON /IPARAM/ nminzo
  INTEGER nminzo
  INTEGER tmp_
  INTEGER set_nminzo
  nminzo = tmp_
  set_nminzo = 0
END FUNCTION set_nminzo

FUNCTION set_ntwo(tmp_)
  COMMON /IPARAM/ ntwo
  INTEGER ntwo
  INTEGER tmp_
  INTEGER set_ntwo
  ntwo = tmp_
  set_ntwo = 0
END FUNCTION set_ntwo

FUNCTION set_imodel(tmp_)
  COMMON /IPARAM/ imodel
  INTEGER imodel
  INTEGER tmp_
  INTEGER set_imodel
  imodel = tmp_
  set_imodel = 0
END FUNCTION set_imodel

FUNCTION set_iprint(tmp_)
  COMMON /IPARAM/ iprint
  INTEGER iprint
  INTEGER tmp_
  INTEGER set_iprint
  iprint = tmp_
  set_iprint = 0
END FUNCTION set_iprint

FUNCTION set_ib3f(tmp_)
  COMMON /IPARAM/ ib3f
  INTEGER ib3f
  INTEGER tmp_
  INTEGER set_ib3f
  ib3f = tmp_
  set_ib3f = 0
END FUNCTION set_ib3f

FUNCTION set_iexch(tmp_)
  COMMON /IPARAM/ iexch
  INTEGER iexch
  INTEGER tmp_
  INTEGER set_iexch
  iexch = tmp_
  set_iexch = 0
END FUNCTION set_iexch

FUNCTION set_tcrit(tmp_)
  COMMON /IPARAM/ tcrit
  DOUBLE PRECISION tcrit
  DOUBLE PRECISION tmp_
  INTEGER set_tcrit
  tcrit = tmp_
  set_tcrit = 0
END FUNCTION set_tcrit

FUNCTION set_tcomp(tmp_)
  COMMON /IPARAM/ tcomp
  DOUBLE PRECISION tcomp
  DOUBLE PRECISION tmp_
  INTEGER set_tcomp
  tcomp = tmp_
  set_tcomp = 0
END FUNCTION set_tcomp

FUNCTION set_qe(tmp_)
  COMMON /IPARAM/ qe
  DOUBLE PRECISION qe
  DOUBLE PRECISION tmp_
  INTEGER set_qe
  qe = tmp_
  set_qe = 0
END FUNCTION set_qe

FUNCTION set_alphal(tmp_)
  COMMON /IPARAM/ alphal
  INTEGER alphal
  INTEGER tmp_
  INTEGER set_alphal
  alphal = tmp_
  set_alphal = 0
END FUNCTION set_alphal

FUNCTION set_alphah(tmp_)
  COMMON /IPARAM/ alphah
  INTEGER alphah
  INTEGER tmp_
  INTEGER set_alphah
  alphah = tmp_
  set_alphah = 0
END FUNCTION set_alphah

FUNCTION set_brakem(tmp_)
  COMMON /IPARAM/ brakem
  DOUBLE PRECISION brakem
  DOUBLE PRECISION tmp_
  INTEGER set_brakem
  brakem = tmp_
  set_brakem = 0
END FUNCTION set_brakem

FUNCTION set_body_1(tmp_)
  COMMON /IPARAM/ body_1
  DOUBLE PRECISION body_1
  DOUBLE PRECISION tmp_
  INTEGER set_body_1
  body_1 = tmp_
  set_body_1 = 0
END FUNCTION set_body_1

FUNCTION set_body_n(tmp_)
  COMMON /IPARAM/ body_n
  DOUBLE PRECISION body_n
  DOUBLE PRECISION tmp_
  INTEGER set_body_n
  body_n = tmp_
  set_body_n = 0
END FUNCTION set_body_n

FUNCTION set_body_fracb(tmp_)
  COMMON /IPARAM/ body_fracb
  DOUBLE PRECISION body_fracb
  DOUBLE PRECISION tmp_
  INTEGER set_body_fracb
  body_fracb = tmp_
  set_body_fracb = 0
END FUNCTION set_body_fracb

FUNCTION set_amin(tmp_)
  COMMON /IPARAM/ amin
  DOUBLE PRECISION amin
  DOUBLE PRECISION tmp_
  INTEGER set_amin
  amin = tmp_
  set_amin = 0
END FUNCTION set_amin

FUNCTION set_amax(tmp_)
  COMMON /IPARAM/ amax
  DOUBLE PRECISION amax
  DOUBLE PRECISION tmp_
  INTEGER set_amax
  amax = tmp_
  set_amax = 0
END FUNCTION set_amax

FUNCTION set_qvir(tmp_)
  COMMON /IPARAM/ qvir
  DOUBLE PRECISION qvir
  DOUBLE PRECISION tmp_
  INTEGER set_qvir
  qvir = tmp_
  set_qvir = 0
END FUNCTION set_qvir

FUNCTION set_rbar(tmp_)
  COMMON /IPARAM/ rbar
  DOUBLE PRECISION rbar
  DOUBLE PRECISION tmp_
  INTEGER set_rbar
  rbar = tmp_
  set_rbar = 0
END FUNCTION set_rbar

FUNCTION set_zmbar(tmp_)
  COMMON /IPARAM/ zmbar
  DOUBLE PRECISION zmbar
  DOUBLE PRECISION tmp_
  INTEGER set_zmbar
  zmbar = tmp_
  set_zmbar = 0
END FUNCTION set_zmbar

FUNCTION set_w_naught(tmp_)
  COMMON /IPARAM/ w_naught
  DOUBLE PRECISION w_naught
  DOUBLE PRECISION tmp_
  INTEGER set_w_naught
  w_naught = tmp_
  set_w_naught = 0
END FUNCTION set_w_naught

FUNCTION set_bmin(tmp_)
  COMMON /IPARAM/ bmin
  DOUBLE PRECISION bmin
  DOUBLE PRECISION tmp_
  INTEGER set_bmin
  bmin = tmp_
  set_bmin = 0
END FUNCTION set_bmin

FUNCTION set_bmax(tmp_)
  COMMON /IPARAM/ bmax
  DOUBLE PRECISION bmax
  DOUBLE PRECISION tmp_
  INTEGER set_bmax
  bmax = tmp_
  set_bmax = 0
END FUNCTION set_bmax

FUNCTION set_tau_naught(tmp_)
  COMMON /IPARAM/ tau_naught
  DOUBLE PRECISION tau_naught
  DOUBLE PRECISION tmp_
  INTEGER set_tau_naught
  tau_naught = tmp_
  set_tau_naught = 0
END FUNCTION set_tau_naught

FUNCTION set_gamma(tmp_)
  COMMON /IPARAM/ gamma
  DOUBLE PRECISION gamma
  DOUBLE PRECISION tmp_
  INTEGER set_gamma
  gamma = tmp_
  set_gamma = 0
END FUNCTION set_gamma

FUNCTION set_xtid(tmp_)
  COMMON /IPARAM/ xtid
  DOUBLE PRECISION xtid
  DOUBLE PRECISION tmp_
  INTEGER set_xtid
  xtid = tmp_
  set_xtid = 0
END FUNCTION set_xtid

FUNCTION set_rplum(tmp_)
  COMMON /IPARAM/ rplum
  DOUBLE PRECISION rplum
  DOUBLE PRECISION tmp_
  INTEGER set_rplum
  rplum = tmp_
  set_rplum = 0
END FUNCTION set_rplum

FUNCTION set_dttp(tmp_)
  COMMON /IPARAM/ dttp
  DOUBLE PRECISION dttp
  DOUBLE PRECISION tmp_
  INTEGER set_dttp
  dttp = tmp_
  set_dttp = 0
END FUNCTION set_dttp

FUNCTION set_dtte(tmp_)
  COMMON /IPARAM/ dtte
  DOUBLE PRECISION dtte
  DOUBLE PRECISION tmp_
  INTEGER set_dtte
  dtte = tmp_
  set_dtte = 0
END FUNCTION set_dtte

FUNCTION set_dtte_naught(tmp_)
  COMMON /IPARAM/ dtte_naught
  DOUBLE PRECISION dtte_naught
  DOUBLE PRECISION tmp_
  INTEGER set_dtte_naught
  dtte_naught = tmp_
  set_dtte_naught = 0
END FUNCTION set_dtte_naught

FUNCTION set_tcrevo(tmp_)
  COMMON /IPARAM/ tcrevo
  DOUBLE PRECISION tcrevo
  DOUBLE PRECISION tmp_
  INTEGER set_tcrevo
  tcrevo = tmp_
  set_tcrevo = 0
END FUNCTION set_tcrevo

FUNCTION set_xtau(tmp_)
  COMMON /IPARAM/ xtau
  DOUBLE PRECISION xtau
  DOUBLE PRECISION tmp_
  INTEGER set_xtau
  xtau = tmp_
  set_xtau = 0
END FUNCTION set_xtau

FUNCTION set_ytau(tmp_)
  COMMON /IPARAM/ ytau
  DOUBLE PRECISION ytau
  DOUBLE PRECISION tmp_
  INTEGER set_ytau
  ytau = tmp_
  set_ytau = 0
END FUNCTION set_ytau

FUNCTION set_ybmin(tmp_)
  COMMON /IPARAM/ ybmin
  DOUBLE PRECISION ybmin
  DOUBLE PRECISION tmp_
  INTEGER set_ybmin
  ybmin = tmp_
  set_ybmin = 0
END FUNCTION set_ybmin

FUNCTION set_zini(tmp_)
  COMMON /IPARAM/ zini
  DOUBLE PRECISION zini
  DOUBLE PRECISION tmp_
  INTEGER set_zini
  zini = tmp_
  set_zini = 0
END FUNCTION set_zini

FUNCTION set_ikroupa(tmp_)
  COMMON /IPARAM/ ikroupa
  INTEGER ikroupa
  INTEGER tmp_
  INTEGER set_ikroupa
  ikroupa = tmp_
  set_ikroupa = 0
END FUNCTION set_ikroupa

FUNCTION set_iflagns(tmp_)
  COMMON /IPARAM/ iflagns
  INTEGER iflagns
  INTEGER tmp_
  INTEGER set_iflagns
  iflagns = tmp_
  set_iflagns = 0
END FUNCTION set_iflagns

FUNCTION set_iflagbh(tmp_)
  COMMON /IPARAM/ iflagbh
  INTEGER iflagbh
  INTEGER tmp_
  INTEGER set_iflagbh
  iflagbh = tmp_
  set_iflagbh = 0
END FUNCTION set_iflagbh

FUNCTION set_nitesc(tmp_)
  COMMON /IPARAM/ nitesc
  INTEGER nitesc
  INTEGER tmp_
  INTEGER set_nitesc
  nitesc = tmp_
  set_nitesc = 0
END FUNCTION set_nitesc

