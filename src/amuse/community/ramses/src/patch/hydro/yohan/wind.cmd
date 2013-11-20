#! /bin/csh -f
#BSUB -J disk
#BSUB -o disk.%J
#BSUB -n 8
#BSUB -c 00:30
#BSUB -r

cd /scratch/scratchdir/teyssier/yohan

set nrestart = `ls | grep backup_| wc -l`  
#set nrestart = 0

cat > start$nrestart.nml <<EOF

&RUN_PARAMS 
hydro=.true.
poisson=.true.
pic=.true.
ncontrol=1
nsubcycle=15*1
nremap=10
nrestart=$nrestart
verbose=.false.
/

&AMR_PARAMS
levelmin=6
levelmax=10
ngridmax=200000
npartmax=100000
boxlen=30.0
/

&INIT_PARAMS
filetype='ascii'
initfile(1)=''
/

&OUTPUT_PARAMS
fbackup=30
foutput=30
noutput=2
tout=0.0,21.0
/

&BOUNDARY_PARAMS
nboundary = 6
bound_type= 2, 2, 2, 2, 2, 2
ibound_min=-1,+1,-1,-1,-1,-1
ibound_max=-1,+1,+1,+1,+1,+1
jbound_min= 0, 0,-1,+1,-1,-1
jbound_max= 0, 0,-1,+1,+1,+1
kbound_min= 0, 0, 0, 0,-1,+1
kbound_max= 0, 0, 0, 0,-1,+1
/

&POISSON_PARAMS
gravity_type=-1
gravity_params=35.,10.,0.1,0.15
epsilon=1d-4
/

&HYDRO_PARAMS
pressure_fix=.true.
gamma=1.666666667
courant_factor=0.8
slope_type=1
scheme='muscl'
riemann='acoustic'
/

&PHYSICS_PARAMS
cooling=.true.
t_star=8.0
n_star=0.1
T2_star=1d4
metal=.true.
eta_sn=0.1
yield=0.1
/

&REFINE_PARAMS 
interpol_var=1
interpol_type=0
mass_sph=5d-5
m_refine=10*20.0
jeans_refine=20*10.0
/

EOF

prun -n8 $HOME/ramses/bin/yohan3d start$nrestart.nml > run$nrestart.log

ls -als


