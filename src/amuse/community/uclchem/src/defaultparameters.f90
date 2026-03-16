!!This file is used to autogenerate the docs. So please ignore the mess!
!!Double !! lines do not show up in docs, single ones do. 
!!If you add a parameter, please take the time to add a useful descriptor comment on the same line
!!and then re-run utils/generate_param_docs.py to update the docs.
!!note the resuting md file needs manually adding to the website.


!---
!id: parameters
!title: Model Parameters
!---
!UCLCHEM will default to these values unless they are overridden by user. Users can override these by adding the variable name as written here in the param_dict argument of any UCLCHEM model function. param_dict is not case sensitive.
!
!## Physical Variables
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
initialTemp=10.0 !Initial gas temperature in Kelvin for all gas parcels in model.
initialDens=1.00d2 !Initial gas density in H nuclei per cm$^{-3}$ for all gas parcels in model.
finalDens=1.00d5 !Final gas density achieved via freefall.
currentTime=0.0 !Time at start of model in years.
finalTime=5.0d6 !Time to stop model in years, if not using `endAtFinalDensity` below.
radfield=1.0 !Interstellar radiation field in Habing
zeta=1.0 !Cosmic ray ionisation rate as multiple of $1.3 10^{-17} s^{-1}$
rout=0.05 !Outer radius of cloud being modelled in pc.
rin=0.0 !Minimum radial distance from cloud centre to consider.
baseAv=2.0 !Extinction at cloud edge, Av of a parcel at rout.
points=1 !Number of gas parcels equally spaced between rin to rout to consider
!
!## Behavioural Controls
!*The following parameters generally turn on or off features of the model. If a parameter is set to `True`, then it is turned on. If it is set to `False`, then it is turned off.*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
freezeFactor=1.0 !Modify freeze out rate of gas parcels by this factor.
endAtFinalDensity=.False. !Choose to end model at final density, otherwise end at final time.
freefall=.False. !Controls whether models density increaes following freefall equation.
freefallFactor=1.0 !Modify freefall rate by factor, usually to slow it.
desorb=.True. !Toggles all non-thermal desoprtion processes on or off.
h2desorb=.True. !Individually toggle non-thermal desorption due to H2 formation.
crdesorb=.True. !Individually toggle non-thermal desorption due to cosmic rays.
uvdesorb=.True. !Individually toggle non-thermal desorption due to uv photons.
thermdesorb=.True. !Toggle continuous thermal desorption.
instantSublimation=.False. !Toggle instantaneous sublimation of the ices at t=0
cosmicRayAttenuation=.False. !Use column density to attenuate cosmic ray ionisation rate following [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
ionModel='L' !L/H model for cosmic ray attenuation [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
improvedH2CRPDissociation=.False. !Use H2 CRP dissociation rate from [Padovani et al. 2018b](https://arxiv.org/abs/1809.04168).
!
!## Input and Output
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
outputFile="output/full.dat" !File to write full output of UCLCHEM. This includes physical parameter values and all abundances at every time step.
columnFile="output/column.dat" !File to write specific species abundances, see outSpecies.
writeStep=1 !Writing to columnFile only happens every writeStep timesteps.
!|abundSaveFile |None| File to store final abundances at the end of the model so future models can use them as the initial abundances. If not provided, no file will be produced.
!|abundLoadFile |None| File from which to load initial abundances for the model, created through `abundSaveFile`. If not provided, the model starts from elemental gas.
!|outSpecies|None| A space separated list of species to output to columnFile. Supplied as a separate list argument to most python functions, see python API docs.
!
!## Initial Abundances
!*Unless otherwise specified, we take all abundances from Jenkins et al. 2009, using the heavily depleted case from Table 4.*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
metallicity=1.0 !Scale the abundances of all elements heavier than He by this factor.
ion=2 !Sets how much elemental C is initially atomic (0= all atomic/1=50:50/2=fully ionized).
fh=0.5 !Total elemental abundance of H is always 1 by definition because abundances are relative to number of H nuclei. Use fh to set how much to initially put in atomic H, the rest goes to H2.
fhe = 0.1 !Total elemental abundance of He.
fc=1.77d-04 !Total elemental abundance of C.
fo  = 3.34d-04 !Total elemental abundance of O.
fn  = 6.18d-05 !Total elemental abundance of N.
fs  = 3.51d-6 !Total elemental abundance of S.
fmg = 2.256d-06 !Total elemental abundance of Mg.
fsi = 1.78d-06 !Total elemental abundance of Si.
fcl = 3.39d-08 !Total elemental abundance of Cl.
fp=7.78d-08 !Total elemental abundance of P.
ffe=2.01d-7!Total elemental abundance of Fe.
ff = 3.6d-08 !fp depleted 1/100 of solar from Asplund 2009.
fd=0.0 ! The following elements are not typically used. We do not recommend any particular value.
fli=0.0 !Total elemental abundance of Li.
fna=0.0 !Total elemental abundance of Na.
fpah=0.0 !Total initial abundance of PAHs.
f15n=0.0 !Total initial abundance of 15N.
f13c=0.0 !Total initial abundance of 13C.
f18O=0.0 !Total initial abundance of 18O.
!!
!! We used to use Asplund et al. 2009,kept here for reference
!! !initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
!! !note fh is fraction of H initially in H atoms. Total H is always 1.
!! !fh=0.5;fhe = 0.1;fc  = 2.6d-04;fo  = 4.6d-04;fn  = 6.1d-05
!! fs  = 1.318d-05;fmg = 3.981d-05;fsi = 1.0d-07;fcl = 3.162d-07;
!! fp=2.57d-09 ; ff = 3.6d-08 !fp depleted 1/100 of solar
!
!## Integration Controls
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
reltol=1d-8 !Relative tolerance for integration, see [integration docs](/docs/trouble-integration) for advice.
abstol_factor=1.0d-14 !Absolute tolerance for integration is calculated by multiplying species abundance by this factor.
abstol_min=1.0d-25 !Minimum value absolute tolerances can take.
MXSTEP=10000 !Maximum steps allowed in integration before warning is thrown.
!
!## Here be Dragons
!*These are not recommended to be changed unless you know what you are doing*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
ebmaxh2=1.21d3 ! Maximum binding energy of species desorbed by H2 formation.
ebmaxcr=1.21d3 ! Maximum binding energy of species desorbed by cosmic ray ionisation.
ebmaxuvcr=1.0d4 ! Maximum binding energy of species desorbed by UV photons.
epsilon=0.01 !Number of molecules desorbed per H2 formation.
uv_yield=0.1 !Number of molecules desorbed per UV photon.
phi=1.0d5 !Number of molecules desorbed per cosmic ray ionisation.
uvcreff=1.0d-3 !Ratio of CR induced UV photons to ISRF UV photons.
omega=0.5 !Dust grain albedo.
!|alpha|{1:0.0,2:0.0}| Set alpha coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how alpha is used for each reaction type.|
!|beta|{1:0.0,2:0.0}| Set beta coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how beta is used for each reaction type.|
!|gama|{1:0.0,2:0.0}| Set gama coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how gama is used for each reaction type.|