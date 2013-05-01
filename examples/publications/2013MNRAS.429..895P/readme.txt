this directory contains the scripts for the runs of

"The formation of planets in circumbinary disks"
Pelupessy & Portegies Zwart 2012, MNRAS
http://adsabs.harvard.edu/doi/10.1093/mnras/sts461

it contains the following files:

./readme.txt - this readme

./disk_script/run.py - run script
./disk_script/binary.py - contains setup and mainloop of run
./disk_script/fast.py - bridge implementation
./disk_script/boxedfi.py - hydrodynamic wrapper of Fi
./disk_script/directsum.py - grav. skeleton code for disk planet interactions

./three_body_script/bin_planet.py - three body run script
./three_body_script/python_interface.py - remote function implementarion
./three_body_script/kepler_16_three_body.py - job farming script

Note that currently these are provided for reference. Some work may be needed
to adapt these to your local AMUSE setup. 

history:

V0.1 - 18/12/2012 - initial release
