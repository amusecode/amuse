#!/bin/bash
#MSUB -r tstranexg        # Nom du job                
#MSUB -n 4           # reservation du nombre de processus
#MSUB -N 1            # Reservation du nombre de noeuds
#MSUB -T 3600         # Limite de temps elapsed du job 
#MSUB -q prod          # Reservation sur la queue
#MSUB -o stdout        # Sortie standard
#MSUB -e stderr        # Sortie d'erreur       
#MSUB -p gen2191

set -x

export DATE=`date +%F_%Hh%M`

cd $SCRATCHDIR/tstranex/ics/boxlen12p5_n256

mpirun /home/cont003/teyssier/tim/mpgrafic/mpgrafic/bin/mpgrafic < $BRIDGE_MSUB_PWD/input_12.5.input > run$DATE.log

