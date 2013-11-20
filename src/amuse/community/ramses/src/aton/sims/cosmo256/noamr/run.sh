#!/bin/bash
#MSUB -r tstranex5          # Nom du job                
#MSUB -n 16                  # Reservation de 4 processus au total
#MSUB -N 8                  # Les processus sont répartis sur 2 noeuds
#MSUB -T 36000                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout          # Sortie standard
#MSUB -e stderr          # Sortie d'erreur       
#MSUB -p gen6667      # Allocation
#MSUB -q gpu                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x

export DATE=`date +%F_%Hh%M`

mpirun $BRIDGE_MSUB_PWD/ramses3d $BRIDGE_MSUB_PWD/rad_aton.nml > log_$DATE.log
