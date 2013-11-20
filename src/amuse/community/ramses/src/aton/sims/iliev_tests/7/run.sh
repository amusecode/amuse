#!/bin/bash
#MSUB -r tstranex2           # Nom du job                
#MSUB -n 16                  # Reservation de 4 processus au total
#MSUB -N 8                  # Les processus sont répartis sur 2 noeuds
#MSUB -T 36000                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_test7aton          # Sortie standard
#MSUB -e stderr_test7aton          # Sortie d'erreur       
#MSUB -p gen2191      # Allocation
#MSUB -q gpu                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x

mpirun $BRIDGE_MSUB_PWD/ramses3d $BRIDGE_MSUB_PWD/test7_aton.nml > log_$DATE.log
