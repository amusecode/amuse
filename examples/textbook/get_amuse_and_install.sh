#! /bin/bash
export GITHASH=0b7d0c5aa8dcaa2fc564aba62dc9a7f686c29382
git clone https://github.com/amusecode/amuse.git
## git checkout ${GITHASH}
git reset --hard ${GITHASH}
mv amuse amuse-${GITHASH}

## export SCRIPT_DIR="${PWD}"
export AMUSE_DIR="${PWD}/amuse-${GITHASH}"

cd ${AMUSE_DIR}
./configure && \
make framework

