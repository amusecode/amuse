#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#


APPVER=4.4.4
APPFILE=starlab.tar.gz
APP_DIR=starlab-${APPVER}
URL=http://www.ids.ias.edu/~starlab/download/starlab.tar.gz

DOWNLOAD_DIR=_downloaded
BUILD_DIR=_build
SOURCE_DIR=_source

#rm -Rf ${DOWNLOAD_DIR}
#rm -Rf ${SOURCE_DIR}
#rm -Rf ${BUILD_DIR}
mkdir ${DOWNLOAD_DIR}
mkdir ${SOURCE_DIR}
mkdir ${BUILD_DIR}

echo "Downloading archive file..."
cd  ${DOWNLOAD_DIR}
if [ -e ${APPFILE} ]; then
	echo "...File already downloaded";
else
	curl -O ${URL} ;
fi

cd ..
echo "...Done"


cd ${SOURCE_DIR}

echo "Unpacking source files.."
rm -Rf ${APPFILE}
tar -xf ../${DOWNLOAD_DIR}/$APPFILE
echo "..Done"

cd ..

echo "Building files.."
cd ${SOURCE_DIR}
cd ${APP_DIR}

#--with-pydebug
./configure --prefix=${PREFIX}  --with-pic --with-qt=no F77=gfortran LDFLAGS=-fPIC CFLAGS=-fPIC FFLAGS=-fPIC
make 
make install
echo "..Done"
