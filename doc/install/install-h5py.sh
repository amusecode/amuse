#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#

PREFIX=${HOME}/_abs

APPVER=1.1.0
APPFILE=h5py-${APPVER}.tar.gz
APP_DIR=h5py-${APPVER}
URL=http://h5py.googlecode.com/files/h5py-1.1.0.tar.gz

DOWNLOAD_DIR=_downloaded
BUILD_DIR=_build
SOURCE_DIR=_source

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
python setup.py \
	configure \
	--hdf5=${PREFIX}

python setup.py build
python setup.py install

echo "..Done"
