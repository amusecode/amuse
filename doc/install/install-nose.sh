#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#

#PREFIX=${HOME}/_abs

APPVER=0.11.1
APPFILE=nose-${APPVER}.tar.gz
APP_DIR=nose-${APPVER}
URL=http://somethingaboutorange.com/mrl/projects/nose/${APPFILE}

DOWNLOAD_DIR=_downloaded
BUILD_DIR=_build
SOURCE_DIR=_source

#rm -Rf ${DOWNLOAD_DIR}
#rm -Rf ${SOURCE_DIR}
#rm -Rf ${BUILD_DIR}
mkdir ${DOWNLOAD_DIR}
mkdir ${SOURCE_DIR}
mkdir ${BUILD_DIR}

echo "Downloading all files"
cd  ${DOWNLOAD_DIR}
curl -O ${URL} 

cd ..
echo "Done"


cd ${SOURCE_DIR}

echo "Unpacking source files.."
rm -Rf ${APP_DIR}
tar -xf ../${DOWNLOAD_DIR}/$APPFILE
echo "..Done"

cd ..

echo "Building files.."
cd ${SOURCE_DIR}
cd ${APP_DIR}

python setup.py build
python setup.py install

echo "..Done"
