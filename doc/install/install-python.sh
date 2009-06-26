#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#

set APPVER=2.5.4
#APPVER=2.6.2
set APPFILE=Python-${APPAPPVER}.tar.bz2
set APP_DIR=Python-${APPVER}
set URL=http://www.python.org/ftp/python/${APPVER}/${APPFILE}

set DOWNLOAD_DIR=_downloaded
set BUILD_DIR=_build
set SOURCE_DIR=_source

mkdir ${DOWNLOAD_DIR}
mkdir ${SOURCE_DIR}
mkdir ${BUILD_DIR}


echo "Downloading archive file..."
cd  ${DOWNLOAD_DIR}
if [ -e ${APPFILE} ] ; then
	echo "...File already downloaded";
else
	curl -O ${URL} ;
fi

cd ..
echo "Done"


cd ${SOURCE_DIR}
rm -Rf ${APP_DIR}
echo "Unpacking source files.."
tar -xf ../${DOWNLOAD_DIR}/$APPFILE
echo "..Done"

cd ..

echo "Building files.."
cd ${BUILD_DIR}
rm -Rf ${APP_DIR}
mkdir ${APP_DIR}
cd ${APP_DIR}

#--with-pydebug
#--with-dyld \
#--enable-framework=${PREFIX}/Framework \
../../${SOURCE_DIR}/${APP_DIR}/configure \
	--prefix=${PREFIX} \
	--program-suffix=.exe 
make 
make install
echo "..Done"
