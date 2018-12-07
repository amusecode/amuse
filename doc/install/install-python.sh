#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#

# Prerequisites directory check
if [ -z ${PREFIX} ]; then
	echo The PREFIX variable is not set, please set it to an user directory
	exit 1
fi
if [ ! -d ${PREFIX} ]; then
	echo ${PREFIX} directory does not exists, please create it first!
	exit 1
fi

# Python version selection (default to Python 2)
if [ "${1}" = "python3" ]
then
	APPVER=3.6.5
else
	APPVER=2.7.15
fi

# Python
APPFILE=Python-${APPVER}.tgz
APP_DIR=Python-${APPVER}
URL=https://www.python.org/ftp/python/${APPVER}/${APPFILE}

# OpenSSL
OPENSSLVERSION="1.1.0g"
OPENSSLFILE=openssl-${OPENSSLVERSION}.tar.gz 
OPENSSLURL=https://www.openssl.org/source/old/1.1.0/${OPENSSLFILE}
OPENSSLDIR=openssl-${OPENSSLVERSION}

# Setting up the directory structure
INSTALL_DIR=${PREFIX}/install
mkdir ${INSTALL_DIR}
cd ${INSTALL_DIR}
DOWNLOAD_DIR=${INSTALL_DIR}/_downloaded
mkdir ${DOWNLOAD_DIR}
BUILD_DIR=${INSTALL_DIR}/_build
mkdir ${BUILD_DIR}
SOURCE_DIR=${INSTALL_DIR}/_source
mkdir ${SOURCE_DIR}

# Download
echo "Downloading"
cd  ${DOWNLOAD_DIR}
if [ -e ${APPFILE} ] ; then
	echo "Python already downloaded";
else
	echo "Downloading Python"
	if which curl >/dev/null; then
		curl -L -O ${URL} ;
	else
		wget ${URL};		
	fi
fi
if [ -e ${OPENSSLFILE} ] ; then
	echo "OpenSSL already downloaded";
else
	echo "Downloading OpenSSL"
	if which curl >/dev/null; then
		curl -L -O ${OPENSSLURL} ;
	else
		wget ${OPENSSLURL};
	fi
fi
cd ..
echo "Done"

# Unpack
echo "Unpacking"
cd ${SOURCE_DIR}
rm -Rf ${APP_DIR}
echo "Unpacking Python..."
tar -xf ${DOWNLOAD_DIR}/${APPFILE}
cd ${SOURCE_DIR}
rm -Rf ${OPENSSLDIR}
echo "Unpacking OpenSSL..."
tar -xf ${DOWNLOAD_DIR}/${OPENSSLFILE} || exit $?
cd ..
echo "Done"

# Build
echo "Building"
MACHINE=`(uname -m) 2>/dev/null`
PLATFORM=`uname`
echo "Building OpenSSL"
cd ${SOURCE_DIR}/${OPENSSLDIR}
if [ ${PLATFORM} == 'Darwin' ]; then
	if [ ${MACHINE} == 'x86_64' ]; then
		./Configure darwin64-x86_64-cc \
			--prefix=${PREFIX}  \
			--openssldir=${PREFIX}/openssl \
			--shared
	else
		./Configure darwin64-i386-cc \
			--prefix=${PREFIX}  \
			--openssldir=${PREFIX}/openssl \
			--shared
	fi
else
	./config \
		--prefix=${PREFIX}  \
		--openssldir=${PREFIX}/openssl \
		--shared
fi
make
make install
echo "Build Python"
cd ${BUILD_DIR}
rm -Rf ${APP_DIR}
mkdir ${APP_DIR}
cd ${APP_DIR}
if [ $PLATFORM == 'Darwin' ] ; then
	${SOURCE_DIR}/${APP_DIR}/configure \
		--with-dyld \
		--prefix=${PREFIX} \
	    --enable-unicode=ucs4\
		--program-suffix=.exe ;
else
	${SOURCE_DIR}/${APP_DIR}/configure \
		--prefix=${PREFIX} \
		--libdir=${PREFIX}/lib \
		--enable-shared \
	    --enable-unicode=ucs4\
		--program-suffix=.exe ;
fi
make 
make install
echo "Done"

echo "Install complete"
