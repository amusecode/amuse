#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#

#APPVER=2.5.4
APPVER=2.6.3
APPFILE=Python-${APPVER}.tar.bz2
APP_DIR=Python-${APPVER}
URL=http://www.python.org/ftp/python/${APPVER}/${APPFILE}

if [ -z ${PREFIX} ]; then
	echo The PREFIX variable is not set, please set it to an user directory
	exit 1
fi

if [ ! -d ${PREFIX} ]; then
	echo ${PREFIX} directory does not exists, please create it first!
	exit 1
 
fi

 

INSTALL_DIR=$PREFIX/install
mkdir $INSTALL_DIR
cd $INSTALL_DIR

DOWNLOAD_DIR=$INSTALL_DIR/_downloaded
BUILD_DIR=$INSTALL_DIR/_build
SOURCE_DIR=$INSTALL_DIR/_source

mkdir ${DOWNLOAD_DIR}
mkdir ${SOURCE_DIR}
mkdir ${BUILD_DIR}


echo "Downloading archive file..."
cd  ${DOWNLOAD_DIR}
if [ -e ${APPFILE} ] ; then
	echo "...File already downloaded";
else
	if which curl >/dev/null; then
		curl -O ${URL} ;
	else
		wget ${URL};		
	fi
fi

cd ..
echo "Done"


cd ${SOURCE_DIR}
rm -Rf ${APP_DIR}
echo "Unpacking source files.."
tar -xf ${DOWNLOAD_DIR}/$APPFILE
echo "..Done"

cd ..

echo "Building files.."
cd ${BUILD_DIR}
rm -Rf ${APP_DIR}
mkdir ${APP_DIR}
cd ${APP_DIR}

UNAME=`uname`
if [ $UNAME == 'Darwin' ] ; then
	${SOURCE_DIR}/${APP_DIR}/configure \
		--with-dyld \
		--enable-framework=${PREFIX}/Framework \
		--enable-universal-sdk \
		--prefix=${PREFIX} \
		--program-suffix=.exe ;
else
	${SOURCE_DIR}/${APP_DIR}/configure \
		--prefix=${PREFIX} \
		--enable-shared \
		--program-suffix=.exe ;
fi

make 
make install
make 
make install
echo "..Done"
