#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#

#APPVER=2.5.4
APPVER=2.7.14
#APPVER=2.6.5
#APPVER=2.7.1
APPFILE=Python-${APPVER}.tar.bz2
APPFILE=Python-${APPVER}.tgz
APP_DIR=Python-${APPVER}
URL=https://www.python.org/ftp/python/${APPVER}/${APPFILE}


OPENSSLVERSION="1.0.2m"
OPENSSLFILE=openssl-${OPENSSLVERSION}.tar.gz 
OPENSSLURL=https://www.openssl.org/source/old/1.0.2/${OPENSSLFILE}
OPENSSLDIR=openssl-${OPENSSLVERSION}

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
		curl -L -O ${URL} ;
	else
		wget ${URL};		
	fi
fi

if [ -e ${OPENSSLFILE} ] ; then
	echo "...openssl file already downloaded";
else
	if which curl >/dev/null; then
		curl -L -O ${OPENSSLURL} ;
	else
		wget ${OPENSSLURL};
	fi
fi

cd ..
echo "Done"


cd ${SOURCE_DIR}
rm -Rf ${APP_DIR}
echo "Unpacking source files.."
tar -xf ${DOWNLOAD_DIR}/${APPFILE}
echo "..Done"

cd ${SOURCE_DIR}
rm -Rf ${OPENSSLDIR}
echo "Unpacking openssl source files.."
tar -xf ${DOWNLOAD_DIR}/${OPENSSLFILE} || exit $?
echo "..Done"
cd ..

echo "Building files.."

cd ${SOURCE_DIR}
cd ${OPENSSLDIR}


MACHINE=`(uname -m) 2>/dev/null`
PLATFORM=`uname`

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

cd ${BUILD_DIR}
rm -Rf ${APP_DIR}
mkdir ${APP_DIR}
cd ${APP_DIR}

UNAME=`uname`
if [ $UNAME == 'Darwin' ] ; then
	${SOURCE_DIR}/${APP_DIR}/configure \
		--with-dyld \
		--prefix=${PREFIX} \
	    --enable-unicode=ucs4\
		--program-suffix=.exe ;
#--enable-framework=${PREFIX}/Framework \
#--enable-universal-sdk \
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
make 
make install
echo "..Done"
