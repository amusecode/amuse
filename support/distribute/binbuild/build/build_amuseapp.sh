#!/bin/bash


BASEDIR=`pwd`
PLATFORM=`uname`
ARCHITECTURE=`uname -m`

TMPDIR=/data1/vanelteren/tmp
SVNURL=http://www.amusecode.org/svn/trunk
SVNOPTS=
NOW=$(date +"%Y%m%d-%H%M%S")
NOWSHORT=$(date +"%Y%m%d")
WORKDIR=amuse_$NOW
#REVISION=`svn --password=svn2amuse --username=svn info -r HEAD $SVNURL | grep Revision | awk -- '{print $2}'`
VERSION=$NOWSHORT

PYTHONMAJOR="2"
PYTHONMINOR="7"
PYTHONRELEASE="2"
PYTHONPRERELEASE=0
PYTHONMAJORMINOR="${PYTHONMAJOR}.${PYTHONMINOR}"
PYTHONVERSION="${PYTHONMAJOR}.${PYTHONMINOR}.${PYTHONRELEASE}"

OPENSSLVERSION="0.9.8t"

if [ ${PYTHONPRERELEASE} == 1 ]; then
    FTPPATH="http://www.python.org/ftp/python/${PYTHONMAJORMINOR}"
else
    FTPPATH="http://www.python.org/ftp/python/${PYTHONVERSION}"
fi

UNICODETYPE="ucs4"

INSTALLDIR="${BASEDIR}/py_install"
SHELLDIR="${BASEDIR}/../shell"

RELEASEDIR=amuse-${VERSION}-${PLATFORM}_${ARCHITECTURE}
DISTFILE=${RELEASEDIR}.tar.gz

rm -f amuse-*-${PLATFORM}_${ARCHITECTURE}.tar.gz

echo "Distfile = ${DISTFILE}"

if [ ! -e "installed" ]; then
    if [ ${PLATFORM} == 'Darwin' ]; then
        # delete previous source
        rm -rf Python-${PYTHONVERSION} || exit $?
        if [ ! -e "Python-${PYTHONVERSION}.tar.bz2" ]; then
            # download
            curl -OL ${FTPPATH}/Python-${PYTHONVERSION}.tar.bz2 || exit $?
        fi
        # extract
        tar -xjvf Python-${PYTHONVERSION}.tar.bz2 || exit $?
        cd Python-${PYTHONVERSION}
    
        if [ ${ARCHITECTURE} == 'i386' ]; then
            export MACOSX_DEPLOYMENT_TARGET=10.5
        else
            export MACOSX_DEPLOYMENT_TARGET=10.6
        fi

        # configure
        
        # Build Python
        ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --disable-framework --disable-universalsdk
        
        # build
        make || exit $?
    
        # install
        make install || exit $?

    else
        # delete previous source
        rm -rf Python-${PYTHONVERSION} || exit $?
        rm -rf openssl-${OPENSSLVERSION} || exit $?
        
        if [ ! -e "openssl-${OPENSSLVERSION}.tar.gz" ]; then
            # download
            wget http://www.openssl.org/source/openssl-${OPENSSLVERSION}.tar.gz || exit $?
        fi
        
        tar zxf openssl-${OPENSSLVERSION}.tar.gz || exit $?
        
        cd openssl-${OPENSSLVERSION}
        
        ./config --prefix=${INSTALLDIR}  --openssldir=${INSTALLDIR}/openssl --shared || exit $?
        
        make || exit $?
        
        make install  || exit $?
        
        cd ${BASEDIR}
    
        if [ ! -e "Python-${PYTHONVERSION}.tar.bz2" ]; then
            # download
            wget ${FTPPATH}/Python-${PYTHONVERSION}.tar.bz2 || exit $?
        fi
        # extract
        tar -xjvf Python-${PYTHONVERSION}.tar.bz2 || exit $?
        
        cd Python-${PYTHONVERSION}
    

        # configure
        ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --enable-shared LD_LIBRARY_PATH=${INSTALLDIR}/lib || exit $?
        
        # patch 
        patch setup.py < ${BASEDIR}/setup_linux.patch || exit $?
        
        # build
        make || exit $?
    
        # install
        make install || exit $?
    
    fi
    # we're done
    cd ${BASEDIR}
    rm -Rf libsinstalled || exit $?
    rm -Rf pipsinstalled || exit $?
    rm -Rf amuseinstalled || exit $?
    touch "installed" || exit $?
fi

 
# run the python build script
export PYTHONHOME="${BASEDIR}/py_install"
export PATH=${PYTHONHOME}/bin:$PATH
export PYTHON=${PYTHONHOME}/bin/python
export LD_LIBRARY_PATH=${PYTHONHOME}/lib
export FC=gfortran
export F77=gfortran
export

if [ ! -e "libsinstalled" ]; then
    ${PYTHON} build_libraries.py install  || exit $?
    ${PYTHONHOME}/bin/easy_install pip
    touch "libsinstalled" || exit $?
fi


#${PYTHONHOME}/bin/pip install ipython || exit $?
if [ ! -e "pipsinstalled"  ]; then
    ${PYTHONHOME}/bin/easy_install pip
    ${PYTHONHOME}/bin/easy_install readline
    ${PYTHONHOME}/bin/easy_install IPython
    ${PYTHONHOME}/bin/pip install Flask || exit $?
    ${PYTHONHOME}/bin/pip install matplotlib || exit $?
    touch "pipsinstalled" || exit $?
fi

if [ ! -e "amuseinstalled" ]; then
    
    export PATH=${BASEDIR}/static_libs/bin:${PYTHONHOME}/bin:$PATH
    export
    
    svn export --password=reviewboard2amuse --username=reviewboard ${SVNURL} amuse-src || exit $?
    
    cd amuse-src

    ./configure --with-fftw=${BASEDIR}/static_libs --with-hdf5=${PYTHONHOME} || exit $?

    ${PYTHON} setup.py install || exit $?

    cd ..
    
    make distclean

    touch "amuseinstalled" || exit $?
fi


if [ ${PLATFORM} == "Darwin" ]; then
    echo 'move refs'
    ${PYTHON} mvpath.py -p ${PYTHONHOME}/lib/
    ${PYTHON} mvref.py -p ${PYTHONHOME}/lib/ -b ${PYTHONHOME}
else
    echo 'move refs'
    ${PYTHON} linux_set_rpath.py --path=${PYTHONHOME}/lib/ --bin-path=${PYTHONHOME}/ || exit $?
fi

rm -Rf ${RELEASEDIR}

cp -R ${INSTALLDIR} ${RELEASEDIR}

cp -R ${SHELLDIR}/* ${RELEASEDIR}

tar -czf ${DISTFILE} ${RELEASEDIR}

rm -Rf ${RELEASEDIR}
