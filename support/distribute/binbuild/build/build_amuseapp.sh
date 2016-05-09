#!/bin/bash


BASEDIR=`pwd`
PLATFORM=`uname`
ARCHITECTURE=`uname -m`

#TMPDIR=/data1/vanelteren/tmp
SVNURL=http://www.amusecode.org/svn/trunk
SVNOPTS=
NOW=$(date +"%Y%m%d-%H%M%S")
NOWSHORT=$(date +"%Y%m%d")
WORKDIR=amuse_$NOW
#REVISION=`svn --password=svn2amuse --username=svn info -r HEAD $SVNURL | grep Revision | awk -- '{print $2}'`

echo "Version to build:" ${VERSION:=${NOWSHORT}}

PYTHONMAJOR="2"
PYTHONMINOR="7"
PYTHONRELEASE="9"
PYTHONPRERELEASE=0
PYTHONMAJORMINOR="${PYTHONMAJOR}.${PYTHONMINOR}"
PYTHONVERSION="${PYTHONMAJOR}.${PYTHONMINOR}.${PYTHONRELEASE}"

OPENSSLVERSION="1.0.1s"

if [ ${PYTHONPRERELEASE} == 1 ]; then
    FTPPATH="http://www.python.org/ftp/python/${PYTHONMAJORMINOR}"
else
    FTPPATH="http://www.python.org/ftp/python/${PYTHONVERSION}"
fi

UNICODETYPE="ucs4"

INSTALLDIR="${BASEDIR}/py_install"
SHELLDIR="${BASEDIR}/../shell"
TUTORIALDIR="${BASEDIR}/../../../../doc/interactive_tutorial"

RELEASEDIR=amuse-${VERSION}-${PLATFORM}_${ARCHITECTURE}
DISTFILE=${RELEASEDIR}.tar.gz

rm -f amuse-*-${PLATFORM}_${ARCHITECTURE}.tar.gz

echo "Distfile = ${DISTFILE}"

if [ ! -e "installed" ]; then
    if [ ${PLATFORM} == 'Darwin' ]; then
        rm -rf openssl-${OPENSSLVERSION} || exit $?
        
        if [ ! -e "openssl-${OPENSSLVERSION}.tar.gz" ]; then
            # download
            curl -OL http://www.openssl.org/source/openssl-${OPENSSLVERSION}.tar.gz || exit $?
        fi
        
        tar zxf openssl-${OPENSSLVERSION}.tar.gz || exit $?
        
        cd openssl-${OPENSSLVERSION}
        
        ./config --prefix=${INSTALLDIR}  --openssldir=${INSTALLDIR}/openssl --shared || exit $?
        if [ ${ARCHITECTURE} == 'x86_64' ]; then
            ./Configure darwin64-x86_64-cc --openssldir=${INSTALLDIR}/openssl --shared  || exit $?
        fi
        
        make || exit $?
        
        make install  || exit $?
        
        cd ${BASEDIR}
        
        tar -xvf certs.tar.gz -C ${INSTALLDIR}/openssl/certs  || exit $?
        
        # delete previous source
        rm -rf Python-${PYTHONVERSION} || exit $?
        if [ ! -e "Python-${PYTHONVERSION}.tgz" ]; then
            # download
            curl -OkL ${FTPPATH}/Python-${PYTHONVERSION}.tgz || exit $?
        fi
        # extract
        tar -xvf Python-${PYTHONVERSION}.tgz || exit $?
        cd Python-${PYTHONVERSION}
    
        if [ ${ARCHITECTURE} == 'i386' ]; then
            export MACOSX_DEPLOYMENT_TARGET=10.5
        else
            export MACOSX_DEPLOYMENT_TARGET=10.6
        fi

        # configure
        
        # Build Python
        ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --disable-framework --disable-universalsdk

        # patch 
        patch setup.py < ${BASEDIR}/setup_linux.patch || exit $?
        
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
        if [ ${ARCHITECTURE} == 'x86_64' ]; then
            ./Configure linux-x86_64 --prefix=${INSTALLDIR} --openssldir=${INSTALLDIR}/openssl --shared  || exit $?
        fi
        make || exit $?
        
        make install  || exit $?
        
        cd ${BASEDIR}
        
        tar -xvf certs.tar.gz -C ${INSTALLDIR}/openssl/certs  || exit $?
    
        if [ ! -e "Python-${PYTHONVERSION}.tgz" ]; then
            # download
            wget ${FTPPATH}/Python-${PYTHONVERSION}.tgz  --no-check-certificate|| exit $?
        fi
        # extract
        tar -xvf Python-${PYTHONVERSION}.tgz || exit $?
        
        cd Python-${PYTHONVERSION}
    

        # configure
        
        if [ ${PLATFORM} == 'Darwin' ]; then
            if [ ${ARCHITECTURE} == 'i386' ]; then
                export MACOSX_DEPLOYMENT_TARGET=10.5
            else
                export MACOSX_DEPLOYMENT_TARGET=10.6
            fi
        
            ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --disable-framework --disable-universalsdk || exit $?
        else
            ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --enable-shared LD_LIBRARY_PATH=${INSTALLDIR}/lib || exit $?
        fi
        
        # patch 
        patch setup.py < ${BASEDIR}/setup_linux.patch || exit $?
        
        # build
        make || exit $?
    
        # install
        make install || exit $?
        
    
    fi
    
    cd ${BASEDIR}
    
    export PYTHONHOME="${BASEDIR}/py_install"
    export PATH=${PYTHONHOME}/bin:$PATH
    export PYTHON=${PYTHONHOME}/bin/python
    export LD_LIBRARY_PATH=${PYTHONHOME}/lib/
    
    ${PYTHONHOME}/bin/python make_cert_links.py ${INSTALLDIR}/openssl/certs || exit $?
    
    # we're done
    rm -Rf libsinstalled || exit $?
    rm -Rf pipsinstalled || exit $?
    rm -Rf ytinstalled || exit $?
    rm -Rf amuseinstalled || exit $?
    touch "installed" || exit $?
fi

 
# run the python build script
export PYTHONHOME="${BASEDIR}/py_install"
export PATH=${PYTHONHOME}/bin:$PATH
export PYTHON=${PYTHONHOME}/bin/python
export LD_LIBRARY_PATH=${PYTHONHOME}/lib/
export PKG_CONFIG_PATH=${PYTHONHOME}/lib/pkgconfig/
export FC=gfortran
export F77=gfortran
if [ ${PLATFORM} == 'Darwin' ]; then
    export CC=cc
    export CXX=c++
fi

if [ ! -e "libsinstalled" ]; then
    ${PYTHON} build_libraries.py install || exit $?
    
    #rebuild python
    
    cd Python-${PYTHONVERSION}
    
    make clean || exit $?
    
    if [ ${PLATFORM} == 'Darwin' ]; then
        if [ ${ARCHITECTURE} == 'i386' ]; then
            export MACOSX_DEPLOYMENT_TARGET=10.5
        else
            export MACOSX_DEPLOYMENT_TARGET=10.6
        fi
    
        ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --disable-framework --disable-universalsdk || exit $?
    else
        ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --enable-shared LD_LIBRARY_PATH=${INSTALLDIR}/lib || exit $?
    fi
            
    make || exit $?
    
    make install || exit $?
    
    cd ..
    
    touch "libsinstalled" || exit $?
    
    ${PYTHONHOME}/bin/easy_install pip
fi


#${PYTHONHOME}/bin/pip install ipython || exit $?
if [ ! -e "pipsinstalled"  ]; then
    ${PYTHONHOME}/bin/easy_install pip
    
    ${PYTHONHOME}/bin/easy_install readline
    
    export PIP_CERT=`python -m pip._vendor.requests.certs`
    
    export PIP_INSTALL_OPTION=--zmq=${PYTHONHOME}
    
    ${PYTHONHOME}/bin/pip install pyzmq || exit $?
    
    export PIP_INSTALL_OPTION=
    
    ${PYTHONHOME}/bin/easy_install tornado
    
    ${PYTHONHOME}/bin/easy_install IPython
    
    ${PYTHONHOME}/bin/pip install Cython || exit $?
    
    ${PYTHONHOME}/bin/pip install Flask || exit $?
    
    mkdir mpl 
    
    py_install/bin/pip install --download mpl 'matplotlib==1.2.1' || exit $?
    
    cd mpl
    
    tar -xvf matplotlib-1.2.1.tar.gz
    
    cd matplotlib-1.2.1
    
    export CFLAGS="-I${PYTHONHOME}/include -I${PYTHONHOME}/include/freetype2"

    export LDFLAGS="-L${PYTHONHOME}/lib"
    
    ${PYTHONHOME}/bin/python setup.py install || exit $?
    
    cd ../../
    
    rm -Rf mpl || exit $?
    
    touch "pipsinstalled" || exit $?
fi

touch "ytinstalled" || exit $?
if [ ! -e "ytinstalled"  ]; then

    rm -Rf yt-hg || exit $?
    
    curl -OL http://bitbucket.org/yt_analysis/yt/get/tip.tar.gz  || exit $?
    
    tar zxf tip.tar.gz || exit $?
    
    cd yt_analysis-yt-*  || exit $?
    
    echo $INSTALLDIR > hdf5.cfg
    
    ${PYTHON} setup.py install || exit $?
    
    cd ..
    
    rm -Rf yt_analysis-yt-*  || exit $?
    
    rm -Rf tip.tar.gz || exit $?
    
    touch "ytinstalled" || exit $?
fi

if [ ! -e "amuseinstalled" ]; then
    
    export PATH=${BASEDIR}/static_libs/bin:${PYTHONHOME}/bin:$PATH
    export BUILD_BINARY=1    
    export DOWNLOAD_CODES=1    
    
    cd ../../../.. # cd to amuse root directory
    
    #make distclean PYTHON=${PYTHON}

    if [ ${PLATFORM} == "Darwin" ]; then
        export CXXCPP="g++ -E"
    fi
    ./configure --with-fftw=${BASEDIR}/static_libs --with-hdf5=${PYTHONHOME} PYTHON=${PYTHON} || exit $?
    
    ${PYTHON} setup.py install || exit $?

    make distclean PYTHON=${PYTHON}

    cd support/distribute/binbuild/build

    touch "amuseinstalled" || exit $?
fi


if [ ${PLATFORM} == "Darwin" ]; then
    echo 'move refs'
   
    cp /usr/local/lib/*.dylib ${PYTHONHOME}/lib 

    chmod u+w ${PYTHONHOME}/lib/lib*.dylib

    ${PYTHON} mvpath.py -p ${PYTHONHOME}/lib/
    ${PYTHON} mvref.py -p ${PYTHONHOME}/lib/ -b ${PYTHONHOME}
    ${PYTHON} mvref.py -p /usr/local/lib/ -b ${PYTHONHOME}/lib -r ./
else
    echo 'move refs'
    ${PYTHON} linux_set_rpath.py --path=${PYTHONHOME}/lib/ --bin-path=${PYTHONHOME}/ || exit $?
fi

rm -Rf ${RELEASEDIR}

cp -R ${INSTALLDIR} ${RELEASEDIR}

cp -R ${SHELLDIR}/* ${RELEASEDIR}

cp -R ${TUTORIALDIR} ${RELEASEDIR}/tutorial

tar -czf ${DISTFILE} ${RELEASEDIR}

rm -Rf ${RELEASEDIR}
