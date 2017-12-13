#!/bin/bash

BUILD=amuse
FIXREFS=yes
for i in "$@"
do
case $i in
    -O|--omuse)
    BUILD=omuse
    ;;
    -R|--nofixrefs)
    FIXREFS=no
    ;;
    -v=*|--version=*)
    VERSION="${i#*=}"
    ;;
    *)
            # unknown option
    ;;
esac
done

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

echo "OS X bit:" ${OSX_BIT:=64bit}

echo "Version to build:" ${VERSION:=${NOWSHORT}}

PYTHONMAJOR="2"
PYTHONMINOR="7"
PYTHONRELEASE="13"
PYTHONPRERELEASE=0
PYTHONMAJORMINOR="${PYTHONMAJOR}.${PYTHONMINOR}"
PYTHONVERSION="${PYTHONMAJOR}.${PYTHONMINOR}.${PYTHONRELEASE}"

OPENSSLVERSION="1.0.2l"

PIPVERSION=8.1.2
READLINEVERSION=6.2.4.1
PYZMQVERSION=15.4.0
TORNADOVERSION=4.4.1
JUPYTERVERSION=1.0.0
IPYTHONVERSION=4.0.0
CYTHONVERSION=0.24.1
FLASKVERSION=0.11.1
PILLOWVERSION=3.3.1
MATPLOTLIBVERSION=1.4.3

if [ ${PYTHONPRERELEASE} == 1 ]; then
    FTPPATH="https://www.python.org/ftp/python/${PYTHONMAJORMINOR}"
else
    FTPPATH="https://www.python.org/ftp/python/${PYTHONVERSION}"
fi

UNICODETYPE="ucs4"

INSTALLDIR="${BASEDIR}/py_install"
SHELLDIR="${BASEDIR}/../shell"
TUTORIALDIR="${BASEDIR}/../../../../doc/interactive_tutorial"
if [ "${BUILD}" == "omuse" ]; then
  TUTORIALDIR="${BASEDIR}/../../../../doc/omuse_tutorial"
fi

if [ ${OSX_BIT} == '32bit' ]; then
	ARCHITECTURE=i386
fi

RELEASEDIR=${BUILD}-${VERSION}-${PLATFORM}_${ARCHITECTURE}
DISTFILE=${RELEASEDIR}.tar.gz

rm -f ${BUILD}-*-${PLATFORM}_${ARCHITECTURE}.tar.gz

echo "Distfile = ${DISTFILE}"

if [ ! -e "installed" ]; then

    if [ ${PLATFORM} != 'Darwin' ]; then
        tar xvzf bzip2-1.0.6_light.tar.gz
        cd bzip2-1.0.6_light
        make install PREFIX="${BASEDIR}/py_install"
        cd ..
    fi

    if [ ${PLATFORM} == 'Darwin' ]; then
        rm -rf openssl-${OPENSSLVERSION} || exit $?
        
        if [ ! -e "openssl-${OPENSSLVERSION}.tar.gz" ]; then
            # download
            curl -OL https://www.openssl.org/source/openssl-${OPENSSLVERSION}.tar.gz || exit $?
        fi
        
        tar zxf openssl-${OPENSSLVERSION}.tar.gz || exit $?
        
        cd openssl-${OPENSSLVERSION}
        
        #./config --prefix=${INSTALLDIR}  --openssldir=${INSTALLDIR}/openssl --shared || exit $?
	if [ ${OSX_BIT} == '64bit' ]; then
        	./Configure darwin64-x86_64-cc --prefix=${INSTALLDIR}  --openssldir=${INSTALLDIR}/openssl --shared || exit $?
	else
        	./Configure darwin-i386-cc --prefix=${INSTALLDIR}  --openssldir=${INSTALLDIR}/openssl --shared || exit $?
	fi
        
        make || exit $?
        
        make install  || exit $?
        
        cd ${BASEDIR}
        
        tar -xvf certs.tar.gz -C ${INSTALLDIR}/openssl/certs  || exit $?
         
	if [ ${OSX_BIT} == '32bit' ]; then
		export CFLAGS="-arch i386"
		export FCFLAGS="-arch i386"
		export CXXFLAGS="-arch i386"
		export LDFLAGS="-arch i386"
    	fi
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
            export MACOSX_DEPLOYMENT_TARGET=10.8
        fi

        # configure
        
        # Build Python
        ./configure --enable-unicode=${UNICODETYPE} --prefix=${INSTALLDIR} --disable-framework --disable-universalsdk CC=/usr/bin/gcc CXX=/usr/bin/g++

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
            wget https://www.openssl.org/source/openssl-${OPENSSLVERSION}.tar.gz || exit $?
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
export PREFIX="${BASEDIR}/py_install"
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
	if [ ${OSX_BIT} == '32bit' ]; then
		export CFLAGS="-arch i386"
		export CXXFLAGS="-arch i386"
		export LDFLAGS="-arch i386"
		export CC="/usr/bin/gcc -arch i386"
		export CXX="/usr/bin/g++ -arch i386"
		export FC="gfortran -arch i386"
		export F77="gfortran -arch i386"
    	fi
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
    
    ${PYTHONHOME}/bin/easy_install pip==${PIPVERSION}
fi


if [ ! -e "pipsinstalled"  ]; then
    ${PYTHONHOME}/bin/easy_install pip==${PIPVERSION}

    export PIP_CERT=`python -m pip._vendor.requests.certs`    

    ${PYTHONHOME}/bin/pip install readline==${READLINEVERSION} || exit $?
        
    export PIP_INSTALL_OPTION=--zmq=${PYTHONHOME}
    
    ${PYTHONHOME}/bin/pip install pyzmq==${PYZMQVERSION} || exit $?
    
    export PIP_INSTALL_OPTION=
    

    ${PYTHONHOME}/bin/pip install tornado || exit $?
        
    #~ ${PYTHONHOME}/bin/pip install ipython[all] || exit $?
    # is this equivalent to..(?)
    ${PYTHONHOME}/bin/pip install ipython==${IPYTHONVERSION}  || exit $?

    ${PYTHONHOME}/bin/pip install jupyter==${JUPYTERVERSION}  || exit $?
    
    ${PYTHONHOME}/bin/pip install Cython==${CYTHONVERSION} || exit $?
    
    ${PYTHONHOME}/bin/pip install Flask==${FLASKVERSION} || exit $?
    
    ${PYTHONHOME}/bin/pip install --global-option="build_ext" --global-option="--disable-jpeg" pillow==${PILLOWVERSION} || exit $?
    
    rm -Rf mpl || exit $?
        
    mkdir mpl 
    
    cd mpl
    
    ${PYTHONHOME}/bin/pip download --no-deps --no-binary matplotlib matplotlib==${MATPLOTLIBVERSION} || exit $?
        
    tar -xvf matplotlib-${MATPLOTLIBVERSION}.tar.gz || exit $?
    
    cd matplotlib-${MATPLOTLIBVERSION} || exit $?
    
    export CFLAGS="-I${PYTHONHOME}/include -I${PYTHONHOME}/include/freetype2"

    export LDFLAGS="-L${PYTHONHOME}/lib"

    # hopefully maxosx backend will be enough
    if [ ${PLATFORM} == 'Darwin' ]; then
      echo "[gui_support]" > setup.cfg
      echo tkagg = False >> setup.cfg
    fi
    
    ${PYTHONHOME}/bin/python setup.py install || exit $?
    
    cd ../../
    
    touch "pipsinstalled" || exit $?
fi

touch "ytinstalled" || exit $?
if [ ! -e "ytinstalled"  ]; then

    rm -Rf yt-hg || exit $?
    
    curl -OL https://bitbucket.org/yt_analysis/yt/get/tip.tar.gz  || exit $?
    
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
    
    make distclean PYTHON=${PYTHON}

    if [ ${PLATFORM} == "Darwin" ]; then
        export CXXCPP="g++ -E"
    fi
    ./configure --with-fftw=${BASEDIR}/static_libs --with-hdf5=${PYTHONHOME} PYTHON=${PYTHON} || exit $?

    export PYTHONPATH=${PYTHONPATH}:`pwd`/src

    # first build omuse codes..
    if [ "${BUILD}" == "omuse" ]; then
      ${PYTHON} setup.py build_codes --codes-dir=src/omuse/community || exit $?
    fi
    ${PYTHON} setup.py install || exit $?

    make distclean PYTHON=${PYTHON}

    cd support/distribute/binbuild/build

    touch "amuseinstalled" || exit $?
fi


if [ "${FIXREFS}" == "yes" ]; then
  if [ ${PLATFORM} == "Darwin" ]; then
      echo 'move refs'
     
      cp /usr/local/lib/*.dylib ${PYTHONHOME}/lib 
  
      chmod u+w ${PYTHONHOME}/lib/lib*.dylib
  
      ${PYTHON} mvpath.py -p ${PYTHONHOME}/lib/
      ${PYTHON} mvref.py -p ${PYTHONHOME}/lib/ -b ${PYTHONHOME}
      ${PYTHON} mvref.py -p /usr/local/lib/ -b ${PYTHONHOME}/lib -r ./
  else
      echo 'move refs'
      
      chmod u+w ${PYTHONHOME}/lib/engines/*.so
      chmod u+w ${PYTHONHOME}/lib/*.so

      ${PYTHON} linux_set_rpath.py --path=${PYTHONHOME}/lib/ --bin-path=${PYTHONHOME}/ || exit $?
  fi
fi

if [ -e local_fixes ]; then
  source local_fixes
fi

rm -Rf ${RELEASEDIR}

cp -R ${INSTALLDIR} ${RELEASEDIR}

cp -R ${SHELLDIR}/* ${RELEASEDIR}

if [ "${BUILD}" == "omuse" ]; then
  mv ${RELEASEDIR}/amuse ${RELEASEDIR}/omuse
  mv ${RELEASEDIR}/amuse-tutorial ${RELEASEDIR}/omuse-tutorial
fi

mkdir ${RELEASEDIR}/tutorial
cp  ${TUTORIALDIR}/*.ipynb ${RELEASEDIR}/tutorial/
cp  ${TUTORIALDIR}/amuserc ${RELEASEDIR}/tutorial/

tar -czf ${DISTFILE} ${RELEASEDIR}

rm -Rf ${RELEASEDIR}
