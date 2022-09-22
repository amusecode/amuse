#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
configuration from config.mk
"""
import os
import warnings


def parse_configmk(filename):
    configfile = open(filename, "r")
    lines = configfile.readlines()
    configfile.close()
    cfgvars = dict()
    if "amuse configuration" not in lines[0]:
        raise Exception(
            "file: {0} is not an amuse configuration file".format(filename)
        )
    for line in lines:
        if "=" in line:
            var, value = line.split("=", 1)
            if value.startswith("@") and value.endswidth("@"):
                warnings.warn(
                    "possible configuration error/ unconfigured variable in"
                    " {0}".format(filename)
                )
            cfgvars[var] = value.strip()
    return cfgvars


try:
    configmk = parse_configmk("config.mk")
except IOError:
    from .support import get_amuse_root_dir
    configmk = parse_configmk(
        os.path.join(get_amuse_root_dir(), "config.mk")
    )


class interpreters():
    python = configmk["PYTHON"]


class compilers():
    cxx = configmk["CXX"]
    cc = configmk["CC"]
    fc = configmk["FC"]

    cxx_flags = configmk["CXXFLAGS"]
    cc_flags = configmk["CFLAGS"]
    fc_flags = configmk["FCFLAGS"]
    ld_flags = configmk["LDFLAGS"]

    found_fftw = configmk["FOUND_FFTW"]
    fftw_flags = configmk["FFTW_FLAGS"]
    fftw_libs = configmk["FFTW_LIBS"]

    found_gsl = configmk["FOUND_GSL"]
    gsl_flags = configmk["GSL_FLAGS"]
    gsl_libs = configmk["GSL_LIBS"]

    gfortran_version = configmk["GFORTRAN_VERSION"]
    ifort_version = configmk["IFORT_VERSION"]

    fc_iso_c_bindings = configmk["FC_ISO_C_AVAILABLE"] == 'yes'

    cython = configmk["CYTHON"]
    pythondev_cflags = configmk["PYTHONDEV_CFLAGS"]
    pythondev_ldflags = configmk["PYTHONDEV_LDFLAGS"]


class mpi():
    is_enabled = configmk["MPI_ENABLED"] == 'yes'
    mpicxx = configmk["MPICXX"]
    mpicc = configmk["MPICC"]
    mpif95 = configmk["MPIFC"]
    mpifc = configmk["MPIFC"]
    mpif90 = configmk["MPIFC"]
    mpiexec = configmk["MPIEXEC"]

    mpi_cflags = configmk["MPI_CFLAGS"]
    mpi_cxxflags = configmk["MPI_CXXFLAGS"]
    mpi_fcflags = configmk["MPI_FCFLAGS"]
    mpi_clibs = configmk["MPI_CLIBS"]
    mpi_cxxlibs = configmk["MPI_CXXLIBS"]
    mpi_fclibs = configmk["MPI_FCLIBS"]


class java():
    is_enabled = configmk["JAVA_ENABLED"] == 'yes'
    java = configmk["JAVA"]
    javac = configmk["JAVAC"]
    jar = configmk["JAR"]
    version = configmk["JAVA_VERSION"]


class cuda():
    is_enabled = configmk["CUDA_ENABLED"] == 'yes'
    compiler = configmk["NVCC"]
    compiler_flags = configmk["NVCC_FLAGS"]
    toolkit_path = configmk["CUDA_TK"]
    sdk_path = "/TOBEFIXED"
    cuda_libs = configmk["CUDA_LIBS"]
    sapporo_version = configmk["SAPPORO_VERSION"]


class openmp():
    is_enabled = configmk["OPENMP_ENABLED"] == 'yes'
    fcflags = configmk["OPENMP_FCFLAGS"]
    cflags = configmk["OPENMP_CFLAGS"]
