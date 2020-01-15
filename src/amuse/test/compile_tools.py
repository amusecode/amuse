import os
import subprocess
import time
import shlex
import numpy
import shutil

from amuse.rfi.tools import create_java
from amuse.rfi.tools import create_c
from amuse.rfi.tools import create_fortran
from amuse.rfi.core import config
from amuse.test.amusetest import get_amuse_root_dir


def get_mpicc_name():
    try:
        from amuse import config
        is_configured = hasattr(config, 'mpi')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.mpi.mpicc
    else:
        return os.environ['MPICC'] if 'MPICC' in os.environ else 'mpicc'


def get_mpicxx_name():
    try:
        from amuse import config
        is_configured = hasattr(config, 'mpi')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.mpi.mpicxx
    else:
        return os.environ['MPICXX'] if 'MPICXX' in os.environ else 'mpicxx'


def get_mpicc_flags():
    try:
        from amuse import config
        is_configured = hasattr(config, 'compilers')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.compilers.cc_flags
    else:
        return ""


def get_mpicxx_flags():
    try:
        from amuse import config
        is_configured = hasattr(config, 'compilers')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.compilers.cxx_flags
    else:
        return ""

def get_mpif90_name():
    try:
        from amuse import config
        is_configured = hasattr(config, 'mpi')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.mpi.mpif95
    else:
        return os.environ['MPIFC'] if 'MPIFC' in os.environ else 'mpif90'


def get_mpif90_arguments():
    name = get_mpif90_name()
    return list(shlex.split(name))


def has_fortran_iso_c_binding():
    try:
        from amuse import config
        is_configured = hasattr(config, 'mpi')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.compilers.fc_iso_c_bindings
    else:
        return False


def get_ld_flags():
    try:
        from amuse import config
        is_configured = hasattr(config, 'compilers')
    except ImportError:
        is_configured = False

    if is_configured:
        return config.compilers.ld_flags
    else:
        return ""


def wait_for_file(filename):
    for dt in [0.01, 0.01, 0.02, 0.05]:
        if os.path.exists(filename):
            return
        time.sleep(dt)


def open_subprocess(arguments, stdin=None):
    process = subprocess.Popen(
        arguments,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    stdout, stderr = process.communicate(input=stdin)
    return process, stdout, stderr


def c_compile(objectname, string, extra_args=[]):
    if os.path.exists(objectname):
        os.remove(objectname)

    mpicc = get_mpicc_name()
    arguments = [mpicc]
    arguments.extend(get_mpicc_flags().split())
    rootdir = get_amuse_root_dir()
    arguments.extend(["-I", rootdir + "/lib/stopcond", "-x", "c", "-c"]
                     + extra_args + ["-o", objectname, "-"])

    process, stderr, stdout = open_subprocess(arguments, stdin=string)
    
    if process.returncode == 0:
        wait_for_file(objectname)

    if process.returncode != 0 or not os.path.exists(objectname):
        raise Exception("Could not compile {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(objectname, stdout, stderr, ' '.join(arguments)))


def cxx_compile(objectname, string, extra_args=[]):
    if os.path.exists(objectname):
        os.remove(objectname)

    mpicxx = get_mpicxx_name()
    arguments = [mpicxx]
    arguments.extend(get_mpicxx_flags().split())
    rootdir = get_amuse_root_dir()
    arguments.extend(["-I", rootdir + "/lib/stopcond", "-x", "c++", "-c"]
                     + extra_args + ["-o", objectname, "-"])

    process, stderr, stdout = open_subprocess(arguments, stdin=string)

    if process.returncode == 0:
        wait_for_file(objectname)

    if process.returncode != 0 or not os.path.exists(objectname):
        raise Exception("Could not compile {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(objectname, stdout, stderr, ' '.join(arguments)))


def c_build(exename, objectnames, extra_args=[]):
    if os.path.exists(exename):
        os.remove(exename)

    mpicxx = get_mpicxx_name()
    arguments = [mpicxx]
    arguments.extend(objectnames)
    arguments.append("-o")
    arguments.append(exename)

    if 'LIBS' in os.environ:
        libs = os.environ['LIBS'].split()
        arguments.extend(libs)

    arguments.extend(get_ld_flags().split())

    arguments.extend(extra_args)

    process, stderr, stdout = open_subprocess(arguments)

    if process.returncode == 0:
        wait_for_file(exename)

    if process.returncode != 0 or not (os.path.exists(exename)
                                       or os.path.exists(exename+'.exe')):
        raise Exception("Could not build {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(exename, stdout, stderr, ' '.join(arguments)))


def build_worker(codestring, path_to_results, specification_class, write_header=True, 
      extra_args=[], needs_mpi = False):
    path = os.path.abspath(path_to_results)
    codefile = os.path.join(path, "code.o")
    headerfile = os.path.join(path, "worker_code.h")
    interfacefile = os.path.join(path, "interface.o")
    exefile = os.path.join(path, "c_worker")

    c_compile(codefile, codestring)

    uc = create_c.GenerateACHeaderStringFromASpecificationClass()
    uc.specification_class = specification_class
    uc.needs_mpi = needs_mpi
    header = uc.result

    if write_header:
      with open(headerfile, "w") as f:
          f.write(header)
      extra_args=extra_args+["-I",path]

    uc = create_c.GenerateACSourcecodeStringFromASpecificationClass()
    uc.specification_class = specification_class
    uc.needs_mpi = needs_mpi
    code = uc.result

    cxx_compile(interfacefile, code if write_header else header+"\n\n"+code, 
                    extra_args=extra_args)
    c_build(exefile, [interfacefile, codefile], extra_args=extra_args)

    return exefile


def c_pythondev_compile(objectname, string):
    root, ext = os.path.splitext(objectname)
    sourcename = root + '.c'
    amuse_root = get_amuse_root_dir()
    if os.path.exists(objectname):
        os.remove(objectname)

    with open(sourcename, "w") as f:
        f.write(string)

    mpicc = get_mpicc_name()
    arguments = [mpicc]
    arguments.extend(get_mpicc_flags().split())
    arguments.extend(["-I",numpy.get_include()])

    arguments.extend(["-I", amuse_root + "/lib/stopcond","-I", amuse_root + "/lib/amuse_mpi",  "-fPIC", "-c",  "-o", objectname, sourcename])
    arguments.extend(shlex.split(config.compilers.pythondev_cflags))

    process, stderr, stdout = open_subprocess(arguments)

    if process.returncode == 0:
        wait_for_file(objectname)

    if process.returncode != 0 or not os.path.exists(objectname):
        raise Exception("Could not compile {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(objectname, stdout, stderr, ' '.join(arguments)))


def c_pythondev_build(exename, objectnames):
    if os.path.exists(exename):
        os.remove(exename)

    mpicxx = get_mpicxx_name()
    arguments = [mpicxx]
    arguments.extend(objectnames)
    arguments.extend(shlex.split(config.compilers.pythondev_ldflags))
    arguments.append("-o")
    arguments.append(exename)

    if 'LIBS' in os.environ:
        libs = os.environ['LIBS'].split()
        arguments.extend(libs)

    arguments.extend(get_ld_flags().split())

    process, stderr, stdout = open_subprocess(arguments)

    if process.returncode == 0:
        wait_for_file(exename)

    if process.returncode != 0 or not (os.path.exists(exename) or os.path.exists(exename+'.exe')):
        raise Exception("Could not build {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(exename, stdout, stderr, ' '.join(arguments)))


def c_pythondev_buildso(soname, objectnames):
    if os.path.exists(soname):
        os.remove(soname)
    amuse_root = get_amuse_root_dir()

    mpicxx = get_mpicxx_name()
    arguments = [mpicxx]
    arguments.extend(objectnames)
    arguments.extend(shlex.split(config.compilers.pythondev_ldflags))
    arguments.append("-shared")
    arguments.append("-o")
    arguments.append(soname)
    arguments.append("-L"+amuse_root+"/lib/amuse_mpi")
    arguments.append("-lamuse_mpi")

    if 'LIBS' in os.environ:
        libs = os.environ['LIBS'].split()
        arguments.extend(libs)

    arguments.extend(get_ld_flags().split())

    process, stderr, stdout = open_subprocess(arguments)

    if process.returncode == 0:
        wait_for_file(soname)

    if process.returncode != 0 or not (os.path.exists(soname)):
        raise Exception("Could not build {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(soname, stdout, stderr, ' '.join(arguments)))


def fortran_pythondev_buildso(soname, objectnames):
    if os.path.exists(soname):
        os.remove(soname)
    amuse_root = get_amuse_root_dir()

    arguments = get_mpif90_arguments()
    arguments.extend(objectnames)
    arguments.extend(shlex.split(config.compilers.pythondev_ldflags))
    arguments.append("-shared")
    arguments.append("-o")
    arguments.append(soname)
    arguments.append("-L"+amuse_root+"/lib/amuse_mpi")
    arguments.append("-lamuse_mpi")

    if 'LIBS' in os.environ:
        libs = os.environ['LIBS'].split()
        arguments.extend(libs)

    arguments.extend(get_ld_flags().split())

    process, stderr, stdout = open_subprocess(arguments)

    if process.returncode == 0:
        wait_for_file(soname)

    if process.returncode != 0 or not (os.path.exists(soname)):
        raise Exception("Could not build {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(soname, stdout, stderr, ' '.join(arguments)))


def fortran_compile(objectname, string, extra_args=[]):
    if os.path.exists(objectname):
        os.remove(objectname)

    root, ext = os.path.splitext(objectname)
    sourcename = root + '.f90'

    with open(sourcename, "w") as f:
        f.write(string)

    arguments = get_mpif90_arguments()
    arguments.extend(["-g",
                     "-I{0}/lib/forsockets".format(get_amuse_root_dir())]
                     + extra_args + ["-c",  "-o", objectname, sourcename])
    arguments.append("-Wall")

    process, stderr, stdout = open_subprocess(arguments)
    process.wait()

    if process.returncode == 0:
        wait_for_file(objectname)

    if process.returncode != 0 or not os.path.exists(objectname):
        raise Exception("Could not compile {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(objectname, stdout, stderr, ' '.join(arguments)))


def fortran_build(exename, objectnames, extra_args=[]):
    if os.path.exists(exename):
        os.remove(exename)

    arguments = get_mpif90_arguments()
    arguments.extend(objectnames)

    arguments.append("-L{0}/lib/forsockets".format(get_amuse_root_dir()))
    arguments.append("-Wall")
    arguments.append("-lforsockets")
    arguments.extend(extra_args)

    arguments.append("-o")
    arguments.append(exename)

    arguments.extend(get_ld_flags().split())

    process, stderr, stdout = open_subprocess(arguments)

    if process.returncode == 0:
        wait_for_file(exename)

    if process.returncode != 0 or not os.path.exists(exename):
        raise Exception("Could not build {0}\nstdout:\n{1}\nstderr:\n{2}\narguments:\n{3}".format(exename, stdout, stderr, ' '.join(arguments)))

def build_java_worker(codestring, path_to_results, specification_class):
    path = os.path.abspath(path_to_results)
    exefile = os.path.join(path,"java_worker")
   
    #generate interface class
    interfacefile = os.path.join(path,"CodeInterface.java")

    uc = create_java.GenerateAJavaInterfaceStringFromASpecificationClass()
    uc.specification_class = specification_class
    interface =  uc.result

    with open(interfacefile, "w") as f:
        f.write(interface)
    
    #generate worker class
    workerfile = os.path.join(path,"Worker.java")

    uc = create_java.GenerateAJavaSourcecodeStringFromASpecificationClass()
    uc.specification_class = specification_class
    worker =  uc.result

    with open(workerfile, "w") as f:
        f.write(worker)

    #write code imlementation to a file

    codefile = os.path.join(path,"Code.java")

    with open(codefile, "w") as f:
        f.write(codestring)

    #compile all code

    if not config.java.is_enabled:
        raise Exception("java not enabled")

    javac = config.java.javac
    if not os.path.exists(javac):
        raise Exception("java compiler not available")
        
    jar = config.java.jar

    jarfile = 'worker.jar'

    tmpdir = os.path.join(path, "_tmp")

    shutil.rmtree(tmpdir, ignore_errors=True)
    os.mkdir(tmpdir)

    returncode = subprocess.call(
        [javac, "-d", tmpdir, "Worker.java", "CodeInterface.java", "Code.java"],
    cwd = path,
    )
        
    if returncode != 0:
        print("Could not compile worker")

#make jar file

    returncode = subprocess.call(
        [jar, '-cf', 'worker.jar', '-C', tmpdir, '.'],
        cwd = path,
    )
        
    if returncode != 0:
        print("Could not compile worker")

    #generate worker script

    uc = create_java.GenerateAJavaWorkerScript()
    uc.specification_class = specification_class
    uc.code_dir = path
    worker_script =  uc.result

    with open(exefile, "w") as f:
        f.write(worker_script)

    os.chmod(exefile, 0o777)

    return exefile

def build_fortran_worker(codestring, path_to_results, interface, needs_mpi = True,
                             extra_fflags=[], extra_ldflags=[]):
    
    path = os.path.abspath(path_to_results)  
    codefile = os.path.join(path,"code.o")
    interfacefile = os.path.join(path,"interface.o")
    exefile = os.path.join(path,"fortran_worker")
    
    fortran_compile(codefile, codestring, extra_args=extra_fflags)
    
    uc = create_fortran.GenerateAFortranSourcecodeStringFromASpecificationClass()
    uc.specification_class = interface
    uc.needs_mpi = needs_mpi
    string =  uc.result
    fortran_compile(interfacefile, string, extra_args=extra_fflags)
    fortran_build(exefile, [interfacefile, codefile], extra_args=extra_ldflags )
    return exefile
