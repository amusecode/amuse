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


def wait_for_file(filename):
    for dt in [0.01, 0.01, 0.02, 0.05]:
        if os.path.exists(filename):
            return
        time.sleep(dt)


def c_compile(objectname, string):
    root, ext = os.path.splitext(objectname)
    sourcename = root + '.c'

    if os.path.exists(objectname):
        os.remove(objectname)

    with open(sourcename, "w") as f:
        f.write(string)

    mpicc = self.get_mpicc_name()
    arguments = [mpicc]
    arguments.extend(self.get_mpicc_flags().split())
    arguments.extend(["-I", "lib/stopcond", "-c",  "-o", objectname, sourcename])

    process = subprocess.Popen(
        arguments,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout, stderr = process.communicate()

    if process.returncode == 0:
        self.wait_for_file(objectname)

    if process.returncode != 0 or not os.path.exists(objectname):
        print "Could not compile {0}, error = {1}".format(objectname, stderr)
        raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))


def cxx_compile(objectname, string):
    root, ext = os.path.splitext(objectname)
    sourcename = root + '.cc'

    if os.path.exists(objectname):
        os.remove(objectname)

    with open(sourcename, "w") as f:
        f.write(string)

    mpicxx = self.get_mpicxx_name()
    arguments = [mpicxx]
    arguments.extend(self.get_mpicxx_flags().split())
    arguments.extend(["-I", "lib/stopcond", "-c",  "-o", objectname, sourcename])

    process = subprocess.Popen(
        arguments,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout, stderr = process.communicate()

    if process.returncode == 0:
        self.wait_for_file(objectname)

    if process.returncode != 0 or not os.path.exists(objectname):
        print "Could not compile {0}, error = {1}".format(objectname, stderr)
        raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))

    print stdout
    print stderr


def c_build(exename, objectnames):
    if os.path.exists(exename):
        os.remove(exename)

    mpicxx = self.get_mpicxx_name()
    arguments = [mpicxx]
    arguments.extend(objectnames)
    arguments.append("-o")
    arguments.append(exename)

    if 'LIBS' in os.environ:
        libs = os.environ['LIBS'].split()
        arguments.extend(libs)

    process = subprocess.Popen(
        arguments,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout, stderr = process.communicate()

    if process.returncode == 0:
        self.wait_for_file(exename)

    if process.returncode != 0 or not (os.path.exists(exename) or os.path.exists(exename+'.exe')):
        print "Could not compile {0}, error = {1}".format(exename, stderr)
        raise Exception("Could not build {0}, error = {1}".format(exename, stderr))

    print stdout
    print stderr


def build_worker(path_to_results, specification_class):
    path = os.path.abspath(path_to_results)
    codefile = os.path.join(path, "code.o")
    headerfile = os.path.join(path, "worker_code.h")
    interfacefile = os.path.join(path, "interface.o")
    exefile = os.path.join(path, "c_worker")

    c_compile(codefile, codestring)

    uc = create_c.GenerateACHeaderStringFromASpecificationClass()
    uc.specification_class = specification_class
    uc.needs_mpi = False
    header = uc.result

    with open(headerfile, "w") as f:
        f.write(header)

    uc = create_c.GenerateACSourcecodeStringFromASpecificationClass()
    uc.specification_class = specification_class
    uc.needs_mpi = False
    code = uc.result

    cxx_compile(interfacefile, code)
    c_build(exefile, [interfacefile, codefile])
