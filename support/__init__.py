supportrc=dict(framework_install=True, package_name="amuse", abort_on_build_failure=True)

def use(arg):
    if arg == "package":
        supportrc["framework_install"]=True
    else:
        if arg not in ["system","installed","environment"]:
            warnings.warn(" assuming framework already installed")
        supportrc["framework_install"]=False

def set_package_name(arg):
    supportrc["package_name"]=arg

def set_abort_on_build_failure(arg):
    supportrc["abort_on_build_failure"]=arg
