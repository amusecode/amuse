supportrc=dict(framework_install=True, package_name="amuse")

def use(arg):
    if arg == "package":
        supportrc["framework_install"]=True
    else:
        if arg not in ["system","installed","environment"]:
            warnings.warn(" assuming framework already installed")
        supportrc["framework_install"]=False

def set_package_name(arg):
    supportrc["package_name"]=arg
