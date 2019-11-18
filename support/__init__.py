supportrc=dict(framework_install=True, package_name="amuse", allow_build_failures='some')

def use(arg):
    if arg == "package":
        supportrc["framework_install"]=True
    else:
        if arg not in ["system","installed","environment"]:
            warnings.warn(" assuming framework already installed")
        supportrc["framework_install"]=False

def set_package_name(arg):
    supportrc["package_name"]=arg

def set_allow_build_failures(arg):
    if arg=="yes" or (arg==True): arg='some'
    if arg=="no" or (arg==False): arg='none'
    supportrc["allow_build_failures"]=arg
