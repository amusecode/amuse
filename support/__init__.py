supportrc=dict(install_mode="framework", package_name="amuse", allow_build_failures='some')

def use(arg):
    if arg in ["framework_install_mode"):
        supportrc["install_mode"]="framework"
    elif arg in ["core_install_mode"]:
        supportrc["install_mode"]="core"
    else:
        if arg not in ["package_install_mode"]:
            warnings.warn("package_install_mode: assuming framework already installed")
        supportrc["install_mode"]="package"

def set_package_name(arg):
    supportrc["package_name"]=arg

def set_allow_build_failures(arg):
    if arg=="yes" or (arg==True): arg='some'
    if arg=="no" or (arg==False): arg='none'
    supportrc["allow_build_failures"]=arg
