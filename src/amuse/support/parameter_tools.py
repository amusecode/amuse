import f90nml
from collections import defaultdict
from omuse.units.quantities import new_quantity, to_quantity, is_quantity


# CodeWithNamelistParameters
# 
# namelist_parameters=dict(
#   parametername  =  dict(group_name="name", short="codename", dtype="int32", default=64, description="description", ptype="nml" [, set_name="name"]),
# )

class CodeWithNamelistParameters(object):
    def __init__(self, namelist_parameters):
        self._namelist_parameters=namelist_parameters
    
    def define_parameters(self,object):
        for name,p in self._namelist_parameters.iteritems():
            if p["ptype"] in ["nml", "nml+normal"]:
                parameter_set_name=p.get("set_name", "parameters_"+p["group_name"])
                object.add_interface_parameter( name, p["description"], p["default"], "before_set_interface_parameter", parameter_set=parameter_set_name)

    def read_namelist_parameters(self, inputfile):

        self._nml_file=inputfile
        self._nml_params = f90nml.read(inputfile)

        for group, d in self._nml_params.iteritems():
            for name, val in d.iteritems():
                if name in self._namelist_parameters:
                    group_name=self._namelist_parameters[name]["group_name"]
                    parameter_set_name=self._namelist_parameters[name].get("set_name", "parameters_"+group_name)
                    parameter_set=getattr(self, parameter_set_name)
                    if is_quantity(self._namelist_parameters[name]["default"]):
                        setattr(parameter_set, name, new_quantity(val, to_quantity(self._namelist_parameters[name]["default"]).unit) )
                    else:
                        setattr(parameter_set, name, val )
                else:
                    print "'%s' of group '%s' not in the namelist_parameters"%(name, group)

    def write_namelist_parameters(self, outputfile, do_patch=False, nml_file=None):
        patch=defaultdict( dict )
        for name, v in self._namelist_parameters.iteritems():
            group_name=v["group_name"]
            group=patch[group_name]
            short=v["short"]
            parameter_set_name=v.get("set_name", "parameters_"+group_name)
            parameter_set=getattr(self, parameter_set_name)
            if getattr(parameter_set, name) is None:  # omit if value is None
                continue
            if is_quantity(self._namelist_parameters[name]["default"]):
                group[short]=to_quantity(getattr(parameter_set, name)).value_in(self._namelist_parameters[name]["default"].unit)
            else:
                group[short]=getattr(parameter_set, name)
        
        if do_patch:
            f90nml.patch(nml_file or self._nml_file,patch,outputfile)
        else:
            f90nml.write(patch, outputfile, force=True)
