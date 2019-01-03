try:
    import f90nml
    HAS_F90NML = True
except:
    HAS_F90NML = False
from collections import defaultdict
from amuse.units.quantities import new_quantity, to_quantity, is_quantity


# CodeWithNamelistParameters
#
# namelist_parameters=(
#   dict(name="name", group_name="name", short="codename", dtype="int32", default=64, description="description", ptype="nml" [, set_name="name"]), ...
# )

class CodeWithNamelistParameters(object):
    def __init__(self, namelist_parameters):
        self._namelist_parameters = dict([((x["short"], x["group_name"]), x) for x in namelist_parameters])

    def define_parameters(self, object):
        for p in self._namelist_parameters.values():
            if p["ptype"] in ["nml", "nml+normal"]:
                parameter_set_name = p.get("set_name", "parameters_" + p["group_name"])
                object.add_interface_parameter(p["name"], p["description"], p["default"], "before_set_interface_parameter", parameter_set=parameter_set_name)

    def read_namelist_parameters(self, inputfile):

        self._nml_file = inputfile
        self._nml_params = f90nml.read(inputfile)

        for group, d in self._nml_params.iteritems():
            for short, val in d.iteritems():
                key = (short, group.upper())
                if key in self._namelist_parameters:
                    group_name = self._namelist_parameters[key]["group_name"]
                    name = self._namelist_parameters[key]["name"]
                    parameter_set_name = self._namelist_parameters[key].get("set_name", "parameters_" + group_name)
                    parameter_set = getattr(self, parameter_set_name)
                    if is_quantity(self._namelist_parameters[key]["default"]):
                        setattr(parameter_set, name, new_quantity(val, to_quantity(self._namelist_parameters[key]["default"]).unit))
                    else:
                        setattr(parameter_set, name, val)
                else:
                    print "'%s' of group '%s' not in the namelist_parameters" % (short, group)

    def write_namelist_parameters(self, outputfile, do_patch=False, nml_file=None):
        patch = defaultdict(dict)
        for p in self._namelist_parameters.values():
            name = p["name"]
            group_name = p["group_name"]
            group = patch[group_name]
            short = p["short"]
            parameter_set_name = p.get("set_name", "parameters_" + group_name)
            parameter_set = getattr(self, parameter_set_name)
            if getattr(parameter_set, name) is None:  # omit if value is None
                continue
            if is_quantity(p["default"]):
                group[short] = to_quantity(getattr(parameter_set, name)).value_in(p["default"].unit)
            else:
                group[short] = getattr(parameter_set, name)

        if do_patch:
            f90nml.patch(nml_file or self._nml_file, patch, outputfile)
        else:
            f90nml.write(patch, outputfile, force=True)
