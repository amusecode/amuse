import numpy

from ConfigParser import ConfigParser

try:
    import f90nml
    HAS_F90NML=True
except:
    HAS_F90NML=False
from collections import defaultdict
from amuse.units.quantities import new_quantity, to_quantity, is_quantity


# CodeWithNamelistParameters
# 
# namelist_parameters=(
#   dict(name="name", group_name="name", short="codename", dtype="int32", default=64, description="description", ptype="nml" [, set_name="name"]), ...
# )

class CodeWithNamelistParameters(object):
    def __init__(self, namelist_parameters):
        self._namelist_parameters=dict([((x["short"],x["group_name"]),x) for x in namelist_parameters])
    
    def define_parameters(self,handler):
        for p in self._namelist_parameters.values():
            if p["ptype"] in ["nml", "nml+normal"]:
                parameter_set_name=p.get("set_name", "parameters_"+p["group_name"])
                handler.add_interface_parameter( p["name"], p["description"], p["default"], "before_set_interface_parameter", parameter_set=parameter_set_name)

    def read_namelist_parameters(self, inputfile):

        self._nml_file=inputfile
        self._nml_params = f90nml.read(inputfile)

        for group, d in self._nml_params.iteritems():
            for short, val in d.iteritems():
                key=(short,group.upper())
                if key in self._namelist_parameters:
                    group_name=self._namelist_parameters[key]["group_name"]
                    name=self._namelist_parameters[key]["name"]
                    parameter_set_name=self._namelist_parameters[key].get("set_name", "parameters_"+group_name)
                    parameter_set=getattr(self, parameter_set_name)
                    if is_quantity(self._namelist_parameters[key]["default"]):
                        setattr(parameter_set, name, new_quantity(val, to_quantity(self._namelist_parameters[key]["default"]).unit) )
                    else:
                        setattr(parameter_set, name, val )
                else:
                    print "'%s' of group '%s' not in the namelist_parameters"%(short, group)

    def write_namelist_parameters(self, outputfile, do_patch=False, nml_file=None):
        patch=defaultdict( dict )
        for p in self._namelist_parameters.values():
            name=p["name"]
            group_name=p["group_name"]
            group=patch[group_name]
            short=p["short"]
            parameter_set_name=p.get("set_name", "parameters_"+group_name)
            parameter_set=getattr(self, parameter_set_name)
            if getattr(parameter_set, name) is None:  # omit if value is None
                continue
            if is_quantity(p["default"]):
                value=to_quantity(getattr(parameter_set, name)).value_in(p["default"].unit)
            else:
                value=getattr(parameter_set, name)            
            if isinstance(value,numpy.ndarray):
                value=list(value)  # necessary until f90nml supports numpy arrays
            group[short]=value
        
        if do_patch:
            f90nml.patch(nml_file or self._nml_file,patch,outputfile)
        else:
            f90nml.write(patch, outputfile, force=True)      

class CodeWithIniFileParameters(object):
    def __init__(self, inifile_parameters=dict()):
        self._inifile_parameters=dict([((x["name"],x["group_name"]),x) for x in inifile_parameters])
        self._optionxform=str
        
    def define_parameters(self, handler):
        _tmp=dict()
        for p in self._inifile_parameters.values():
            if p["ptype"] in ["ini", "ini+normal"]:
                parameter_set_name=p.get("set_name", p["group_name"])
                if parameter_set_name not in _tmp:
                    _tmp[parameter_set_name]=[ x.name for x in handler.definitions[parameter_set_name] ]
                if not p["name"] in _tmp[parameter_set_name]:  
                    handler.add_interface_parameter( p["name"], p["description"], p["default"], 
                                      "before_set_interface_parameter", parameter_set=parameter_set_name)
                    
        self.set_parameters()

    def set_parameters(self):
        for p in self._inifile_parameters.values():
              parameter_set_name=p.get("set_name", None) or p["group_name"]
              parameter_set=getattr(self, parameter_set_name)
              name=p["name"]
              value=p.get("value", None) or p["default"]
              setattr(parameter_set, name, value)
              
    def read_inifile_parameters(self, configfile):
        self._configfile=configfile
        parser=ConfigParser()
        parser.optionxform=self._optionxform
        parser.read(configfile)
        for section in parser.sections():
            group=section
            for option in parser.options(section):
                key=(option,group)
                if key in self._inifile_parameters:
                    ptype=self._inifile_parameters[key]["ptype"]
                    dtype=self._inifile_parameters[key]["dtype"]
                    value=self.interpret_value(parser.get(group, option), dtype=dtype)
                    if is_quantity(self._inifile_parameters[key]["default"]):
                        value= new_quantity(val, to_quantity(self._inifile_parameters[key]["default"]).unit)
                    self._inifile_parameters[key]["value"]=value
                else:
                    value=self.interpret_value(parser.get(group, option))
                    self._inifile_parameters[key]=dict(
                        group_name=group,
                        name=option,
                        set_name=group,
                        default=value,
                        value=value,
                        short=option,
                        ptype="ini",
                        dtype="unknown",
                        description="unknown parameter read from %s"%configfile
                        )

    def _convert(self, value, dtype):
        if dtype is "bool":
            if value.lower() in ["0", "false", "off","no"]:
                return False
            else:
                return True
        if dtype in ["str", None]:
            return value
        return numpy.fromstring(value, sep=",")[0]

    def interpret_value(self,value, dtype=None):
        if value.find(',')>=0:
            return [self._convert(x, dtype) for x in value.split(",")]
        return self._convert(value, dtype)

    def output_format_value(self,value):
        if isinstance(value, list):
          return ','.join(value)
        else:
          return value


    def write_inifile_parameters(self, outputfile):
        parser=ConfigParser()
        parser.optionxform=self._optionxform

        for p in self._inifile_parameters.values():
            name=p["name"]
            group_name=p["group_name"]
            short=p["short"]
            parameter_set_name=p.get("set_name", group_name)
            parameter_set=getattr(self, parameter_set_name)
            if is_quantity(p["default"]):
                value=to_quantity(getattr(parameter_set, name)).value_in(p["default"].unit)
            else:
                value=self.output_format_value(getattr(parameter_set, name))            
            if not parser.has_section(group_name):
                parser.add_section(group_name)
            parser.set(group_name,short,value)
        
        f=open(outputfile, "w")
        parser.write(f)
        f.close()

