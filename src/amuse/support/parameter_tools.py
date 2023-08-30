import numpy

from configparser import ConfigParser
from collections import defaultdict, OrderedDict
from io import StringIO

try:
    import f90nml
    HAS_F90NML=True
except:
    HAS_F90NML=False

from amuse.units.quantities import new_quantity, to_quantity, is_quantity

# parameters can be supplied as:
# 
# parameters=(
#   dict(name="name", group_name="name", short="codename", dtype="int32", default=64, description="description", ptype="nml" [, set_name="name"]), ...
# )

dtype_str={ str: "str", 
            bool: "bool",
            int: "int32",
            float: "float64",
            complex: "complex",
            list: "list",
            dict: "dict",
            OrderedDict: "OrderedDict",
          }

def parameter_list_py_code(parameters, label="parameters"):
    header="""from omuse.units import units
from collections import OrderedDict

{label} = (

""".format(label=label)

    template='  dict(group_name={group}, name={name}, short={short}, dtype={dtype}, default={default}, description={description}, ptype={ptype}),\n'

    footer=""")
    """
    
    by_group=defaultdict(dict)
    for key in parameters:
        short,group=key
        by_group[group][short]=parameters[key]

    body=[]
    for group in by_group:
        _body=""
        for short in by_group[group]:
          _body+=template.format(name=by_group[group][short]["name"].__repr__(),
                                group=by_group[group][short]["group_name"].__repr__(),
                                short=by_group[group][short]["short"].__repr__(),
                                dtype=by_group[group][short]["dtype"].__repr__(),
                                default=by_group[group][short]["default"].__repr__(),
                                description=by_group[group][short]["description"].__repr__(),
                                ptype=by_group[group][short]["ptype"].__repr__(),)
        
        body.append(_body)
    body="#\n".join(body)

    return header+body+footer


class _CodeWithFileParameters(object):
    _ptypes=None
    def _write_file(self, inputfile, **kwargs):
        raise Exception("not implemented")
    def _read_file(self, inputfile, rawvals, **kwargs):
        raise Exception("not implemented")

    def define_parameters(self, handler):
        _tmp=dict()
        for p in self._parameters.values():
            if p["ptype"] not in self._ptypes:
                continue
            parameter_set_name=p.get("set_name", None) or self._prefix+p["group_name"].replace(" ","_")
            if parameter_set_name not in _tmp:
                _tmp[parameter_set_name]=[ x.name for x in handler.definitions[parameter_set_name] ]
            if not p["name"] in _tmp[parameter_set_name]:  
                handler.add_interface_parameter( p["name"], p["description"], p["default"], 
                                  "before_set_interface_parameter", parameter_set=parameter_set_name)

        self.set_parameters()

    def set_parameters(self):
        for p in self._parameters.values():
            if p["ptype"] not in self._ptypes:
                continue
            parameter_set_name=p.get("set_name", None) or self._prefix+p["group_name"].replace(" ","_")
            parameter_set=getattr(self, parameter_set_name)
            name=p["name"]
            value=p.get("value", p["default"])
            setattr(parameter_set, name, value)

    def interpret_value(self,value, dtype=None):
        raise Exception("not implemented")

    def read_parameters(self, inputfile, add_missing_parameters=False):
        self._file=inputfile

        _nml_params = f90nml.read(inputfile)

        rawvals, comments = self._read_file(inputfile)

        for key, rawval in rawvals.items():      
                if key in self._parameters:
                    group_name=self._parameters[key]["group_name"]
                    name=self._parameters[key]["name"]
                    dtype=self._parameters[key]["dtype"]
                    val=self.interpret_value( rawval, dtype=dtype)
                    if is_quantity(self._parameters[key]["default"]):
                        self._parameters[key]["value"]=new_quantity(val, to_quantity(self._parameters[key]["default"]).unit)
                    else:
                        self._parameters[key]["value"]=val 
                else:
                    if not add_missing_parameters:
                        print("'{0}' of group '{1}' not in the parameters list".format(*key))
                    else:
                        value=rawval
                        description=comments.get(key, "unknown parameter read from {0}".format(inputfile))
                        self._parameters[key]=dict(
                            group_name=key[1],
                            name=key[0],
                            short_name=key[0],
                            default=value,
                            value=value,
                            short=key[0],
                            ptype=self._ptypes[0],
                            dtype=dtype_str[type(value)],
                            description=description
                            )                        

    def write_parameters(self, outputfile, **options):

        rawvals=dict()

        for key, p in self._parameters.items():
            name=p["name"]
            group_name=p["group_name"]
            short=p["short"]
            parameter_set_name=p.get("set_name", None) or self._prefix+p["group_name"].replace(" ","_")
            parameter_set=getattr(self, parameter_set_name)
            if is_quantity(p["default"]):
                value=to_quantity(getattr(parameter_set, name)).value_in(p["default"].unit)
            else:
                value=getattr(parameter_set, name)
            
            rawvals[key]=self.output_format_value(value)
            
        self._write_file(outputfile, rawvals, **options)

class CodeWithNamelistParameters(_CodeWithFileParameters):
    """
    Mix-in class to 1) namelist file support to code interfaces and 2) automatically generate
    parameter sets from descriptions or namelist files. 
    
    This class takes a list of parameter descriptions (optional) and has functions to
    read and write namelist files. Every namelist section corresponds to a different 
    parameter set.
    """
    _ptypes=["nml", "nml+normal"]

    def __init__(self, _parameters, prefix="parameters_"):
        if not HAS_F90NML:
            raise Exception("f90nml package not available")
        self._parameters=dict([((x["short"].lower(),x["group_name"]),x) for x in _parameters])
        self._prefix=prefix
        self._file=None

    def _read_file(self, inputfile):
        _nml_params = f90nml.read(inputfile).todict()
        rawvals=OrderedDict()

        for group, d in _nml_params.items():
            for short, val in d.items():
                key=(short.lower(),group.upper())
                rawvals[key]=val

        return rawvals, dict()

    def _write_file(self, outputfile, rawvals, do_patch=False, nml_file=None):
        patch=OrderedDict()

        for key,rawval in rawvals.items():
            if rawval is None:  # omit if value is None
                continue
            if isinstance(rawval,numpy.ndarray):
                rawval=list(rawval)  # necessary until f90nml supports numpy arrays
            if key[1] not in patch:
              patch[key[1]]=OrderedDict()
            patch[key[1]][key[0]]=rawval
        
        if do_patch:
            _tmp=f90nml.read(nml_file or self._nml_file)
            _tmp.update(patch)
            f90nml.write(_tmp, outputfile, force=True)      
            # workaround because this can produce errors (f90nml 1.1.2):
            #~ f90nml.patch(nml_file or self._nml_file,f90nml.Namelist(patch),outputfile)
        else:
            f90nml.write(patch, outputfile, force=True)      

    def write_namelist_parameters(self, outputfile, do_patch=False, nml_file=None):
        return self.write_parameters(outputfile, do_patch=do_patch, nml_file=nml_file)

    def read_namelist_parameters(self, inputfile, add_missing_parameters=False):
        return self.read_parameters(inputfile,add_missing_parameters)

    def output_format_value(self,value):
        return value

    def interpret_value(self,value, dtype=None):
        return value   # dtype, arrays should be handled by f90nml 


class CodeWithIniFileParameters(_CodeWithFileParameters):
    """
    Mix-in class to 1) INI-like file support to code interfaces and 2) automatically generate
    parameter sets from descriptions or Ini files. 
    
    This class takes a list of parameter descriptions (optional) and has functions to
    read and write INI files. Every section corresponds to a different parameter set.
    """
    _ptypes=["ini", "ini+normal"]
    def __init__(self, _parameters=None, prefix="ini_"):
        if _parameters is None:
          _parameters=dict()
        self._parameters=dict([((x["name"],x["group_name"]),x) for x in _parameters])
        self._optionxform=str
        self._prefix=prefix
        self._file=None
                      
    def _read_file(self, inputfile):
        f=open(inputfile,"r")
        values=StringIO()
        comments=StringIO()
        for line in f.readlines():
            if "=" in line:
              key, val=line.split("=", 1)
              if "#" in val:
                val,comment=val.split("#",1)
                comments.write("=".join([key,comment])+"\n")
              values.write("=".join([key,val])+"\n")
            else:
              values.write(line+"\n")
              comments.write(line+"\n")

        values.seek(0)
        comments.seek(0)

        rawvals=self.parse_fp(values)
        comments=self.parse_fp(comments)

        return rawvals, comments


    def parse_fp(self, fp):
        parser=ConfigParser()
        parser.optionxform=self._optionxform
        parser.readfp(fp)

        rawvals=dict()
        for section in parser.sections():
            group=section
            for option in parser.options(section):
                key=(option,group)

                rawvals[key]=parser.get(group, option)
                
        return rawvals

    def _convert(self, value, dtype):
        if dtype=="bool":
            if value.lower() in ["0", "false", "off","no"]:
                return False
            else:
                return True
        if dtype in ["str", None]:
            return value
        return numpy.fromstring(value, sep=",")[0]

    def interpret_value(self,value, dtype=None):
        if value=="":
            return value
        if value.find(',')>=0:
            return [self._convert(x, dtype) for x in value.split(",")]
        return self._convert(value, dtype)

    def _write_file(self, outputfile, rawvals):
        parser=ConfigParser()
        parser.optionxform=self._optionxform

        for key, rawval in rawvals.items():
            section=key[1]
            short=key[0]
            
            if not parser.has_section(section):
                parser.add_section(section)

            if isinstance(rawval, list):
                rawval=','.join(rawval)
            parser.set(section,short,self.output_format_value(rawval))

        f=open(outputfile, "w")
        parser.write(f)
        f.close()

    def output_format_value(self,value):
        if isinstance(value, list):
          return ','.join([str(v) for v in value])
        else:
          return str(value)
        
    def write_inifile_parameters(self, outputfile):
        return self.write_parameters(outputfile)

    def read_inifile_parameters(self, inputfile, add_missing_parameters=False):
        return self.read_parameters(inputfile,add_missing_parameters)

