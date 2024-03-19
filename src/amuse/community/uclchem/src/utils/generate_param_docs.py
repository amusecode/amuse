"""Simple script to parse defaultparameters.f90 and generate a markdown
file that can be copied over to the website docs. Run it with
python generate_param_docs.py default_parameter_fortran_file output_markdown_file
"""
from sys import argv

param_file=argv[1]
output_file=argv[2]
with open(param_file,"r") as f:
    with open(output_file,"w") as output:
        for line in f.readlines():
            if line.startswith("!"):
                if not line.startswith("!!"):
                    output.write(line.replace("!",""))
            else:
                if "=" in line:
                    new_line=line.split("=")
                    line="|"+new_line[0]+"|"
                    new_line=new_line[1].split("!")
                    line=line+new_line[0]+"|"+new_line[1].strip()+"|\n"
                    output.write(line)