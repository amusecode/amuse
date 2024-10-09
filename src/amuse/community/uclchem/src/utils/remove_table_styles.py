"""Simple script to parse the markdown files created from Tutorial notebooks and
remove style tags. For some reason nbconvert to markdown produces style tags
that cannot be compiled by mdx. So we just remove them.
"""
from sys import argv
for mdx_file in argv[1:]:
    lines=[]
    ignore_lines=False
    with open(mdx_file,"r") as f:
        for line in f.readlines():
            #if a line has style in it, change the flag
            #this way, <style> makes the code start ignoring lines until </style>
            if "style" in line:
                ignore_lines= not ignore_lines
                #exception are styles embeded in other tags. Just keep the flag and don't ignore later lines
                if "<tr" in line:
                    lines.append("<tr>\n")
                    ignore_lines=False
            #if no style, just add the line
            else:
                if not ignore_lines:
                    lines.append(line)
    with open(mdx_file,"w") as f:
        f.writelines(lines)

