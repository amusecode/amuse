# relative import hack
# https://stackoverflow.com/questions/16981921/relative-imports-in-python-3
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
