import amuse

import sys

print(sys.path)

def test():
  amuse.config.compilers.fc_iso_c_bindings

if __name__=="__main__":
  test()
