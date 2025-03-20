from amuse.support.import_helper import load_code


MESAInterface = load_code("mesa_r15140", "MESAInterface")
MESA = load_code("mesa_r15140", "MESA")

Mesa = MESA
