import sys
# from interface import GENEC_STAR_PARAMETERS as INPUT
from interface import GENEC_STAR_PROPERTIES as INPUT
from amuse.community import NO_UNIT

def main():
    for star_property in INPUT:
        dtype = INPUT[star_property][0]
        if dtype in ["i", "int32"]:
            f_dtype = "integer"
        elif dtype in ["d", "float64"]:
            f_dtype = "real(kindreal)"
        elif dtype in ["b", "bool"]:
            f_dtype = "logical"
        elif dtype in ["s", "string"]:
            f_dtype = "character(256)"
        else:
            print(f"wrong dtype: {dtype}")
            sys.exit()
        if INPUT[star_property][1] != '':
            dtype = f"'{INPUT[star_property][0]}' | units.{INPUT[star_property][1]}"
        else:
            dtype = f"'{INPUT[star_property][0]}'"
        with open("interface_temp.py", 'a') as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{star_property}(index_of_the_particle='i'):\n"
                f"        returns ({star_property}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{star_property}(index_of_the_particle='i', {star_property}={dtype}):\n"
                f"        returns ()\n"
                f"\n"
            )

        with open("interface_temp.f90", 'a') as f_out:
            f_out.write(
                f"integer function get_{star_property}(index_of_the_particle, {star_property})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    {f_dtype}, intent(out):: {star_property}\n"
                f"    {star_property} = GenecStar%{star_property}\n"
                f"    get_{star_property} = 0\n"
                f"end function get_{star_property}\n"
                f"\n"
                f"integer function set_{star_property}(index_of_the_particle, {star_property})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    {f_dtype}, intent(in):: {star_property}\n"
                f"    GenecStar%{star_property} = {star_property}\n"
                f"    set_{star_property} = 0\n"
                f"end function set_{star_property}\n"
                f"\n"
            )


if __name__ == "__main__":
    main()
