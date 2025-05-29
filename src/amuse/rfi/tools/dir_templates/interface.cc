#include "{code}_worker.h"

extern int echo(int input);

int echo_int(int input, int * output) {{
    *output = echo(input);
    return 0;
}}

