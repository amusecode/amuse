
// par_to_use.C: Turn the list of parameters used by pgetopt
//		 into a "Usage" statement.

// Name shortened to satisfy some library managers.

#include "stdinc.h"

void params_to_usage(ostream& s, char* name, const char *param_string)
{
    int i = 0;

    while (*param_string) {
	if (++i == 1) s << "Usage: " << name;

	if (*param_string == ':')
	    s << " #";
	else if (*param_string == '.')
	    s << " [#]";
	else {
	    if (i > 1) s << "]";
	    s << " [-" << *param_string;
	}

	param_string++;
    }

    if (i > 0) s << "]" << endl;
}
