#include "stdinc.h"

//// Test Starlab run-time help function
////
//// Options:
////       --help     print this message
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.


#ifndef TOOLBOX

local char* stredit(const char* s,
		    char c1, char c2)	// duplicate in kira_init.C
{
    if (!s) return NULL;

    char* s1 = new char[strlen(s)+1];
    if (!s1) return NULL;

    strcpy(s1, s);

    char* ss = s1;
    while (*ss != '\0') {
	if (*ss == c1) *ss = c2;
	ss++;
    }

    return s1;
}

void get_runtime_help(const char* source_file,
		      const char* date, const char *time, int level)
{
    // Extract help information from the specified source file.
    // source_file was the location of the file in question at the
    // time the program was compiled, but we should allow for the
    // possibilities that the name has changed or the file has
    // become inaccessible.

    cerr << endl
	 << "    Starlab version " << VERSION << endl;

    char* sd = stredit(date, '_', ' ');		// unnecessary now
    char* st = stredit(time, '_', ' ');		// unnecessary now
    if (sd && st) {
	cerr << "    program created on " << sd << " at " << st << endl;
	delete sd;
	delete st;
    }

    // Look for the source file.  First try the name verbatim, then
    // try prepending the current STARLAB_PATH environment string.

    char src[1024];
    bool exists = false;

    strcpy(src, source_file);
    ifstream file(src);

    if (file) {

	exists = true;
	file.close();

    } else {

	// See if we can construct a file name from the environment.

	if (getenv("STARLAB_PATH")) {
	    sprintf(src, "%s/%s", getenv("STARLAB_PATH"), source_file);
	    ifstream file(src);
	    if (file) {
		exists = true;
		file.close();
	    }
	}
    }

    cerr << "    source file ";
    if (!exists) {
	cerr << "unavailable" << endl;
	exit(0);
    } else
	cerr << "    " << src << endl << endl;

    // First, find and print out level-1 help lines.  Reformatting by
    // Steve (10/04) to let help2man turn the output into man pages.

    char cmd[1024];
    strcpy(cmd, "grep '^////' ");
    strcat(cmd, src);
    //    strcat(cmd, " | sed s%////%\"   \"%");
    strcat(cmd, " | sed s%////\\ %% | sed s%////%%");

    // Assume that help lines begin with "//// " (level 1)
    // or "//++ " (level 2).

//    PRL(cmd);
    system(cmd);
    cout<< endl;

    // Now check for level-2 help.

    if (level > 1) {
	strcpy(cmd, "grep '^//++' ");
	strcat(cmd, src);
	strcat(cmd, " | sed s%//++%\"    + \"%");

	system(cmd);
	cerr << endl;
    }

    exit(0);
}

void check_runtime_help(int argc, char** argv,
			const char* source_file,
			const char* date, const char *time)
{
    int help_level = 0;

    for (int i = 1; i < argc; i++) {
	if (strstr(argv[i], "-help")) help_level++;
	if (strstr(argv[i], "-HELP")) help_level = 2;
    }

    if (help_level) get_runtime_help(source_file, date, time, help_level);
}

#else

main(int argc, char** argv)
{
    check_help();
    extern char *poptarg;
    int c;
    const char *param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1) {}
}

#endif
