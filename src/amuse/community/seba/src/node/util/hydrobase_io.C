
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// hydrobase_io.C
//

#include "hydrobase.h"
#include "util_io.h"

istream & hydrobase::scan_hydro_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while(get_line(s,input_line), !matchbracket(END_HYDRO, input_line)) {
	char keyword[MAX_INPUT_LINE_LENGTH];
	if (getequals(input_line, keyword))		// demand "="
	    add_story_line(hydro_story, input_line);
    }
    return s;
}

ostream& hydrobase::print_hydro_story(ostream& s)
{
    put_story_header(s, HYDRO_ID);

    if (hydro_story)
        put_story_contents(s, *hydro_story);

    put_story_footer(s, HYDRO_ID);
    
    return s;
}
