//
// star_io.C
//

#include "star.h"
#include "util_io.h"

istream & star::scan_star_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
//    char compare_string[MAX_INPUT_LINE_LENGTH];

//    if (is_star_in_binary())
//       sscanf(compare_string,")Star");
//    else
//       sscanf(compare_string,")Binary");

    cerr << "\n\nUsing star_io version of scan_star_story()...\n\n";

    while(get_line(s,input_line), strcmp(END_STAR, input_line)){
//    while(get_line(s,input_line), strcmp(compare_string,input_line)){
	char keyword[MAX_INPUT_LINE_LENGTH];
	if (getequals(input_line, keyword))		// demand "="
	    add_story_line(star_story, input_line);
    }
    return s;
}

ostream& star::print_star_story(ostream& s,
				int short_output)	// default = 0
{
//    char star_string[MAX_INPUT_LINE_LENGTH];

//    if (is_star_in_binary())
//       sscanf(star_string,")Star");
//    else
//       sscanf(star_string,")Binary");

    put_story_header(s, STAR_ID);

    if (star_story)
        put_story_contents(s, *star_story);

    put_story_footer(s, STAR_ID);
    
    return s;
}







