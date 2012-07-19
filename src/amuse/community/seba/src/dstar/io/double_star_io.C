//
// base_element_io.C
//

#include "double_star.h"
#include "util_io.h"

//#include "hdyn.h"

istream & double_star::scan_star_story(istream& s) {

    char input_line[MAX_INPUT_LINE_LENGTH];

    while(get_line(s,input_line), strcmp(END_STAR, input_line)){
	char keyword[MAX_INPUT_LINE_LENGTH];
	getequals(input_line, keyword);		// demand "=", do nothing??
    }

    return s;
}

#if 0
        real number;
    	if(0){   // trick to keep the if() statement, for the else below
        }else if(!strcmp("Type",keyword)){
            char type_string[MAX_INPUT_LINE_LENGTH];
            sscanf(input_line,"%*s%*s%s",type_string);
            set_bin_type(extract_binary_type_string(type_string));
        }else if(!strcmp("semi",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_semi(semi);
        }else if(!strcmp("ecc",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_eccentricity(number);
	}else{
	    add_story_line(star_story, input_line);
	}
    }
#endif

ostream& double_star::print_star_story(ostream& s,
				       int short_output) { // default = 0

    put_story_header(s, STAR_ID);
//    put_story_header(s, "Binary");
//cerr<<"addresses: "<<this<<" "<<primary<<" "<<secondary<<endl;
//cerr<<get_primary()<<" "<<get_secondary()<<endl;
//cerr<<primary->get_root()<<" "<<secondary->get_root()<<endl;
//cerr<<get_primary()->get_root()<<" "<<get_secondary()->get_root()<<endl;
//cerr<<get_secondary()->get_companion()<<" "<<get_primary()->get_companion()<<endl;


//    put_string(s, "  Type  =  ", type_string(get_bin_type()));
//    put_real_number(s, "  semi  =  ", get_semi());
//    put_real_number(s, "  ecc   =  ", get_eccentricity());
    
    if (star_story)
        put_story_contents(s, *star_story);

    put_story_footer(s, STAR_ID);

    return s;
}

void extract_line_text(binary_type& type, real& semi, real& ecc, story& s) {

        char keyword[MAX_INPUT_LINE_LENGTH];
        char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];
        char line [MAX_INPUT_LINE_LENGTH];
        char type_string[MAX_INPUT_LINE_LENGTH];

        strcpy(line, s.get_text());
        sscanf(line,"%s%s",keyword,should_be_equal_sign);
        if(strcmp("=",should_be_equal_sign)){
            cerr << "Expected '=', but got '"<< should_be_equal_sign <<"'\n";
            exit(1);
        }

        real number;
        if(0){   // trick to keep the if() statement, for the else below
       }else if(!strcmp("Type",keyword)){
            char str_tpe[MAX_INPUT_LINE_LENGTH];
            sscanf(line,"%*s%*s%s",type_string);
            type = extract_binary_type_string(type_string);
        }else if(!strcmp("semi",keyword)){
            sscanf(line,"%*s%*s%lf",&semi);
        }else if(!strcmp("ecc",keyword)){
            sscanf(line,"%*s%*s%lf",&ecc);
        }
    }

void extract_story_chapter(binary_type& type, real& sma, real& ecc, story& s) {

    if (!s.get_chapter_flag()) {
        cerr << "extract_story_chapter: not a story\n";
        exit(1);
    }
    for (story * d = s.get_first_daughter_node(); d != NULL;
                 d = d->get_next_story_node()) {
        if (d->get_chapter_flag())
            extract_story_chapter(type, sma, ecc, *d);
        else
            extract_line_text(type, sma, ecc, *d);
    }
 }


