
/// @file util_io.h  Dyn I/O functions.

#define MAX_INPUT_LINE_LENGTH 255

int get_line(istream & s, char * line);
int check_input_line(istream &s, const char* reference_string);
int check_and_skip_input_line(istream &s, const char* reference_string);
int get_data_line(istream & s,char * input_line);
void set_vector_from_input_line(vec & v, char * input_line);
void set_vector_from_string(vec & v, char *str);
xreal get_xreal_from_input_line(char * input_line);
int matchbracket(const char *bracket, const char *line); // ("(Particle", line)
						         // matches "(P" or "(Particle"
const char *getequals(const char *line, char *keyword); // demand "keyword = value"

bool use_short_story_keywords( bool useshort );	// controls put_story_{head,foot}er

real read_unformatted_real( istream & s );
void read_unformatted_vector( istream & s, vec & v );
real read_unformatted32_real( istream & s );
void read_unformatted32_vector( istream & s, vec & v );

void write_unformatted_real( ostream & s, real r );
void write_unformatted_vector( ostream & s, vec & v );
void write_unformatted32_real( ostream & s, real r );
void write_unformatted32_vector( ostream & s, vec & v );

// Kludges to make GNU/HP output work:

void put_story_header(ostream & s, const char * id);
void put_story_footer(ostream & s, const char * id);
void put_real_number(ostream & s, const char *, xreal x);
void put_real_number(ostream & s, const char *, real x);
void put_real_vector(ostream & s, const char *, vec v);
void put_integer(ostream & s, const char *, int i);
void put_string(ostream & s, const char *, const char * str);

 
