#ifndef  STARLAB_STORY_H
#  define  STARLAB_STORY_H

/// @file story.h Story class definitions.

#include  "starlab_vector.h"

//----------------------------------------------------------------------
//
// Define strings to start and end various story elements.

#if 1

// Old version:

#define PARTICLE_ID	"Particle"
#define START_PARTICLE	"(Particle"
#define END_PARTICLE	")Particle"

#define DYNAMICS_ID	"Dynamics"
#define START_DYNAMICS	"(Dynamics"
#define END_DYNAMICS	")Dynamics"

#define LOG_ID		"Log"
#define START_LOG	"(Log"
#define END_LOG		")Log"

#define HYDRO_ID	"Hydro"
#define START_HYDRO	"(Hydro"
#define END_HYDRO	")Hydro"

#define STAR_ID		"Star"
#define START_STAR	"(Star"
#define END_STAR	")Star"

#else

// New short (too short?) version:

#define PARTICLE_ID	"P"
#define START_PARTICLE	"(P"
#define END_PARTICLE	")P"

#define DYNAMICS_ID	"D"
#define START_DYNAMICS	"(D"
#define END_DYNAMICS	")D"

#define LOG_ID		"L"
#define START_LOG	"(L"
#define END_LOG		")L"

#define HYDRO_ID	"H"
#define START_HYDRO	"(H"
#define END_HYDRO	")H"

#define STAR_ID		"S"
#define START_STAR	"(S"
#define END_STAR	")S"

#endif

// Note that old and new versions are *not* compatible (Steve, 10/00).
//
//----------------------------------------------------------------------

#define MAX_STORY_LINE_LENGTH 255
#define chapter_begin_char  '('
#define chapter_end_char    ')'

/// \a story: A tree-structured list of character strings.
class story
    {
    private:

        story* next_story_node;
        story* first_daughter_node;
	story* last_daughter_node;
	char* text;
	int  chapter_flag;   // 0: this story is a line text
                             // 1: this story is a chapter

    public:

	/// Default is a bare story holding a single text line.

	/// A general story consists of a linked list of chapters,
	/// each containing a list of text lines.

	story(int flag = 0)  // default: a bare story, to hold a text line only
	    {
	    next_story_node = first_daughter_node = last_daughter_node = NULL;
	    text = NULL;
	    chapter_flag = flag;
	    }

	~story();

	/// Next line in the story, possibly in a new chapter.

	story * get_next_story_node()  {return next_story_node;}

	/// First line in this chapter.

	story * get_first_daughter_node()  {return first_daughter_node;}

	/// Last line in this chapter.

	story * get_last_daughter_node()  {return last_daughter_node;}

	/// Contents of the text line in this story node.

	char * get_text()  {return text;}

	/// Is this a chapter (1) or a text line(0)?

	int  get_chapter_flag()   {return chapter_flag;}

	/// Add to a story.

	void set_next_story_node(story * s)  {next_story_node = s;}

	/// Set first text line.

	void set_first_daughter_node(story * s)  {first_daughter_node = s;}

	/// Add/set last text line.

	void set_last_daughter_node(story * s)  {last_daughter_node = s;}

	/// Add text to a story node.

        void set_text(const char * a_string)
	    {
	    if(text != NULL)
	        delete [] text;
	    text = new char[strlen(a_string)+1];
	    strcpy(text, a_string);		  // Note that a_string
						  // isn't deleted here.
	    }

//	friend istream& operator>>(istream & , story * & );
//        friend ostream& operator<<(ostream & , story * );
    };

story* mk_story_line();
story* mk_story_line(const char *);
story* mk_story_chapter();
story* mk_story_chapter(const char *);

story* get_story(istream &);
story* get_story(istream &, const char *);
void put_story(ostream &, story &);
void put_story_contents(ostream &, story &, const char *prefix = NULL);
void put_simple_story_contents(ostream& str, story& s, const char *prefix = NULL);
void put_simple_story_contents(FILE *fp, story& s, const char *prefix = NULL);
void add_story_line(story *, const char *);
void rm_daughter_story(story * s, story * d);

story * find_qmatch(story *, const char *);
int  rmq(story *, const char *);

int  getiq(story *, const char *, bool verbose=false);
unsigned long getulq(story *, const char *, bool verbose=false);
unsigned long long getullq(story *, const char *, bool verbose=false);
real  getrq(story *, const char *, bool verbose=false);
char *getsq(story *, const char *, bool verbose=false);
vec  getvq(story *, const char *, bool verbose=false);
void getra(story *, const char *, real *, int, bool verbose=false);
void getia(story *, const char *, int *, int, bool verbose=false);
void getia(story *, const char *, unsigned long *, int, bool verbose=false);
void getia(story *, const char *, unsigned long long *, int, bool verbose=false);

void putiq(story *, const char *, int);
void putiq(story *, const char *, unsigned long);
void putiq(story *, const char *, unsigned long long);
void putrq(story *, const char *, real, int precision = STD_PRECISION);
void putra(story *, const char *, real *, int);
void putia(story *, const char *, int *, int);
void putia(story *, const char *, unsigned long *, int);
void putia(story *, const char *, unsigned long long *, int);
void putsq(story *, const char *, const char *);
void putvq(story *, const char *, vec &, int precision = STD_PRECISION);

void dump_story(story* s, int indent = 0);

#endif
 
