
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Convert an extended version of ASCII "dumbp" format
////
////            (id1, mass1, pos1, vel1, var11, var21, var31,...
////             id2, mass2, pos2, vel2, var12, var22, var32,...
////             id3, mass3, pos3, vel3, var13, var23, var33,...
////             etc.)
////
//// data into a Starlab snapshot (flat tree).  This is the inverse
//// function to dumbp.  It differs from readp in that up to 128 extra
//// variables are stored in the dyn story as "var1 = xxx" etc.  The
//// labels actually used may be set on the command line.
////
//// Usage:  readp2 [OPTIONS]
////
//// Options:
////              -c    add a comment to the output snapshot [false]
////              -i    number the particles sequentially [don't number]
////              -v    specify name for extra variable (e.g. -v 2 temp)
////                    (note that numbering starts at 1)
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//	     Steve McMillan, June 2008

#include "dyn.h"

#ifdef TOOLBOX

#define MAX_LINE 1024
#define MAX_DATA  128

// Unabashed ~C code!

local bool is_space(char *s)
{
  if (!s) return false;
  return (*s == '\t' || *s == ' ');		// space characters: '\t', ' '
}

local bool is_term(char *s)
{
  if (!s) return true;				// terminating characters:
  return (*s < ' ' && *s != '\t');		// anything < ' ' except '\t'
}

char *current_word(char *s)
{
  if (!s) return NULL;
  while (is_space(s)) s++;			// skip whitespace
  if (is_term(s)) return NULL;
  return s;
}
 
char *next_word(char *s)
{
  if (!s) return NULL;
  if (!is_space(s))
    while (!is_space(s) && !is_term(s)) s++;	// go to end of current word
  return current_word(s);
}

int get_data(int &i, real &m, vec &pos, vec &vel, real data[])
{
  char input[MAX_LINE], *s = input;
  int count = 0;
  s = fgets(input, MAX_LINE, stdin);		// NULL terminated, w/newline
  if (s) {
    s = current_word(s);
    while (s) {
      real x;
      int j;
      if (count == 0)
	j = sscanf(s, "%d", &i);
      else
	j = sscanf(s, "%lf", &x);
      count += j;
      if (j > 0) {
	if (count == 1)
	  ;
	else if (count == 2)
	  m = x;
	else if (count <= 5)
	  pos[count-3] = x;
	else if (count <= 8)
	  vel[count-6] = x;
	else
	  data[count-9] = x;
      } else
	break;
      s = next_word(s);
    }
  }
  return count;					// return # of variables read
}

int main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    const char *param_string = "c:iv::";

    char *comment;
    bool c_flag = false;
    bool i_flag = false;

    char *varname[MAX_DATA], temp[32];
    int k;

    for (k = 0; k < MAX_DATA; k++) {
      sprintf(temp, "var%d", k+1);
      varname[k] = (char *)malloc(strlen(temp)+1);
      strcpy(varname[k], temp);
    }

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.3 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'v': k = atoi(poparr[0]);
		      if (k >= 1 && k <= MAX_DATA)
			strcpy(varname[k-1], poparr[1]);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    // Circumvent the normal input function, since we want to store
    // any extra data found on each line.  Build the tree from scratch.

    dyn *b, *bb, *bp = NULL;
    b = new dyn();
    b->set_root(b);

    int i, n;
    real mass;
    vec pos, vel;
    real vardata[MAX_DATA];

    while ((n = get_data(i, mass, pos, vel, vardata)) > 0) {
      bb = new dyn();
      bb->set_parent(b);
      if (!bp)
	b->set_oldest_daughter(bb);
      else {
	bb->set_elder_sister(bp);
	bp->set_younger_sister(bb);
      }
      bp = bb;

      if (i_flag)
	bb->set_label(1);
      else
	bb->set_label(i);
      bb->set_mass(mass);
      bb->set_pos(pos);
      bb->set_vel(vel);

      // The returned value n is the total number of variables read
      // from the line.  The first 8 are "known": ID, mass, pos, vel.
      // The 9th and above are extra.

      for (k = 0; k < n-8; k++)
	putrq(bb->get_dyn_story(), varname[k], vardata[k]);
    }

    if (c_flag) b->set_col_output(true);
    if (c_flag) b->log_comment(comment);

    b->log_history(argc, argv);

    dyn::set_col_output(false);		// force dyn output (default is col)
    put_dyn(b);
}

#endif
