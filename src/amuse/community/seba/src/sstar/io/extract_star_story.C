
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  extract_star_story.C: creates a star part for each body
 *.............................................................................
 *    version 1.0:  Apr 1993   Piet Hut, Steve McMillan, Jun Makino
 *.............................................................................
 *  This file includes the following non-local functions:
 *    addstar1 - creates a star part for each node
 *.............................................................................
 */
#include "main_sequence.h"
#include "util_io.h"

#ifndef TOOLBOX

void dump_line_text(ostream& str, story& s)
    {
#ifndef BAD_GNU_IO
    str << s.get_text() << "\n";
#else
    fprintf(stdout, "%s\n", s.get_text());
#endif
    }

void extract_line_text(char* type, real& t_cur,  real& t_rel,
		       real& m_rel, real& m_env,
		       real& m_core, real& co_core,
		       real& t_eff, real& l_eff,
		       real& p_rot, real& b_fld,
		       story& s)
    {

        char keyword[MAX_INPUT_LINE_LENGTH];
        char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];
        char line [MAX_INPUT_LINE_LENGTH];
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
            sscanf(line,"%*s%*s%s",type);
        }else if(!strcmp("T_cur",keyword)){
            sscanf(line,"%*s%*s%lf",&t_cur);
        }else if(!strcmp("T_rel",keyword)){
            sscanf(line,"%*s%*s%lf",&t_rel);
        }else if(!strcmp("M_rel",keyword)){
            sscanf(line,"%*s%*s%lf",&m_rel);
        }else if(!strcmp("M_env",keyword)){
            sscanf(line,"%*s%*s%lf",&m_env);
        }else if(!strcmp("M_core",keyword)){
            sscanf(line,"%*s%*s%lf",&m_core);
        }else if(!strcmp("M_COcore",keyword)){
            sscanf(line,"%*s%*s%lf",&m_core);
        }else if(!strcmp("T_eff",keyword)){
            sscanf(line,"%*s%*s%lf",&t_eff);
        }else if(!strcmp("L_eff",keyword)){
            sscanf(line,"%*s%*s%lf",&l_eff);
        }else if(!strcmp("P_rot",keyword)){       // Experimental extra
	    sscanf(line,"%*s%*s%lf",&p_rot);     // information for
        }else if(!strcmp("B_fld",keyword)){       // Only for neutron stars
  	    sscanf(line,"%*s%*s%lf",&b_fld);
        }

    }


//! currently this does not print metallicity Feb 19th 2009 ST
void extract_story_chapter(char* type, real& z, real& t_cur, real& t_rel,
                           real& m_rel, real& m_env, real& m_core,
			   real& co_core,
			   real& t_eff, real& l_eff,
			   real& p_rot, real& b_fld,
			   story& s) {


    if (!s.get_chapter_flag())
        {
        cerr << "extract_story_chapter: not a story\n";
        exit(1);
        }


    t_rel=-1;

    for (story * d = s.get_first_daughter_node(); d != NULL;
                                                  d = d->get_next_story_node())
      {
	
	
	if (d->get_chapter_flag())
	  extract_story_chapter(type, z, t_cur, t_rel,
				m_rel, m_env, m_core, co_core,
				t_eff, l_eff,
				p_rot, b_fld, *d);
	else
	  extract_line_text(type, t_cur, t_rel,
			    m_rel, m_env, m_core, co_core,
			    t_eff, l_eff,
			    p_rot, b_fld, *d);
      }

    if (t_rel<=0)
      t_rel = t_cur;
}

/*-----------------------------------------------------------------------------
 *  addstar1  -- for all particles, add a star part using "new star()".
 *-----------------------------------------------------------------------------
 */
void  addstar1(node * b, real t_rel, stellar_type type,
	       real mf, real rf, real tf) {		// defaults = 1
      node * bi;

      if (b->get_oldest_daughter() != NULL)
  	 for (bi=b->get_oldest_daughter(); bi != NULL;
	                                   bi=bi->get_younger_sister())
	    addstar1(bi, t_rel, type, mf, rf, tf);
      else {
         int id = b->get_index();
	 real z = cnsts.parameters(solar_metalicity); //added by ST on Feb 19th 2009
         real t_cur=0, t_rel=0, m_rel=1, m_env=0, m_core=0.01, co_core=0;
	 real t_eff=0, l_eff=0;
	 real p_rot=0, b_fld=0;
         real m_tot;
         char type_string[MAX_INPUT_LINE_LENGTH];
	 starbase * old_starbase = b->get_starbase();
         story * s = old_starbase->get_star_story();

         extract_story_chapter(type_string, z, t_cur, t_rel,
                               m_rel, m_env, m_core, co_core,
			       t_eff, l_eff,
			       p_rot, b_fld, *s);
         m_tot = m_env+m_core;
         stellar_type local_type = extract_stellar_type_string(type_string);

         if (local_type!=NAS) {		// is the star properly defined.
            type = local_type;
            single_star* new_star = new_single_star(type, id, z, t_cur, t_rel,
			 m_rel, m_tot, m_core, co_core, p_rot, b_fld,  b);

         }
         else {		// No star story present, at least no proper
			// definition for a single star.
            real t_cur=0,t_rel=0, r_eff=0, m_rel=1, m_core=0.01;
            real m_tot, p_rot=0, b_fld=0;
            starbase * old_starbase = b->get_starbase();
            id = b->get_index();
            m_rel = m_tot = b->get_mass()/mf;

	    // Treat by dynamics pre-requested black holes
	    if(getiq(b->get_log_story(), "black_hole")==1) {
	      single_star* new_star = new_single_star(Black_Hole, id, z,
						      t_cur, t_rel,
						      m_rel, m_tot, m_tot,
						      m_tot,
						      p_rot, b_fld, b);
	    }
	    else {
	      single_star* new_star = new_single_star(type, id, z,
						    t_cur, t_rel,
						    m_rel, m_tot, m_core,
						    co_core,
						    p_rot, b_fld, b);

//            single_star* new_star = new_single_star(type, id, t_cur, t_rel,
//                         m_rel, m_tot, m_core, co_core, p_rot, b_fld, b);

	  }
	}
     }
  }

#  define  FALSE  0
#  define  TRUE   1

#else

/*-----------------------------------------------------------------------------
 *  main  --
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    int  c;
    bool  t_flag = FALSE;
    bool  M_flag = FALSE;
    bool  R_flag = FALSE;
    bool  T_flag = FALSE;
    bool  s_flag = FALSE;
    bool  c_flag = FALSE;
    real  mf = 1;
    real  rf = 1;
    real  tf = 1;
    real  t_rel = 0;
    stellar_type type = Main_Sequence;

    char  *comment;
    extern char *poptarg;

    while ((c = pgetopt(argc, argv, "M:R:T:t:s:c:",
		    "$Revision: 1.5 $", _SRC_)) != -1)
	switch(c)
	    {
            case 'M': M_flag = TRUE;
                      mf = atof(poptarg);
                      break;
            case 'R': R_flag = TRUE;
                      rf = atof(poptarg);
                      break;
            case 'T': T_flag = TRUE;
                      tf = atof(poptarg);
                      break;
	    case 't': T_flag = TRUE;
		      t_rel = atof(poptarg);
		      break;
            case 's': s_flag = TRUE;
                      type = (stellar_type)atoi(poptarg);
                      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': cerr <<
		      "usage: extract_star_story [-t #] [-T #] [-c \"..\"]\n";
		      exit(1);
	    }

    node *b;

    while (b = get_node()) {
        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

	addstar1(b, t_rel, type, mf, rf, tf);

	put_node(b);
	delete b;
    }
}

#endif
