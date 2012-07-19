/*
 *  Bison parser specifying Starlab's dyn format grammar.
 *  Copyright (C) 2003  StarCluster team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


%{

  #include <vector>
  #include <stack>
  #include "dyn.h"

  extern char* yytext;
  enum section_ { P, D, H, S };
  stack<dyn*> p;
  stack<section_> section;
  unsigned long lineno = 1;
  vector<double>* v_ptr;
  struct { char* s; size_t n; } story_line;

  void yyerror(const string err) {
    cerr << err << " on line #" << lineno << " at `" << yytext << "'\n";
    exit(1);
  }

  int yylex(void);

%}

%union {
  double real;
  char* string;
  dyn* dyn_ptr;
  vector<double>* realvec_ptr;
}

%token PARTICLE LOG DYNAMICS HYDRO STAR
%token <string> KEYWORD LOG_STORY STRING
%token <real> NUMBER
%type <dyn_ptr> particle particles
%type <realvec_ptr> real_vec
// this grammar has two shift/reduce conflicts, the empty particle and particles
%expect 2


%%

particle:	{ /* $$ = $1 */ } | '(' PARTICLE {
  section.push(P);
  p.push(static_cast<dyn*>(new_dyn(new_hydrobase, new_starbase, true)));
}
		  key_val_pairs log dynamics hydro star particles ')' PARTICLE {
  if (p.size() > 1) $$ = p.top(), p.pop(); else return 0;
  section.pop();
};

particles:	  { $$ = 0; } | particles particle {
  $2->set_parent(p.top());
  if ($1) $1->set_younger_sister($2); else p.top()->set_oldest_daughter($2);
  $$ = $2;
};

log:		  '(' LOG log_story ')' LOG;

log_story:	| log_story LOG_STORY { if (*$2) p.top()->log_comment($2); free($2); };

dynamics:	  '(' DYNAMICS { section.push(D); }
		  key_val_pairs ')' DYNAMICS { section.pop(); };

hydro:		  '(' HYDRO { section.push(H); }
		  key_val_pairs ')' HYDRO { section.pop(); };

star:		  '(' STAR  { section.push(S); }
		  key_val_pairs ')' STAR { section.pop(); };

key_val_pairs:	| key_val_pairs key_val_pair;

key_val_pair:	  KEYWORD '=' real_vec {
  char line[strlen($1)+2+story_line.n+1];
  line[0] = '\0';
  strcpy(stpcpy(stpcpy(line, $1), " ="), story_line.s);
  switch (section.top()) {
  case P:
    if      (!strcmp($1, "i")) p.top()->set_index(int((*$3)[0]));
    else if (!strcmp($1, "N"));
    else add_story_line(p.top()->get_dyn_story(), line);
    break;
  case D:
    if      (!strcmp($1, "system_time")) p.top()->set_system_time((*$3)[0]);
    else if (!strcmp($1, "m")) p.top()->set_mass((*$3)[0]);
    else if (!strcmp($1, "r"))
      p.top()->set_pos(vec((*$3)[0],(*$3)[1],(*$3)[2]));
    else if (!strcmp($1, "v"))
      p.top()->set_vel(vec((*$3)[0],(*$3)[1],(*$3)[2]));
    else if (!strcmp($1, "a"))
      p.top()->set_acc(vec((*$3)[0],(*$3)[1],(*$3)[2]));
    else add_story_line(p.top()->get_dyn_story(), line);
    break;
  case H: add_story_line(p.top()->get_hydro_story(), line); break;
  case S: add_story_line(p.top()->get_star_story(), line); break;
  default: yyerror("inside unkown section");
  }
  free($1);
  $3->clear();
  *story_line.s = story_line.n = 0;
}
		| KEYWORD '=' STRING {
  char line[strlen($1)+3+strlen($3)+1];
  line[0] = '\0';
  stpcpy(stpcpy(stpcpy(line, $1), " = "), $3);
  switch (section.top()) {
  case P:
    if (!strcmp($1, "name")) p.top()->set_name($3);
    else add_story_line(p.top()->get_dyn_story(), line);
    break;
  case D: add_story_line(p.top()->get_dyn_story(), line); break;
  case H: add_story_line(p.top()->get_hydro_story(), line); break;
  case S: add_story_line(p.top()->get_star_story(), line); break;
  default: yyerror("inside unkown section");
  }
  free($1), free($3);
};

real_vec:	{ $$ = v_ptr; } | real_vec NUMBER { $1->push_back($2); };


%%


dyn* fget_dyn(FILE* fp) {

#ifdef STARLAB_USE_FDYN

  extern FILE* yyin;
  v_ptr = new vector<double>;
  yyin = fp;
  yyparse();
  dyn* root = 0;
  if (!p.empty()) root = p.top(), p.pop();
  delete v_ptr;
  return root;

#endif

}
