#include "usm.h"

void usm::add_shell(mass_shell& shell) {
  shell_ll *new_shell = new shell_ll;
  new_shell->next = shells;
  new_shell->shell = shell;
  shells = new_shell;
  num_shells++;
}

void usm::build_hashtable() {
  shell_ll *ll_cell;
  if (num_shells == 0) {
    cerr << " usm::build_hashtable - no shells found " << endl;
    exit(-1);
  }

  if (shells_ptr != NULL) delete[] shells_ptr;
  shells_ptr = new mass_shell*[num_shells];

  ll_cell = shells;
  
  int i = 1;
  do {
    mass_shell &shell = ll_cell->shell;
    int j = num_shells - i;
    ll_cell = ll_cell->next;
    i++;
    shell.id = j;
    shells_ptr[j] = &shell;
  } while (ll_cell != NULL);
}

void usm::write(FILE *fout) {
  fprintf(fout, "(usm\n");

  writei(num_shells, "num_shells");
  writef(star_mass, "mass");
  writef(star_radius, "radius");
  writef(star_age, "age");

  for (int i = 0; i < num_shells; i++) shells_ptr[i]->write(fout);

  fprintf(fout, ")usm\n");
}


void usm::read(FILE *fin, int self_check) {
  char line[LINESIZE], variable[128], equal[10], value1[128], value2[128], value3[128];

  if (self_check == 1) {
    fgets(line, LINESIZE, fin);
    sscanf(line, "%s", variable);
    if (strcmp(variable, "(usm") != 0) {
      cerr << " --- Error in reading usm structure --- " << endl;
      cerr << "It seems for me that the file I'm trying to read is wrong one ... I am giving up ... bye bye " << endl;
      sqrt(-1.0); exit(-1);
    }
  }

  int nshells = 0;
  num_shells = 0;
  int failed = 0;
  while(1) {
    fgets(line, LINESIZE, fin);
    sscanf(line, "%s %s  %s %s %s", variable, equal,  value1, value2, value3);
    
    readi(nshells, "num_shells");
    readf(star_mass, "star_mass");
    readf(star_radius, "star_radius");
    readf(star_mass, "mass");
    readf(star_radius, "radius");
    readf(star_age, "star_age");
    readf(star_age, "age");

    check("(mass_shell") {
      mass_shell shell;
      shell.read(fin, 0);
      if (num_shells > 2) {
	
// 	fprintf(stderr, "num_shells= %d: ", num_shells);
// 	fprintf(stderr, "mnew= %lg  mold= %lg\n",
// 		shell.mass, shells->shell.mass);
	if (shell.mass > shells->shell.mass)
	  add_shell(shell);
	else
	  failed++;

      }	else {
	add_shell(shell);
      }
    }

    check(")usm") break;
  }

  build_hashtable();

  if (nshells != num_shells+failed) {
    cerr << "usm::read - incosistency found ... quit " << endl;
    PRL(nshells);
    PRL(num_shells);
    sqrt(-1.0); exit(-1);
  }
  
  if (failed > 0) {
    fprintf(stderr, "  *** WARNING *** \n");
    fprintf(stderr, "few mass shells do not have monotically increasing mass\n");
    fprintf(stderr, "nshells= %d  num_read= %d  n_failed= %d\n", nshells, num_shells, failed);
  }
}

/* this routine build a void stellar model */
void usm::build_void_model() {
  if (shells != NULL) {
    cerr << "usm::build_void_model() - cannot build a void model, because it is occupied with some other model \n";
    cerr << " if a void model is required please use: \n";
    cerr << "    usm void_model; void_model.build_void_model(); \n";
  }
}
