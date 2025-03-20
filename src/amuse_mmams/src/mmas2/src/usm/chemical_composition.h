#ifndef _CHEMICAL_COMPOSITION_H_
#define _CHEMICAL_COMPOSITION_H_

#include "include/stdinc.h"

class chemical_composition {
protected:
public:
  real H1;
  real He4;
  real O16;
  real N14;
  real C12;
  real Ne20;
  real Mg24;
  real Si28;
  real Fe56;

  chemical_composition(real Z = 0.02) {
    H1 = 0;
    He4 = 0;
    O16 = 0;
    N14 = 0;
    C12 = 0;
    Ne20 = 0;
    Mg24 = 0;
    Si28 = 0;
    Fe56 = 0;
  }
  ~chemical_composition() {
  }

  void operator=(chemical_composition &che) {
    H1 = che.H1;
    He4 = che.He4;
    O16 = che.O16;
    N14 = che.N14;
    C12 = che.C12;
    Ne20 = che.Ne20;
    Mg24 = che.Mg24;
    Si28 = che.Si28;
    Fe56 = che.Fe56;
  }

  void read(FILE *fin, int self_check) {
    char line[LINESIZE], variable[128], equal[10], value1[128], value2[128], value3[128];

    if (self_check == 1) {
      fgets(line, LINESIZE, fin);
      sscanf(line, "%s", variable);
      PRL(line);
      if (strcmp(variable, "(chemical_composition") != 0) {
	cerr << " I tried to read, but it seems for me that it is not a chemical_composition data " << endl;
	cerr << "       ..... I am sorry, but I have to give up. bye bye ... " << endl;
	PRL(sqrt(-1.0));
	exit(-1);
      }
    }

    while(1) {
      fgets(line, LINESIZE, fin);
      sscanf(line, "%s %s  %s %s %s", variable, equal,  value1, value2, value3);
      
      readf(H1, "H1");     readf(H1, "h1");
      readf(He4, "He4");   readf(He4, "he4");
      readf(O16, "O16");   readf(O16, "o16");
      readf(N14, "N14");      readf(N14, "n14");
      readf(C12, "C12");      readf(C12, "c12");
      readf(Ne20, "Ne20");      readf(Ne20, "ne20");
      readf(Mg24, "Mg24");      readf(Mg24, "mg24");
      readf(Si28, "Si28");      readf(Si28, "si28");
      readf(Fe56, "Fe56");       readf(Fe56, "fe56");


      check(")chemical_composition") break;
    }
    
  }

  void write(FILE *fout) {
    fprintf(fout, "(chemical_composition\n");

    writef(H1, "H1");
    writef(He4, "He4");
    writef(O16, "O16");
    writef(N14, "N14");
    writef(C12, "C12");
    writef(Ne20, "Ne20");
    writef(Mg24, "Mg24");
    writef(Si28, "Si28");
    writef(Fe56, "Fe56");
    
    fprintf(fout, ")chemical_composition\n");
  }
};

#endif //  _CHEMICAL_COMPOSITION_H_
