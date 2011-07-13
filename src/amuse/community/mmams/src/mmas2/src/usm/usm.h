#ifndef _USM_
#define _USM_

#include "mass_shell.h"

class usm {
private:
  int num_shells;
  
  struct shell_ll {
    mass_shell shell;
    shell_ll *next;
  } *shells;
  mass_shell **shells_ptr;

public:
  
  real star_mass;         // in solar mass
  real star_radius;       // in solar radius
  real star_age;          // in Myr;

  usm() {
    num_shells = 0;
    shells = NULL;
    shells_ptr = NULL;
  }

  ~usm() {
    if (num_shells > 0) {
      if (shells_ptr != NULL) delete[] shells_ptr;
      
      shell_ll *ll_cell = shells;
      shell_ll *next_shell;
      do {
	next_shell = ll_cell->next;
	delete ll_cell;
	ll_cell = next_shell;
      } while (ll_cell != NULL);
    }
  }

  usm& operator=(const usm &in_model) {
      build_void_model();
      for (int i = 0; i < in_model.get_num_shells(); i++)
          add_shell(in_model.get_shell(i));
      
      if (num_shells > 0) build_hashtable();
      star_mass = in_model.star_mass;
      star_radius = in_model.star_radius;
      star_age = in_model.star_age;
      return *this;
  }
  
  int get_num_shells() const {return num_shells;}
  mass_shell& get_shell(int i) const {return *shells_ptr[i];}
  mass_shell& get_last_shell() const {return shells->shell;}

  void add_shell(mass_shell &shell);
  void build_hashtable();
  void build_void_model();

  void write(FILE*);
  void read(FILE*, int);
};

#endif //  _USM_
