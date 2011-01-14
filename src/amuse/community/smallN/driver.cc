#include "interface.h"

#define N 3

main()
{
  // Test the multiple interface.  Store two binaries, then have them
  // interact with a third star.

  // Set up the interaction.

  int id[N] = {0, 0, 5};
  real mass[N] = {2, 2, 1};
  real pos[3][3] = {-1,-1,-1, 1,-1, 1, 0, 4, 0};
  real vel[3][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};

  id[0] = add_binary(1, 2, 1, 1, 0.5, 0.5);
  id[1] = add_binary(3, 4, 1, 1, 1, 0.25);

  clear_multiple();
  for (int k = 0; k < N; k++)
    add_to_interaction(id[k], mass[k], pos[k], vel[k]);

  // Integrate to completion and update all internal arrays.

  cerr << "initial configuration:" << endl;
  report_multiples(1);

  cerr << endl << "calling integrate_multiple" << endl;
  int n_new = integrate_multiple();

  cerr << endl << "final configuration:" << endl;

  // Retrieve the results.

  cerr << "top-level objects:" << endl;
  for (int k = 0; k < n_new; k++) {
    int i, isn;
    real m, x[3], v[3];
    get_particle(k, i, isn, m, x, v);
    vec pos(x[0],x[1],x[2]);
    PRC(i); PRC(m); PRC(isn); PRL(pos);
    if (is_multiple(i)) {
      PRI(2); PRC(get_mass(i)); PRL(get_radius(i));
    }
  }
  report_multiples(1);
}
