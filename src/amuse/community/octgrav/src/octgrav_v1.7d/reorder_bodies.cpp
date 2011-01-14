#include "octgrav.h"

#define BITS_PER_DIMENSION 18
typedef long long peanokey;

struct peano_struct {
  peanokey key;
  int  value;
};

static int quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};

static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static int rotx_table[8]  = { 3, 0, 0, 2, 2, 0, 0, 1 };
static int roty_table[8]  = { 0, 1, 1, 2, 2, 3, 3, 0 };

static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };

// static int flag_quadrants_inverse = 1;
// static char quadrants_inverse_x[24][8];
// static char quadrants_inverse_y[24][8];
// static char quadrants_inverse_z[24][8];

/*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
 *  with x,y,z in the range between 0 and 2^bits-1.
 */

peanokey peano_hilbert_key(int x, int y, int z, int bits) {
  int i, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;
  peanokey key;
  
  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;

  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;
      
      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }

  return key;
}


int compare_peanokey(const void *a, const void *b) {
  if(((struct peano_struct *) a)->key < (((struct peano_struct *) b)->key))
    return -1;
  
  if(((struct peano_struct *) a)->key > (((struct peano_struct *) b)->key))
    return +1;
  
  return 0;
}

void octgrav::reorder_bodies() {

  double t1 = get_time();
  int n_cells = cell_list.size();
  float4 r_min = cell_list[0]->get_pos();
  float4 r_max = r_min;
  for (int i = 0; i < n_cells; i++) {
    float4 pos = cell_list[i]->get_pos();
    r_min.x = min(r_min.x, pos.x);
    r_min.y = min(r_min.y, pos.y);
    r_min.z = min(r_min.z, pos.z);
    
    r_max.x = max(r_max.x, pos.x);
    r_max.y = max(r_max.y, pos.y);
    r_max.z = max(r_max.z, pos.z);
  }
  float size = max(r_max.z - r_min.z,
		   max(r_max.y - r_min.y, r_max.x - r_min.x));
 
  float domain_fac = 1.0 / size * (((peanokey)1) << (BITS_PER_DIMENSION));
  peano_struct *keys = (peano_struct*)malloc(n_cells * sizeof(peano_struct));

  for (int i = 0; i < n_cells; i++) {
    keys[i].value = i;
    float4 pos = cell_list[i]->get_pos();
    int x = (int)((pos.x - r_min.x) * domain_fac);
    int y = (int)((pos.y - r_min.y) * domain_fac);
    int z = (int)((pos.z - r_min.z) * domain_fac);
    keys[i].key   = peano_hilbert_key(x, y, z, BITS_PER_DIMENSION);
  }

  qsort(keys, n_cells, sizeof(peano_struct), compare_peanokey);

  vector<oct_tree<N_PER_CELL>*> sorted_cells(n_cells);
  for (int i = 0; i < n_cells; i++) {
    int j = keys[i].value;
    sorted_cells[i] = cell_list[j];
  }
  fprintf(stderr, " \n ---- %lg sec\n", get_time() - t1);

  
  int offset = 0;
  int pc = 0;
  for (int i = 0; i < n_cells; i++) {
    cell_list[i] = sorted_cells[i];
    int n1 = leaf_list.size();
    cell_list[i]->generate_leaf_list(leaf_list);
    int n2 = leaf_list.size();
    cell_list[i]->offset = offset;
//     fprintf(stderr, "n1= %d  n2= %d\n", n1, n2);
    for (int j = n1; j < n2; j++) {
      leaf_list[j]->offset = offset;
      vector<int> &bodies_id = leaf_list[j]->get_bodies_id();
      int k = leaf_list[j]->get_n_bodies();
      offset += k;
      while (k-- > 0) {
	int id = bodies_id[k];
	bodies_id_orig[pc]  = id;
	hst.bodies_pos[pc]  = bodies_pos[id];
	hst.bodies_grav[pc] = (float4){INT_AS_FLOAT(k), 0,0,0};
	pc++;
      }
    }
//     int offset1 = cell_list[i]->offset + cell_list[i]->get_total_bodies();
//     fprintf(stderr, "cell_offset= %d leaf_offset= %d\n",
// 	    offset1, offset);
  }

//   fprintf(stderr, "pc= %d  n_bodies =%d\n", pc, n_bodies);
//   double tot_mass = 0;
//   for (int i = 0; i < pc; i++) {
//     tot_mass += hst.bodies_pos[i].w;
//   }
//   fprintf(stderr, "tot_mass= %lg\n", tot_mass);
  free(keys);

}
