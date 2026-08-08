/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/mesh/voronoi/voronoi.h
 * \date        05/2018
 * \brief       Header for Voronoi mesh-construcion
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 29.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef HAVE_H_VORONOI
#define HAVE_H_VORONOI

#include <gmp.h>

#define STACKSIZE_TETRA 10000
#define MIN_ALLOC_NUMBER 1000
#define ALLOC_INCREASE_FACTOR 1.1
#define ALLOC_DECREASE_FACTOR 0.7
#define MAX_VORONOI_ITERATIONS 500

#define GENTLE_DEREFINE_FACTOR 1.2

#define USEDBITS 52

#if USEDBITS > 31
typedef signed long long int IntegerMapType;
void MY_mpz_set_si(mpz_t dest, signed long long int val);
void MY_mpz_mul_si(mpz_t prod, mpz_t mult, signed long long int val);
void MY_mpz_sub_ui(mpz_t prod, mpz_t mult, unsigned long long int val);
#else /* #if USEDBITS > 31 */
typedef signed long int IntegerMapType;
#define MY_mpz_set_si mpz_set_si
#define MY_mpz_mul_si mpz_mul_si
#define MY_mpz_sub_ui mpz_sub_ui
#endif /* #if USEDBITS > 31 #else */

#define DOUBLE_to_VORONOIINT(y) ((IntegerMapType)(((*((long long *)&y)) & 0xFFFFFFFFFFFFFllu) >> (52 - USEDBITS)))

/*    Prerequisites for this function:
 *    sizeof(double)==sizeof(unsigned long long)
 *    doubles must be stored according to IEEE 754
 */
static inline IntegerMapType double_to_voronoiint(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu) >> (52 - USEDBITS);
}

static inline double mask_voronoi_int(double x)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d   = x;
  u.ull = u.ull & (~((1llu << (52 - USEDBITS)) - 1));
  return u.d;
}

#ifndef TWODIMS

#define EDGE_0 1  /* points 0-1 */
#define EDGE_1 2  /* points 0-2 */
#define EDGE_2 4  /* points 0-3 */
#define EDGE_3 8  /* points 1-2 */
#define EDGE_4 16 /* points 1-3 */
#define EDGE_5 32 /* points 2-3 */
#define EDGE_ALL 63

#else /* #ifndef TWODIMS */

#define EDGE_0 1 /* points 1-2 */
#define EDGE_1 2 /* points 0-2 */
#define EDGE_2 4 /* points 0-1 */
#define EDGE_ALL 7

#endif /* #ifndef TWODIMS #else */

#define HSML_INCREASE_FACTOR 1.3

#ifdef TWODIMS /* will only be compiled in 2D case */
#define DIMS 2
#else /* #ifdef TWODIMS */
#define DIMS 3
#endif /*#ifdef TWODIMS #else */

typedef struct
{
  double x, y, z;  // The 3-space position of the point
  MyIDType ID;
  int task;   // The MPI task owning this cell
  int index;  // The hydro quantity index of the cell
  int originalindex, timebin;
  unsigned int image_flags;

#ifndef OPTIMIZE_MEMORY_USAGE
  double xx, yy, zz;
  IntegerMapType ix, iy, iz;
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */

#ifdef DOUBLE_STENCIL
  MyFloat Hsml;
  int first_connection;
  int last_connection;
  char flag_primary_triangle;
#endif /* #ifdef DOUBLE_STENCIL */
} point;

typedef struct tetra_data
{
  int p[DIMS + 1];           /*!< oriented tetrahedron points */
  int t[DIMS + 1];           /*!< adjacent tetrahedrons, always opposite to corresponding point */
  unsigned char s[DIMS + 1]; /*!< gives the index of the point in the adjacent tetrahedron that
                                lies opposite to the common face */

  /* Note: if t[0] == -1, the tetrahedron has been deleted */
} tetra;

typedef struct tetra_center_data
{
#ifndef OPTIMIZE_MEMORY_USAGE
  double cx, cy, cz; /*!< describes circumcircle center */
#else                /* #ifndef OPTIMIZE_MEMORY_USAGE */
  MyFloat cx, cy, cz;
#endif               /*#ifndef OPTIMIZE_MEMORY_USAGE */
} tetra_center;

typedef struct tri_data
{
  double p[DIMS + 1][DIMS];
  int owner;
} triangle;

extern unsigned char *Edge_visited;

extern struct list_export_data
{
  unsigned int image_bits;
  int origin, index;
  int nextexport;
} * ListExports;

extern int Ninlist, MaxNinlist;

extern struct area_list_data
{
  int task, index;
  double darea;
} * AreaList;

extern int Narea, MaxNarea;

extern int NumGasInMesh;
extern int *List_InMesh;

extern struct list_P_data
{
  int firstexport, currentexport;

} * List_P;

typedef struct connection_data
{
  int task;
  int index;
  int image_flags;
  int next;

  int dp_index; /*!< this seems to be needed always the way voronoi_makeimage is implemented at the moment */
  int vf_index; /*!< index to the corresponding face */
#if defined(TETRA_INDEX_IN_FACE)
  int dt_index;
#endif /* #if defined(TETRA_INDEX_IN_FACE)*/
  MyIDType ID;
} connection;

/*! This structure contains the points where a line segment intersects
 *  the tetrahedron faces and the internal voronoi faces. Is returned
 *  by calc_voronoi_intersections().
 */
typedef struct intersection_list_data
{
  double s;       /*!< the distance from the entry point (fraction of whole segment) */
  point p;        /*!< the intersection point */
  int indA, indB; /*!< the indices of the tetra points (0-4) defining the face */
} intersection_list;

extern int CountInSphereTests, CountInSphereTestsExact;
extern int CountConvexEdgeTest, CountConvexEdgeTestExact;
extern int CountFlips, Count_1_to_3_Flips2d, Count_2_to_4_Flips2d;
extern int Count_1_to_4_Flips, Count_2_to_3_Flips, Count_3_to_2_Flips, Count_4_to_4_Flips;
extern int Count_EdgeSplits, Count_FaceSplits;
extern int Count_InTetra, Count_InTetraExact;
extern int Largest_N_DP_Buffer;

extern int Ninlist, MaxNinlist;

typedef struct individual_alloc_data
{
  double AllocFacNdp;
  double AllocFacNdt;
  double AllocFacNvf;
  double AllocFacNinlist;
  double AllocFacN_DP_Buffer;
  double AllocFacNflux;
  double AllocFacNradinflux;
  double AllocFacNvc;
} mesh_alloc_facs;

typedef struct tessellation_data
{
  int Ndp;    /*!< number of delaunay points */
  int MaxNdp; /*!< maximum number of delaunay points */
  point *DP;  /*!< delaunay points */

  int Ndt;
  int MaxNdt;        /*!< number of delaunary tetrahedra */
  tetra *DT;         /*!< Delaunay tetrahedra */
  tetra_center *DTC; /*!< circumcenters of delaunay tetrahedra */
  char *DTF;

  int Nvf;    /*!< number of Voronoi faces */
  int MaxNvf; /*!< maximum number of Voronoi faces */
  face *VF;   /*!< Voronoi faces */

  mesh_alloc_facs Indi;
} tessellation;

extern tessellation Mesh, DeRefMesh;

extern int DPinfinity;

extern int Nvc;    /* number of connections */
extern int MaxNvc; /* maximum number of connections */
extern int Largest_Nvc;
extern connection *DC; /* Connections */
extern int FirstUnusedConnection;

extern double CentralOffsetX, CentralOffsetY, CentralOffsetZ, ConversionFac;

int derefine_add_point_and_split_tri(int q, triangle *trilist, int n, int max_n, double vol);
void derefine_refine_process_edge(tessellation *T, double *vol, int tt, int nr);
void derefine_refine_compute_volumes(double *vol);
int derefine_refine_get_triangles(tessellation *T, int tt, int nr, point *dtip, triangle *trilist, int ntri, int max_n_tri);
void create_mesh(void);
void mesh_setup_exchange(void);
void free_mesh(void);
void free_mesh_structures_not_needed_for_derefinement_refinement(void);
void free_all_remaining_mesh_structures(void);
void apply_area_list(void);
int area_list_data_compare(const void *a, const void *b);
void write_voronoi_mesh(tessellation *T, char *fname, int writeTask, int lastTask);
void initialize_and_create_first_tetra(tessellation *T);
void compute_voronoi_faces_and_volumes(void);
void get_line_segments(int sphp_index, int dp_index, double *segments, unsigned int *nof_elements, unsigned int max_elements);
double cross_section_plane_cell(int sphp_index, int dp_index, double *center, double *n);
void intersections_plane_cell(int sphp_index, int dp_index, double *center, double *n, double *polygon, unsigned int *nof_elements);
void intersection_plane_grid(double *center, double *n, const char *filename);
void process_edge_faces_and_volumes(tessellation *T, int tt, int nr);
int insert_point(tessellation *T, int pp, int ttstart);
void make_an_edge_split(tessellation *T, int tt0, int edge_nr, int count, int pp, int *ttlist);
void make_a_face_split(tessellation *T, int tt0, int face_nr, int pp, int tt1, int tt2, int qq1, int qq2);
double calculate_tetra_volume(point *p0, point *p1, point *p2, point *p3);
void make_a_4_to_4_flip(tessellation *T, int tt, int tip_index, int edge_nr);
double get_tri_volume(int i, triangle *trilist);
void make_a_1_to_4_flip(tessellation *T, int pp, int tt0, int tt1, int tt2, int tt3);
void make_a_3_to_2_flip(tessellation *T, int tt0, int tt1, int tt2, int tip, int edge, int bottom);
void make_a_2_to_3_flip(tessellation *T, int tt0, int tip, int tt1, int bottom, int qq, int tt2);
int get_tetra(tessellation *T, point *p, int *moves, int ttstart, int *flag, int *edgeface_nr);
int InTetra(tessellation *T, int tt, point *pp, int *edgeface_nr, int *nexttetra);
double InSphere(point *p0, point *p1, point *p2, point *p3, point *p);
void update_circumcircle(tessellation *T, int tt);
int test_tetra_orientation(point *p0, point *p1, point *p2, point *p3);
int voronoi_ghost_search_alternative(tessellation *T);
void compute_circumcircles(tessellation *T);
int compute_max_delaunay_radius(void);
void check_for_min_distance(tessellation *T);
void check_links(tessellation *T);
void check_orientations(tessellation *T);
void check_tetras(tessellation *T, int npoints);
int voronoi_get_local_particles(void);
int convex_edge_test(tessellation *T, int tt, int tip, int *edgenr);
void calculate_gradients(void);
void limit_gradient(double *d, double phi, double min_phi, double max_phi, MySingle *dphi);
void exchange_primitive_variables(void);
void exchange_primitive_variables_and_gradients(void);
int compare_primexch(const void *a, const void *b);

/* 2D voronoi routines */
void check_edge_and_flip_if_needed(tessellation *T, int ip, int it);
int get_triangle(tessellation *T, int pp, int *moves, int *degenerate_flag, int ttstart);
double InCircle(point *p0, point *p1, point *p2, point *p);
void make_a_1_to_3_flip(tessellation *T, int pp, int tt0, int tt1, int tt2);
double test_triangle_orientation(tessellation *T, int pp0, int pp1, int pp2);
void make_a_2_to_4_flip(tessellation *T, int pp, int tt0, int tt1, int tt2, int tt3, int i0, int j0);
void dump_points(tessellation *T);
void set_integers_for_pointer(point *p);

#if !defined(ONEDIMS)
#ifndef OPTIMIZE_MEMORY_USAGE
static inline void set_integers_for_point(tessellation *T, int pp)
{
  point *p = &T->DP[pp];
  set_integers_for_pointer(p);
}
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
static inline void get_integers_for_point(point *p, IntegerMapType ixyz[], double xyz[])
{
  xyz[0] = (p->x - CentralOffsetX) * ConversionFac + 1.0;
  xyz[1] = (p->y - CentralOffsetY) * ConversionFac + 1.0;
  xyz[2] = (p->z - CentralOffsetZ) * ConversionFac + 1.0;

  ixyz[0] = double_to_voronoiint(xyz[0]);
  ixyz[1] = double_to_voronoiint(xyz[1]);
  ixyz[2] = double_to_voronoiint(xyz[2]);

  xyz[0] = mask_voronoi_int(xyz[0]);
  xyz[1] = mask_voronoi_int(xyz[1]);
  xyz[2] = mask_voronoi_int(xyz[2]);
}
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

#else  /* #if !defined(ONEDIMS) */
void set_integers_for_point(tessellation *T, int pp);
#endif /* #if !defined(ONEDIMS) #else */

/* quick function to compare a point to the infinity point */
static inline int isInfinity(point *p) { return p->x == MAX_DOUBLE_NUMBER; }

int solve_linear_equations(double *m, double *res);
void check_triangles(tessellation *T, int npoints);
int InCircle_Quick(tessellation *T, int pp0, int pp1, int pp2, int pp);
int InCircle_Errorbound(tessellation *T, int pp0, int pp1, int pp2, int pp);
int InCircle_Exact(tessellation *T, int pp0, int pp1, int pp2, int pp);
int Orient2d_Exact(tessellation *T, int pp0, int pp1, int pp2);
int Orient2d_Quick(tessellation *T, int pp0, int pp1, int pp2);
int FindTriangle(tessellation *T, int tt, int pp, int *degnerate_flag, int *nexttetra);
int InSphere_Exact(point *p0, point *p1, point *p2, point *p3, point *p);
int InSphere_Quick(point *p0, point *p1, point *p2, point *p3, point *p);
int InSphere_Errorbound(point *p0, point *p1, point *p2, point *p3, point *p);
int Orient3d_Quick(point *p0, point *p1, point *p2, point *p3);
int Orient3d(point *p0, point *p1, point *p2, point *p3);
int Orient3d_Exact(point *p0, point *p1, point *p2, point *p3);
int count_undecided_tetras(tessellation *T);
int ngb_treefind_ghost_search(tessellation *T, MyDouble searchcenter[3], MyDouble refpos[3], MyFloat hsml, MyFloat maxdist, int target,
                              int origin, int mode, int thread_id, int numnodes, int *firstnode);
int voronoi_ghost_search_evaluate(tessellation *T, int target, int mode, int q, int thread_id);
int voronoi_ghost_search(tessellation *T);
double distance_to_border(int cell);

#endif /* HAVE_H_VORONOI */
