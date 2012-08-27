#ifndef PROTOTYPES_H
#define PROTOTYPES_H 
#include "copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions from the /src directory.      */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "athena.h"
#include "defs.h"
#include "config.h"

/* Include prototypes from /src sub-directories */
#ifdef FFT_ENABLED
#include "fftsrc/prototypes.h"
#endif

#include "gravity/prototypes.h"
#include "integrators/prototypes.h"
#include "microphysics/prototypes.h"
#include "particles/prototypes.h"
#include "reconstruction/prototypes.h"
#include "rsolvers/prototypes.h"

/*----------------------------------------------------------------------------*/
/* main.c */
int athena_main(int argc, char *argv[]);

/*----------------------------------------------------------------------------*/
/* ath_array.c */
void*   calloc_1d_array(                      size_t nc, size_t size);
void**  calloc_2d_array(           size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);
void free_1d_array(void *array);
void free_2d_array(void *array);
void free_3d_array(void *array);

/*----------------------------------------------------------------------------*/
/* ath_log.c */
void ath_log_set_level(const int out, const int err);
void ath_log_open(const char *basename, const int lazy, const char *mode);
void ath_log_close(void);
FILE *athout_fp(void);
FILE *atherr_fp(void);
void ath_flush_out(void);
void ath_flush_err(void);
int ath_perr(const int level, const char *fmt, ...);
int ath_pout(const int level, const char *fmt, ...);

/*----------------------------------------------------------------------------*/
/* ath_files.c */
char *ath_fname(const char *path, const char *basename,
                const char *levstr, const char *domstr,
                const int dlen, const int idump, 
                const char *id, const char *ext);

/*----------------------------------------------------------------------------*/
/* ath_signal.c */
void ath_sig_init(void);
int  ath_sig_act(int *piquit);

/*----------------------------------------------------------------------------*/
/* baton.c */
void baton_start(const int Nb, const int tag);
void baton_stop(const int Nb, const int tag);

/*----------------------------------------------------------------------------*/
/* bvals_mhd.c  */
void bvals_mhd_init(MeshS *pM);
void bvals_mhd_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc);
void bvals_mhd(DomainS *pDomain);

/*----------------------------------------------------------------------------*/
/* bvals_shear.c  */
#ifdef SHEARING_BOX
void ShearingSheet_ix1(DomainS *pD);
void ShearingSheet_ox1(DomainS *pD);
void RemapEy_ix1(DomainS *pD, Real ***emfy, Real **remapEyiib);
void RemapEy_ox1(DomainS *pD, Real ***emfy, Real **remapEyoib);
void bvals_shear_init(MeshS *pM);
void bvals_shear_destruct(void);
#ifdef FARGO
void Fargo(DomainS *pD);
#endif
#endif /* SHEARING_BOX */

#if defined (FARGO) && defined (CYLINDRICAL)
void bvals_shear_init(MeshS *pM);
void bvals_shear_destruct(void);
void Fargo(DomainS *pD);
#endif 

/*----------------------------------------------------------------------------*/
/* cc_pos.c */
void cc_pos(const GridS *pG, const int i, const int j,const int k,
            Real *px1, Real *px2, Real *px3);
void fc_pos(const GridS *pG, const int i, const int j,const int k,
            Real *px1, Real *px2, Real *px3);
#ifdef CYLINDRICAL
Real x1vc(const GridS *pG, const int i);
#endif
#ifdef PARTICLES
int celli(const GridS *pGrid, const Real x, const Real dx1_1, int *i, Real *a);
Real x1cc(const GridS *pGrid, const int i);
int cellj(const GridS *pGrid, const Real y, const Real dx2_1, int *j, Real *b);
Real x2cc(const GridS *pGrid, const int j);
int cellk(const GridS *pGrid, const Real z, const Real dx3_1, int *k, Real *c);
Real x3cc(const GridS *pGrid, const int k);
#endif

/*----------------------------------------------------------------------------*/
/* convert_var.c */
PrimS Cons_to_Prim(const ConsS *pU);
ConsS Prim_to_Cons(const PrimS *pW);
Prim1DS Cons1D_to_Prim1D(const Cons1DS *pU, const Real *pBx);
Cons1DS Prim1D_to_Cons1D(const Prim1DS *pW, const Real *pBx);
#ifndef SPECIAL_RELATIVITY
Real cfast(const Cons1DS *U, const Real *Bx);
#endif
#ifdef SPECIAL_RELATIVITY
PrimS check_Prim(const ConsS *pU);
#ifdef MHD
PrimS fix_vsq (const ConsS *pU);
PrimS entropy_fix (const ConsS *pU, const Real *ent);
Prim1DS check_Prim1D (const Cons1DS *pU, const Real *pBx);
#endif /* MHD */
#endif /* SPECIAL_RELATIVITY */


/*----------------------------------------------------------------------------*/
/* init_grid.c */
void init_grid(MeshS *pM);

/*----------------------------------------------------------------------------*/
/* init_mesh.c */
void init_mesh(MeshS *pM);
void get_myGridIndex(DomainS *pD, const int my_id, int *pi, int *pj, int *pk);

/*----------------------------------------------------------------------------*/
/* new_dt.c */
void new_dt(MeshS *pM);

/*----------------------------------------------------------------------------*/
/* output.c - and related files */
void init_output(MeshS *pM);
void data_output(MeshS *pM, const int flag);
void add_rst_out(OutputS *new_out);
void data_output_destruct(void);
void dump_history_enroll(const ConsFun_t pfun, const char *label);
Real ***OutData3(GridS *pGrid, OutputS *pOut, int *Nx1, int *Nx2, int *Nx3);
Real  **OutData2(GridS *pGrid, OutputS *pOut, int *Nx1, int *Nx2);
Real   *OutData1(GridS *pGrid, OutputS *pOut, int *Nx1);

void output_pdf  (MeshS *pM, OutputS *pOut);
void output_pgm  (MeshS *pM, OutputS *pOut);
void output_ppm  (MeshS *pM, OutputS *pOut);
void output_vtk  (MeshS *pM, OutputS *pOut);
void output_tab  (MeshS *pM, OutputS *pOut);

void dump_binary  (MeshS *pM, OutputS *pOut);
void dump_history (MeshS *pM, OutputS *pOut);
void dump_tab_cons(MeshS *pM, OutputS *pOut);
void dump_tab_prim(MeshS *pM, OutputS *pOut);
void dump_vtk     (MeshS *pM, OutputS *pOut);

/*----------------------------------------------------------------------------*/
/* par.c */
void   par_open(char *filename);
void   par_cmdline(int argc, char *argv[]);
int    par_exist(char *block, char *name);

char  *par_gets(char *block, char *name);
int    par_geti(char *block, char *name);
double par_getd(char *block, char *name);

char  *par_gets_def(char *block, char *name, char   *def);
int    par_geti_def(char *block, char *name, int    def);
double par_getd_def(char *block, char *name, double def);

void   par_sets(char *block, char *name, char *sval, char *comment);
void   par_seti(char *block, char *name, char *fmt, int ival, char *comment);
void   par_setd(char *block, char *name, char *fmt, double dval, char *comment);

void   par_dump(int mode, FILE *fp);
void   par_close(void);

#ifdef MPI_PARALLEL
void par_dist_mpi(const int mytid, MPI_Comm comm);
#endif

/*----------------------------------------------------------------------------*/
/* prob/PROBLEM.c ; linked to problem.c */
void problem(DomainS *pD);
void Userwork_in_loop(MeshS *pM);
void Userwork_after_loop(MeshS *pM);
void problem_read_restart(MeshS *pM, FILE *fp);
void problem_write_restart(MeshS *pM, FILE *fp);
ConsFun_t get_usr_expr(const char *expr);
VOutFun_t get_usr_out_fun(const char *name);
#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name);
void gasvshift(const Real x1, const Real x2, const Real x3, Real *u1, Real *u2, Real *u3);
void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3, const Real v1, const Real v2, const Real v3);
#endif

/*----------------------------------------------------------------------------*/
/* restart.c  */
void dump_restart(MeshS *pM, OutputS *pout);
void restart_grids(char *res_file, MeshS *pM);

/*----------------------------------------------------------------------------*/
/* show_config.c */
void show_config(void);
void show_config_par(void);

/*----------------------------------------------------------------------------*/
/* smr.c */
void RestrictCorrect(MeshS *pM);
void Prolongate(MeshS *pM);
void SMR_init(MeshS *pM);

/*----------------------------------------------------------------------------*/
/* utils.c */
char *ath_strdup(const char *in);
int ath_gcd(int a, int b);
int ath_big_endian(void);
void ath_bswap(void *vdat, int sizeof_len, int cnt);
void ath_error(char *fmt, ...);
void minmax1(Real   *data, int nx1,                   Real *dmin, Real *dmax);
void minmax2(Real  **data, int nx2, int nx1,          Real *dmin, Real *dmax);
void minmax3(Real ***data, int nx3, int nx2, int nx1, Real *dmin, Real *dmax);
void do_nothing_bc(GridS *pG);
Real compute_div_b(GridS *pG);
int sign_change(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *a, Real *b);
int bisection(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *root);
Real trapzd(Real (*func)(Real), const Real a, const Real b, const int n, const Real s);
Real qsimp(Real (*func)(Real), const Real a, const Real b);
Real avg1d(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real avg2d(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real avg3d(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real avgXZ(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real vecpot2b1i(Real (*A2)(Real,Real,Real), Real (*A3)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k);
Real vecpot2b2i(Real (*A1)(Real,Real,Real), Real (*A3)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k);
Real vecpot2b3i(Real (*A1)(Real,Real,Real), Real (*A2)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k);
#ifdef PARTICLES
void InverseMatrix(Real **a, int n, Real **b);
void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c);
#endif
#endif /* PROTOTYPES_H */
