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
 * \file        src/main/proto.h
 * \date        05/2018
 * \brief       Function declarations.
 * \details     No particular order.
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 29.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef PROTO_H
#define PROTO_H

#include "../gravity/forcetree.h"
#include "../main/allvars.h"
#include "../utils/timer.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef IMPOSE_PINNING
#ifndef __USE_GNU
#define __USE_GNU
#endif /* #ifndef __USE_GNU */
#include <sched.h>
#endif /* #ifdef IMPOSE_PINNING */

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif /* #ifdef HAVE_HDF5 */

#if defined(COOLING)
#include "../cooling/cooling_proto.h"
#endif /* #if defined(COOLING) */

void sfr_init();
void sfr_create_star_particles(void);
void ngb_finish_rangebounds_update(int nchanged, int *nodelist);
void ngb_update_rangebounds(int i, int *nchanged, int *nodelist);
int ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport,
                          int *nsend_local);
int ngb_treebuild(int npart);
void ngb_treeallocate(void);
void ngb_treefree(void);
int ngb_treefind_export_node_threads(int no, int target, int thread_id, int image_flag);
int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int mode, int thread_id, int numnodes,
                                  int *firstnode);

void drift_node(struct NgbNODE *current, integertime time1);
void drift_all_particles(void);
double get_desired_softening_from_mass(double mass);
void log_restart_debug(void);
int get_thread_num(void);
void report_pinning(void);
void detect_topology(void);
void pin_to_core_set(void);
void get_core_set(void);
int derefine_should_this_cell_be_merged(int i, int flag);

void gravity_external(void);
void gravity(int timebin, int fullflag);
int my_ffsll(peanokey i);
void set_cosmo_factors_for_current_time(void);
void calc_exact_gravity_for_particle_type(void);
void calculate_non_standard_physics_with_valid_gravity_tree(void);
void calculate_non_standard_physics_with_valid_gravity_tree_always(void);
int get_softeningtype_for_hydro_cell(int i);
void gravity_forcetest_testforcelaw(void);
void *myfree_query_last_block(void);

void subdivide_evenly(int N, int pieces, int index, int *first, int *count);
void force_evaluate_direct(int target, int result_idx, int nimport);
void gravity_direct(int timebin);
double dabs(double a);
double dmax(double a, double b);
double dmin(double a, double b);
double max_array(double *a, int num_elements);
int imax(int a, int b);
int imin(int a, int b);
double mysort(void *base, size_t nel, size_t width, int (*compar)(const void *, const void *));

int myflush(FILE *fstream);
int flush_everything(void);
void gravity_force_finalize(int timebin);
void permutate_chunks_in_list(int ncount, int *list);
double get_default_softening_of_particletype(int type);
double get_random_number_aux(void);
void sumup_large_ints_comm(int n, int *src, long long *res, MPI_Comm comm);
void ngb_update_velocities(void);
void hello(void);
void find_long_range_step_constraint(void);

void ngb_treemodifylength(int delta_NgbMaxPart);
void domain_resize_storage(int count_get, int count_get_sph, int option_flag);
void init_individual_softenings(void);
void do_derefinements_and_refinements();
void mark_active_timebins(void);
void voronoi_test(void);
void execute_resubmit_command(void);
void output_compile_time_options(void);
void init_io_fields();
void produce_dump(void);

void create_snapshot_if_desired(void);
void output_log_messages(void);
void mpi_report_committable_memory(void);
long long report_comittable_memory(long long *MemTotal, long long *Committed_AS, long long *SwapTotal, long long *SwapFree);
int check_for_interruption_of_run(void);
void set_non_standard_physics_for_current_time(void);
void calculate_non_standard_physics_prior_mesh_construction(void);
void calculate_non_standard_physics_end_of_step(void);
void compute_statistics(void);
void face_limit_fluxes(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R,
                       struct fluxes *flux, double dt, double *count, double *count_reduced);

double get_sound_speed(int p);
void set_pressure_of_cell(int i);
void gradient_init(MyFloat *addr, MyFloat *addr_exch, MySingle *addr_grad, int type);
void limit_vel_gradient(double *d, MySingle *grad_vx, MySingle *grad_vy, MySingle *grad_vz, double csnd);
void subfind_density_hsml_guess(void);
void peano_hilbert_key_inverse(peanokey key, int bits, peano1D *x, peano1D *y, peano1D *z);
void find_nearest_meshpoint_global(mesh_search_data *searchdata, int n, int hsmlguess, int verbose);
void reorder_DP(void);
void peano_hilbert_order_DP(void);
void validate_vertex_velocities(void);

double get_cell_radius(int i);
double nearest_x(double d);
double nearest_y(double d);
double nearest_z(double d);
int voronoi_get_connected_particles(tessellation *T);
void voronoi_init_connectivity(tessellation *T);
void voronoi_update_connectivity(tessellation *T);
int compare_foreign_connection(const void *a, const void *b);
void voronoi_remove_connection(int i);
int pmforce_is_particle_high_res(int type, MyDouble *pos);

void cooling_only(void);
void report_VmRSS(void);
void tree_based_timesteps_setsoundspeeds(void);
void voronoi_update_ghost_velvertex(void);
int should_this_cell_be_split(int i);
int do_refinements(void);
int should_this_cell_be_merged(int i, int flag);
int do_derefinements(void);
void move_collisionless_particle(int new_i, int old_i);
void dump_memory_table(void);

void report_detailed_memory_usage_of_largest_task(void);
void calculate_vertex_velocity_divergence(void);
void make_list_of_active_particles(void);
void find_gravity_timesteps_and_do_gravity_step_first_half(void);
void do_gravity_step_second_half(void);
void voronoi_1D_reorder_gas(void);
int voronoi_1D_compare_key(const void *a, const void *b);
void voronoi_1D_order(void);
void pm2d_init_periodic(void);
void pm2d_init_periodic_allocate(void);

void pm2d_init_periodic_free(void);
void pm2d_force_periodic(int mode);
int pm2d_periodic_compare_sortindex(const void *a, const void *b);
void pm2d_mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
int timestep_evaluate(int target, int mode, int threadid);
void tree_based_timesteps(void);
int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
int MPI_hypercube_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcount, int *displs,
                             MPI_Datatype recvtype, MPI_Comm comm);
double parallel_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));

double parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *), MPI_Comm comm);
int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);
void drift_particle(int i, integertime time1);
void put_symbol(char *string, double t0, double t1, char c);
void write_cpu_log(void);
void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr, int clear_flag, char *origin);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line, char *origin);
void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);
void mymalloc_init(void);
void calculate_maxid(void);
void determine_compute_nodes(void);
double INLINE_FUNC hubble_function(double a);
void fof_fof(int num);
double fof_find_groups(MyIDType *vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask);
void fof_compile_catalogue(void);
void fof_save_groups(int num);

double fof_periodic(double x);
double fof_periodic_wrap(double x);
double fof_find_nearest_dmparticle(MyIDType *vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask);
void fof_compute_group_properties(int gr, int start, int len);
int fof_compare_FOF_PList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b);
int fof_compare_FOF_GList_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b);
int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b);

int fof_compare_Group_GrNr(const void *a, const void *b);
int fof_compare_Group_MinIDTask(const void *a, const void *b);
int fof_compare_Group_MinID(const void *a, const void *b);
int fof_compare_ID_list_GrNrID(const void *a, const void *b);
int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b);
int fof_compare_Group_Len(const void *a, const void *b);
int fof_compare_aux_sort_Type(const void *a, const void *b);
int fof_compare_aux_sort_GrNr(const void *a, const void *b);
int fof_compare_aux_sort_OriginTask_OriginIndex(const void *a, const void *b);
int fof_compare_aux_sort_FileOrder(const void *a, const void *b);

int fof_compare_local_sort_data_targetindex(const void *a, const void *b);
void fof_subfind_exchange(MPI_Comm Communicator);
void fof_prepare_output_order(void);
void fof_compute_group_properties(int gr, int start, int len);
void fof_exchange_group_data(void);
void fof_finish_group_properties(void);
double fof_get_comoving_linking_length(void);
void fof_assign_group_numbers(void);
void fof_reorder_PS(int *Id, int Nstart, int N);
void fof_subfind_write_file(char *fname, int writeTask, int lastTask);

void fof_subfind_prepare_ID_list(void);
int subfind_compare_procassign_GrNr(const void *a, const void *b);
double subfind_so_potegy(double *egypot);
void subfind_distlinklist_get_two_heads(long long ngb_index1, long long ngb_index2, long long *head, long long *head_attach);
void fof_check_for_full_nodes_recursive(int no);
int fof_return_a_particle_in_cell_recursive(int no);
void subfind_loctree_copyExtent(void);
int subfind_distlinklist_get_tail_set_tail_increaselen(long long index, long long *tail, long long newtail);
void subfind_reorder_according_to_submp(void);
int subfind_compare_submp_OldIndex(const void *a, const void *b);

int subfind_compare_submp_GrNr_DM_Density(const void *a, const void *b);
double subfind_exchange(void);
void subfind_coll_domain_decomposition(void);
void subfind_coll_domain_combine_topleaves_to_domains(int ncpu, int ndomain);
void subfind_coll_domain_free(void);
void subfind_coll_domain_allocate(void);
int subfind_coll_domain_determineTopTree(void);
void subfind(int num);
double subfind_density(int mode);
double subfind_overdensity(void);

void subfind_save_final(int num);
void subfind_process_group_collectively(int nsubgroups_cat);
void subfind_coll_findExtent(void);
void subfind_reorder_PS(int *Id, int Nstart, int N);
void subfind_reorder_P(int *Id, int Nstart, int N);
void subfind_distribute_particles(MPI_Comm Communicator);
void subfind_coll_domain_walktoptree(int no);
int subfind_compare_densities(const void *a, const void *b);
int subfind_compare_binding_energy(const void *a, const void *b);
int subfind_compare_dist_rotcurve(const void *a, const void *b);

int subfind_compare_coll_candidates_rank(const void *a, const void *b);
int subfind_compare_coll_candidates_boundlength(const void *a, const void *b);
int subfind_compare_coll_candidates_nsubs(const void *a, const void *b);
int subfind_compare_coll_candidates_subnr(const void *a, const void *b);
void subfind_col_find_coll_candidates(int totgrouplen);
void subfind_unbind_independent_ones(int count);
void subfind_distribute_groups(void);
void subfind_potential_compute(int num, struct unbind_data *d, int phase, double weakly_bound_limit);
int subfind_col_unbind(struct unbind_data *d, int num, int *num_non_gas);
void subfind_find_linkngb(void);

int subfind_loctree_treebuild(int npart, struct unbind_data **mp);
void subfind_loctree_update_node_recursive(int no, int sib, int father);
double subfind_loctree_treeevaluate_potential(int target);
void subfind_loctree_copyExtent(void);
double subfind_locngb_treefind(MyDouble xyz[3], int desngb, double hguess);
void subfind_loctree_findExtent(int npart, struct unbind_data *mp);
int subfind_locngb_treefind_variable(MyDouble searchcenter[3], double hguess);
size_t subfind_loctree_treeallocate(int maxnodes, int maxpart);
void subfind_loctree_treefree(void);
void subfind_find_nearesttwo(void);

int subfind_process_group_serial(int gr, int offset, int nsubgroups_cat);
int subfind_unbind(struct unbind_data *ud, int len, int *len_non_gas);
int subfind_locngb_compare_key(const void *a, const void *b);
int subfind_compare_serial_candidates_subnr(const void *a, const void *b);
int subfind_compare_serial_candidates_rank(const void *a, const void *b);
int subfind_compare_dens(const void *a, const void *b);
int subfind_compare_serial_candidates_boundlength(const void *a, const void *b);
int subfind_compare_dist_rotcurve(const void *a, const void *b);
int subfind_compare_binding_energy(const void *a, const void *b);
int subfind_compare_densities(const void *a, const void *b);

int subfind_compare_ID_list(const void *a, const void *b);
int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b);
void subfind_poll_for_requests(void);
long long subfind_distlinklist_setrank_and_get_next(long long index, long long *rank);
long long subfind_distlinklist_get_rank(long long index);
void subfind_distlinklist_set_next(long long index, long long next);
void subfind_distlinklist_add_particle(long long index);
void subfind_distlinklist_add_bound_particles(long long index, int nsub);
void subfind_distlinklist_mark_particle(long long index, int target, int submark);
long long subfind_distlinklist_get_next(long long index);

long long subfind_distlinklist_get_head(long long index);
void subfind_distlinklist_set_headandnext(long long index, long long head, long long next);
void subfind_distlinklist_set_tailandlen(long long index, long long tail, int len);
void subfind_distlinklist_get_tailandlen(long long index, long long *tail, int *len);
void subfind_distlinklist_set_all(long long index, long long head, long long tail, int len, long long next);
long long subfind_distlinklist_set_head_get_next(long long index, long long head);
int subfind_compare_dist_rotcurve(const void *a, const void *b);
void subfind_coll_treeallocate(int maxpart, int maxindex);
void subfind_coll_treefree(void);
void subfind_coll_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z);

void subfind_coll_exchange_topleafdata(void);
void subfind_coll_update_node_recursive(int no, int sib, int father, int *last);
void subfind_coll_insert_pseudo_particles(void);
int subfind_coll_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, unsigned long long xc, unsigned long long yc,
                                    unsigned long long zc, unsigned long long ilen);
int subfind_coll_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char levels);
int subfind_coll_treebuild_construct(int npart, struct unbind_data *mp);
int subfind_coll_treebuild(int npart, struct unbind_data *mp);
double subfind_get_particle_balance(void);
int subfind_fof_compare_ID(const void *a, const void *b);
void write_file(char *fname, int readTask, int lastTask, int subbox_flag);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);
int get_values_per_blockelement(enum iofields blocknr);
int get_datatype_in_block(enum iofields blocknr, int mode);
void get_dataset_name(enum iofields blocknr, char *buf);
int blockpresent(enum iofields blocknr, int write);
void fill_write_buffer(void *buffer, enum iofields blocknr, int *pindex, int pc, int type, int subbox_flag);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);
int get_particles_in_block(enum iofields blocknr, int *typelist);
int get_bytes_per_blockelement(enum iofields blocknr, int mode);
void read_file(const char *fname, int filenr, int readTask, int lastTask, int);

void get_Tab_IO_Label(enum iofields blocknr, char *label);
void long_range_init_regionsize(void);
int find_files(const char *fname);
double get_random_number(void);
int peano_compare_key(const void *a, const void *b);
void mysort_domain(void *b, size_t n, size_t s);
void mysort_peano(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
int density_isactive(int n);
size_t sizemax(size_t a, size_t b);
void my_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno);

void reconstruct_timebins(void);
peanokey peano_hilbert_key(peano1D x, peano1D y, peano1D z, int bits);
void enable_core_dumps_and_fpu_exceptions(void);
void find_next_sync_point(void);
void set_units_sfr(void);
void gravity_forcetest(void);
void allocate_memory(void);
void begrun0(void);
void begrun1(void);
void begrun2(void);

int init(void);
void loadrestart(void);
void reread_params_after_loading_restart(void);
void check_omega(void);
void close_logfiles(void);
void compute_grav_accelerations(int timebin, int fullflag);
void compute_global_quantities_of_system(void);
void cooling_and_starformation(void);
void density(void);
void do_box_wrapping(void);

void domain_Decomposition(void);
double enclosed_mass(double R);
void endrun(void);
void energy_statistics(void);
void ewald_corr(double dx, double dy, double dz, double *fper);
void ewald_force(double x, double y, double z, double force[3]);
int my_fls(int x);
void ewald_init(void);
double ewald_psi(double x, double y, double z);
double ewald_pot_corr(double dx, double dy, double dz);

integertime find_next_outputtime(integertime time);
void minimum_large_ints(int n, long long *src, long long *res);
double get_starformation_rate(int i);
double calc_egyeff(int i, double gasdens, double *ne, double *x, double *tsfr, double *factorEVP);
void gravity_tree(int timebin);
void init_clouds(void);
void integrate_sfr(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void open_logfiles(void);

void peano_hilbert_order(void);
void predict(double time);
void read_ic(const char *fname, int);
void read_header_attributes(FILE *fd);
MyIDType determine_ids_offset(void);
int read_outputlist(char *fname);
void read_parameter_file(char *fname);
void check_parameters();
void reorder_gas(int *Id);
void reorder_particles(int *Id);

void restart(int mod);
void run(void);
void savepositions(int num, int subbox_flag);
void mpi_printf(const char *fmt, ...);
void mpi_fprintf(FILE *stream, const char *fmt, ...);
void mpi_printf_each(const char *fmt, ...);
FILE *open_file(char *);
double second(void);
void set_softenings(void);
void set_units(void);

void setup_smoothinglengths(void);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);
void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
double get_hydrokick_factor(integertime time0, integertime time1);
double get_gravkick_factor(integertime time0, integertime time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);

double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1);
double measure_time(void);
void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(int mode, int *typelist);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);

int pmforce_nonperiodic(int grnr);
void readjust_timebase(double TimeMax_old, double TimeMax_new);
void pm_setup_nonperiodic_kernel(void);
void init_gradients();
void init_scalars();
void print_particle_info(int i);
void print_state_info(struct state *st);
void print_state_face_info(struct state_face *st);
void face_set_scalar_states_and_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
void face_turn_momentum_flux(struct fluxes *flux, struct geometry *geom);

void face_clear_fluxes(struct fluxes *flux);
int face_check_responsibility_of_this_task(tessellation *T, int p1, int p2, struct state *st_L, struct state *st_R);
int face_get_normals(tessellation *T, int i, struct geometry *geom);
int face_get_state(tessellation *T, int p, int i, struct state *st);
void face_boundary_check(point *p, double *velx, double *vely, double *velz);
void face_boundary_check_vertex(tessellation *T, int p, MyFloat *velx, MyFloat *vely, MyFloat *velz);
double face_timestep(struct state *state_L, struct state *state_R, double *hubble_a, double *atime);
void state_convert_to_local_frame(struct state *st, double *vel_face, double hubble_a, double atime);
void face_do_time_extrapolation(struct state *delta, struct state *st, double atime);
void face_do_spatial_extrapolation(struct state *delta, struct state *st, struct state *st_other);

void face_do_spatial_extrapolation_single_quantity(double *delta, double st, double st_other, MySingle *grad, double *dx, double *r);
void face_add_extrapolations(struct state *st_face, struct state *delta_time, struct state *delta_space, struct fvs_stat *stat);
void face_add_extrapolation(struct state *st_face, struct state *delta, struct fvs_stat *stat);
void face_turn_velocities(struct state *st, struct geometry *geom);
void solve_advection(struct state *st_L, struct state *st_R, struct state_face *st_face, struct geometry *geom, double *vel_face);
void face_turnback_velocities(struct state_face *st_face, struct geometry *geom);
void face_get_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux, struct geometry *geom,
                     double *vel_face);
void face_add_fluxes_advection(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face);
double godunov_flux_3d(struct state *st_L, struct state *st_R, struct state_face *st_face);
void sample_solution_vacuum_left_3d(double S, struct state *st_R, struct state_face *st_face);

void sample_solution_vacuum_right_3d(double S, struct state *st_L, struct state_face *st_face);
void sample_solution_vacuum_generate_3d(double S, struct state *st_L, struct state *st_R, struct state_face *st_face);
void get_mach_numbers(struct state *st_L, struct state *st_R, double Press);
void sample_solution_3d(double S, struct state *st_L, struct state *st_R, double Press, double Vel, struct state_face *st_face);
int riemann(struct state *st_L, struct state *st_R, double *Press, double *Vel);
void pressure_function(double P, struct state *st, double *F, double *FD);
double guess_for_pressure(struct state *st_L, struct state *st_R);
void riemann_isotherm(struct state *st_L, struct state *st_R, double *Rho, double *Vel, double csnd);
void isothermal_function(double rhostar, double rho, double *F, double *FD);
void sample_solution_isothermal3d(double S, struct state *st_L, struct state *st_R, double Rho, double Vel, struct state_face *st_face,
                                  double csnd);

void apply_flux_list(void);
int flux_list_data_compare(const void *a, const void *b);
void set_vertex_velocities(void);
int scalar_init(MyFloat *addr, MyFloat *addr_mass, int type);
void compute_interface_fluxes(tessellation *T);
void update_primitive_variables(void);
void set_pressure_of_cell_internal(struct particle_data *P, struct sph_particle_data *SphP, int i);
void do_validity_checks(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);
void update_primitive_variables_single(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);

void update_internal_energy(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);
void mpi_exchange_buffers(void *send_buf, int *send_count, int *send_offset, void *recv_buf, int *recv_count, int *recv_offset,
                          int item_size, int commtag, int include_self);
int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);
void *sort_based_on_mesh_search(mesh_search_data *search, void *data, int n_items, int item_size);
void *sort_based_on_field(void *data, int field_offset, int n_items, int item_size);
void mpi_distribute_items_from_search(mesh_search_data *search, void *data, int *n_items, int *max_n, int item_size, int commtag,
                                      int task_offset, int cell_offset);
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size, int commtag);
void tile_ics(void);
void reallocate_memory_maxpart(void);
void reallocate_memory_maxpartsph(void);

void share_particle_number_in_file(const char *fname, int filenr, int readTask, int lastTask, int readTypes);
int dump_memory_table_buffer(char *p);
void calc_memory_checksum(void *base, size_t bytes);
void allreduce_sparse_double_sum(double *loc, double *glob, int N);
void allreduce_sparse_imin(int *loc, int *glob, int N);
void myMPI_Alltoallv(void *sendb, size_t *sendcounts, size_t *sdispls, void *recvb, size_t *recvcounts, size_t *rdispls, int len,
                     int big_flag, MPI_Comm comm);
int myMPI_Sendrecv(void *sendb, size_t sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvb, size_t recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
size_t roundup_to_multiple_of_cacheline_size(size_t n);
void init_cpu_log(void);

void write_error(int check, size_t nwritten, size_t nmemb);
size_t smax(size_t a, size_t b);
void init_field(enum iofields field, const char *label, const char *datasetname, enum types_in_memory type_in_memory,
                enum types_in_file type_in_file_output, enum types_in_file type_in_file_input, int values_per_block, enum arrays array,
                void *pointer_to_field, void (*io_func)(int, int, void *, int), int typelist_bitmask);
void init_units(enum iofields field, double a, double h, double L, double M, double V, double c);
void init_snapshot_type(enum iofields field, enum sn_type type);

void swap_Nbyte(char *data, int n, int m);
void swap_header(void);

#if defined(COOLING)
void cool_cell(int i);
#endif /* #if defined(COOLING) */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
void special_particle_create_list();
void special_particle_update_list();
#endif /* #ifdef  EXACT_GRAVITY_FOR_PARTICLE_TYPE */

#ifdef HAVE_HDF5

hid_t my_H5Fcreate(const char *fname, unsigned flags, hid_t fcpl_id, hid_t fapl_id);
hid_t my_H5Gcreate(hid_t loc_id, const char *groupname, size_t size_hint);
hid_t my_H5Dcreate(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id);
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
hid_t my_H5Screate_simple(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims);
herr_t my_H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf,
                   const char *datasetname);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
hid_t my_H5Fopen(const char *fname, unsigned int flags, hid_t fapl_id);
hid_t my_H5Dopen(hid_t file_id, const char *datasetname);

hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname);
herr_t my_H5Dread(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void *buf,
                  const char *datasetname);
hid_t my_H5Gopen(hid_t loc_id, const char *groupname);
hid_t my_H5Aopen_name(hid_t loc_id, const char *attr_name);
herr_t my_H5Aread(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Dclose(hid_t dataset_id, const char *datasetname);
herr_t my_H5Gclose(hid_t group_id, const char *groupname);
herr_t my_H5Fclose(hid_t file_id, const char *fname);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);

hid_t my_H5Tcopy(hid_t type_id);
herr_t my_H5Tclose(hid_t type_id);
herr_t my_H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count,
                              const hsize_t *block);
size_t my_H5Tget_size(hid_t datatype_id);
herr_t my_H5Tset_size(hid_t datatype_id, size_t size);
herr_t my_H5Sset_extent_simple(hid_t space_id, int rank, const hsize_t *current_size, const hsize_t *maximum_size,
                               const char *attr_name);
hid_t my_H5Dget_space(hid_t dataset_id, const char *datasetname);

#ifdef HDF5_FILTERS
htri_t my_H5Pall_filters_avail(hid_t plist_id);
hid_t my_H5Pcreate(hid_t class_id);
herr_t my_H5Pclose(hid_t plist);
herr_t my_H5Pset_chunk(hid_t plist, int ndims, const hsize_t *dim);
herr_t my_H5Pset_shuffle(hid_t plist_id);
herr_t my_H5Pset_deflate(hid_t plist_id, uint level);
herr_t my_H5Pset_fletcher32(hid_t plist_id);
#endif /* #ifdef HDF5_FILTERS */

#endif /* #ifdef HAVE_HDF5 */

#ifdef HOST_MEMORY_REPORTING
void check_maxmemsize_setting(void);
#endif /* #ifdef HOST_MEMORY_REPORTING */

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
int get_softening_type_from_mass(double mass);
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

#ifdef MHD
void do_mhd_source_terms_first_half(void);
void do_mhd_source_terms_second_half(void);
#endif /* #ifdef MHD */

#ifdef ONEDIMS_SPHERICAL
void gravity_monopole_1d_spherical();
#endif /* #ifdef ONEDIMS_SPHERICAL */

#if defined(PMGRID)
void my_slab_based_fft(fft_plan *plan, void *data, void *workspace, int forward);
void my_slab_based_fft_c2c(fft_plan *plan, void *data, void *workspace, int forward);
void my_slab_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
void my_slab_transposeA(fft_plan *plan, fft_real *field, fft_real *scratch);
void my_slab_transposeB(fft_plan *plan, fft_real *field, fft_real *scratch);
void my_column_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
void my_column_based_fft_init_c2c(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
void my_column_based_fft(fft_plan *plan, void *data, void *workspace, int forward);
void my_column_based_fft_c2c(fft_plan *plan, void *data, void *workspace, int forward);
void my_fft_swap23(fft_plan *plan, fft_real *data, fft_real *out);

void my_fft_swap13(fft_plan *plan, fft_real *data, fft_real *out);
void my_fft_swap23back(fft_plan *plan, fft_real *data, fft_real *out);
void my_fft_swap13back(fft_plan *plan, fft_real *data, fft_real *out);
#endif /* #if defined(PMGRID) */

#ifdef RIEMANN_HLLC
double godunov_flux_3d_hllc(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
#endif /* #ifdef RIEMANN_HLLC */

#if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD)
void flux_convert_to_lab_frame(struct state *st_L, struct state *st_R, double *vel_face, struct fluxes *flux);
#endif /* #if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD) */

#ifdef RIEMANN_HLLD
double godunov_flux_3d_hlld(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux);
#endif /* #ifdef RIEMANN_HLLD */

#ifdef SUBFIND_EXTENDED_PROPERTIES
void subfind_fof_calc_am_collective(int snapnr, int ngroups_cat);
int subfind_fof_calc_am_serial(int gr, int Offs, int snapnr, int ngroups_cat);
void subfind_add_grp_props_calc_fof_angular_momentum(int num, int ngroups_cat);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

#ifdef USE_SFR
void convert_cell_into_star(int i, double birthtime);
void spawn_star_from_cell(int igas, double birthtime, int istar, MyDouble mass_of_star);
void make_star(int idx, int i, double prob, MyDouble mass_of_star, double *sum_mass_stars);
#endif /* #ifdef USE_SFR */

#endif /* #ifndef PROTO_H */
