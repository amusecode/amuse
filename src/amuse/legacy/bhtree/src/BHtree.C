// 
// BHtree.C
// Version 1999/1/1 --- cleaned up and some comments  added.
//
// The tree construction/handling package designed for TREE
// implementation of SPH/NBODY programs
//
//
// This program uses a new tree construction method, based on
// Morton ordering and top-down tree construction.
//
// non-local functions (not complete yet...)
//
// setup_tree()
//    allocate the memory for tree if necessary
//    create the tree structure
// set_cm_quantities_for_default_tree()
//    calculate mass and position for tree nodes
//
// set_hmax_for_sph() (SPH use only)
// check_and_set_nbl() (SPH use only)
//
// calculate_gravity_using_tree() (NO GRAPE)
//
// evaluate_gravity_using_tree_and_list()(GRAPE)

extern "C" double cpusec();

#include "stdinc.h"
#include "BHtree.h"

#include  <stdlib.h>
#include  <math.h>
#include  <iostream>
#include  <cstdio>

// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>
long supported_conditions = COLLISION_DETECTION_BITMAP | TIMEOUT_DETECTION_BITMAP;
// AMUSE STOPPING CONDITIONS SUPPORT

void dump_octal(BHlong x)
{
    char st[256];
    sprintf(st," %lo ",x);
    cerr <<  st ;
}
    

BHlong conv_to_morton(int ix, int keybits)
{
    BHlong dum = 0;
    //    char st[256];
    //    cerr << "conv_to_morton "; PR(ix); PRL(keybits);
    //    sprintf(st,"%lo",ix);
    //    cerr << "ix = " << st << endl;
    int i, j;
    for(i = j= 0; i<keybits; i++,j+=3){
	if (ix & (1<<i)){
	    dum |= ((BHlong) 1)<<j;
	}
    }
    //sprintf(st,"%lo",dum);
    //    cerr << "dum = " << st << endl;
    return dum;

}
    
inline BHlong  construct_key(const vec & pos, real rscale, int ioffset, int keybits)
{
    long ix[3];
    vec p = pos;
    for(int i = 0; i<3; i++){
	ix[i] = (long) (p[i]*rscale+ioffset);
    }
    return (conv_to_morton(ix[0],keybits)<<2 )
	|(conv_to_morton(ix[1],keybits)<<1 )
	|(conv_to_morton(ix[2],keybits));
}
    

void bhparticle::set_key(real rscale, int ioffset, int keybits)
{
    key =  construct_key(rp->get_pos(), rscale, ioffset, keybits);
    //    PRL(key);
    
}

int compare_key(bhparticle * p1, bhparticle * p2)
{
    long comp = ((long) p1->get_key()) - ((long) p2->get_key());
    if (comp > 0L){
	return 1;
    }else if (comp == 0L){
	return 0;
    }else{
	return -1;
    }
}

void sort_bh_array( bhparticle * r, int lo, int up )
{
    int i, j;
    bhparticle tempr;
    while ( up>lo ) {
	i = lo;
	j = up;
	tempr = r[lo];
	/*** Split file in two ***/
	while ( i<j ) {
	    for ( ; r[j].key > tempr.key; j-- );
	    for ( r[i]=r[j]; i<j && r[i].key<=tempr.key; i++ );
	    r[j] = r[i];
	}
	r[i] = tempr;
	/*** Sort recursively, the smallest first ***/
	if ( i-lo < up-i ) { sort_bh_array(r,lo,i-1);  lo = i+1; }
	else    { sort_bh_array(r,i+1,up);  up = i-1; }
    }
}
void check_bh_array( bhparticle * r, int size )
{
    for(int i = 0; i<size-1;i++){
	if(r[i].get_key() > r[i+1].get_key()){
	    PR(i); PR(r[i].get_key()); PRL(r[i+1].get_key());
	    cerr << "Sort failed ... \n";
	    exit (1);
	}
    }
}

real initialize_key(int nbody,
		    real_particle * rp,
		    int & nkeysize,
		    bhparticle * &bhp)
{
    if (nbody > nkeysize || bhp == NULL){
	if (bhp != NULL){
	    delete [] bhp;
	}
	nkeysize = nbody+100;
	bhp = new bhparticle[nkeysize];
#ifdef REUSE_PREVIOS_DATA
	// With present quicksort routine, the pre-sorted data
	// actuallt DEGRADE its performance. So DO NOT ACTIVATE
	// THIS PART --- JM 1998/12/22
	for(int i = 0; i<nbody; i++){
	    bhparticle * p = bhp + i;
	    p->set_rp(rp+i);
	}
#endif	
    }
    real rmax = 1;
    for(int i = 0; i<nbody; i++){
	vec p = (rp+i)->get_pos();
	for (int k = 0; k<3; k++){
	    if (fabs(p[k])>=rmax) rmax *= 2;
	}
    }
    real rscale = 1.0/rmax*default_ix_offset;
	
    for(int i = 0; i<nbody; i++){
	bhparticle * p = bhp + i;
#ifndef REUSE_PREVIOS_DATA	
	p->set_rp(rp+i);
#endif	
	p->set_key(rscale, default_ix_offset, default_key_length);
	//	PR(i); PRL(p->get_key());
    }
//    cerr << "Call quicksort, cpu = " <<cpusec() << endl;
    sort_bh_array(bhp,0,nbody-1);
    //    qsort(bhp, nbody, sizeof(bhparticle), compare_key);
    // The private sort routine is for some unknow reason
    // much faster than qsort of the system for large N
//    cerr << "Exit quicksort, cpu = " <<cpusec() << endl;
    for(int i = 0; i<nbody; i++){
	// bhparticle * p = bhp + i;
	// PR(i); PR(p->get_key()); PRL(p->get_rp()->get_index());
    }
    return rmax;
}

void bhnode::assign_root(vec root_pos, real length, bhparticle * bp, int np)
{
    pos = root_pos;
    l = length;
    bpfirst = bp;
    nparticle = np;
}


    


void bhnode::create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
				   BHlong current_key,
				   int current_level,
				   int n_critical)
{
//    cerr << "create tree called "; PRC(nparticle); PRL(n_critical);
//    PRL(heap_remainder);
    if (heap_remainder <= 0){
	cerr << "create_tree: no more free node... exit\n";
	exit(1);
    }
    if (nparticle <= n_critical) return;
    if (current_level == 0) return;
    //    cerr << "Enter recursion\n";
    //    dump();
    BHlong keyscale = ((BHlong) 1)<<((current_level-1)*3);
    bhparticle * bptmp = bpfirst;
    int npremain = nparticle;
    for(int i=0; i<8;i++)child[i] = NULL;
    isleaf = 1;
    for(int i=0; i<8;i++){
	BHlong new_key = current_key + keyscale * i;
	vec new_pos = pos + vec( ((i&4)*0.5-1)*l/4,
				       ((i&2)    -1)*l/4,
				       ((i&1)*2  -1)*l/4);
	
	if(bptmp->get_key() - new_key <keyscale){
	    // current bptmp is actually in the current subnode
	    // search for the end location
	    int p0 = 0;
	    int p1 = npremain-1;
	    if ((bptmp+p1)->get_key() - new_key >=keyscale){
		while (p1 - p0 > 1){
		    int pnew = (p0+p1)/2;
		    if ((bptmp+pnew)->get_key() - new_key <keyscale){
			p0 = pnew;
		    }else{
			p1 = pnew;
		    }
		}
		p1 = p0;
	    }
	    p1 ++;
	    isleaf = 0;
	    child[i] = heap_top;
	    heap_top ++;
	    heap_remainder -- ;
	    child[i]->bpfirst = bptmp;
	    child[i]->pos = new_pos;
	    child[i]->l = l*0.5;
	    child[i]->nparticle = p1;
	    child[i]->isleaf = 1;
	    child[i]->create_tree_recursive(heap_top, heap_remainder,
					    new_key, current_level-1, n_critical);
	    bptmp += p1;
	    npremain -= p1;
	    if (npremain <= 0) return;
				  
	}
    }

    //dump();
}


void spc(int indent)
{
    for(int i=0;i<indent;i++)cerr << " ";
}

void bhnode::dump(int indent = 0)
{
    int i;
    spc(indent); cerr << "node pos " << pos ;
#ifdef SPH    
    cerr << " h " << hmax_for_sph;
#endif
    cerr << endl;
    spc(indent); cerr << "node cm  " << cmpos << " m " << cmmass ;
    if (isleaf){
	cerr << " IS LEAF" ;PRL(nparticle);
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    for(int j=0;j<indent+2;j++)cerr << " ";
	    real_particle * p = (bp+i)->get_rp();
	    PR(p->get_index()); PRL(p->get_pos());
	}
    }else{
	cerr << " IS _not_ LEAF ";PRL(nparticle);
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->dump(indent + 2);
	    }
	}
    }
}

// inbox: returns 0 if the particle is in the box;

int inbox(vec  & cpos, // center of the box
	  vec  & pos,  // position of the particle
	  real l)         // length of one side of the box
    
{
    for(int  i = 0; i< ndim; i++){
	if (fabs(pos[i]-cpos[i]) > l*0.5) return 1;
    }
    return 0;
}
	
	
int bhnode::sanity_check()
{
    int i;
    int iret = 0;
    if (isleaf){
	// this is the lowest level node. Things to check:
	// all particles are in the cell
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    vec ppos = p->get_pos();
	    if(inbox(pos,ppos,l)){
		cerr << "Error, particle out of box ... \n";
		dump();
		return 1;
	    }
	}
    }else{

	// This is the non-leaf node. Check the position and side
	// length of the child cells and then check recursively..
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		int err = 0;
	        err = child[i]->sanity_check();
		if (l*0.5 != child[i]->get_length()) err += 2;
		vec relpos = pos-child[i]->get_pos();
		for (int k = 0 ; k<ndim;k++){
		    if (fabs(relpos[k]) !=l*0.25)err += 4;
		}
		if (err){
		    cerr << "Child " << i << " Error type = " << err << endl;
		    dump();
		}
		iret += err;
	    }
	}
    }
    return iret;
}

#ifdef SPH	
void  bhnode::set_hmax_for_sph()
{
    int i;
    hmax_for_sph = 0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real hp = (bp+i)->get_rp()->get_h();
	    if(hmax_for_sph < hp) hmax_for_sph = hp;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_hmax_for_sph();
		if (hmax_for_sph < child[i]->hmax_for_sph)
		    hmax_for_sph = child[i]->hmax_for_sph;
	    }
	}
    }
}
#endif

void  bhnode::set_cm_quantities()
{
    int i;
    cmpos = 0.0;
    cmmass = 0.0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    cmpos += mchild*(bp+i)->get_rp()->get_pos();
	    cmmass += mchild;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_cm_quantities();
		real mchild = child[i]->cmmass;
		cmpos += mchild*child[i]->cmpos;
		cmmass += mchild;
	    }
	}
    }
    cmpos /= cmmass;
}

real separation_squared(bhnode * p1,bhnode * p2)
{
    real r2 = 0;
    real xmin = (p1->get_length()+ p2->get_length())*0.5;
    vec dx = p1->get_pos() - p2->get_pos();
    for (int k = 0; k<ndim; k++){
	real adx = fabs(dx[k]);
	if (adx > xmin){
	    adx -= xmin;
	}else{
	    adx = 0;
	}
	r2 += adx*adx;
    }
    return r2;
}

real separation_squared(bhnode * p1, vec & pos2)
{
    real r2 = 0;
    real xmin = p1->get_length()*0.5;
    vec dx = p1->get_pos() - pos2;
    for (int k = 0; k<ndim; k++){
	real adx = fabs(dx[k]);
	if (adx > xmin){
	    adx -= xmin;
	}else{
	    adx = 0;
	}
	r2 += adx*adx;
    }
    return r2;
}

int  are_overlapped(bhnode * p1,bhnode * p2)
{
    real xmin = (p1->get_length()+ p2->get_length())*0.4999999999999999;
    vec dx = p1->get_pos() - p2->get_pos();
    for (int k = 0; k<ndim; k++){
	if(fabs(dx[k]) > xmin) return 0;
    }
    return 1;
}
#ifdef SPH
int check_and_set_nbl(bhnode * p1,bhnode * p2);
int check_and_set_nbl(bhnode * p1,bhnode * p2)
{
    int iret = 0;
    if(p1 == NULL) return iret;
    if(p2 == NULL) return iret;
    real rcrit = p1->hmax_for_sph;
    if (rcrit <  p2->hmax_for_sph)rcrit = p2->hmax_for_sph;
    rcrit *= 2;
    real rcrit2 = rcrit*rcrit;
    if(separation_squared(p1,p2) > rcrit2) return iret;
    if (p1->isleaf == 0 || p2->isleaf == 0){
	//
	// either node is not leaf. Go down
	//
	if (p1->isleaf || (p2->isleaf == 0) && p2->l > p1->l){
	    //
	    // p1 is already leaf, or both are node but p2 is larger.
	    // go down p2 by 
	    for(int i = 0; i<8;i++)
		iret |= check_and_set_nbl(p1, p2->child[i]);
	}else{
	    //
	    // now, we can go down p1 ...
	    for(int i = 0; i<8;i++)
		iret |= check_and_set_nbl(p1->child[i], p2);
	}
    }else{
	//
	// both are leaves; can directly evaluate particles
	bhparticle * bpi = p1->bpfirst;
	bhparticle * bpj = p2->bpfirst;
	for(int i = 0; i < p1->nparticle; i++)
	    for(int j = 0; j < p2->nparticle; j++){
		if((bpi+i)->get_key() <(bpj+j)->get_key() ){
		    //
		    // here, comparison of key guarantee that a pair
		    // is evaluated only once
	  
		    iret |=check_and_set_nbl(*((bpi+i)->get_rp()),
					     *((bpj+j)->get_rp()));
		}
	    }
    }
    return iret;
}

#endif
static bhparticle * bp = NULL;
static int bhpsize = 0;
static int bnsize = 0;
static bhnode * bn;
void set_cm_quantities_for_default_tree()
{
    bn->set_cm_quantities();
}


void real_system::setup_tree()
{
    real rsize = initialize_key(n,get_particle_pointer(),bhpsize,bp);
    if (bnsize < bhpsize/2){
	if (bnsize != 0){
	    delete [] bn;
	}
	bnsize = (int)(bhpsize*0.4+100);
	bn = new bhnode[bnsize];
    }
    for(int j = 0; j<bnsize;j++) (bn+j)->clear();
    bn->assign_root(vec(0.0), rsize*2, bp, n);
    bhnode * btmp = bn+1;
    int heap_remainder = bnsize-1;
    BHlong key = 0;
    
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, ncrit_for_tree);
    //cello: last argument was 12 for some reason 
                  
    
    set_cm_quantities_for_default_tree();
//    PR(bnsize);    PRL(heap_remainder);
    //    PRL(bn->sanity_check());
}

	
#ifdef SPH	
int sph_system::set_nnb_using_tree()
{
    setup_tree();
    real_particle * psph = get_particle_pointer();
    apply_vf(real_particle::clear_nnb);
    bn->set_hmax_for_sph();
    int iret = check_and_set_nbl(bn, bn);
    apply_vf(real_particle::sort_nblist);
    return iret;
}
#endif

void accumulate_force_from_point(vec dx, real r2, real eps2, 
				 vec & acc,
				 real & phi,
				 real jmass)
{
    double r2inv = 1/(r2+eps2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    phi -= jmass*rinv;
    acc += jmass*r3inv*dx;
}

static real total_interactions;
static int tree_walks;
static int nisum;
void clear_tree_counters()
{
    total_interactions = 0;
    tree_walks = 0;
    nisum = 0;
}
void print_tree_counters()
{
    real avg = total_interactions/nisum;
    PRC(nisum); PRC(tree_walks); PRC(total_interactions); PRL(avg);
    cout <<"tree_walks = " << tree_walks << " ntaverage = " << avg << endl;
}

void calculate_force_from_interaction_list(const vec & pos,
					   real eps2, 
					    vec & acc,
					    real & phi,
					    vec * poslist,
					    real * masslist,
					    int list_length)
{
    acc = 0.0;
    phi = 0.0;
    for(int i = 0; i<list_length; i++){
	vec dx = *(poslist+i)-pos;
	real r2 = dx*dx;
	accumulate_force_from_point(dx, r2, eps2, acc, phi,*(masslist+i));
    }
}

int my_index;		// horrible hack to provide tree traversal
int nn_index;		// data to the calling function below!
real_particle *nn_ptr;
real r_nn_2;

void bhnode::accumulate_force_from_tree(vec & ipos, real eps2, real theta2,
					vec & acc,
					real & phi)
{
    vec dx = cmpos - ipos;
    real r2 = dx*dx;
    if (r2*theta2 > l*l){
	// node and position is well separated;
	accumulate_force_from_point(dx, r2, eps2, acc, phi, cmmass);
    }else{
	int i;
	if (isleaf){
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		vec dx = (bp+i)->get_rp()->get_pos()-ipos;
		real r2 = dx*dx;

		// Added by Steve (8/07).

		if (r2 < r_nn_2) {
		    int iindex = (bp+i)->get_rp()->get_index();
		    if (iindex != my_index) {
		      r_nn_2 = r2;
		      nn_index = iindex;
		      nn_ptr = (bp+i)->get_rp();
		    }
		}

		accumulate_force_from_point(dx, r2, eps2, acc, phi,
					    (bp+i)->get_rp()->get_mass());
		total_interactions += 1;
	    }
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		child[i]->accumulate_force_from_tree(ipos, eps2, theta2,
						     acc, phi);
		}
	    }
	}
    }
}

void bhnode::add_to_interaction_list(bhnode & dest_node, real theta2,
				     vec * pos_list,
				     real * mass_list,
				     int & nlist,
				     int list_max,
				     int & first_leaf)
{
#if 0    
    real r2 = separation_squared(this, &dest_node);
    if (r2*theta2 > l*l){
#else
    if(!are_overlapped(this,&dest_node) && (separation_squared(&dest_node,cmpos)*theta2 > l*l)){
#endif
	// node and position is well separated;
	*(pos_list+nlist) = cmpos;
	*(mass_list+nlist) = cmmass;
	nlist ++;
	if (nlist > list_max){
	    cerr << "List length exceeded\n";
	    exit(1);
	}
	
    }else{
	int i;
	if (isleaf || (this == (&dest_node))){
	    if (this == (&dest_node)){
		// adding the particles in the node itself
		first_leaf = nlist;
	    }
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		*(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
		*(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
		nlist ++;
		if (nlist > list_max){
		    cerr << "List length exceeded\n";
		    exit(1);
		}
	    }
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		    child[i]->add_to_interaction_list(dest_node, theta2,
						      pos_list, mass_list,
						      nlist, list_max,
						      first_leaf);
		}
	    }
	}
    }
}



#ifdef HARP3

extern "C" void h3open_();
extern "C" void h3close_();
extern "C" void accel_by_harp3_separate_noopen_(int * ni, vec * xi,
						int * nj, vec * xj,
						real *m,
						vec *  a,
						real *p,
						real * eps2);

void calculate_force_from_interaction_list_using_grape4(vec * pos_list, real * mass_list,
							int list_length, int first_leaf, int ni,
							real eps2,
							vec * acc_list, real * phi_list)
{
    static call_count = 0;
    static h3_open_state = 0;
    if (h3_open_state == 0){
	h3open_();
	h3_open_state = 1;
    }
    //    PR(ni);PRL(list_length);
    nisum += ni;
    tree_walks += 1;
    total_interactions += ((real)ni)*list_length;
    accel_by_harp3_separate_noopen_(&ni,pos_list+first_leaf, &list_length,pos_list, mass_list,
				    acc_list, phi_list, &eps2);
    call_count += ni;
    if (call_count > 500000){
	cerr << "Close and release GRAPE-4\n";
	h3close_();
	h3_open_state = 0;
	call_count = 0;
    }
}
#endif

void bhnode::evaluate_gravity_using_tree_and_list(bhnode & source_node,
						  real theta2,
						  real eps2,
						  int ncrit)
{
    const int list_max = 40000;
    static real mass_list[list_max];
    static vec pos_list[list_max];
    real epsinv = 1.0/sqrt(eps2);
#ifdef HARP3    
    static vec * acc_list = NULL;
    static real * phi_list = NULL;
    if (acc_list == NULL){
	acc_list = new vec[ncrit + 100];
	phi_list = new real[ncrit + 100];
    }
#endif
    //    PR(pos); PR(nparticle); PRL(isleaf);
    if((nparticle > ncrit) && (isleaf==0)){
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->evaluate_gravity_using_tree_and_list(source_node,
							       theta2,
							       eps2,
							       ncrit);
	    }
	}
    }else{
	//
	// node is below critical ... first create list
	//
	int list_length = 0;
	int first_leaf = -1;
	source_node.add_to_interaction_list(*this,  theta2,
					    pos_list,
					    mass_list,
					    list_length,
					    list_max,
					    first_leaf);
	if (first_leaf == -1){
	    cerr << "evaluate_gravity: impossible error \n";
	    cerr << "failed to find the node in the tree \n";
	    exit(1);
	}
	bhparticle * bp = bpfirst;
#ifndef HARP3	
	for(int i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    vec acc;
	    real phi;
	    calculate_force_from_interaction_list(pos_list[i+first_leaf],eps2, acc, phi,
					  pos_list,mass_list,list_length);
	    p->set_acc_gravity(acc);
	    p->set_phi_gravity(phi + p->get_mass()*epsinv);
	}
#else
	calculate_force_from_interaction_list_using_grape4(pos_list, mass_list,list_length, first_leaf,
							   nparticle, eps2, acc_list, phi_list);
	for(int i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    p->set_acc_gravity(acc_list[i]);
	    p->set_phi_gravity(phi_list[i] + p->get_mass()*epsinv);
	}
#endif	
    }
}

void evaluate_gravity_using_default_tree_and_list(real theta2,
					  real eps2,
					  int ncrit)
{
    bn->evaluate_gravity_using_tree_and_list(*bn, theta2,eps2, ncrit);
}

int id_collision_1, id_collision_2;	// used extern in ../muse_dynamics.C
real r_collision_12;

void real_particle::calculate_gravity_using_tree(real eps2, real theta2)
{
    acc_gravity = 0;
    phi_gravity = mass/sqrt(eps2);

    my_index = index;	// global: defined above!
    nn_index = -1;
    r_nn_2 = 1.e30;

    bn->accumulate_force_from_tree(pos, eps2, theta2,
				   acc_gravity, phi_gravity);
    nisum += 1;


    
// AMUSE STOPPING CONDITIONS SUPPORT (AVE May 2010)
    if(nn_index >= 0 && (COLLISION_DETECTION_BITMAP & enabled_conditions) ) {
        real rad1 = radius;
        real rad2 = nn_ptr->get_radius();
        if (r_nn_2 <= pow(rad1+rad2, 2))
        {
            int stopping_index  = next_index_for_stopping_condition();
            set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
            set_stopping_condition_particle_index(stopping_index, 0, my_index);
            set_stopping_condition_particle_index(stopping_index, 1, nn_index);
        }
    }
// AMUSE STOPPING CONDITIONS SUPPORT

}


//AVE Mar 2010
vec real_system::calculate_gravity_at_point(vec pos, real eps2, real theta2)
{
    vec acc_gravity = 0;
    real phi_gravity = mass/sqrt(eps2);

    nn_index = -1;
    r_nn_2 = 0.0; // turn-off collision detection

    bn->accumulate_force_from_tree(pos, eps2, theta2, acc_gravity, phi_gravity);
    
    return acc_gravity;
}

//AVE Mar 2010
real real_system::calculate_potential_at_point(vec pos, real eps2, real theta2)
{
    vec acc_gravity = 0;
    real phi_gravity = 0.0;

    nn_index = -1;
    r_nn_2 = 0.0; // turn-off collision detection

    bn->accumulate_force_from_tree(pos, eps2, theta2, acc_gravity, phi_gravity);
      
    return phi_gravity;
}


#ifdef TESTXXX
//
// do the serious test of
// construction of tree
// consistency of tree
// validity of the neighbour list (by comparing with the result
// of direct calculation

int main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    pb.initialize_h_and_nbl(pow(1.0/n,0.33333));
    static sph_system pbcopy = pb;
    copy_sph_particles(&pb, &pbcopy);
    real_particle * psph = pb.get_particle_pointer();
    real_particle * psphcopy = pbcopy.get_particle_pointer();
    pb.set_nnb_using_tree();
    cerr << "Dumping copy ... \n";
    cerr << "checking NB \n";
    int error = 0;
    for(int i = 0; i<n; i++){
	(psph+i)->sort_nblist();
	int err = 0;
	if((psph+i)->get_nnb() != (psphcopy+i)->get_nnb()){
	    cerr << "Neighbour count differs for "; PRL(i);
	    err = 1;
	    
	}
	if (err == 0){
	    for(int j = 0; (j< (psph+i)->get_nnb()) && (err == 0); j++){
		if ((psph+i)->get_neighbor(j)->get_index()!=
		    (psphcopy+i)->get_neighbor(j)->get_index()) err = 1;
	    }
	}
	if(err){
	    (psph+i)->dump();
	    (psphcopy+i)->dump();
	    error ++;
	}
    }
    PRL(error);
}
#endif

#ifdef TEST
//
// do the serious test of
// construction of tree
// consistency of tree
// validity of the neighbour list (by comparing with the result
// of direct calculation

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    bhparticle * bp = NULL;

    bn = new bhnode[n];
    for(int i = 0; i<1; i++){
	real rsize = initialize_key(n,pb.get_particle_pointer(),nkey,bp);
	for(int j = 0; j<n;j++) (bn+j)->clear();
	bn->assign_root(vec(0.0), rsize*2, bp, n);
	bhnode * btmp = bn+1;
	int heap_remainder = n-1;
	BHlong key = 0;
        bn->create_tree_recursive(btmp,  heap_remainder,key, default_key_length, 4);
    }
    PRL(bn->sanity_check());
    pb.initialize_h_and_nbl(pow(1.0/n,0.33333));
    bn->set_hmax_for_sph();
    //    bn->dump();
    bn->set_cm_quantities();
    //    bn->dump();
    static sph_system pbcopy = pb;
    copy_sph_particles(&pb, &pbcopy);
    real_particle * psph = pb.get_particle_pointer();
    real_particle * psphcopy = pbcopy.get_particle_pointer();
    for(int i = 0; i<n; i++){
	(psph+i)->clear_nnb();
    }
    PRL(check_and_set_nbl(bn, bn));
    cerr << "Dumping copy ... \n";
    cerr << "checking NB \n";
    int error = 0;
    for(int i = 0; i<n; i++){
	(psph+i)->sort_nblist();
	int err = 0;
	if((psph+i)->get_nnb() != (psphcopy+i)->get_nnb()){
	    cerr << "Neighbour count differs for "; PRL(i);
	    err = 1;
	    
	}
	if (err == 0){
	    for(int j = 0; (j< (psph+i)->get_nnb()) && (err == 0); j++){
		if ((psph+i)->get_neighbor(j)->get_index()!=
		    (psphcopy+i)->get_neighbor(j)->get_index()) err = 1;
	    }
	}
	if(err){
	    (psph+i)->dump();
	    (psphcopy+i)->dump();
	    error ++;
	}
    }
    PRL(error);
    pb.use_self_gravity = 1;
    pb.eps2_for_gravity = 0.01;
#define COMPARISON_WITH_DIRECT    
#ifdef COMPARISON_WITH_DIRECT    
    pb.calculate_uncorrected_gravity_direct();
    copy_sph_particles(&pb, &pbcopy);
    psphcopy = pbcopy.get_particle_pointer();
    cerr << "Direct force \n";
    for(int i = 0; i<n; i++){
	real phi = (psphcopy+i)->get_phi_gravity();
	vec acc  = (psphcopy+i)->get_acc_gravity();
	PR(i); PR(phi); PRL(acc);
    }
#endif


    
    cerr << "Tree   force \n";
    for(int j = 0; j<10; j++){
	PRL(j);
	pb.apply_vf(real_particle::clear_acc_phi_gravity);
	for(int i = 0; i<n; i++){
	    (psph+i)->calculate_gravity_using_tree(pb.eps2_for_gravity, 0.4);
	}
    }
    pb.apply_vf(real_particle::clear_acc_phi_gravity);
    bn->evaluate_gravity_using_tree_and_list(*bn,0.4,pb.eps2_for_gravity,1);
#ifdef COMPARISON_WITH_DIRECT    
    real perrmax = 0;
    real ferrmax = 0;
    for(int i = 0; i<n; i++){
	real phi = (psph+i)->get_phi_gravity();
	real phierr = (psphcopy+i)->get_phi_gravity()-phi;
	vec acc  = (psph+i)->get_acc_gravity();
	vec accerr  = (psphcopy+i)->get_acc_gravity()-acc;
	PR(i); PR(phi); PRC(acc); PRC(phierr); PRL(accerr);
	real prelerr = fabs(phierr/phi);
	real frelerr = abs(accerr)/abs(acc);
	if(perrmax < prelerr) perrmax = prelerr;
	if(ferrmax < frelerr) ferrmax = frelerr;
    }
    PR(perrmax);    PRL(ferrmax);
#else
    for(int i = 0; i<n; i++){
	real phi = (psph+i)->get_phi_gravity();
	vec acc  = (psph+i)->get_acc_gravity();
	PR(i); PR(phi); PRL(acc); 
    }
	
#endif    
    
}
#endif

	

#ifdef TESTXX
//
// Sample test for timing purpose...

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    bhparticle * bp = NULL;

    bn = new bhnode[n];
    for(int i = 0; i<10; i++){
	real rsize = initialize_key(n,pb.get_particle_pointer(),nkey,bp);
	for(int j = 0; j<n;j++) (bn+j)->clear();
	bn->assign_root(vec(0.0), rsize*2, bp, n);
	bhnode * btmp = bn+1;
	int heap_remainder = n-1;
	BHlong key = 0;
        bn->create_tree_recursive(btmp,
				  heap_remainder,key,
				  default_key_length, 8 );
	PRL(heap_remainder);
    }
    PRL(bn->sanity_check());
    real_particle * psph = pb.get_particle_pointer();
    real h0 = pow(1.0/n,0.33333);
    for(int i = 0; i<10; i++){
	
	pb.apply_vf(real_particle::set_h, h0);
	pb.apply_vf(real_particle::clear_nnb);
	bn->set_hmax_for_sph();
	//    bn->dump();
	PRL(check_and_set_nbl(bn, bn));
	pb.apply_vf(real_particle::sort_nblist);
		
    }
}
#endif






