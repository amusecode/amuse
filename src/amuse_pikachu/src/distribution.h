#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include"sort.h"
#include"const.h"
#include"mpi_interface.h"
#include"BHtree.h"

extern void gdr_open();
extern void gdr_close();
extern void halt_program();


inline int determine_sample_freq(const int &Nprt_global,
				 const int &maxsample){
    int Nreal = Nprt_global;
    int sample_freq = (Nreal+maxsample-1)/maxsample;
    mpi_bcast_T(&sample_freq, 1);
    return sample_freq;
}

inline int determine_sample_freq(const int &Nprt_global){
    int Nreal = Nprt_global;
    int maxsample = (int)(NSAMPLE_MAX*0.8); // 0.8 is safety factor
    int sample_freq = (Nreal+maxsample-1)/maxsample;
    mpi_bcast_T(&sample_freq, 1);
    return sample_freq;
}

inline void create_division(const int n, int &nx, int &ny, int &nz){
    static double over3 = 1.0/3.0;
    int n0, n1;
    n0 = (int)pow(n+0.1, over3);
    while(n%n0)n0--;
    nx = n0;
    n1 = n/nx;
    n0 = (int)sqrt(n1+0.1);
    while(n1%n0)n0++;
    ny = n0; nz = n1/n0;
    int ntmp;
    if (nz > ny){
	ntmp = nz; nz = ny; ny = ntmp;
    }
    if (ny > nx){
	ntmp = nx; nx = ny; ny = ntmp;
    }
    if (nz > ny){
	ntmp = nz; nz = ny; ny = ntmp;
    }
    if (nx*ny*nz != n){
	cerr << "create_division: Intenal Error " << n << " " << nx
	     << " " << ny << " " << nz <<endl;
    }
}

inline void sort_coord_array( Vector3 *pos, int lo, int up, int cid){
  int i, j;
  Vector3 postmp;
  while( up>lo ){
    i = lo;
    j = up;
    postmp = pos[lo];
    /*** Split file in two ***/
    while ( i<j ) {
      for ( ; pos[j][cid] > postmp[cid]; j-- );
      for ( pos[i]=pos[j]; i<j && pos[i][cid]<=postmp[cid]; i++ );
      pos[j] = pos[i];
    }
    pos[i] = postmp;
    /*** Sort recursively, the smallest first ***/
    if ( i-lo < up-i ) { sort_coord_array(pos,lo,i-1,cid);  lo = i+1; }
    else    { sort_coord_array(pos,i+1,up,cid);  up = i-1; }
  }
}

inline void calculate_boxdim(const int &Nsampletot, Vector3 pos[], 
		      const int &cid, const int &istart, 
		      const int &iend, const double &rmax,
                      double &xlow, double &xhigh){
  if(istart == 0) {xlow = -rmax;}
  else{xlow = (pos[istart][cid]+ pos[istart-1][cid])*0.5;}
  if(iend == Nsampletot-1) {xhigh = rmax;}
  else{xhigh = (pos[iend][cid]+ pos[iend+1][cid])*0.5;}
}

inline void collect_sample_particles(Particle prt[], 
				     const int &Nprt_local, 
				     const int &sample_freq, 
				     int &Nsample_global, 
				     Vector3 sample_coord[]){ 
  int Nsample_local = 0;
  int Nproc = mpi_get_size();
  static int Nsample_array[NPROC_MAX];
  static int Nsample_disp[NPROC_MAX+1];
  for(int i=0; i<Nprt_local; i+=sample_freq){
    sample_coord[Nsample_local] = prt[i].pos;
    Nsample_local++;
  }
  MPI_Allgather(&Nsample_local, 1, MPI_INT, Nsample_array, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allreduce(&Nsample_local, &Nsample_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  if(mpi_get_rank() == 0){
    cerr<<"Nsample_global="<<Nsample_global<<endl;
  }
  Nsample_disp[0] = 0;
  for(int i=0; i<Nproc; i++){
    Nsample_disp[i+1] = Nsample_disp[i] + Nsample_array[i];
  }
  Vector3 sample_coord_recv[NSAMPLE_MAX];
  mpi_gatherv_T(sample_coord, Nsample_local, sample_coord_recv, Nsample_array, Nsample_disp);
  for(int i=0; i<NSAMPLE_MAX; i++){
    sample_coord[i] = sample_coord_recv[i];
  }
}


// ramx is half lengh of a box.
// don't initialize
inline void calc_rmax(Particle prt[],  
		      double &rmax, 
		      const int &Nprt_local){
  double rmax_local = rmax;
  for(int i=0; i<Nprt_local; i++){
    for(int k=0; k<3; k++){
      if(fabs(prt[i].pos[k]) > rmax_local){
	rmax_local = fabs(prt[i].pos[k]);
      }
    }
  }
  rmax = mpi_max_double(&rmax_local)*1.00000000001;
}

inline void determine_division(Vector3 sample_coord[], 
			       const int &Nsample_global, 
			       const int npdim[], 
			       const double &rmax, 
			       Vector3 xlow[], 
			       Vector3 xhigh[]){
  int Nproc = npdim[0]*npdim[1]*npdim[2];
  if(mpi_get_rank() == 0){
    sort_coord_array(sample_coord, 0, Nsample_global-1, 0);
    int istart[NPROC_MAX]; // istart th particles is first particle in the ith box
    int iend[NPROC_MAX]; // iend th particles is last particle in the ith box
    for(int i=0;i<Nproc; i++){
      istart[i] = (i*Nsample_global)/Nproc;
      if(i>0) iend[i-1]=istart[i]-1;
      //if(i>=0) iend[i-1]=istart[i]-1;
    }
    iend[Nproc-1] = Nsample_global-1;

    for(int ix=0; ix<npdim[0]; ix++){
      double x0=0.0;
      double x1=0.0;
      int ix0 = ix*npdim[1]*npdim[2]; // first box id in given x region
      int ix1 = (ix+1)*npdim[1]*npdim[2]; // first box id in next given x region
      int Nx = Nsample_global;
      calculate_boxdim(Nx,  sample_coord,  0, 
		       istart[ix0], iend[ix1-1],  rmax,  x0,  x1);
      for(int i=ix0; i<ix1; i++){
	xlow[i][0]=x0;
	xhigh[i][0]=x1;
      }
    }

    for(int ix=0; ix<npdim[0]; ix++){
      int ix0 = ix*npdim[1]*npdim[2]; // first box id in given x region
      int ix1 = (ix+1)*npdim[1]*npdim[2]; // first box id in next given x region
      int Ny = iend[ix1-1] - istart[ix0] + 1;
      sort_coord_array(sample_coord, istart[ix0],iend[ix1-1], 1);
      for(int iy=0; iy<npdim[1]; iy++){
	double y0, y1;
	int iy0 = ix0+(iy*npdim[2]); // first box id in given x region and given y region.
	int iy1 = ix0+((iy+1)*npdim[2]); // first box id in given x region and next given y region.
	calculate_boxdim(Ny, sample_coord+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], rmax, y0,y1);
	for(int i=iy0; i<iy1; i++){
	  xlow[i][1]=y0;
	  xhigh[i][1]=y1;
	}
      }
    }
    for(int ix=0; ix<npdim[0]; ix++){
      int ix0 = ix*npdim[1]*npdim[2];
      for(int iy=0; iy<npdim[1]; iy++){
	int iy0 = ix0+(iy*npdim[2]);
	int iy1 = ix0+((iy+1)*npdim[2]);
	int Nz = iend[iy1-1] - istart[iy0] + 1;
	sort_coord_array(sample_coord, istart[iy0], iend[iy1-1], 2);
	for(int iz=0; iz<npdim[2]; iz++){
	  double z0, z1;
	  int iz0 = iy0+iz;
	  calculate_boxdim(Nz, sample_coord+istart[iy0], 2, istart[iz0]-istart[iy0],
			   iend[iz0]-istart[iy0], rmax, z0,z1);
	  xlow[iz0][2]=z0;
	  xhigh[iz0][2]=z1;
	}
      }
    }
  }
  mpi_bcast_T(xlow, Nproc);
  mpi_bcast_T(xhigh, Nproc);
}


inline int isinbox(Vector3 pos, Vector3 xlow, Vector3 xhigh){
  int inbox = 1;
  for(int k = 0; k<3;k++){
    if((pos[k] <xlow[k])||(pos[k] >xhigh[k])){
      inbox = 0;
      break;
    }
  }
  return inbox;
}


inline int mpi_exchange_particle(int box_dest, 
				 Particle *prt, 
				 int firstloc, 
				 const int &Nprt_block, // particle to be sent
				 int box_source, 
				 int &i_local){ 
  static MPI_Status status;
  // last particles sent initially (e.g. initially, 0 node send to Nproc-1 node)
  // received particles contained backward (e.g. initially, node contain particles into Nprt_loca_max-1 address)
  int myrank = mpi_get_rank();
  // exchange available buffer size...
  // initially i_local = Nprt_local_max,  firstloc is first address of sent particles, Nprt_block is the number of sent particles
  int mybuffsize = i_local - (firstloc+Nprt_block);
  int buffsize;

  // my buff size sent to source node
  MPI_Sendrecv(&mybuffsize, 1, MPI_INT, box_source, myrank*10+9,
	       &buffsize, 1, MPI_INT, box_dest, box_dest*10+9, MPI_COMM_WORLD, &status);
  int Nsend = Nprt_block;
  if(buffsize < Nsend){
    return 1;
  }
  //first send&get the number of particles to send&get
  int Nrecv;
  MPI_Sendrecv(&Nsend, 1, MPI_INT, box_dest, myrank*10,
	       &Nrecv, 1, MPI_INT, box_source, box_source*10, MPI_COMM_WORLD, &status);
  i_local -= Nrecv;
  int size_of_particle = sizeof(Particle);
  MPI_Sendrecv(prt+firstloc,  size_of_particle*Nsend, MPI_BYTE, box_dest, myrank*10+1,
	       prt+i_local, size_of_particle*Nrecv, MPI_BYTE, box_source, box_source*10+1, MPI_COMM_WORLD, &status);
  return 0;
} 

inline int exchange_particles(Particle prt[],
			      int &Nprt_local, 
			      const int &Nprt_local_max, 
			      Vector3 xlow[], 
			      Vector3 xhigh[],
			      int DST_RANK[],
			      int SRC_RANK[]){
    // Initially, Nprt_local particles are divide into Nproc
    // Particles in 1st block don't go anywhere, those in 2nd block go to myrank+1 th node, those in 3rd block go to myrank+2 th node ....
    int myrank = mpi_get_rank();
    int Nproc = mpi_get_size();
    int i_local = 0;
    int totalsent = 0;
    Particle prt_tmp;
    int firstloc[NPROC_MAX]; // first particle address in each block.
    int Nprt_block[NPROC_MAX]; // the number of particles in each block.
    // Loop over particles and determine which particle wants to go where
    // In ith box, particles in first block (prt[0]-prt[NFS_LOC+NBH_LOC])
    for(int ib=0; ib<Nproc; ib++){
        int ibox = ib == 0 ? MYRANK : DST_RANK[ib-1];
        firstloc[ibox] = i_local;
        for(int i=i_local; i<Nprt_local; i++){
            if(isinbox(prt[i].pos, xlow[ibox], xhigh[ibox])){
                prt_tmp = prt[i_local];
                prt[i_local] = prt[i];
                prt[i] = prt_tmp;
                i_local++;
            }
        }
        Nprt_block[ibox] = i_local-firstloc[ibox];
    }
    totalsent = mpi_sum_int_np(Nprt_local - Nprt_block[myrank]);
    if (myrank == 0){
        cerr << "total number of exchanged particles (is integer?)= " << totalsent/2 << endl;
    }
    i_local = Nprt_local_max;
    for(int ib=0; ib<Nproc-1; ib++){
        int box_dest = DST_RANK[ib];
        int box_source = SRC_RANK[ib]; 
        if( mpi_exchange_particle(box_dest, prt, firstloc[box_dest], Nprt_block[box_dest], box_source, i_local) ){
            cerr<<"local particles is overflowed (in function exchange_particles_with_overflow_check)"<<endl;
            halt_program();
        }
    }
    int is, id;
    int idfirst = Nprt_block[myrank];
    // packing particlese which come from other node
    // first idfirst Particles are original particles in mynode.
    // next particles come from other node.
    for(id=idfirst, is=i_local; is<Nprt_local_max; is++,id++){
        prt[id] = prt[is];
    }
    Nprt_local = id;
    
    return 0;
}


inline void distribute_particle(Particle prt_local[],
				Particle BH_glb[],
				const int &NBH_global_org,
				int &NFS_local,
				int &NBH_local,
				const int &Nprt_local_max,
				const int npdim[3],
				const int &sample_freq,
				Vector3 xlow[],
				Vector3 xhigh[],
				int DST_RANK[],
				int SRC_RANK[],
				double &rmax){


    int Nsample_global = 0;
    int Nprt_local = NFS_local + NBH_local;
    static Vector3 sample_coord[NSAMPLE_MAX];

    double toffset = mpi_wtime();
    collect_sample_particles(prt_local,  Nprt_local,
                             sample_freq,   Nsample_global, 
                             sample_coord);
    double tcollect = mpi_wtime() - toffset;

    toffset = mpi_wtime();
    rmax = 0.0;
    calc_rmax(prt_local,  rmax,  Nprt_local);
    determine_division(sample_coord,  Nsample_global,  npdim,  rmax, xlow, xhigh);
    double tdivision = mpi_wtime() - toffset;

    toffset = mpi_wtime();
    exchange_particles(prt_local,  Nprt_local,  Nprt_local_max, xlow, xhigh, DST_RANK, SRC_RANK);
    double texprt = mpi_wtime() - toffset;


    toffset = mpi_wtime();
    NBH_local = 0;
    Particle prt_tmp;
    for(int i=0; i<Nprt_local; i++){
        prt_local[i].node_org = MYRANK;
        if(prt_local[i].index < NBH_global_org){
            prt_tmp = prt_local[NBH_local];
            prt_local[NBH_local] = prt_local[i];
            prt_local[i] = prt_tmp;
            NBH_local++;
        }
    }
    double tBH0 = mpi_wtime() - toffset;

    toffset = mpi_wtime();
    for(int i=0; i<NBH_global_org; i++){
        for(int box=0; box<NPROC; box++){
            if(isinbox(BH_glb[i].pos, xlow[box], xhigh[box]) ){
                BH_glb[i].node_org = box;
                break;
            }
        }
    }
    NFS_local = Nprt_local - NBH_local;
    double tBH1 = mpi_wtime() - toffset;


    if(MYRANK == 0){
	cerr<<"tcollect="<<tcollect<<endl;
	cerr<<"tdivision="<<tdivision<<endl;
	cerr<<"texprt="<<texprt<<endl;
	cerr<<"tBH0="<<tBH0<<endl;
	cerr<<"tBH1="<<tBH1<<endl;
    }

}

#endif //DISTRIBUTION_H
