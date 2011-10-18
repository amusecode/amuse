#include "sapporohostclass.h"

struct real3 {
  double x, y, z;
};

int main(int argc, char *argv[]) {
  int n = 1024;
  if (argc > 1) n = atoi(argv[1]);
  cerr << " n = " << n << endl;

  int nx = 1;
  int n1 = 0;
  int n2 = n;

  double (*pos)[3] = new double[n][3];
  double (*vel)[3] = new double[n][3];
  double (*acc)[3] = new double[n][3];
  double (*jrk)[3] = new double[n][3];
  double (*snp)[3] = new double[n][3];
  double (*crk)[3] = new double[n][3];
  
  double *pot  = new double[n];
  double *mass = new double[n];
  int    *nnb  = new int[n];
  double *h2  = new double[n];
  int    *nngb  = new int[n];
  int    *ngb_list = new int[n];
  int    *id   = new int[n];

  double tm = 0;
  for (int i = 0; i < n; i++) {
    pos[i][0] = drand48();
    pos[i][1] = drand48();
    pos[i][2] = drand48();

    vel[i][0] = drand48();
    vel[i][1] = drand48();
    vel[i][2] = drand48();
    
    snp[i][0] = drand48();
    snp[i][1] = drand48();
    snp[i][2] = drand48();
    
    crk[i][0] = drand48();
    crk[i][1] = drand48();
    crk[i][2] = drand48();    
 
    h2[i] = 1.0 * pow(12.0/n, 1.0/3);
    h2[i] = h2[i]*h2[i];

    mass[i] = 1.0/n * drand48();
    tm += mass[i];
    id[i] = i + 1;
  }

  for (int i = 0; i < n; i++) {
    mass[i] *= 1.0/tm;
  }
  
  sapporo grav;
  
  int cluster_id;
//   int sapporo::open(std::string kernelFile, int *devices, int nprocs = 1, int order = FOURTH)
//   grav.open("CUDA/kernels4th.ptx", cluster_id, 1, 1);
  int devices[] = {1,0,2,3};
  grav.open("CUDA/kernels6thDP.ptx",devices , 1, 2);
  
  int ipmax = grav.get_n_pipes();
  double *i_nene = new double[ipmax];

  double null3[3] = {0,0,0};
  for (int i = 0; i < n; i++) {
    grav.set_j_particle(i, id[i],
			0, 0,
			mass[i], null3, null3, null3,
			vel[i], pos[i], null3, null3, 0);
  }
  
 
  grav.set_time(0);
 
  double eps2 = 0;
  n1 = 0;
  n2 = 2*ipmax;
  for (int i = n1; i < n2; i += ipmax) {
    int npart = min(n2 - i, ipmax);
    
//     n = n - 130000;
    grav.startGravCalc(n, npart,
			id + i,
			pos+i, vel+i,
			acc+i, jrk+i, pot+i, eps2, h2, NULL);
                        

    grav.getGravResults(n, npart,
			id+i, pos+i, vel+i,
			eps2, h2,
			acc+i, jrk+i,snp+i, crk+i, pot+i, nnb+i, NULL, true);
                        
//     for(int j=0; j < ipmax; j++){
//     printf("calchalf2: %d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n", j,
//            acc[j][0], acc[j][1], acc[j][2],
//            pot[j], jrk[j][0],
//            jrk[j][1], jrk[j][2], nnb[j]);
//     }

//       grav.read_ngb_list(cluster_id);
//       for (int i1 = i; i1 < i + npart; i1++) {
//         grav.get_ngb_list(cluster_id,
//                           i1 - i,
//                           n,
//                           nngb[i1],
//                           ngb_list);
//         fprintf(stderr," ipipe= %d: n_ngb= %d\n", i1 - i, nngb[i1]);
//       }
                        
  }

//  cerr << "After last half \n";  exit(0);

  for (int kk = 0; kk < 1; kk++) {
  fprintf(stderr, "Computing forces on GPU\n");
  double t1 = get_time();
  int n_evaluate = 0;

  for (int cycle = 0; cycle < 100; cycle++) {
    int an = 1024;
    int n1 = int(an*drand48());
    int n2 = int(an*drand48());
    if (n2 < n1) {
      int tmp = n2;
      n2 = n1;
      n1 = tmp;
    }
//     n1 = 0;
//     n2 = n;
    n_evaluate += n2 - n1;
    fprintf(stderr, "cycle= %d dn= %d .... ", cycle, n2-n1);
    for (int i = n1; i < n2; i++) {
      grav.set_j_particle(i, id[i],
			  0, 0,
			  mass[i], null3, null3, null3,
			  vel[i], pos[i], null3, null3, 0);
    }
//     fprintf(stderr, "n1= %d n2= %d \n", n1, n2);
    grav.set_time(0);
    double tx = get_time();
    for (int i = n1; i < n2; i += ipmax) {
      int npart = min(n2 - i, ipmax);
      
//       int ntest = n-2;
   int ntest = n;   
      grav.startGravCalc(ntest, npart,
			  id + i,
			  pos+i, vel+i,
			  acc+i, jrk+i, pot+i, eps2, h2, NULL);

      grav.getGravResults(ntest, npart,
			 id+i, pos+i, vel+i,
			 eps2, h2,
			 acc+i, jrk+i, snp+i, crk+i, pot+i, NULL, NULL, false);
                         
/*      int j = i;
      for(int j=i; j < npart; j++){      
       printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", j,
           acc[j][0], acc[j][1], acc[j][2],
           pot[j], jrk[j][0],
           jrk[j][1], jrk[j][2]);
      }
       exit(0);  */   
                         
                         
      grav.read_ngb_list(cluster_id);
      for (int i1 = i; i1 < i + npart; i1++) {
	grav.get_ngb_list(cluster_id,
			  i1 - i,
			  n,
			  nngb[i1],
			  ngb_list);
//   	fprintf(stderr," ipipe= %d: n_ngb= %d\n", i1 - i, nngb[i1]);
      }
    }
    
//     exit(0);
    
    double dt = get_time() - tx;
    double dt1 = get_time() - t1;
    double gflop1 = 60.0e-9*(n2 - n1)*n/dt;
    double gflop2 = 60.0e-9*n_evaluate*n/dt1;
    fprintf(stderr,"  GFLOP/s: %g %g \n", gflop1, gflop2);
  }
  double dt = get_time() - t1;
  fprintf(stderr, " Completed in %lf seconds [ni_total= %d  %lg GFLOP/s] \n", dt, n_evaluate, 60.0e-9*n_evaluate*n/dt);
  grav.close();
  }

  double (*accx)[3] = new double[n][3];
  double (*jrkx)[3] = new double[n][3];
  double *potx = new double[n];

  fprintf(stderr, "Computing forces on CPU\n");
  double t1 = get_time();

  for (int i = nx*n1; i < nx*n2; i++) {
    if (i%100 == 0) fprintf(stderr, "i= %d [%d]\n", i, nx);
    double posi[3] = {pos[i][0], pos[i][1], pos[i][2]};
    double veli[3] = {vel[i][0], vel[i][1], vel[i][2]};

    double acct[4] = {0,0,0,0};
    double jrkt[3] = {0,0,0};

    for (int j = 0; j < n; j++) {
      double dr[3] = {posi[0] - pos[j][0],
		      posi[1] - pos[j][1],
		      posi[2] - pos[j][2]};

      double dv[3] = {veli[0] - vel[j][0],
		      veli[1] - vel[j][1],
		      veli[2] - vel[j][2]};
      double ds2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      
      if (ds2 == 0) continue;

      double s = sqrt(ds2 + eps2);
      double c = 1.0/s/s/s;
      
      acct[0] -= mass[j] * c * dr[0];
      acct[1] -= mass[j] * c * dr[1];
      acct[2] -= mass[j] * c * dr[2];
	acct[3] -= mass[j] / s;
      
      // jerk
      double s5 = c/s/s;
      double rv = 3.0 * s5 * (dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2]);

      jrkt[0] -= mass[j] * (c * dv[0] - rv * dr[0]);
      jrkt[1] -= mass[j] * (c * dv[1] - rv * dr[1]);
      jrkt[2] -= mass[j] * (c * dv[2] - rv * dr[2]);
    }
    accx[i][0] = acct[0];
    accx[i][1] = acct[1];
    accx[i][2] = acct[2];
    potx[i]    = acct[3];
    
    jrkx[i][0] = jrkt[0];
    jrkx[i][1] = jrkt[1];
    jrkx[i][2] = jrkt[2];
  }
  cerr << " done in " << get_time() - t1 << " sec \n";

  
  double force[3]  = {0,0,0};
  double  torque[3] = {0,0,0};
  for (int i = n1; i < n2; i++) {
	    
    force[0] += mass[i] * acc[i][0];
    force[1] += mass[i] * acc[i][1];
    force[2] += mass[i] * acc[i][2];
    
    torque[0] += mass[i] * (acc[i][1]*pos[i][2] - acc[i][2]*pos[i][1]);
    torque[1] += mass[i] * (acc[i][2]*pos[i][0] - acc[i][0]*pos[i][2]);
    torque[2] += mass[i] * (acc[i][0]*pos[i][1] - acc[i][1]*pos[i][0]);
    
    if (nx == 0) continue;
    real3 da = (real3){acc[i][0] - accx[i][0],
		       acc[i][1] - accx[i][1],
		       acc[i][2] - accx[i][2]};

    real3 dj = (real3){jrk[i][0] - jrkx[i][0],
		       jrk[i][1] - jrkx[i][1],
		       jrk[i][2] - jrkx[i][2]};
    
    double dpot = pot[i] - potx[i];

    fprintf(stdout, "%lg %lg %lg\n",
	    sqrt((da.x*da.x + da.y*da.y+da.z*da.z)/(accx[i][0]*accx[i][0] + accx[i][1]*accx[i][1]+ accx[i][2]*accx[i][2])),
	    sqrt((dj.x*dj.x + dj.y*dj.y+dj.z*dj.z)/(jrkx[i][0]*jrkx[i][0] + jrkx[i][1]*jrkx[i][1]+ jrkx[i][2]*jrkx[i][2])),
	    fabs(dpot)/potx[i]);

//     fprintf(stdout, "%g %g %g  %g %g %g  %g\n",
// 	    acc[i][0], acc[i][1], acc[i][2],
// 	    jrk[i][0], jrk[i][1], jrk[i][2], pot[i]);
//     double a = pot[i];
//     double b = potx[i];
//     fprintf(stdout, " %g %g %g \n",
// 	    a, b, (a-b)/b);
//     fprintf(stdout, "%g %g %g  %g %g %g  %g\n",
// 	    accx[i][0], accx[i][1], accx[i][2],
// 	    jrkx[i][0], jrkx[i][1], jrkx[i][2], potx[i]);

  }

  fprintf(stderr, "total force  = [%lg %lg %lg] \n",
	  force[0], force[1], force[2]);
  fprintf(stderr, "total torque = [%lg %lg %lg] \n",
	  torque[0], torque[1], torque[2]);

 
  cerr << "done!\n";
}
