#include "main.h"

// random skip factor, must be larger than the number of random draws per particle
#define SKIP  100000

main(argc,argv)
int argc;
char **argv;
{
        int i, j, k, nobj=10000;
        long long lseed=0;
        int seed= -123;
        float q=0.9, rtrunc=5.0;
        float rhoran, massapp;
        float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
        float phi, cph, sph, cth, sth, vR, vp, vz;
        float E, Lz, rad, rad2;
        float f0, frand, fmax, psi;
        float Fhalo(), halodens_(), pot_();
        float ran1();
        void ran_seed();
        float dr, rhomax1;
        float t, mass;
        float u1, v1, u1max, v1max;
        float xcm, ycm, zcm, vxcm, vycm, vzcm;
        float zip = 0.0, psi0;
        float rhomax=0.15, rhomin, rhotst;
        float stream=0.5, massfrac=1.0;
        int icofm=1;
        char harmfile[60];
        FILE *outfile;
        strcpy(harmfile,"dbh.dat");
        fquery("Enter streaming fraction",&stream);
        iquery("Enter the number of particles",&nobj);
      	llquery("Enter negative integer seed",&lseed);
        iquery("Center particles (0=no, 1=yes)",&icofm);
        cquery("Enter harmonics file",harmfile);

    readmassrad(); /* reads in mass and truncation radius of components */
        readharmfile_(harmfile,&gparam);
        A = gparam.A; B = gparam.B; C = gparam.C; 
        v0 = gparam.v0; q = gparam.q; psi0 = gparam.psi0;
        v02 = v0*v0;haloflag=gparam.ihaloflag;
        
        fprintf(stderr,"flags: %d %d %d\n",
                gparam.idiskflag,gparam.ibulgeflag,gparam.ihaloflag);

        rtrunc = haloedge;  mass = halomass/nobj;
        fprintf(stderr,"haloedge: %g\n",haloedge);
        u1max = rtrunc;
        v1max = 0.5*M_PI;

        r = (phase *) calloc(nobj,sizeof(phase));

/* Find maximum of rho*r^2 */
        dr = rtrunc/100;
    rhomax1 = 0.0;
/* check the R-axis first */
        for(i=0; i<100; i++) {
                float z, rcyl, rhocur;
        rcyl = i*dr;
        z = 0;
                rhocur = halodens_(&rcyl,&z);
        rhocur *= (rcyl*rcyl);
                if( rhocur > rhomax1  )
                        rhomax1 = rhocur;
        }
    fprintf(stderr,"rhomax on R-axis is %g\n",rhomax1);
/* then check the z-axis */
        for(i=0; i<100; i++) {
                float z, rcyl, rhocur;
        rcyl = 0;
        z = i*dr;
                rhocur = halodens_(&rcyl,&z);
        rhocur *= (z*z);
                if( rhocur > rhomax1  )
                        rhomax1 = rhocur;
        }
        fprintf(stderr,"rhomax1 = %g\n",rhomax1);
        rhomax = 1.5*rhomax1;
        rhomin = 1.0e-10*rhomax;  /* define density contour cutoff */
        fprintf(stderr,"rhomin, rhomax: %g %g\n",rhomin,rhomax); 

        fprintf(stderr,"Calculating halo positions and velocities\n");
        if(lseed<0) lseed=-lseed;
        ran_seed(SKIP*lseed);
        for(i=0; i<nobj;) {

                u1 = u1max*ran1(&seed);
                v1 = 2.0*v1max*(ran1(&seed) - 0.5);
                R = u1;
                z = R*tan(v1);
                
                if(R*R+z*z > 4*haloedge*haloedge)
                        continue;

                rhotst = halodens_(&R,&z);
                rhotst = rhotst*(R*R + z*z);

                j++;
                if( rhotst < rhomin )
                        continue;
                if(z> 100000.) printf("v1, rhotst: %g %g\n",z,rhotst);


                rhoran = (rhomax - rhomin)*ran1(&seed);
                if( rhoran > rhotst )
                        continue;
                phi = 2.0*M_PI*ran1(&seed);
                x = R*cos(phi);
                y = R*sin(phi);

                psi = pot_(&R,&z);
                psi = -psi;

                vmax2 = 2*psi;
                vmax = sqrt(vmax2);
                fmax = 1.1*Fhalo(psi,0.0);
                f0 = 0.0; frand = 1.0; /* dummy starters */
                j = 0;
                while( frand > f0 ) {

                                /* select a random speed <3x the circular speed */

                        v2 = 1.1*2.0*psi;
                        while( v2 > vmax2 ) {
                                vR = 2*vmax*(ran1(&seed) - 0.5);
                                vp = 2*vmax*(ran1(&seed) - 0.5);
                                vz = 2*vmax*(ran1(&seed) - 0.5);
                                v2 = vR*vR + vp*vp + vz*vz;
                        }
                        E = psi - 0.5*v2;
                        Lz = R*vp;
                        f0 = Fhalo(E,Lz);
                        frand = fmax*ran1(&seed);
                        j++;
                }

                if( ran1(&seed) < stream )
                        vp = fabs(vp);
                else
                        vp = -fabs(vp);

                cph = x/R; sph = y/R;
                vx = vR*cph - vp*sph;
                vy = vR*sph + vp*cph;
        /*  vz stays the same */

                r[i].x = (float) x;
                r[i].y = (float) y;
                r[i].z = (float) z;
                r[i].vx = (float)vx;
                r[i].vy = (float)vy;
                r[i].vz = (float)vz;
                i++;ran_seed(SKIP*(lseed+i));
                if( i % 1000 == 0 ) fprintf(stderr,".");
        }
        fprintf(stderr,"\n");

        if( icofm ) {
                xcm = ycm =zcm = vxcm =vycm =vzcm = 0;
                for(i=0; i<nobj; i++) {
                        xcm += r[i].x;
                        ycm += r[i].y;
                        zcm += r[i].z;
                        vxcm += r[i].vx;
                        vycm += r[i].vy;
                        vzcm += r[i].vz;
                }
                xcm /= nobj; ycm /=nobj; zcm /= nobj;
                vxcm /= nobj; vycm /=nobj; vzcm /= nobj;

                for(i=0; i<nobj; i++) {
                        r[i].x -= xcm;
                        r[i].y -= ycm;
                        r[i].z -= zcm;
                        r[i].vx -= vxcm;
                        r[i].vy -= vycm;
                        r[i].vz -= vzcm;
                }
        }
    t = 0.0;
#ifdef ASCII
	fprintf(stdout,"%d %f\n",nobj,t);
	for(i=0; i<nobj; i++) {
	  fprintf(stdout,"% 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e\n", mass, r[i].x, r[i].y, r[i].z, r[i].vx, r[i].vy, r[i].vz);
	}
#else
	for(i=0; i<nobj; i++) {
	  fwrite(&mass,sizeof(float),1,stdout);
      fwrite(r+i,sizeof(phase),1,stdout);
    }
#endif        
        exit(0);
}
/* This defines the distribution function */




float Fhalo1(E, Lz)
float E, Lz;
{
        float f0;

        if( E < 0 )
                return 0.0;
        f0 = ((A*Lz*Lz + B)*exp(2.0*E/v02) + C)*(exp(2.0*E/v02) - 1.0);

        return f0;
}



float Fhalo(E,lz)
float E, lz;
{
  float fhalo3_();
  if(haloflag == 1) return Fhalo1(E,lz);
  if(haloflag == 3) {E=-E;return fhalo3_(&E);}
  printf("error\n");
  exit(1);  
  }
