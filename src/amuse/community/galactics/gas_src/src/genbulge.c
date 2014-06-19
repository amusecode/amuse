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
        int icofm=1;
        float q=0.9, rtrunc=5.0;
        float rhoran, massapp;
        float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
        float phi, cph, sph, cth, sth, vR, vp, vz;
        float E, Lz, rad, rad2;
        float f0, frand, fmax, psi;
        float Fbulge();
        float ran1();
        void ran_seed();
        float t, mass;
        float dr;
        float u1, v1, u1max, v1max;
        float xcm, ycm, zcm, vxcm, vycm, vzcm;
        float zip = 0.0, psi0;
        float rhomax=0.15, rhomin, rhotst;
    float rhomax1;
        float stream=0.5;
        float bulgedens_(), pot_();
        char harmfile[40];
        FILE *outfile;


        fprintf(stderr,"The Bulge\n");
        (void) strcpy(harmfile,"dbh.dat");
    fquery("Enter streaming fraction",&stream);
        iquery("Enter the number of particles",&nobj);
      	llquery("Enter negative integer seed",&lseed);
    iquery("Center particles (0=no, 1=yes)",&icofm);
        cquery("Enter harmonics file",harmfile);

        readmassrad();
        readharmfile_(harmfile,&gparam);
        psicutb = gparam.psicut; sigb2 = gparam.sigbulge; rho1 = gparam.rho1;
//        sigb2 = sigb*sigb;
        psi0=gparam.psi0;
        fbulgeconst = rho1/pow(2.0*M_PI*sigb2,1.5);
        bulgeflag=gparam.ibulgeflag;
        mass = bulgemass/nobj; 
        rtrunc = bulgeedge;
        u1max = rtrunc;
        v1max = 0.5*M_PI;

        R=0;z=0;
//        psi0=pot_(&R,&z);
        fprintf(stderr,"psi0 %f %f %f\n",psi0,sigb2,psicutb);
        fflush(stderr);

        r = (phase *) calloc(nobj,sizeof(phase));

/* Find maximum of rho*r^2 */
        dr = rtrunc/100;
    rhomax1 = 0.0;
        for(i=0; i<100; i++) {
                float z, rcyl, rhocur;
        rcyl = i*dr;
        z = 0;
                rhocur = bulgedens_(&rcyl,&z);
        rhocur *= (rcyl*rcyl);
                if( rhocur > rhomax1  )
                        rhomax1 = rhocur;
        }
        fprintf(stderr,"rhomax1 = %g\n",rhomax1);
        rhomax = 1.5*rhomax1;
        rhomin = 1.0e-10*rhomax;  /* define density contour cutoff */

        fprintf(stderr,"Calculating bulge particles positions and velocities\n");
        if(lseed<0) lseed=-lseed;
        ran_seed(SKIP*lseed);
        for(i=0; i<nobj;) {

                u1 = u1max*ran1(&seed);
                v1 = 2.0*v1max*(ran1(&seed) - 0.5);
                R = u1;
                z = R*tan(v1);

                rhotst = bulgedens_(&R,&z);
                rhotst = rhotst*(R*R + z*z);

                j++;
                if( rhotst < rhomin )
                        continue;

                rhoran = (rhomax - rhomin)*ran1(&seed);
                if( rhoran > rhotst )
                        continue;
                phi = 2.0*M_PI*ran1(&seed);
                x = R*cos(phi);
                y = R*sin(phi);

                psi = pot_(&R,&z);
//                vmax2 = -2.0*(psi-psi0-psicutb);
                vmax2=-2*psi;
                vmax = sqrt(vmax2);
                fmax = Fbulge(psi,psi0);
                f0 = 0.0; frand = 1.0; /* dummy starters */
                j = 0;
                while( frand > f0 ) {
                /* select a random speed < the escape speed */

                        v2 = 1.1*vmax2;
                        while( v2 > vmax2 ) {
                                vx = 2*vmax*(ran1(&seed) - 0.5);
                                vy = 2*vmax*(ran1(&seed) - 0.5);
                                vz = 2*vmax*(ran1(&seed) - 0.5);
                                v2 = vx*vx + vy*vy + vz*vz;
                        }
                        E = 0.5*v2 + psi;
                        f0 = Fbulge(E,psi0);
                        frand = fmax*ran1(&seed);
                        j++;
                }

/* Streaming of the bulge */
                vR = (vx*x + vy*y)/R;
                vp = (-vx*y + vy*x)/R;
        if( ran1(&seed) < stream )
            vp = fabs(vp);
        else
            vp = -fabs(vp);
                vx = (vR*x - vp*y)/R;
                vy = (vR*y + vp*x)/R;

                r[i].x = (float) x;
                r[i].y = (float) y;
                r[i].z = (float) z;
                r[i].vx = (float)vx;
                r[i].vy = (float)vy;
                r[i].vz = (float)vz;
                i++;ran_seed(SKIP*(lseed + i));
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
    return 0;
    
}

/* This defines the distribution function */

float Fbulge1(E,psi0)
float E,psi0;
{
        float f0;

        if( E-psi0 > psicutb )
                return 0.0;

        f0 = exp(-(E-psi0-psicutb)/sigb2)-1;

        return f0;
}

float Fbulge(E,psi0)
float E, psi0;
{
  float fbulge3_();  
  if(bulgeflag == 1) return Fbulge1(E,psi0);
  if(bulgeflag == 3) return fbulge3_(&E,&psi0);
  printf("error2\n");
  exit(1);  
  }

