/* -------------------------------- plumix ------------------------------- */

/*  program for generating of a mass segregated star cluster specified     */
/*  by its mass and/or number of stars and a power-law or single-mass      */
/*  mass function. Output is written to files 'n6.pos' and 'n6.check'.     */
/*  The first one is a list of stellar masses, positions and velocities    */
/*  in the format suitable for Sverre Aarseth's NBODYx code; the latter    */
/*  one is a human-readable dump of most of the quantities used in the     */
/*  code. All quantities are in 'NBODY' units; total energy of the cluster */
/*  is implicitely set to 0.5.                                             */
/*                                                                         */
/*  For details see 'A new method to create initially mass segregated star */
/*  clusters in virial equilibrium', by Subr, Kroupa & Baumgardt, 2008     */
/* ----------------------------------------------------------------------- */

#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <rand.h>

// some physical constants in CGS:
#define _G      6.673e-8
#define _pc     3.0856e18
#define _M_sun  1.989e33
#define _R_sun  6.96265e10
#define _AU     1.49597870e13
#define _GMsun  1.3272597e26

#define RSCALE  0.58904862254809  // = 3pi / 16 = 3pi / (16 * |U_tot|)

#define BUFF_STEP 1024

/* Quantities related to individual particles are split into two structures
   for the sake of optimal utilization on the processor cache memory.
   The first one holds variables intensively used in the inner loop
   for calculating potentials. The second one holds the rest.  */
struct t_star1
{
	double mass, U_tmp, Ui;
	double r[3];
};

struct t_star2
{
	double v[3];
	double U, U_sub, E;
	double UT_add, UT;
	double rad, M_sub;
	long   ntry;
};

struct t_star1 *star1;
struct t_star2 *star2;

void position(int id, double U_sub, double UT_sub, double S)
{
	long   i;
	double rx, ry, rz, r1;
	double rfac;
	
	rfac = pow(star2[id].M_sub, 2. * S) / (1. - S);
	
	star2[id].ntry = 0;
	
	do
	 {
		star2[id].U_sub = 0.;
		star2[id].ntry++;

		r1 = rand2() * INV_RAND_MAX2;
		r1 = exp(-2. * log(r1) / 3.);
		r1 = RSCALE * sqrt(1. / (r1 - 1.));

		r1 *= rfac;
		
		isorand(star1[id].r);
		star1[id].r[0] *= r1;
		star1[id].r[1] *= r1;
		star1[id].r[2] *= r1;
		star2[id].rad = r1;

		for (i = 0; i < id; i++)
		 {
			rx = star1[i].r[0] - star1[id].r[0];
			ry = star1[i].r[1] - star1[id].r[1];
			rz = star1[i].r[2] - star1[id].r[2];
			r1 = sqrt(rx*rx + ry*ry + rz*rz);
			star1[i].U_tmp = star1[id].mass / r1;
			star2[id].U_sub += star1[i].mass / r1;
		 };
		r1 = U_sub + star1[id].mass * star2[id].U_sub;
	 }
	while (fabs(r1 - UT_sub) > 0.1 * (r1 + UT_sub) / sqrt(1. + id));
//	while (fabs(r1 - policajt) > 0.5 * (r1 + policajt) / sqrt(1. + id));

	for (i = 0; i < id; i++) star2[i].U += star1[i].U_tmp;
};

void find_beta(double *avg, double *beta)
{
	int i;
	double a1, a2, I1, I2, J1, J2;
	double alo, ah, Ilo, Ih, Jlo, Jh;
	double tmpavg;

	if (*avg > 0.75)
	 {
		*beta = -0.5;
		*avg = 1.;
		return;
	 };
	
	Ilo = I1 = 3. * M_PI / 16.;
	Ih  = I2 = 0.2;
	Jlo = J1 = M_PI / 4;
	Jh  = J2 = 1. / 3.;
	alo = a1 = -0.5;
	ah  = a2 = 0.;
	
	tmpavg = I2 / J2;
	i = 0;
	while (tmpavg > *avg)
	 {
		if (!(i & 1))
		 {
			a1 += 1.;
			I1 /= 1. + 5. / (2. * a1);
			J1 /= 1. + 3. / (2. * a1);
			Ilo = I2; Ih = I1;
			Jlo = J2; Jh = J1;
			alo = a2; ah = a1;
		 }
		else
		 {
			a2 += 1.;
			I2 /= 1. + 5. / (2. * a2);
			J2 /= 1. + 3. / (2. * a2);
			Ilo = I1; Ih = I2;
			Jlo = J1; Jh = J2;
			alo = a1; ah = a2;
		 };
		tmpavg = Ih / Jh;
		i++;
	 };
	 
	*beta = alo + (ah - alo) * (*avg - Ilo / Jlo) / (tmpavg - Ilo / Jlo);
	
	if (*beta > 0.) *avg = exp(*beta * log(*beta / (1. + *beta))) / (1. + *beta);
	else *avg = 2. * exp((2.* *beta + 1.) * log((2.* *beta + 1.) / (2.* *beta + 3.))) / (2.* *beta + 3.);
};

#define SWAP(a,b) { temp = star1[a].mass; star1[a].mass = star1[b].mass; star1[b].mass = temp; };

void quick(int start, int stop)
{
	int i, j;
	double temp, median;

	median = 0.5 * (star1[start].mass + star1[stop].mass);

	i = start;
	j = stop;
	while (i <= j)
	 {
		while ((i <= stop) && (star1[i].mass >= median)) i++;
		while ((j >= start) && (star1[j].mass <= median)) j--;
		if (j > i)
		 {
			SWAP(i,j);
			i++;
			j--;
		 }
		else if (i == j)
		 {
			i++;
			j--;
		 };
	 };
	
	if (start + 1 < j) quick(start, j);
	else if ((start < j) && (star1[start].mass < star1[j].mass)) SWAP(start, j);

	if (i + 1 < stop) quick(i, stop);
	else if ((i < stop) && (star1[i].mass < star1[stop].mass)) SWAP(i, stop)
};

void usage()
{
	printf("\nUSAGE:\n plumix [-S s] [-hv]\n\n"
		"OPTIONS:\n"
		" -S  set the value of the index of mass segregation; default value: 0\n"
		" -v  print some kind of 'progress bar' during the program work\n"
		" -h  print this help and exit\n\n");
};

/* int main(int argv, char **argc) */
int plumix_(double *body, int *n, double *xtemp, double *vtemp, double *w0)
{
	int    option;
	int    verbose;
	int    buff_size;
	long   i, j, icheck;
	double r1, r2, tmp;
	double U_tot, U_tot1, K_tot;
	double UT_sub, M_sub;
	double maxim, beta;
	FILE   *pos, *check;
	char   buffer[1024];
/* star cluster parameters */
	int    N_star;
	double M_tot;
	double S;
/* CMS */
	double cmx, cmy, cmz;
	double cmdx, cmdy, cmdz;
	
	verbose = 1;
/* 	S = 0.; */
	S = *w0;
/* 	while ((option = getopt(argv, argc, "hS:v")) != -1) switch (option) */
/* 	 { */
/* 		case 'S': S = atof(optarg); break; */
/* 		case 'v': verbose = 1; break; */
/* 		case 'h': usage(); return 0; */
/* 		case ':': */
/* 		case '?': usage(); return 1; */
/* 	 }; */
	 
	printf("in plumix\n");
	printf("%f %i %f\n",body[0],*n,S);
	time(&i); // inicializace generatoru nahodnych cisel
	init_rand(i);
	rand2();
	
	star1 = NULL;
	star2 = NULL;
	N_star = buff_size = 0;
	M_tot = 0.;

/* 	while (fgets(buffer, 1024, stdin)) */
	while (N_star < *n)
	 {
/* 		if (*buffer == '#') continue; */
		N_star++;
		if (N_star > buff_size)
		 {
			buff_size += BUFF_STEP;
			star1 = (struct t_star1 *)realloc(star1, buff_size * sizeof(struct t_star1));
			star2 = (struct t_star2 *)realloc(star2, buff_size * sizeof(struct t_star2));
		 };
/*  		star1[N_star-1].mass = atof(buffer);  */
  		star1[N_star-1].mass = body[N_star-1];
 		M_tot += star1[N_star-1].mass; 
	 };
	printf("%i %f %f\n",N_star,star1[0].mass,star1[N_star-1].mass);
/* 	return 0; */
	star1 = (struct t_star1 *)realloc(star1, N_star * sizeof(struct t_star1));
	star2 = (struct t_star2 *)realloc(star2, N_star * sizeof(struct t_star2));

	quick(0, N_star-1);

	if (verbose) printf("initialising..."); fflush(stdout);
	M_sub = U_tot1 = 0.;
	for (i = 0; i < N_star; i++)
	 {
//		printf("%8.4f\n", star1[i].mass);
		star2[i].U = star1[i].U_tmp = star2[i].UT_add = 0.;
		star1[i].mass /= M_tot;
		M_sub += star1[i].mass;
		star2[i].M_sub = M_sub;
		if (S > 0.) star1[i].Ui = star1[i].mass * pow(M_sub, -S);
		else star1[i].Ui = star1[i].mass;
		for (j = 0; j < i; j++)
		 {
			tmp = star1[i].Ui * star1[j].Ui;
			star2[i].UT_add += tmp;
			star1[i].U_tmp += tmp;
			star1[j].U_tmp += tmp;
		 };
		U_tot1 += star2[i].UT_add;
	 };
	
	UT_sub = U_tot = K_tot = 0.;
	cmx = cmy = cmz = 0.;
	cmdx = cmdy = cmdz = 0.;
	 
	icheck = N_star / 10;

	if (verbose) printf(" done\nworking"); fflush(stdout);
	
/* vyrobime pozice a napocitame energie */
	for (i = 0; i < N_star; i++)
	 {
		star2[i].UT_add *= 0.5 / U_tot1;
		star2[i].UT = star1[i].U_tmp * 0.5 / U_tot1;
		UT_sub += star2[i].UT_add;
		position(i, U_tot, UT_sub, S);
		if (verbose && !((i - 1) % icheck)) { printf("%i\n",i-1); fflush(stdout); }
		U_tot += star1[i].mass * star2[i].U_sub;
	 };
	 
	if (verbose) printf(" done\n");

	for (i = 0; i < N_star; i++)
	 {
		maxim = 0.25 * star1[i].mass / star2[i].UT;
		find_beta(&maxim, &beta);
		if (beta > 0.)
		 {
			do
			 {
				r1 = maxim * rand2() * INV_RAND_MAX2;
				r2 = rand2() * INV_RAND_MAX2;
				r2 *= r2;
				tmp = 1. - r2;
			 }
			while (r1 > r2 * exp(beta * log(tmp)));
		 }
		else
		 {
			do
			 {
				r1 = maxim * rand2() * INV_RAND_MAX2;
				r2 = sin(M_PI_2 * rand2() * INV_RAND_MAX2);
				r2 *= r2;
				tmp = sqrt(1. - r2);
			 }
			while (r1 > r2 * exp((2.* beta + 1.) * log(tmp)));
		 };
	
		star2[i].E = r2 * (star2[i].U + star2[i].U_sub);
		
		tmp = sqrt(2. * star2[i].E);
		isorand(star2[i].v);
		star2[i].v[0] *= tmp;
		star2[i].v[1] *= tmp;
		star2[i].v[2] *= tmp;
		K_tot += star1[i].mass * star2[i].E;
		
/* centre of mass and its momentum */
		cmx += star1[i].mass * star1[i].r[0];
		cmy += star1[i].mass * star1[i].r[1];
		cmz += star1[i].mass * star1[i].r[2];
		cmdx += star1[i].mass * star2[i].v[0];
		cmdy += star1[i].mass * star2[i].v[1];
		cmdz += star1[i].mass * star2[i].v[2];
	 };
	
	icheck = 0;

/* CMS correction */
	for (i = 0; i < N_star; i++)
	 {
		star1[i].r[0] -= cmx;
		star1[i].r[1] -= cmy;
		star1[i].r[2] -= cmz;
		star2[i].v[0] -= cmdx;
		star2[i].v[1] -= cmdy;
		star2[i].v[2] -= cmdz;
		star2[i].rad  = sqrt(star1[i].r[0]*star1[i].r[0] + star1[i].r[1]*star1[i].r[1] + star1[i].r[2]*star1[i].r[2]);
		icheck += star2[i].ntry;
	 };
	 
/* write masses, positions and velocities to 'n6.pos' in a format convenient for NBODYx code;
   some other diagnostics stuff is stored in 'n6.check' */
	pos = fopen("n6.pos", "w");
	check = fopen("n6.check", "w");
	
/* 	fprintf(check, "# %s", argc[0]); */
/* 	for (i = 1; i < argv; i++) fprintf(check, " %s", argc[i]); */
	fprintf(check, "\n#\n# N    = %d\n", N_star);
	fprintf(check, "# <m>  = %f\n", M_tot / N_star);
	fprintf(check, "# M_c  = %.2f\n", M_tot);
	fprintf(check, "# U_t  = %f\n", U_tot);
	fprintf(check, "# K_t  = %f\n", K_tot);
	fprintf(check, "# V_t  = %f\n", K_tot / U_tot);
	fprintf(check, "# NT   = %ld\n", icheck);
	fprintf(check, "# CM   = (%7.4f , %7.4f , %7.4f)\n", cmx, cmy, cmz);
	fprintf(check, "# CMD  = (%7.4f , %7.4f , %7.4f)\n", cmdx, cmdy, cmdz);
	
	U_tot = 0.;
	for (i = 0; i < N_star; i++)
	 {
		fprintf(pos, "%.12f  %18.14f  %18.14f  %18.14f  %16.12f  %16.12f  %16.12f\n",
			star1[i].mass, star1[i].r[0], star1[i].r[1], star1[i].r[2], star2[i].v[0], star2[i].v[1], star2[i].v[2]);
		
		U_tot += star1[i].mass * star2[i].U_sub;
		fprintf(check, "%7.4f  %e  %e  %e  %e  %e  %e  %5ld\n",
			star1[i].mass * M_tot, star2[i].rad, star2[i].M_sub, star2[i].E, star2[i].U, star2[i].U_sub, U_tot, star2[i].ntry);
		xtemp[3*(i)] = star1[i].r[0];
		xtemp[3*(i)+1] = star1[i].r[1];
		xtemp[3*(i)+2] = star1[i].r[2];

		vtemp[3*(i)] = star2[i].v[0];
		vtemp[3*(i)+1] = star2[i].v[1];
		vtemp[3*(i)+2] = star2[i].v[2];
	 };
	
	fclose(pos);
	fclose(check);
	
	free(star1);
	free(star2);
	
	return 0;
};
