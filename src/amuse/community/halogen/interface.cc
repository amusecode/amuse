#include "interface.h"
#define MAXLEN_PATH 250

// Default parameters:
DOUBLE randomseed = 42.0;

// Other globals:
INT outputgridr, outputgriddf;
DOUBLE t0, t1, t2, t3, t4, t5, t6, t7;
PARTICLE *bh;
SI *halo;
GI *gi;
CHAR FILENAME[STRINGSIZE + MAXLEN_PATH + 20], INPUTNAME[STRINGSIZE + MAXLEN_PATH];
FILE *file;

bool particles_generated = false;
bool particles_generated_previously = false;
CHAR output_path[MAXLEN_PATH];      /*!< output directory of the code */
CHAR output_basename[STRINGSIZE];   /*!< basename of output files */



/*
 * Interface code for the stripped-down version of HALOGEN for the AMUSE project
 */

int initialize_code(){
    initialise_fixed_structures(&bh, &halo, &gi);
    set_standard_values_for_parameters(&outputgridr, &outputgriddf, gi, halo);
    halo->p = NULL;
    strcpy(output_path, "./");
    strcpy(output_basename, "halogen");
    strcpy(INPUTNAME, output_path);
    strcat(INPUTNAME, output_basename);
    initialise_black_hole(bh);
    initialise_parameters(halo);
    return 0;
}

int cleanup_code(){
    free(bh);
    free(halo->sp);
    free(halo->griddf);
    free(halo->p);
    free(halo);
    free(gi->stuff);
    free(gi->gridr);
    free(gi);
    return 0;
}

int commit_parameters(){
    fprintf(stderr,"Checking parameters, calculating halo properties and initialising grid in r... \n");
    srand(randomseed);

    /*
    ** Check main input parameters
    */

    if (strcmp(INPUTNAME, "none") == 0) {
        fprintf(stderr,"You have not set a name for the output model.\n");
        usage();
	}
    if ((NGRIDR-1) % (NGRIDDF-1) != 0) {
        fprintf(stderr,"Bad choice of NGRIDR and NGRIDDF!\n");
        fprintf(stderr,"These numbers have to fulfill the condition (NGRIDR-1) mod (NGRIDDF-1) == 0.\n");
        usage();
	}

    check_main_parameters(halo);
    calculate_parameters(gi,halo);    
    initialise_all_grids(gi, bh, halo, outputgridr, outputgriddf, t0, &t1, FILENAME, INPUTNAME);
    initialise_structure(gi,halo);
    t2 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\n",t2-t1);
    return 0;
}

int recommit_parameters(){
    free(halo->p);
    return commit_parameters();
}

int generate_particles(){
    INT i, j, index;
    
    // Set particle positions
    fprintf(stderr,"Setting particle positions... \n");
    set_positions(halo);
    
    // Set particle velocities
    t3 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle velocities... \n",t3-t2);
    set_velocities(gi,halo);
    
    // Set remaining attributes
    t4 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting remaining particle attributes... \n",t4-t3);
    set_attributes(gi,halo);
    
    // Calculate a few things and do center of mass correction
    t5 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nCalculating a few things, doing mass scaling and correct center of mass position and velocity... \n",t5-t4);
    double_particles(halo);
    calculate_stuff(gi,bh,halo);

    /*
    ** Write Output
    */

    t6 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nWriting Output... \n",t6-t5);

    sprintf(FILENAME,"%s.IC.ascii",INPUTNAME);
    file = fopen(FILENAME,"w");

    index = 1;
    if (bh->mass != 0) {
        assert(fprintf(file,OFI1" ",index) > 0);
        assert(fprintf(file,OFD5" ",bh->mass) > 0);
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",bh->r[j+1]) > 0);
        }
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",bh->v[j+1]) > 0);
        }
        assert(fprintf(file,"\n") > 0);
        index++;
    }
    for (i = 0; i < halo->N; i++) {
        assert(fprintf(file,OFI1" ",index) > 0);
        assert(fprintf(file,OFD5" ",halo->p[i].mass) > 0);
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",halo->p[i].r[j+1]) > 0);
            }
        for (j = 0; j < 3; j++) {
            assert(fprintf(file,OFD6" ",halo->p[i].v[j+1]) > 0);
            }
	assert(fprintf(file,"\n") > 0);
        index++;
	}
    fclose(file);
    
    /* 
    ** Print some output in file
    */

    sprintf(FILENAME,"%s.out",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    fprintf(file,"HALOGEN4MUSE by Marcel Zemp / Version 30. August 2007\n\n");
    fprintf(file,"Command line\n\n");
    //for (i = 0; i < argc; i++) fprintf(file,"%s ",argv[i]);
    fprintf(file,"N/A (executed by AMUSE)");
    fprintf(file,"\n\n");
    fprintf(file,"Model properties\n\n");
    fprintf(file,"alpha = "OFD1"\n",halo->sp->alpha);
    fprintf(file,"beta  = "OFD1"\n",halo->sp->beta);
    fprintf(file,"gamma = "OFD1"\n",halo->sp->gamma);
    fprintf(file,"rho0  = "OFD3" MU LU^-3\n",halo->sp->rho0);
    fprintf(file,"\n");
    if (bh->mass > 0) {
	fprintf(file,"MBH = "OFD3" MU\n",bh->mass);
	fprintf(file,"\n");
	}
    fprintf(file,"rs      = "OFD3" LU\n",halo->sp->rs);
    fprintf(file,"rhalf   = "OFD3" LU\n",halo->sp->rhalf);
    fprintf(file,"rvir    = "OFD3" LU\n",halo->sp->rvir);
    fprintf(file,"rinner  = "OFD3" LU\n",gi->rinner);
    fprintf(file,"router  = "OFD3" LU\n",gi->router);
    if (halo->sp->rcutoff != SBI) {
        fprintf(file,"rcutoff = "OFD3" LU\n",halo->sp->rcutoff);
        fprintf(file,"rdecay  = "OFD3" LU\n",halo->sp->rdecay);
        }
    fprintf(file,"\n");
    fprintf(file,"M(rs)      = "OFD3" MU\n",Menc(halo->sp->rs,gi));
    fprintf(file,"M(rhalf)   = "OFD3" MU\n",Menc(halo->sp->rhalf,gi));
    fprintf(file,"M(rvir)    = "OFD3" MU\n",Menc(halo->sp->rvir,gi));
    fprintf(file,"M(rinner)  = "OFD3" MU\n",Menc(gi->rinner,gi));
    fprintf(file,"M(router)  = "OFD3" MU\n",Menc(gi->router,gi));
    if (halo->sp->rcutoff != SBI) {
	fprintf(file,"M(rcutoff) = "OFD3" MU\n",Menc(halo->sp->rcutoff,gi));
	}
    fprintf(file,"\n");
    fprintf(file,"Sampling properties\n\n");
    fprintf(file,"|Cr| = "OFD3" MU       Cr = ("OFD4", "OFD4", "OFD4") MU\n",gi->stuff->Cr[0],gi->stuff->Cr[1],gi->stuff->Cr[2],gi->stuff->Cr[3]);
    fprintf(file,"|Cv| = "OFD3" MU TU^-1 Cv = ("OFD4", "OFD4", "OFD4") MU TU^-1\n",gi->stuff->Cv[0],gi->stuff->Cv[1],gi->stuff->Cv[2],gi->stuff->Cv[3]);
    fprintf(file,"\n");
    fprintf(file,"Etot = "OFD4" MU LU^2 TU^-2\n",gi->stuff->Etot);
    fprintf(file,"Ekin = "OFD4" MU LU^2 TU^-2\n",gi->stuff->Ekin);
    fprintf(file,"Epot = "OFD4" MU LU^2 TU^-2\n",gi->stuff->Epot);
    fprintf(file,"Rvir = |2*Ekin/Epot| = %g\n",fabs(2*gi->stuff->Ekin/gi->stuff->Epot));
    fprintf(file,"\n");
    fprintf(file,"Ntot                = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->N,gi->stuff->N);
    fprintf(file,"rimp                = "OFD3" LU\n",halo->rimp);
    fprintf(file,"r1                  = "OFD3" LU\n",halo->r1);
    fprintf(file,"r100                = "OFD3" LU\n",halo->r100);
    fprintf(file,"Mtheo               = "OFD3" MU\n",bh->mass + halo->sp->M);
    fprintf(file,"Msamp               = "OFD3" MU\n",gi->stuff->Mp);
    fprintf(file,"(Msamp-Mtheo)/Mtheo = "OFD3"\n",gi->stuff->Mp/(bh->mass + halo->sp->M)-1.0);
    fprintf(file,"Random seed         = "OFD3"\n",randomseed);
    fprintf(file,"\n");
    fprintf(file,"Times for individual steps\n\n");
    fprintf(file,"Calculation of halo properties and initialisation of grid in r: "OFD1" seconds.\n",t1-t0);
    fprintf(file,"Initialisation of grid for distribution function: "OFD1" seconds.\n",t2-t1);
    fprintf(file,"Setting particle positions: "OFD1" seconds\n",t3-t2);
    fprintf(file,"Setting particle velocities: "OFD1" seconds\n",t4-t3);
    fprintf(file,"Setting remaining particle attributes: "OFD1" seconds\n",t5-t4);
    fprintf(file,"Calculating a few things and correct center of mass: "OFD1" seconds\n",t6-t5);
    t7 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(file,"Writing output: "OFD1" seconds\n",t7-t6);
    fprintf(file,"Total time: "OFD1" seconds\n",t7-t0);
    fclose(file);
   
    fprintf(stderr,"Done in "OFD1" seconds\nTotal time needed was "OFD1" seconds\n",t7-t6,t7-t0);
    
    particles_generated = true;
    return 0;
}

int clear_generated_particles(){
    if (particles_generated){
        particles_generated_previously = true;
        particles_generated = false;
    }
    return 0;
}

int get_number_of_particles_updated(int *number_of_particles_updated){
    if (!particles_generated)
        *number_of_particles_updated = 0;
    else {
        if (bh->mass > 0)
            *number_of_particles_updated = halo->N0 + 1;
        else
            *number_of_particles_updated = halo->N0;
        particles_generated_previously = true;
        particles_generated = false;
    }
    return 0;
}

int get_mass(int index_of_the_particle, double *mass){
    if (index_of_the_particle >= halo->N0){
        if (index_of_the_particle == halo->N0 && bh->mass > 0){
            *mass = bh->mass;
            return 0;
        } else return -1;
    }
    *mass = halo->p[index_of_the_particle].mass;
    return 0;
}

int get_position(int index_of_the_particle, double *x, double *y, double *z){
    if (index_of_the_particle >= halo->N0){
        if (index_of_the_particle == halo->N0 && bh->mass > 0){
            *x = bh->r[1];
            *y = bh->r[2];
            *z = bh->r[3];
            return 0;
        } else return -1;
    }
    *x = halo->p[index_of_the_particle].r[1];
    *y = halo->p[index_of_the_particle].r[2];
    *z = halo->p[index_of_the_particle].r[3];
    return 0;
}

int get_velocity(int index_of_the_particle, double *vx, double *vy, double *vz){
    if (index_of_the_particle >= halo->N0){
        if (index_of_the_particle == halo->N0 && bh->mass > 0){
            *vx = bh->v[1];
            *vy = bh->v[2];
            *vz = bh->v[3];
            return 0;
        } else return -1;
    }
    *vx = halo->p[index_of_the_particle].v[1];
    *vy = halo->p[index_of_the_particle].v[2];
    *vz = halo->p[index_of_the_particle].v[3];
    return 0;
}


//  Parameters

int get_model_alpha(double *model_alpha){
    *model_alpha = halo->sp->alpha;
    return 0;
}
int set_model_alpha(double model_alpha){
    halo->sp->alpha = model_alpha;
    return 0;
}

int get_model_beta(double *model_beta){
    *model_beta = halo->sp->beta;
    return 0;
}
int set_model_beta(double model_beta){
    halo->sp->beta = model_beta;
    return 0;
}

int get_model_gamma(double *model_gamma){
    *model_gamma = halo->sp->gamma;
    return 0;
}
int set_model_gamma(double model_gamma){
    halo->sp->gamma = model_gamma;
    return 0;
}

int get_total_mass(double *total_mass){
    *total_mass = halo->sp->M;
    return 0;
}
int set_total_mass(double total_mass){
    halo->sp->M = total_mass;
    return 0;
}

int get_scale_radius(double *scale_radius){
    *scale_radius = halo->sp->rs;
    return 0;
}
int set_scale_radius(double scale_radius){
    halo->sp->rs = scale_radius;
    return 0;
}

int get_cutoff_radius(double *cutoff_radius){
    *cutoff_radius = halo->sp->rcutoff;
    return 0;
}
int set_cutoff_radius(double cutoff_radius){
    halo->sp->rcutoff = cutoff_radius;
    return 0;
}

int get_target_number_of_particles(int *target_number_of_particles){
    *target_number_of_particles = halo->N0;
    return 0;
}
int set_target_number_of_particles(int target_number_of_particles){
    halo->N0 = target_number_of_particles;
    return 0;
}

int get_black_hole_mass(double *black_hole_mass){
    *black_hole_mass = bh->mass;
    return 0;
}
int set_black_hole_mass(double black_hole_mass){
    bh->mass = black_hole_mass;
    return 0;
}

int get_random_seed(double *random_seed_out){
    *random_seed_out = randomseed;
    return 0;
}
int set_random_seed(double random_seed_in){
    randomseed = random_seed_in;
    return 0;
}

int get_do_exact_virial_radius_flag(int *do_exact_virial_radius_flag){
    *do_exact_virial_radius_flag = gi->dorvirexact;
    return 0;
}
int set_do_exact_virial_radius_flag(int do_exact_virial_radius_flag){
    if (do_exact_virial_radius_flag) {
        gi->dorvirexact = 1;
    } else {
        gi->dorvirexact = 0;
    }
    return 0;
}

int get_outputgridr_flag(int *outputgridr_flag){
    *outputgridr_flag = outputgridr;
    return 0;
}
int set_outputgridr_flag(int outputgridr_flag){
    if (outputgridr_flag) {
        outputgridr = 1;
    } else {
        outputgridr = 0;
    }
    return 0;
}
int get_outputgriddf_flag(int *outputgriddf_flag){
    *outputgriddf_flag = outputgriddf;
    return 0;
}
int set_outputgriddf_flag(int outputgriddf_flag){
    if (outputgriddf_flag) {
        outputgriddf = 1;
    } else {
        outputgriddf = 0;
    }
    return 0;
}

int get_output_basename(char **output_basename_out){
    *output_basename_out = output_basename;
    return 0;
}
int set_output_basename(char *output_basename_in){
    if (strlen(output_basename_in) >= STRINGSIZE)
        return -1;
    strcpy(output_basename, output_basename_in);
    strcpy(INPUTNAME, output_path);
    strcat(INPUTNAME, output_basename);
    return 0;
}
int get_output_path(char **output_path_out){
    *output_path_out = output_path;
    return 0;
}
int set_output_path(char *output_path_in){
    int length = strlen(output_path_in);
    if (length >= MAXLEN_PATH)
        return -1;
    strcpy(output_path, output_path_in);
    if (length > 0 && output_path[length - 1] != '/')
        strcat(output_path, "/");
    strcpy(INPUTNAME, output_path);
    strcat(INPUTNAME, output_basename);
    return 0;
}
