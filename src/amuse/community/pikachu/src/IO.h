#ifndef IO_H
#define IO_H

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include"mpi.h"
#include"mpi_interface.h"
#include"particle.h"
#include"sort.h"
#include"global.h"

inline void read0_scatter(Particle prt[],
			  int &NFS_global,
			  int &NFS_local,
			  double &Tsys, 
			  char sinput[],
			  const int &NBH_global,
			  const int &NBH_local){
    ifstream finput;
    if(MYRANK == 0){
	finput.open(sinput);
	if(!finput){ 
	    cerr<<"input file can't open "<<endl; 
	    halt_program();
	}
	cerr<<"Input File: "<<sinput<<endl;
	finput>>NFS_global;
	int dim;
	finput>>dim;
	cerr<<"dim="<<dim<<endl;
	finput>>Tsys;
	cerr<<"Tsys="<<Tsys<<endl;
    }
    MPI_Bcast(&NFS_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Tsys, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int Nproc = NPROC;
    int NFS_global_over_Nproc = NFS_global/Nproc;

    int *NFS_local_array = new int[Nproc];
    int *NFS_local_disp = new int[Nproc+1];
    NFS_local_disp[0] = 0;
    for(int i=0; i<Nproc; i++){
	NFS_local_array[i] = NFS_global_over_Nproc;
	if(i < NFS_global % Nproc){
	    NFS_local_array[i]++;
	}
	NFS_local_disp[i+1] = NFS_local_disp[i] + NFS_local_array[i];
    }

    double *mass_array = new double[NFS_local_array[0]+1];
    Vector3 *pos_array = new Vector3[NFS_local_array[0]+1];
    Vector3 *vel_array = pos_array;

    int tag = 0;
    MPI_Status status;
    for(int i=0; i<Nproc; i++){
	for(int j=NFS_local_disp[i]; j<NFS_local_disp[i+1]; j++){
	    if(MYRANK==0 ){
		finput>>mass_array[j-NFS_local_disp[i]];
	    }
	}
	if(i != 0){
	    if(MYRANK==0){
		MPI_Send(mass_array, NFS_local_array[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
	    }
	    else if(MYRANK==i){
		MPI_Recv(mass_array, NFS_local_array[i], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	    }
	}
	if(MYRANK==i){
	    for(int k=0; k<NFS_local_array[i]; k++){
		prt[k+NBH_local].mass = mass_array[k];
		prt[k+NBH_local].index = NFS_local_disp[i] + k + NBH_global;
		prt[k+NBH_local].type = star;
	    }
	}
    }

    int vec_size = sizeof(Vector3);
    for(int i=0; i<Nproc; i++){
	for(int j=NFS_local_disp[i]; j<NFS_local_disp[i+1]; j++){
	    if(MYRANK==0){
		finput>>pos_array[j - NFS_local_disp[i]];
	    }
	}
	if(i != 0){
	    if(MYRANK==0){
		MPI_Send(pos_array, NFS_local_array[i]*vec_size, MPI_BYTE, i, tag, MPI_COMM_WORLD);
	    }
	    else if(MYRANK==i){
		MPI_Recv(pos_array, NFS_local_array[i]*vec_size, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status);
	    }
	}
	if(MYRANK==i){
	    for(int k=0; k<NFS_local_array[i]; k++){
		prt[k+NBH_local].pos = pos_array[k];
	    }
	}
    }

    for(int i=0; i<Nproc; i++){
	for(int j=NFS_local_disp[i]; j<NFS_local_disp[i+1]; j++){
	    if(MYRANK==0){
		finput>>vel_array[j - NFS_local_disp[i]];
	    }
	}
	if(i != 0){
	    if(MYRANK==0){
		MPI_Send(vel_array, NFS_local_array[i]*vec_size, MPI_BYTE, i, tag, MPI_COMM_WORLD);
	    }
	    else if(MYRANK==i){
		MPI_Recv(vel_array, NFS_local_array[i]*vec_size, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status);
	    }
	}
	if(MYRANK==i){
	    for(int k=0; k<NFS_local_array[i]; k++){
		prt[k+NBH_local].vel = vel_array[k];
	    }
	}
    }
    NFS_local = NFS_local_array[MYRANK];
    //for(int i=0; i<NFS_local+NBH_local; i++){ prt[i].time = Tsys;}

    delete [] mass_array;
    delete [] pos_array;
    delete [] NFS_local_array;
    delete [] NFS_local_disp;

}


inline void read1_scatter(Particle prt[],
                          Particle BH_glb[],
                          Particle prt_dead_glb[],
			  double &Tsys, 
			  char sinput[]){

    ifstream finput;
    int Ntot;
    double Egr;
    if(MYRANK == 0){
        finput.open(sinput);
        if(!finput){ 
            cerr<<"input file can't open "<<endl; 
            halt_program();
        }
        cerr<<"Input File: "<<sinput<<endl;
        finput.read((char*)&Ntot, sizeof(int));
        finput.read((char*)&NBH_GLB_ORG, sizeof(int));
        NFS_GLB_ORG = Ntot - NBH_GLB_ORG;
        cerr<<"NFS_GLB_ORG="<<NFS_GLB_ORG<<endl;
        cerr<<"NBH_GLB_ORG="<<NBH_GLB_ORG<<endl;
        finput.read((char*)&Tsys, sizeof(double));
        cerr<<"Tsys="<<Tsys<<endl;
        finput.read((char*)&Egr, sizeof(double));
        cerr<<"Egr="<<Egr<<endl;
    }
    MPI_Bcast(&NFS_GLB_ORG, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NBH_GLB_ORG, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Tsys, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Egr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int Nprt_global_org = NFS_GLB_ORG + NBH_GLB_ORG;
    int Nproc = NPROC;
    int Nprt_global_over_Nproc = Nprt_global_org/Nproc;

    int *Nprt_local_array = new int[Nproc];
    int *Nprt_local_disp = new int[Nproc+1];
    Nprt_local_disp[0] = 0;
    for(int i=0; i<Nproc; i++){
        Nprt_local_array[i] = Nprt_global_over_Nproc;
        if(i < Nprt_global_org % Nproc){
            Nprt_local_array[i]++;
        }
        Nprt_local_disp[i+1] = Nprt_local_disp[i] + Nprt_local_array[i];
    }

    double *mass_array = new double[Nprt_local_array[0]+1];
    Vector3 *pos_array = new Vector3[Nprt_local_array[0]+1];
    Vector3 *vel_array = pos_array;

    int tag = 0;
    MPI_Status status;
    for(int i=0; i<Nproc; i++){
        for(int j=Nprt_local_disp[i]; j<Nprt_local_disp[i+1]; j++){
            if(MYRANK==0 ){
                finput.read((char*)(&mass_array[j-Nprt_local_disp[i]]), sizeof(double));
            }
        }
        if(i != 0){
            if(MYRANK==0){
                MPI_Send(mass_array, Nprt_local_array[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            }
            else if(MYRANK==i){
                MPI_Recv(mass_array, Nprt_local_array[i], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
            }
        }
        if(MYRANK==i){
            for(int k=0; k<Nprt_local_array[i]; k++){
                //cerr<<"Nprt_local_array[i]="<<Nprt_local_array[i]<<endl;
                //cerr<<"mass_array[k]="<<mass_array[k]<<endl;
                prt[k].mass = mass_array[k];
                prt[k].index = Nprt_local_disp[i] + k;
                prt[k].type = star;
            }
        }
    }

    int vec_size = sizeof(Vector3);
    for(int i=0; i<Nproc; i++){
        for(int j=Nprt_local_disp[i]; j<Nprt_local_disp[i+1]; j++){
            if(MYRANK==0){
                //finput>>pos_array[j - NFS_local_disp[i]];
                finput.read((char*)&pos_array[j-Nprt_local_disp[i]], sizeof(Vector3));
            }
        }
        if(i != 0){
            if(MYRANK==0){
                MPI_Send(pos_array, Nprt_local_array[i]*vec_size, MPI_BYTE, i, tag, MPI_COMM_WORLD);
            }
            else if(MYRANK==i){
                MPI_Recv(pos_array, Nprt_local_array[i]*vec_size, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status);
            }
        }
        if(MYRANK==i){
            for(int k=0; k<Nprt_local_array[i]; k++){
                prt[k].pos = pos_array[k];
            }
        }
    }

    for(int i=0; i<Nproc; i++){
        for(int j=Nprt_local_disp[i]; j<Nprt_local_disp[i+1]; j++){
            if(MYRANK==0){
                finput.read((char*)&vel_array[j-Nprt_local_disp[i]], sizeof(Vector3));
            }
        }
        if(i != 0){
            if(MYRANK==0){
                MPI_Send(vel_array, Nprt_local_array[i]*vec_size, MPI_BYTE, i, tag, MPI_COMM_WORLD);
            }
            else if(MYRANK==i){
                MPI_Recv(vel_array, Nprt_local_array[i]*vec_size, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status);
            }
        }
        if(MYRANK==i){
            for(int k=0; k<Nprt_local_array[i]; k++){
                prt[k].vel = vel_array[k];
            }
        }
    }


    int Nprt_local = Nprt_local_array[MYRANK];
/*
// NO CHECK
//merge case
    int Nprt_local = NFS_local_array[MYRANK];
    int Ntail = Nprt_local - 1;
    int Ndead_local = 0;
    for(int i=0; i<=Ntail;){
        if(prt[i].mass == 0.0){
            prt_dead_glb[Ndead_local] = prt[i];
            Ndead_local++;
            prt[i] = prt[Ntail];
            Ntail--;
        }
        i++;
    }
*/

    int nomerge = 1;
    if(nomerge){
        NFS_GLB = NFS_GLB_ORG;
        NBH_GLB = NBH_GLB_ORG;
        NDEAD_GLB = 0;
        if(MYRANK == 0 && Nprt_local <= NBH_GLB_ORG){
            halt_program();
        }
        if(MYRANK == 0){
            NBH_LOC = NBH_GLB;
            for(int nbh=0; nbh<NBH_GLB; nbh++){
                BH_glb[nbh] = prt[nbh];
            }
        }
        else{
            NBH_LOC = 0;
        }
        NFS_LOC = Nprt_local - NBH_LOC;
    }

    MPI_Bcast(BH_glb, NBH_GLB_ORG*sizeof(Particle), MPI_BYTE, 0, MPI_COMM_WORLD);

    delete [] mass_array;
    delete [] pos_array;
    delete [] Nprt_local_array;
    delete [] Nprt_local_disp;
}

inline void readparam(char paramfile[],  char initial_file[],  char output_dir[],
		      Particle BH_global[],  int &NBH,
		      double &eps2_FS_FS,  double &eps2_FS_BH,  double &eps2_BH_BH, 
		      double &rcut_out_FS_FS,  double &rcut_out_FS_BH,  double &rcut_out_BH_BH,
		      double &rsearch_FS_FS,  double &rsearch_FS_BH,  double &rsearch_BH_BH,
		      double &Tend,  double &dt_global,
		      int &snpid,  double &snpinterval,  int &read_flag,
		      double &theta2,  int &Nleaf,  int &Ncrit,	 int &quad_flag,
		      double &eta_s,  double &eta_FS, double &eta_BH, 
		      double &vel_light){
    ifstream finput;
    if(MYRANK == 0){
	finput.open(paramfile);
	cerr<<"Paramater File: "<<paramfile<<endl;
	static const int REAL    =  1;
	static const int STRING  =  2;
	static const int INT     =  3;
	static const int VEC     =  4;
	static const int VEC7 = 5;
	static const int MAX_TAG =  100;
  
	char buf[STRINGSIZE*9], 
	    buf1[STRINGSIZE],
	    buf2[STRINGSIZE],
	    buf3[STRINGSIZE],
	    buf4[STRINGSIZE],
	    buf5[STRINGSIZE],
	    buf6[STRINGSIZE],
	    buf7[STRINGSIZE],
	    buf8[STRINGSIZE],
	    buf9[STRINGSIZE];

	int  tag_num=0;
	int  id[MAX_TAG];
	void *addr[MAX_TAG];
	char tag[MAX_TAG][64];

	strcpy(tag[tag_num],"initial_file");
	addr[tag_num]=initial_file;
	id[tag_num++]=STRING;

	strcpy(tag[tag_num],"output_dir");
	addr[tag_num]=output_dir;
	id[tag_num++]=STRING;

	read_flag=0;
	strcpy(tag[tag_num],"read_flag");
	addr[tag_num]=&read_flag;
	id[tag_num++]=INT;

	strcpy(tag[tag_num],"dt_global");
	addr[tag_num]=&dt_global;
	id[tag_num++]=REAL;

	double inv_dt_global=0.0;
	strcpy(tag[tag_num],"inv_dt_global");
	addr[tag_num]=&inv_dt_global;
	id[tag_num++]=REAL;

	double inv_eps_FS_FS=0;
	strcpy(tag[tag_num],"inv_eps_FS_FS");
	addr[tag_num]=&inv_eps_FS_FS;
	id[tag_num++]=REAL;

	double inv_eps_FS_BH=0;
	strcpy(tag[tag_num],"inv_eps_FS_BH");
	addr[tag_num]=&inv_eps_FS_BH;
	id[tag_num++]=REAL;

	double inv_eps_BH_BH=0;
	strcpy(tag[tag_num],"inv_eps_BH_BH");
	addr[tag_num]=&inv_eps_BH_BH;
	id[tag_num++]=REAL;

	double inv_rcut_out_FS_FS = 0.0;
	strcpy(tag[tag_num],"inv_rcut_out_FS_FS");
	addr[tag_num]=&inv_rcut_out_FS_FS;
	id[tag_num++]=REAL;

	double inv_rcut_out_FS_BH = 0.0;
	strcpy(tag[tag_num],"inv_rcut_out_FS_BH");
	addr[tag_num]=&inv_rcut_out_FS_BH;
	id[tag_num++]=REAL;

	double inv_rcut_out_BH_BH = 0.0;
	strcpy(tag[tag_num],"inv_rcut_out_BH_BH");
	addr[tag_num]=&inv_rcut_out_BH_BH;
	id[tag_num++]=REAL;


	double inv_rsearch_FS_FS = 0.0;
	strcpy(tag[tag_num],"inv_rsearch_FS_FS");
	addr[tag_num]=&inv_rsearch_FS_FS;
	id[tag_num++]=REAL;

	double inv_rsearch_FS_BH = 0.0;
	strcpy(tag[tag_num],"inv_rsearch_FS_BH");
	addr[tag_num]=&inv_rsearch_FS_BH;
	id[tag_num++]=REAL;

	double inv_rsearch_BH_BH = 0.0;
	strcpy(tag[tag_num],"inv_rsearch_BH_BH");
	addr[tag_num]=&inv_rsearch_BH_BH;
	id[tag_num++]=REAL;


	double vel_disp = 0.0;
	strcpy(tag[tag_num],"vel_disp");
	addr[tag_num]=&vel_disp;
	id[tag_num++]=REAL;

	double search_factor = 0.0;
	strcpy(tag[tag_num],"search_factor");
	addr[tag_num]=&search_factor;
	id[tag_num++]=REAL;



	strcpy(tag[tag_num],"Tend");
	addr[tag_num]=&Tend;
	id[tag_num++]=REAL;

	double theta=0.0;
	strcpy(tag[tag_num],"theta");
	addr[tag_num]=&theta;
	id[tag_num++]=REAL;

	strcpy(tag[tag_num],"Nleaf");
	addr[tag_num]=&Nleaf;
	id[tag_num++]=INT;

	strcpy(tag[tag_num],"Ncrit");
	addr[tag_num]=&Ncrit;
	id[tag_num++]=INT;

	double inv_snpinterval = 0.0;
	strcpy(tag[tag_num],"inv_snpinterval");
	addr[tag_num]=&inv_snpinterval;
	id[tag_num++]=REAL;

	strcpy(tag[tag_num],"snpinterval");
	addr[tag_num]=&snpinterval;
	id[tag_num++]=REAL;




	strcpy(tag[tag_num],"NBH");
	addr[tag_num]=&NBH;
	id[tag_num++]=INT;

	strcpy(tag[tag_num],"snpid");
	addr[tag_num]=&snpid;
	id[tag_num++]=INT;

	strcpy(tag[tag_num],"BHmxv");
	addr[tag_num]=NULL;
	id[tag_num++]=VEC7;

	strcpy(tag[tag_num],"eta_s");
	addr[tag_num]=&eta_s;
	id[tag_num++]=REAL;

	strcpy(tag[tag_num],"eta_FS");
	addr[tag_num]=&eta_FS;
	id[tag_num++]=REAL;
	
	strcpy(tag[tag_num],"eta_BH");
	addr[tag_num]=&eta_BH;
	id[tag_num++]=REAL;
	
	strcpy(tag[tag_num],"quad_flag");
	addr[tag_num]=&quad_flag;
	id[tag_num++]=INT;

	strcpy(tag[tag_num],"vel_light");
	addr[tag_num]=&vel_light;
	id[tag_num++]=REAL;

	int BHid=0;

	while(!finput.eof()){
	    finput.getline(buf,STRINGSIZE*9);
	    if(sscanf(buf,"%s %s %s %s %s %s %s %s %s",buf1,buf2,buf3,buf4,buf5,buf6,buf7,buf8,buf9)<2)continue;
	    if(buf1[0]=='%' || buf1[0]=='#')continue;
	    cerr<<"buf1: "<<buf1<<endl;
	    cerr<<"buf2: "<<buf2<<endl;
	    for(int i=0; i<tag_num; i++){
		if(strcmp(buf1,tag[i])==0){
		    switch (id[i]){
		    case REAL:
			*((double *) addr[i]) = atof(buf2);
			break;
		    case STRING:
			strcpy((char*)addr[i], buf2);
			break;
		    case INT:
			*((int *) addr[i]) = atoi(buf2);
			break;
		    case VEC:
			*((Vector3 *) addr[i]) = atov(buf2, buf3, buf4);
			break;
		    case VEC7:
			BH_global[BHid].mass = atof(buf2);
			BH_global[BHid].pos = atov(buf3, buf4, buf5);
			BH_global[BHid].vel = atov(buf6, buf7, buf8);
			BH_global[BHid].index = BHid;
			BH_global[BHid].type = blackhole;
			BHid++;
			break;
		    default:
			cerr<<"invarid flag is used"<<endl;
			halt_program();
			//exit(-1);
		    }
		    cerr<<buf2<<endl;
		}
	    }
	}
	if(inv_dt_global != 0.0){
	    dt_global = 1.0/inv_dt_global;
	}
	if(inv_snpinterval != 0.0){
	    snpinterval = 1.0/inv_snpinterval;
	}

	theta2=theta*theta;

	if(inv_rcut_out_FS_FS > 0.0){rcut_out_FS_FS = 1.0 / inv_rcut_out_FS_FS;}
	else{ rcut_out_FS_FS = 0.0;}
	if(inv_rcut_out_FS_BH > 0.0){rcut_out_FS_BH = 1.0 / inv_rcut_out_FS_BH;}
	else{ rcut_out_FS_BH = 0.0;}
	if(inv_rcut_out_BH_BH > 0.0){rcut_out_BH_BH = 1.0 / inv_rcut_out_BH_BH;}
	else{rcut_out_BH_BH = 0.0;}
	

	if(search_factor <= 0.0){
	    if(inv_rsearch_FS_FS > 0.0){rsearch_FS_FS = 1.0 / inv_rsearch_FS_FS;}
	    else{rsearch_FS_FS = 0.0;}
	    if(inv_rsearch_FS_BH > 0.0){rsearch_FS_BH = 1.0 / inv_rsearch_FS_BH;}
	    else{rsearch_FS_BH = 0.0;}
	    if(inv_rsearch_BH_BH > 0.0){rsearch_BH_BH = 1.0 / inv_rsearch_BH_BH;}
	    else{rsearch_BH_BH = 0.0;}
	}
	else{
	    rsearch_FS_FS = rcut_out_FS_FS + search_factor*vel_disp*dt_global;
	    rsearch_FS_BH = rcut_out_FS_BH + search_factor*vel_disp*dt_global;
	    rsearch_BH_BH = rcut_out_BH_BH + search_factor*vel_disp*dt_global;
	}

	
	if(inv_eps_FS_FS > 0.0){eps2_FS_FS = 1.0/(inv_eps_FS_FS*inv_eps_FS_FS);}
	else{eps2_FS_FS = 0.0;}
	if(inv_eps_FS_BH > 0.0){eps2_FS_BH = 1.0/(inv_eps_FS_BH*inv_eps_FS_BH);}
	else{eps2_FS_BH = 0.0;}
	if(inv_eps_BH_BH > 0.0){eps2_BH_BH = 1.0/(inv_eps_BH_BH*inv_eps_BH_BH);}
	else{eps2_BH_BH = 0.0;}
	
	if(snpinterval < dt_global){
	    halt_program();
	}

	if(rsearch_FS_FS < rcut_out_FS_FS || rsearch_FS_BH < rcut_out_FS_BH || rsearch_BH_BH < rcut_out_BH_BH){
	    halt_program();
	}

    }

    //cerr<<"read_flag="<<read_flag<<endl;

    MPI_Bcast(&NBH, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps2_FS_FS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps2_FS_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps2_BH_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Tend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&theta2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snpinterval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rcut_out_FS_FS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rcut_out_FS_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rcut_out_BH_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsearch_FS_FS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsearch_FS_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsearch_BH_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nleaf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ncrit, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snpid, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&read_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(BH_global, NBH*sizeof(Particle), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eta_s, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eta_FS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eta_BH, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&quad_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&vel_light, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}


struct Double_List{
  double val;
  int idx;
};
struct Vector3_List{
  Vector3 val;
  int idx;
};

inline void write0_gather(Particle prt[],
			  Particle prt_dead_glb[],
			  const double &Tsys,
			  const double &Egr, 
			  char *dirname, 
			  int &snpid){

  static const int Ndev = 10000;
  static Double_List mass_list_send[Ndev];
  static Double_List mass_list_recv[Ndev];
  static int adr_from_idx[Ndev];

  ofstream fout;
  char sout[STRINGSIZE];
  sprintf(sout,"%s/snap%5d.dat", dirname, snpid);
  for(int i=0;i<STRINGSIZE;i++)if(sout[i]==' ')sout[i]='0';
  int Nprt_glb_org = NFS_GLB_ORG + NBH_GLB_ORG;
  if(MYRANK==0){
    cerr<<"creating snapshot..."<<sout<<endl;
    fout.open(sout);
    fout.write((const char*)&Nprt_glb_org, sizeof(int));
    fout.write((const char*)&NBH_GLB_ORG, sizeof(int));
    fout.write((const char*)&Tsys, sizeof(double));
    fout.write((const char*)&Egr, sizeof(double));
    cerr<<"Nprt_glb_org="<<Nprt_glb_org<<endl;
    cerr<<"NBH_GLB_ORG="<<NBH_GLB_ORG<<endl;
    cerr<<"Tsys="<<Tsys<<endl;
    cerr<<"Egr="<<Egr<<endl;
  }

  //cerr<<"NFS_LOC+NBH_LOC="<<NFS_LOC+NBH_LOC<<endl;
  int *idx_array = new int[NFS_LOC+NBH_LOC];
  int *adr_array = new int[NFS_LOC+NBH_LOC];
  for(int i=0; i<NFS_LOC+NBH_LOC; i++){
    adr_array[i] = i;
    idx_array[i] = prt[i].index;
  }

  //cerr<<"-d: MYRANK="<<MYRANK<<endl;
  if(NFS_LOC+NBH_LOC > 0){
    Qsort_index(idx_array, adr_array, 0, NFS_LOC+NBH_LOC-1);
  }

  /*
  for(int i=0; i<NFS_LOC+NBH_LOC; i++){
    cerr<<"i="<<i<<endl;
    cerr<<"adr_array[i]="<<adr_array[i]<<endl;
    cerr<<"prt[adr_array[i]].index="<<prt[adr_array[i]].index<<endl;
  }
  */


  int *idx_dead_array = new int[NDEAD_GLB];
  int *adr_dead_array = new int[NDEAD_GLB];
  for(int i=0; i<NDEAD_GLB; i++){
    adr_dead_array[i] = i;
    idx_dead_array[i] = prt_dead_glb[i].index;
  }
  if(NDEAD_GLB > 0){
    Qsort_index(idx_dead_array, adr_dead_array, 0, NDEAD_GLB-1);
  }


  //cerr<<"-a: MYRANK="<<MYRANK<<endl;

  int i_loc = 0;
  int i_dead_loc = 0;
  for(int n=0; n<NFS_GLB_ORG+NBH_GLB_ORG; n += Ndev){
    int list_len_loc = 0;
    int list_len_glb = 0;
    if(i_loc < NFS_LOC+NBH_LOC){
      while(prt[adr_array[i_loc]].index < n+Ndev){
	mass_list_send[list_len_loc].val = prt[adr_array[i_loc]].mass;
	mass_list_send[list_len_loc].idx = prt[adr_array[i_loc]].index;
	list_len_loc++;
	i_loc++;
	if(i_loc >= NFS_LOC+NBH_LOC){break;}
      }
    }
    /*
    while(i_loc < NFS_LOC+NBH_LOC && prt[adr_array[i_loc]].index < n+Ndev){
      //cerr<<"i_loc="<<i_loc<<endl;
      //cerr<<"adr_array[i_loc]="<<adr_array[i_loc]<<endl;
      //cerr<<"prt[adr_array[i_loc]].index="<<prt[adr_array[i_loc]].index<<endl;
      mass_list_send[list_len_loc].val = prt[adr_array[i_loc]].mass;
      mass_list_send[list_len_loc].idx = prt[adr_array[i_loc]].index;
      list_len_loc++;
      i_loc++;
      if(i_loc >= NFS_LOC+NBH_LOC){break;}
    }
    */

    if(NDEAD_GLB > 0 && MYRANK == 0 && i_dead_loc < NDEAD_GLB){
      while(prt_dead_glb[adr_dead_array[i_dead_loc]].index < n+Ndev){
	mass_list_send[list_len_loc].val = prt_dead_glb[adr_dead_array[i_dead_loc]].mass;
	mass_list_send[list_len_loc].idx = prt_dead_glb[adr_dead_array[i_dead_loc]].index;
	list_len_loc++;
	i_dead_loc++;
	if(i_dead_loc >= NDEAD_GLB){break;}
      }
    }


    //cerr<<"-0: MYRANK="<<MYRANK<<endl;
    //cerr<<"list_len_loc="<<list_len_loc<<endl;
    //cerr<<"i_loc="<<i_loc<<endl;
    mpi_allgatherv_T(mass_list_send, list_len_loc, mass_list_recv, list_len_glb, NPROC);
    for(int i=0; i<list_len_glb; i++){
      adr_from_idx[mass_list_recv[i].idx-n] = i;
    }
    if(MYRANK == 0){
      //cerr<<"list_len_glb="<<list_len_glb<<endl;
      for(int i=0; i<list_len_glb; i++){
	//cerr<<"i="<<i<<endl;
	//cerr<<"adr_from_idx[i]="<<adr_from_idx[i]<<endl;
	//cerr<<"mass_list_recv[adr_from_idx[i]].idx="<<mass_list_recv[adr_from_idx[i]].idx<<endl;
	fout.write((const char*)&mass_list_recv[adr_from_idx[i]].val, sizeof(double));
	//cout<<mass_list_recv[adr_from_idx[i]].val<<endl;
      }
    }
  }

  //cerr<<"a: MYRANK="<<MYRANK<<endl;


  static Vector3_List pos_list_send[Ndev];
  static Vector3_List pos_list_recv[Ndev];
  i_loc = 0;
  i_dead_loc = 0;
  for(int n=0; n<NFS_GLB_ORG+NBH_GLB_ORG; n += Ndev){
    int list_len_loc = 0;
    int list_len_glb = 0;
    if(i_loc < NFS_LOC+NBH_LOC){
      while(prt[adr_array[i_loc]].index < n+Ndev){
	pos_list_send[list_len_loc].val = prt[adr_array[i_loc]].pos;
	pos_list_send[list_len_loc].idx = prt[adr_array[i_loc]].index;
	list_len_loc++;
	i_loc++;
	if(i_loc >= NFS_LOC+NBH_LOC){break;}
      }
    }
    if(NDEAD_GLB > 0 && MYRANK == 0 && i_dead_loc < NDEAD_GLB){
      while(prt_dead_glb[adr_dead_array[i_dead_loc]].index < n+Ndev){
	pos_list_send[list_len_loc].val = prt_dead_glb[adr_dead_array[i_dead_loc]].pos;
	pos_list_send[list_len_loc].idx = prt_dead_glb[adr_dead_array[i_dead_loc]].index;
	list_len_loc++;
	i_dead_loc++;
	if(i_dead_loc >= NDEAD_GLB){break;}
      }
    }
    mpi_allgatherv_T(pos_list_send, list_len_loc, pos_list_recv, list_len_glb, NPROC);
    for(int i=0; i<list_len_glb; i++){
      adr_from_idx[pos_list_recv[i].idx-n] = i;
    }
    if(MYRANK == 0){
      for(int i=0; i<list_len_glb; i++){
	fout.write((const char*)&pos_list_recv[adr_from_idx[i]].val, sizeof(Vector3));
	//cout<<pos_list_recv[adr_from_idx[i]].val<<endl;
      }
    }
  }

  //cerr<<"b: MYRANK="<<MYRANK<<endl;

  static Vector3_List *vel_list_send = pos_list_send;
  static Vector3_List *vel_list_recv = pos_list_recv;
  i_loc = 0;
  i_dead_loc = 0;
  for(int n=0; n<NFS_GLB_ORG+NBH_GLB_ORG; n += Ndev){
    int list_len_loc = 0;
    int list_len_glb = 0;
    if(i_loc < NFS_LOC+NBH_LOC){
      while(prt[adr_array[i_loc]].index < n+Ndev){
	vel_list_send[list_len_loc].val = prt[adr_array[i_loc]].vel;
	vel_list_send[list_len_loc].idx = prt[adr_array[i_loc]].index;
	list_len_loc++;
	i_loc++;
	if(i_loc >= NFS_LOC+NBH_LOC){break;}
      }
    }
    if(NDEAD_GLB > 0 && MYRANK == 0 && i_dead_loc < NDEAD_GLB){
      while(prt_dead_glb[adr_dead_array[i_dead_loc]].index < n+Ndev){
	vel_list_send[list_len_loc].val = prt_dead_glb[adr_dead_array[i_dead_loc]].vel;
	vel_list_send[list_len_loc].idx = prt_dead_glb[adr_dead_array[i_dead_loc]].index;
	list_len_loc++;
	i_dead_loc++;
	if(i_dead_loc >= NDEAD_GLB){break;}
      }
    }
    mpi_allgatherv_T(vel_list_send, list_len_loc, vel_list_recv, list_len_glb, NPROC);
    for(int i=0; i<list_len_glb; i++){
      adr_from_idx[vel_list_recv[i].idx-n] = i;
    }
    if(MYRANK == 0){
      for(int i=0; i<list_len_glb; i++){
	fout.write((const char*)&vel_list_recv[adr_from_idx[i]].val, sizeof(Vector3));
	//cout<<vel_list_recv[adr_from_idx[i]].val<<endl;
      }
    }
  }

  //cerr<<"c: MYRANK="<<MYRANK<<endl;

  static Double_List *pot_list_send = mass_list_send;
  static Double_List *pot_list_recv = mass_list_recv;
  i_loc = 0;
  i_dead_loc = 0;
  for(int n=0; n<NFS_GLB_ORG+NBH_GLB_ORG; n += Ndev){
    int list_len_loc = 0;
    int list_len_glb = 0;
    if(i_loc < NFS_LOC+NBH_LOC){
      while(prt[adr_array[i_loc]].index < n+Ndev){
	pot_list_send[list_len_loc].val = prt[adr_array[i_loc]].pot;
	pot_list_send[list_len_loc].idx = prt[adr_array[i_loc]].index;
	list_len_loc++;
	i_loc++;
	if(i_loc >= NFS_LOC+NBH_LOC){break;}
      }
    }
    if(NDEAD_GLB > 0 && MYRANK == 0 && i_dead_loc < NDEAD_GLB){
      while(prt_dead_glb[adr_dead_array[i_dead_loc]].index < n+Ndev){
	pot_list_send[list_len_loc].val = prt_dead_glb[adr_dead_array[i_dead_loc]].pot;
	pot_list_send[list_len_loc].idx = prt_dead_glb[adr_dead_array[i_dead_loc]].index;
	list_len_loc++;
	i_dead_loc++;
	if(i_dead_loc >= NDEAD_GLB){break;}
      }
    }
    mpi_allgatherv_T(pot_list_send, list_len_loc, pot_list_recv, list_len_glb, NPROC);
    for(int i=0; i<list_len_glb; i++){
      adr_from_idx[pot_list_recv[i].idx-n] = i;
    }
    if(MYRANK == 0){
      for(int i=0; i<list_len_glb; i++){
	fout.write((const char*)&pot_list_recv[adr_from_idx[i]].val, sizeof(double));
      }
    }
  }

  delete[] idx_array;
  delete[] adr_array;

  delete[] idx_dead_array;
  delete[] adr_dead_array;

  if(MYRANK==0){
    fout.close();
    cerr<<"finish writing snap shot: "<<sout<<endl;
  }

}




#endif //IO_H
