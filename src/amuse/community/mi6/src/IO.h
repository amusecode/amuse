#include"Particle.h"
#include"external_field.h"

inline void read0(Particle *&prt, 
		  int &Ntot, 
		  const int &NBH, 
		  int &NFS, 
		  double &Tsys, 
		  char *sinput){
  ifstream finput;
  finput.open(sinput);
  cerr<<"Input File: "<<sinput<<endl;
  finput>>NFS;
  Ntot = NFS + NBH;
  cerr<<"Ntot="<<Ntot<<endl;
  int dim;
  finput>>dim;
  cerr<<"dim="<<dim<<endl;
  finput>>Tsys;
  cerr<<"Tsys="<<Tsys<<endl;
  prt = new Particle[Ntot];
  cerr<<"creat new particle"<<endl;
  for(int i=NBH; i<Ntot; i++){finput>>prt[i].mass;}
  for(int i=NBH; i<Ntot; i++){finput>>prt[i].pos;}
  for(int i=NBH; i<Ntot; i++){finput>>prt[i].vel;}
  finput.close();
  cerr<<"Close inputfile: "<<sinput<<endl;
}  
  
inline void readpara0(Particle *&prt, 
		      int &Ntot, 
		      int &NSMBH, 
		      int &NIMBH, 
		      int &NBH, 
		      int &NFS, 
		      double &eps2_fs, 
		      double &eps2_bh, 
		      double &eps2_fs_smbh, 
		      double &eps2_fs_imbh, 
		      double &eta_s, 
		      double &eta_fs, 
		      double &eta_smbh, 
		      double &eta_imbh,
		      double &Tsys, 
		      double &Tmerge, 
		      double &Tend, 
		      char *sinput){
  /*
////////// input parameta file /////////
1: name of nemo ascii file
2: number of SMBHs
3: number of IMBHs
4: eps2_fs
5: eps2_bh
6: eps2_fs_smbh
7: eps2_fs_imbh
8: eta start
9: eta FS
10: eta SMBH
11: eta IMBH
12: end time
13: BH0 mass BH1 mass BH2 mass ...
14: BH0 pos BH1 pos BH2 pos ...
15: BH0 vel BH1 vel BH2 vel ...
  */
  ifstream finput;
  cerr<<"sinput: "<<sinput<<endl;
  finput.open(sinput);
  char sinput2[1024];
  finput>>sinput2;
  cerr<<"sinput2: "<<sinput2<<endl;
  finput>>NSMBH;
  cerr<<"NSMBH="<<NSMBH<<endl;
  finput>>NIMBH;
  cerr<<"NIMBH="<<NIMBH<<endl;
  NBH = NSMBH + NIMBH;
  cerr<<"NBH="<<NIMBH<<endl;
  finput>>eps2_fs;
  cerr<<"eps2_fs="<<eps2_fs<<endl;
  finput>>eps2_bh; 
  cerr<<"eps2_bh="<<eps2_bh<<endl;
  finput>>eps2_fs_smbh; 
  cerr<<"eps2_fs_smbh="<<eps2_fs_smbh<<endl;
  finput>>eps2_fs_imbh; 
  cerr<<"eps2_fs_imbh="<<eps2_fs_imbh<<endl;
  finput>>eta_s; 
  cerr<<"eta_s="<<eta_s<<endl;
  finput>>eta_fs;
  cerr<<"eta_fs="<<eta_fs<<endl;
  finput>>eta_smbh;
  cerr<<"eta_smbh="<<eta_smbh<<endl;
  finput>>eta_imbh;
  cerr<<"eta_imbh="<<eta_imbh<<endl;
  finput>>Tend;
  cerr<<"Tend="<<Tend<<endl;
  read0(prt, Ntot, NBH, NFS, Tsys, sinput2);
  Tmerge=Tsys;
  cerr<<"Tsys="<<Tsys<<endl;
  cerr<<"Tmerge="<<Tmerge<<endl;
  for(int i=0; i<NBH; i++){ finput>>prt[i].mass;}
  for(int i=0; i<NBH; i++){ finput>>prt[i].pos;}
  for(int i=0; i<NBH; i++){ finput>>prt[i].vel;}
  if(NBH > 0){
    cerr<<prt[0].mass<<endl;
    cerr<<prt[0].pos<<endl;
    cerr<<prt[0].vel<<endl;
  }
  for(int i=0; i<Ntot; i++){
    prt[i].index = i;
    //prt[i].index = i+NBH;
    prt[i].address = i;
  }
  finput.close();
}

// EXFLAG 0:no extr field,   1:point mass at the center
inline void write0(Particle *prt, 
		   const int &Ntot, 
		   const int &NBH, 
		   const int &Nmerge, 
		   const double &Tsys, 
		   const double &Tmerge, 
		   const double &Egr,
		   char *dirname, 
		   int &snpid, 
		   const int &EXFLAG){

  int ERROR=0;
  while(ERROR==0){
    ofstream fout;
    char sout[1024];
    sprintf(sout,"%s/snap%5d.dat",dirname,snpid);
    for(int i=0;i<1024;i++)if(sout[i]==' ')sout[i]='0';
    cerr<<"creating SNAP file..."<<sout<<endl;
    cerr<<"Tsys="<<Tsys<<endl;
    fout.open(sout);
    int TMP_Ntot=Ntot;
    int TMP_NBH=NBH;
    if(EXFLAG==1){
      TMP_Ntot += 1; 
      TMP_NBH += 1;
    }
    fout.write((const char*)&TMP_Ntot,sizeof(int));
    fout.write((const char*)&TMP_NBH,sizeof(int));
    fout.write((const char*)&Tsys,sizeof(double));
    fout.write((const char*)&Egr,sizeof(double));

    double SMBHmass = 0.0;
    Vector3 SMBHpos = 0.0;
    Vector3 SMBHvel = 0.0;
    get_SMBH(SMBHmass, SMBHpos, SMBHvel);
    double SMBHphi = 0.0;
    int SMBHindex = 0;
    if(EXFLAG==1){fout.write((const char*)&SMBHmass,sizeof(double));}
    for(int i=0; i<Ntot; i++){fout.write((const char*)&prt[i].mass,sizeof(double));}

    if(EXFLAG==1){fout.write((const char*)&SMBHpos,sizeof(Vector3));}
    for(int i=0; i<Ntot; i++){fout.write((const char*)&prt[i].pos,sizeof(Vector3));}

    if(EXFLAG==1){fout.write((const char*)&SMBHvel,sizeof(Vector3));}
    for(int i=0; i<Ntot; i++){fout.write((const char*)&prt[i].vel,sizeof(Vector3));}

    if(EXFLAG==1){fout.write((const char*)&SMBHphi,sizeof(double));}
    for(int i=0; i<Ntot; i++){fout.write((const char*)&prt[i].phi,sizeof(double));}

    if(EXFLAG==1){fout.write((const char*)&SMBHindex,sizeof(int));}
    for(int i=0; i<Ntot; i++){
      int index_tmp = prt[i].index + 1;
      fout.write((const char*)&index_tmp,sizeof(int));
    }
    fout.write((const char*)&Nmerge,sizeof(int));
    fout.close();
    cerr<<"finish writing : "<<sout<<endl;
    ERROR=1;
    /*
    ifstream fout2;
    fout2.open(sout);
    fout2.seekg(0,ios::end);
    long size = fout2.tellg();
    const static long SIZE0 = size;
    fout2.close();
    cerr<<sout<<": T="<<setprecision(15)<<Tsys<<", LOG size: "<<size<<endl;
    if(SIZE0==size){
      cerr<<"finish writing : "<<sout<<endl;
      ERROR=1;
    }
    else{
      cerr<<"snap size ERROR file name: "<<sout<<endl;
    }
    */

  }
  snpid++;
}

inline void read1(Particle *&prt, 
		  int &Ntot, 
		  int &NBH, 
		  int &NFS, 
		  int &Nmerge, 
		  double &Tsys, 
		  char *sinput, 
		  double &Egr, 
		  const int &EXFLAG){
  cerr<<"reading SNAP file..."<<sinput<<endl;
  ifstream fin;
  fin.open(sinput);
  fin.read((char*)&Ntot,sizeof(int));
  fin.read((char*)&NBH,sizeof(int));

  if(EXFLAG==1){
    Ntot -= 1; 
    NBH -= 1;
  }
  cerr<<"Ntot="<<Ntot<<endl;
  cerr<<"NBH="<<NBH<<endl;
  NFS=Ntot-NBH;
  fin.read((char*)&Tsys,sizeof(double));
  cerr<<"Tsys="<<Tsys<<endl;
  fin.read((char*)&Egr,sizeof(double));
  cerr<<"Egr="<<Egr<<endl;
  Egr=0.0;
  cerr<<"creat new particle"<<endl;
  prt = new Particle[Ntot];
  double SMBHmass = 0.0;
  Vector3 SMBHpos = 0.0;
  Vector3 SMBHvel = 0.0;
  double SMBHphi = 0.0;
  int SMBHindex = 0;
  if(EXFLAG==1){fin.read((char*)&SMBHmass,sizeof(double));}
  for(int i=0; i<Ntot; i++){fin.read((char*)&prt[i].mass,sizeof(double));}
  if(EXFLAG==1){fin.read((char*)&SMBHpos,sizeof(Vector3));}
  for(int i=0; i<Ntot; i++){fin.read((char*)&prt[i].pos,sizeof(Vector3));}
  if(EXFLAG==1){fin.read((char*)&SMBHvel,sizeof(Vector3));}
  for(int i=0; i<Ntot; i++){fin.read((char*)&prt[i].vel,sizeof(Vector3));}
  if(EXFLAG==1){fin.read((char*)&SMBHphi,sizeof(double));}
  for(int i=0; i<Ntot; i++){fin.read((char*)&prt[i].phi,sizeof(double));}
  if(EXFLAG==1){fin.read((char*)&SMBHindex,sizeof(int));}
  if(EXFLAG==1){
    set_SMBH(SMBHmass, SMBHpos, SMBHvel);
  }
  int itmp = 0;
  for(int i=0; i<Ntot; i++){
    int indextmp=0;
    fin.read((char*)&indextmp,sizeof(int));
    prt[i].index = indextmp - 1;
    if(prt[i].mass != 0.0){
      prt[i].address = itmp;
      itmp++;
    }
  }
  fin.read((char*)&Nmerge,sizeof(int));
  cerr<<"N_merge="<<Nmerge<<endl;
  for(int i=0; i<Ntot; i++){prt[i].time = Tsys;}
  fin.close();

}  


inline void readpara1(Particle *&prt, 
		      int &Ntot, 
		      int &NSMBH, 
		      int &NIMBH, 
		      int &NBH, 
		      int &NFS, 
		      int &Nmerge,
		      double &eps2_fs, 
		      double &eps2_bh, 
		      double &eps2_fs_smbh, 
		      double &eps2_fs_imbh, 
		      double &eta_s, 
		      double &eta_fs, 
		      double &eta_smbh, 
		      double &eta_imbh,
		      double &Egr, 
		      double &Tsys, 
		      double &Tmerge, 
		      double &Tend, 
		      char *sinput, 
		      const int &EXFLAG){

  /*
////////// input parameta file /////////
1: name of nemo ascii file
2: number of SMBHs  // added snap shot
3: number of IMBHs  // added snap shot
4: eps2_fs
5: eps2_bh
6: eps2_fs_smbh
7: eps2_fs_imbh
8: eta start
9: eta FS
10: eta SMBH
11: eta IMBH
12: end time
13: BH0 mass BH1 mass BH2 mass ...
14: BH0 pos BH1 pos BH2 pos ...
15: BH0 vel BH1 vel BH2 vel ...
  */

  ifstream fin;
  cerr<<"sinput: "<<sinput<<endl;
  fin.open(sinput);
  char sinput2[1024];
  fin>>sinput2;
  cerr<<"sinput2: "<<sinput2<<endl;

  int NSMBHadd;
  fin>>NSMBHadd;
  cerr<<"NSMBHadd="<<NSMBHadd<<endl;
  int NIMBHadd;
  fin>>NIMBHadd;
  cerr<<"NIMBHadd="<<NIMBHadd<<endl;
  int NBHadd;
  NBHadd = NSMBHadd + NIMBHadd;
  cerr<<"NBHadd="<<NBHadd<<endl;

  fin>>eps2_fs; fin>>eps2_bh; fin>>eps2_fs_smbh; fin>>eps2_fs_imbh;
  cerr<<"eps2_fs="<<eps2_fs<<endl;
  cerr<<"eps2_bh="<<eps2_bh<<endl;
  cerr<<"eps2_fs_smbh="<<eps2_fs_smbh<<endl;
  cerr<<"eps2_fs_imbh="<<eps2_fs_imbh<<endl;
  fin>>eta_s; fin>>eta_fs; fin>>eta_smbh; fin>>eta_imbh;
  fin>>Tend;
  cerr<<"Tend="<<Tend<<endl;
  Particle *prttmp;
  read1(prttmp, Ntot, NBH, NFS, Nmerge, Tsys, sinput2, Egr, EXFLAG);
  Tmerge=Tsys;

  if(EXFLAG==1){NSMBH=0; NIMBH=NBH;}
  Ntot += NBHadd;
  NBH += NBHadd;
  NIMBH = NBH;
  prt = new Particle[Ntot];
  for(int i=0; i<NBHadd; i++){ fin>>prt[i].mass;}
  for(int i=0; i<NBHadd; i++){ fin>>prt[i].pos;}
  for(int i=0; i<NBHadd; i++){ fin>>prt[i].vel;}
  for(int i=0; i<Ntot-NBHadd; i++){prt[i+NBHadd]=prttmp[i];}
  for(int i=0; i<Ntot; i++){prt[i].index=i;}

  /*
  cerr<<"prttmp[0].mass="<<prttmp[0].mass<<endl;
  cerr<<"prttmp[0].pos="<<prttmp[0].pos<<endl;
  cerr<<"prttmp[0].vel="<<prttmp[0].vel<<endl;
  cerr<<"prttmp[1].mass="<<prttmp[1].mass<<endl;
  cerr<<"prttmp[1].pos="<<prttmp[1].pos<<endl;
  cerr<<"prttmp[1].vel="<<prttmp[1].vel<<endl;
  cerr<<"prttmp[2].mass="<<prttmp[2].mass<<endl;
  cerr<<"prttmp[2].pos="<<prttmp[2].pos<<endl;
  cerr<<"prttmp[2].vel="<<prttmp[2].vel<<endl;

  cerr<<"prt[0].mass="<<prt[0].mass<<endl;
  cerr<<"prt[0].pos="<<prt[0].pos<<endl;
  cerr<<"prt[0].vel="<<prt[0].vel<<endl;
  cerr<<"prt[1].mass="<<prt[1].mass<<endl;
  cerr<<"prt[1].pos="<<prt[1].pos<<endl;
  cerr<<"prt[1].vel="<<prt[1].vel<<endl;
  cerr<<"prt[2].mass="<<prt[2].mass<<endl;
  cerr<<"prt[2].pos="<<prt[2].pos<<endl;
  cerr<<"prt[2].vel="<<prt[2].vel<<endl;
  */
  fin.close();
}

