

void main
{
  
  return;
}


int ReadRamses(KD kd,FILE *fp)
{
  int dummy,npart,i,j,ncpu,ndim;
  float *temp;

  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&ncpu, sizeof(ncpu), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&ndim, sizeof(ndim), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&npart, sizeof(ndim), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);

  kd->nActive = npart;
  kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
  
  printf("N part = %d \n",kd->nActive);
  printf("N cpus = %d \n",ncpu);
  printf("N dimensions = %d \n",ndim);

  printf("Memory Allocated\n");


  printf("Reading Positions\n");

  temp=(float *)malloc(kd->nActive*sizeof(float));

  for(i=0;i<=ndim-1;i++)
    {
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&temp,sizeof(float),npart,fp);
      fread(&dummy, sizeof(dummy), 1, fp);
      
      for(j=0;j<=npart-1;j++)
	{
	  kd->p[j].r[i] = temp[j];
	}
    }
  printf("Positions Done\n");
  
  // Skipping velocities
  for(i=0;i<=ndim-1;i++)
    {
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&temp,sizeof(float),npart,fp);
      fread(&dummy, sizeof(dummy), 1, fp);
    }

  //Reading Mass
  for(i=0;i<=0;i++)
    {
      fread(&dummy, sizeof(dummy), 1, fp);
      fread(&temp,sizeof(float),npart,fp);
      fread(&dummy, sizeof(dummy), 1, fp);
    }


  kd->fMass = temp[1];	
  printf("Mass Part. =%f\n", kd->fMass);
  printf("Mass Part. =%f\n", temp[100]);

  return kd->nActive;	

}
