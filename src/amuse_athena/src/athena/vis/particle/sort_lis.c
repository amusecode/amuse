/*==============================================================================
 * FILE: sort_lis.c
 *
 * PURPOSE: Sort the particles in the binary output by particle id for 
 *   visualization and analysis. The output file will use the same file name.
 *
 * COMPILE USING: gcc -Wall -W -o sort_lis sort_lis.c -lm
 *
 * USAGE: ./sort_lis -d <dir> -i <basename-in> -s <post-name> -f <# range(f1:f2)>
 *
 * WRITTEN BY: Xuening Bai, September 2009
 *============================================================================*/

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static void quicksort(int *cpuid, long *pid, long *order, long ns, long ne);
static void sort_error(const char *fmt, ...);
static void usage(const char *arg);

static void* calloc_1d_array(size_t nc, size_t size);
static void** calloc_2d_array(size_t nr, size_t nc, size_t size);
static void free_2d_array(void *array);

/* ========================================================================== */

int main(int argc, char* argv[])
{
  /* argument variables */
  int f1=0,f2=0,fi=1;
  char *defdir = ".";
  char *inbase = NULL, *postname = NULL;
  char *indir = defdir;
  /* fild variables */
  FILE *fid;
  char fname[100];
  /* data variables */
  int i,*cpuid,*property,ntype;
  long p,n,*pid,*order;
  float header[20],**data,time[2],*typeinfo=NULL;

  /* Read Arguments */
  for (i=1; i<argc; i++) {
/* If argv[i] is a 2 character string of the form "-?" then: */
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'd':                                /* -d <indir> */
        indir = argv[++i];
        break;
      case 'i':                                /* -i <basename>   */
        inbase = argv[++i];
        break;
      case 's':                                /* -s <post-name> */
        postname = argv[++i];
        break;
      case 'f':                                /* -f <# range(f1:f2:fi)>*/
        sscanf(argv[++i],"%d:%d:%d",&f1,&f2,&fi);
        if (f2 == 0) f2 = f1;
        break;
      case 'h':                                /* -h */
        usage(argv[0]);
        break;
      default:
        usage(argv[0]);
        break;
      }
    }
  }

  /* Checkpoints */
  if (inbase == NULL)
    sort_error("Please specify input file basename using -i option!\n");

  if (postname == NULL)
    sort_error("Please specify posterior file name using -s option!\n");

  if ((f1>f2) || (f2<0) || (fi<=0))
    sort_error("Wrong number sequence in the -f option!\n");

  /* ====================================================================== */

  for (i=f1; i<=f2; i+=fi)
  {
    fprintf(stderr,"Processing file number %d...\n",i);

    /* Step 1: Read the file */
    sprintf(fname,"%s/%s.%04d.%s.lis",indir,inbase,i,postname);

    fid = fopen(fname,"rb");
    if (fid == NULL)
        sort_error("Fail to open output file %s!\n",fname);

    /* read header */
    fread(header,sizeof(float),12,fid);
    fread(&ntype,sizeof(int),1,fid);

    if (i == f1)
      typeinfo = (float*)calloc_1d_array(ntype,sizeof(float));

    fread(typeinfo,sizeof(float),ntype,fid);
    fread(&time,sizeof(float),2,fid);
    fread(&n,sizeof(long),1,fid);

    data  = (float**)calloc_2d_array(n,7,sizeof(float));
    property = (int*)calloc_1d_array(n,sizeof(int));
    pid   = (long*)calloc_1d_array(n,sizeof(long));
    cpuid = (int*)calloc_1d_array(n,sizeof(int));
    order = (long*)calloc_1d_array(n,sizeof(long));

    /* read data */
    for (p=0; p<n; p++)
    {
      fread(data[p],sizeof(float),7,fid);
      fread(&(property[p]),sizeof(int),1,fid);
      fread(&(pid[p]),sizeof(long),1,fid);
      fread(&(cpuid[p]),sizeof(int),1,fid);
    }

    fclose(fid);

    /* Step 2: sort the particles */
    for (p=0; p<n; p++)
      order[p] = p;

    quicksort(cpuid, pid, order, 0, n-1);

    /* Step 3: write back the ordered particle list */
    fid = fopen(fname,"wb");

    /* write header */
    fwrite(header,sizeof(float),12,fid);
    fwrite(&ntype,sizeof(int),1,fid);
    fwrite(typeinfo,sizeof(float),ntype,fid);
    fwrite(time,sizeof(float),2,fid);
    fwrite(&n,sizeof(long),1,fid);

    /* write data */
    for (p=0; p<n; p++)
    {
      fwrite(data[order[p]],sizeof(float),7,fid);
      fwrite(&property[p],sizeof(int),1,fid);
      fwrite(&pid[order[p]],sizeof(long),1,fid);
      fwrite(&cpuid[order[p]],sizeof(int),1,fid);
    }

    free_2d_array(data);          free(property);
    free(pid);    free(cpuid);    free(order);

    fclose(fid);
  }

  if (typeinfo != NULL)
    free(typeinfo);

  return 0;
}


/* ========================================================================== */

static void quicksort(int *cpuid, long *pid, long *order, long ns, long ne)
{
  long i, a, b, ind, pivot;

  if (ne <= ns) return;

  /* location of the pivot at half chain length */
  pivot = (long)((ns+ne+1)/2);

  /* move the pivot to the start */
  ind = order[pivot];
  order[pivot] = order[ns];
  order[ns] = ind;

  /* initial configuration */
  pivot = ns;
  i = ns + 1;

  /* move the particles that are "smaller" than the pivot before it */
  while (i <= ne)
  {
    a = order[i];   b = order[pivot];

    if ((cpuid[a] < cpuid[b]) || ((cpuid[a] == cpuid[b]) && (pid[a] < pid[b])))
    {/* the ith particle is smaller, move it before the pivot */
      order[pivot] = a;
      pivot++;
      order[i] = order[pivot];
      order[pivot] = b;
    }
    i++;
  }

  /* recursively call this routine to complete sorting */
  quicksort(cpuid,pid,order,ns,pivot-1);
  quicksort(cpuid,pid, order,pivot+1,ne);

  return;
}



/* Write an error message and terminate the simulation with an error status. */
static void sort_error(const char *fmt, ...){
  va_list ap;

  va_start(ap, fmt);         /* ap starts after the fmt parameter */
  vfprintf(stderr, fmt, ap); /* print the error message to stderr */
  va_end(ap);                /* end stdargs (clean up the va_list ap) */

  fflush(stderr);            /* flush it NOW */
  exit(1);                   /* clean up and exit */
}

static void usage(const char *arg)
{
  fprintf(stderr,"\nUsage: %s [options] [block] ...\n", arg);
  fprintf(stderr,"\nOptions:\n");
  fprintf(stderr,"  -d <directory>  name of the input directory\n");
  fprintf(stderr,"                  Default: current directory\n");
  fprintf(stderr,"  -i <name>       basename of input file\n");
  fprintf(stderr,"  -s <name>       posterior name of input file\n");
  fprintf(stderr,"  -f f1:f2:fi     file number range and interval\n");
  fprintf(stderr,"                  Default: <0:0:1>\n");
  fprintf(stderr,"  -h              this help\n");

  fprintf(stderr,"\nExample:\n");
  fprintf(stderr,"%s -d mydir -i streaming2d -s ds -f 0:500\n\n", arg);

  exit(0);
}

static void* calloc_1d_array(size_t nc, size_t size)
{
  void *array;

  if ((array = (void *)calloc(nc,size)) == NULL) {
    sort_error("[calloc_1d] failed to allocate memory (%d of size %d)\n",
              (int)nc,(int)size);
    return NULL;
  }
  return array;
}

static void** calloc_2d_array(size_t nr, size_t nc, size_t size)
{
  void **array;
  size_t i;

  if((array = (void **)calloc(nr,sizeof(void*))) == NULL){
    sort_error("[calloc_2d] failed to allocate memory for %d pointers\n",(int)nr);
    return NULL;
  }

  if((array[0] = (void *)calloc(nr*nc,size)) == NULL){
    sort_error("[calloc_2d] failed to allocate memory (%d X %d of size %d)\n",
              (int)nr,(int)nc,(int)size);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<nr; i++){
    array[i] = (void *)((unsigned char *)array[0] + i*nc*size);
  }

  return array;
}

static void free_2d_array(void *array)
{
  void **ta = (void **)array;

  free(ta[0]);
  free(array);
}

