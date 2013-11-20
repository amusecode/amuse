/*
  Copyright (c) 2006-2011 IDRIS/CNRS
  Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
  Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.
*/

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/statvfs.h> 

#define BUFSIZE 8388608 /* 8 Mo */
#define MAXLINE 256

#ifndef NOUNDERSCORE
#define title title_
#define transferfile transferfile_
#define gethname gethname_
#define getfreespace getfreespace_
#define createdirs createdirs_
#define makedir makedir_
#endif

extern int errno;

extern void title(const int *,char *);

#ifdef IOPROC
void transferfile(const char *src,const char *dest,int *size)
{
  int fi,fo;
  int noct;
  char *buff;

  buff = malloc ((size_t) BUFSIZE);
  if(buff==NULL) {
    perror("buffer allocation failed:");
    return;
  }

  *size=-1;
  if((fi=open(src,O_RDONLY)) == -1) {
    perror(src);
    free(buff);
    return; /* Dangerous? */
  }
  if((fo=creat(dest,0644)) == -1) {
    perror(dest);
    close(fi);
    free(buff);
    return; /* Dangerous? */
  }

  *size=0;
  while(noct=read(fi,buff,sizeof(buff))){
    *size+=noct;
    if(write(fo,buff,noct)==-1){
      perror(dest);
      *size=-1;
      break;
    }
  }

  close(fi);close(fo);

  free(buff);

  /* Do not delete scratch file if error */
  if(*size==-1) return;

  /* Delete src file */
  if(unlink(src) == -1) perror(src);
}


void createdirs(const char *scratch,const char *perm,
               const int *my_iogroup,const int *ncpu_iogroup,const int *first_id)
{
  int i;
  char basename[MAXLINE],dirname[MAXLINE];
  char iochar[6],cpuchar[6];

  /* This is necessary because Fortran does not terminate strings */
  iochar [5] = '\0';
  cpuchar[5] = '\0';

  title(my_iogroup,iochar);

  /* Scratch */
  strcpy(basename,scratch);
  strcat(basename,"ionode_");
  strcat(basename,iochar);
  if(mkdir(basename,0750)==-1 && errno!=EEXIST) perror(basename);
  for(i=*first_id;i<*first_id+*ncpu_iogroup-1;i++){
    strcpy(dirname,basename);
    strcat(dirname,"/process_");
    title(&i,cpuchar);
    strcat(dirname,cpuchar);
    if(mkdir(dirname,0750)==-1 && errno!=EEXIST) perror(dirname);
  }

  /* Perm */
  if(strcmp(scratch,perm)!=0){
    strcpy(basename,perm);
    strcat(basename,"ionode_");
    strcat(basename,iochar);
    if(mkdir(basename,0750)==-1 && errno!=EEXIST) perror(basename);
    for(i=*first_id;i<*first_id+*ncpu_iogroup-1;i++){
      strcpy(dirname,basename);
      strcat(dirname,"/process_");
      title(&i,cpuchar);
      strcat(dirname,cpuchar);
      if(mkdir(dirname,0750)==-1 && errno!=EEXIST) perror(dirname);
    }
  }
}

void getfreespace(const char *dir,double *freespace)
{
  struct statvfs stats;

  if (statvfs(dir,&stats)!=0) {
    perror(dir);
    *freespace=-1.;
    return;
  }

  *freespace  = stats.f_bavail; /* = free blocks */
  *freespace *= stats.f_frsize; /* *= fragment size (the block size does not give the right answer!) */
}

void gethname(char *host)
{
  size_t len;

  if(gethostname(host,MAXLINE)!=0) perror("Gethostname error:");

  len = strlen(host);

  host[len]=' '; /* Fortran does not manage well the \0 character */
}
#endif

void makedir(const char *dir)
{
  if(mkdir(dir,0750)==-1 && errno!=EEXIST) perror(dir);
}
