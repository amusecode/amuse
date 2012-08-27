#include "copyright.h"
/*============================================================================*/
/*! \file ath_files.c
 *  \brief Function for creating descriptive output filenames.
 *
 * PURPOSE: Function for creating descriptive output filenames with form:
 *  -     [path]<basename>[-lev#][-dom#][.idump][.id].<ext>
 *
 *   where 
 *     -   path     = optional path
 *     -   basename = basename of file (usually problem name, e.g. "Sod")
 *     -   lev#     = level number of dump, only included when level > 0
 *     -   dom#     = domain number of dump, only included when domain != 0
 *     -   dlen     = number of digits to use for numeric extension (1..10)
 *     -   idump    = optional dump number (0,1,2,.....)
 *                    if(dlen > 0 and idump >= 0) we use the dump number
 *                    <idump> uses C-format descriptor "%0<dlen>d"
 *     -   id       = optional additional identifier set in <input> block
 *     -   ext      = file extension, e.g. ".tab", ".bin", ".dx", ".vtk"
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ath_fname()							      */
/*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn char *ath_fname(const char *path, const char *basename,
 *                      const char *levstr, const char *domstr, 
 *                      const int dlen, const int idump,
 *                      const char *id, const char *ext)
 *  \brief Creates descriptive output filenames.
 *
 *   Creates filenames with form 
 *    - [path]<basename>[-lev#][-dom#][.idump][.id].<ext>
 *
 *   Used by most of the dump_*() and output_*() functions.
 */
char *ath_fname(const char *path, const char *basename,
                const char *levstr, const char *domstr,
                const int dlen, const int idump,
                const char *id, const char *ext)
{
  char fmt[80];
  size_t size,slen=0;
  char *cp, *fname;

/* 2 = "." following the basename + NULL terminator */
  size = 2 + strlen(basename) + strlen(ext);
  if(path != NULL) size += 1 + strlen(path);        /* add 1 for possible '/' */
  if(levstr != NULL) size += 1 + strlen(levstr);    /* add 1 for the "-" */
  if(domstr != NULL) size += 1 + strlen(domstr);    /* add 1 for the "-" */
  if(dlen > 0) size += (dlen > 10 ? dlen : 10) + 1; /* add 1 for the "." */
  if(id != NULL) size += 1 + strlen(id);            /* add 1 for the "." */

  if((fname = (char*)malloc(size*sizeof(char))) == NULL){
    printf("[ath_fname]: malloc returned a NULL pointer\n");
    return NULL;
  }

/* Build the filename. Start with the optional path */

  cp = fname;
  if(path != NULL){
    slen = strlen(path);
    strcpy(cp, path);
    cp = &(cp[slen]); /* point cp at the '\0' terminator */

    /* Append a '/' if necessary */
    if(cp[-1] != '/'){
      cp[0] = '/';
      cp[1] = '\0';
      cp = &(cp[1]);
    }
  }

/* Append the basename of the file */

  slen = strlen(basename);
  strcpy(cp,basename);
  cp = &(cp[slen]); /* point cp at the '\0' terminator */

/* Append the optional level number, preceded by '-' */

  if (levstr != NULL){
    slen = strlen(levstr);
    cp[0] = '-';
    cp[1] = '\0';
    cp = &(cp[1]);
    strcpy(cp,levstr);
    cp = &(cp[slen]); /* point cp at the '\0' terminator */
  }

/* Append the optional domain number, preceded by '-' */

  if (domstr != NULL){
    slen = strlen(domstr);
    cp[0] = '-';
    cp[1] = '\0';
    cp = &(cp[1]);
    strcpy(cp,domstr);
    cp = &(cp[slen]); /* point cp at the '\0' terminator */
  }

/* Append the optional dump number, id string, and extension */

/* id not NULL */
  if(id){
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,".%%0%dd.%%s.%%s",dlen);
      sprintf(cp,fmt,idump,id,ext);
    }
    else sprintf(cp,".%s.%s",id,ext);
  }
/* id NULL */
  else{
    if(dlen > 0 && idump >= 0){
      sprintf(fmt,".%%0%dd.%%s",dlen);
      sprintf(cp,fmt,idump,ext);
    }
    else sprintf(cp,".%s",ext);
  }

  return fname;
}
