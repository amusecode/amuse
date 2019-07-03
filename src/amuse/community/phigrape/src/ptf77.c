/***********************************************************************
 * ptf77.c -- FORTRAN interface to some pt.c stuff
 *
 * Author: Mark Hays <hays@math.arizona.edu>
 */

#include "pt.h"

#include <stdio.h>

static pt_pipeline_t pipeline;
static int inuse=0;

/* See the bottom of pt.h for the definitions of
 * the F77PIPLINExxxx macros.
 */

void F77PIPELINEINIT (pt_startroutine_t setupproc,
		      pt_startroutine_t stageproc)
{
  if (inuse) {
    fprintf(stderr,"pipeline_init: already in use by FORTRAN\n");
    exit(1);
  }
  pt_pipeline_init(&pipeline,NULL,setupproc,stageproc);
  inuse=1;
}

void F77PIPELINEDONE ()
{
  if (!inuse) {
    fprintf(stderr,
	  "pipeline_done: must call pipeline_init() before pipeline_done()\n");
    exit(1);
  }
  pt_pipeline_destroy(&pipeline);
  inuse=0;
}

void F77PIPELINEEXEC ()
{
  if (!inuse) {
    fprintf(stderr,"pipeline_execute: must call pipeline_init() first\n");
    exit(1);
  }
  pt_pipeline_execute(&pipeline);
}

/* EOF ptf77.c */

