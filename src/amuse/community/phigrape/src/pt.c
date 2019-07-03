/***********************************************************************
 * pt.c -- thread utility routines
 *
 * Author: Mark Hays <hays@math.arizona.edu>
 */

#include "pt.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>	/* for EBUSY */

/*************************************************
 * attempt to lock a mutex
 */
int pt_mutex_trylock(pthread_mutex_t *m,char *msg)
{
  int res;

  /* returns EBUSY if mutex is already locked,
   * and EINVAL if the ptr is bad (on RedHat5.0)
   *
   * might this return EAGAIN on some systems??
   * i can't find any docs on this one's retval!
   *
   */
  if ((res=pthread_mutex_trylock(m)) != EBUSY) {
    fprintf(stderr,"%s: %d: %s\n",msg,res,strerror(res));
    exit(1);
  }
  return(res ? 1 : 0);
}

/*************************************************
 * run nthreads threads in the routine start
 */
void _pt_fork(int nthreads,
	      pt_startroutine_t start,
	      pt_addr_t arg,
	      pt_addr_t *exitcodes)
{
  int i;
  pt_arg_t *args;

  if (nthreads<1) {
    fprintf(stderr,"pt_fork: nthreads<1\n"); exit(1);
  }
  if ((args=(pt_arg_t *) malloc(nthreads*sizeof(pt_arg_t)))==NULL) {
    fprintf(stderr,"pt_fork: malloc failed!\n"); exit(1);
  }
  for (i=0; i<nthreads; i++) {
    args[i].nthreads=nthreads; args[i].myid=i; args[i].data=arg;
  }
  for (i=0; i<nthreads; i++) {
    pt_create(&args[i].self,start,args+i,"pt_fork: pt_create");
  }
  for (i=0; i<nthreads; i++) {
    pt_wait(&args[i].self,
	    exitcodes==NULL?NULL:exitcodes+i,
	    "pt_fork: pt_wait");
  }
  free(args);
}

/*************************************************
 * initialize a gate
 */
void pt_gate_init(pt_gate_t *gate,int nthreads)
{
  gate->ngate=0; gate->nthreads=nthreads;
  pt_mutex_init(  &gate->mutex, "gate: init mutex");
  pt_mutex_init(  &gate->block, "gate: init block");
  pt_cond_init (&gate->condvar, "gate: init condvar");
  pt_cond_init (   &gate->last, "gate: init last");
}

/*************************************************
 * destroy a gate variable
 */
void pt_gate_destroy(pt_gate_t *gate)
{
  gate->ngate=gate->nthreads=0;
  pt_mutex_destroy(  &gate->mutex, "gate: destroy mutex");
  pt_mutex_destroy(  &gate->block, "gate: destroy block");
  pt_cond_destroy (&gate->condvar, "gate: destroy condvar");
  pt_cond_destroy (   &gate->last, "gate: destroy last");
}

/*************************************************
 * enter the gate
 */
void pt_gate_sync(pt_gate_t *gate)
{
  if (gate->nthreads<2) return;           /* trivial case            */
  pt_mutex_lock(&gate->block,             /* lock the block -- new   */
		"gate: lock block");            /*   threads sleep here    */
  pt_mutex_lock(&gate->mutex,             /* lock the mutex          */
		"gate: lock mutex");
  if (++(gate->ngate) < gate->nthreads) { /* are we the last one in? */
    pt_mutex_unlock(&gate->block,         /* no, unlock block and    */
		    "gate: unlock block 1");
    pt_cond_wait(&gate->condvar,          /*   go to sleep           */
		 &gate->mutex,
		 "gate: wait condvar");
  } else {                                /* yes, we're last         */
    pt_cond_broadcast(&gate->condvar,     /* wake everyone up and    */
		      "gate: bcast condvar");
    pt_cond_wait(&gate->last,&gate->mutex,/* go to sleep til they're */
		 "gate: wait last");            /* all awake... then       */
    pt_mutex_unlock(&gate->block,         /* release the block       */
		    "gate: unlock block 2");
  }
  if (--(gate->ngate)==1) {               /* next to last one out?   */
    pt_cond_broadcast(&gate->last,        /* yes, wake up last one   */
		      "gate: bcast last");
  }
  pt_mutex_unlock(&gate->mutex,           /* release the mutex       */
		  "gate: unlock mutex");
}

/*************************************************
 * Pipeline stage: the idea:
 *
 *   main thread   I/O thread
 *           \      /      \
 *            \    /        \
 *            gate1         |
 *           /     \        |
 *          /       \       |
 *       setup       |     work
 *          \       /       |
 *           \     /        |
 *            gate2         |
 *           /    \         /
 *          /      \_______/
 *         |
 *   main continues
 */

/*************************************************
 * couple of convenient macros
 */
#define GATE1(pipeline) pt_gate_sync(&((pipeline)->gate1))
#define GATE2(pipeline) pt_gate_sync(&((pipeline)->gate2))
#define STAGE(pipeline) (*((pipeline)->stageproc))((pipeline)->gdata)
#define SETUP(pipeline) \
  { pt_startroutine_t fp; \
    \
    if ((fp=(pipeline)->setupproc)!=NULL) (*fp)(pipeline->gdata); \
  }

/*************************************************
 * slave thread executes this
 */
static void _pt_pipeline_slave_code(pt_pipeline_t *pipeline)
{
  while (1) {
    GATE1(pipeline);
    if (pipeline->terminate) break;
    GATE2(pipeline);
    STAGE(pipeline);
  }
  pt_exit(NULL);
}

/*************************************************
 * init the info struct and start up the slave
 */
void _pt_pipeline_init(pt_pipeline_t *pipeline,
		       pt_addr_t gdata,
		       pt_startroutine_t setupproc,
		       pt_startroutine_t stageproc)
{
  pt_gate_init(&(pipeline->gate1),2);
  pt_gate_init(&(pipeline->gate2),2);
  pipeline->terminate=0;
  pipeline->gdata=gdata;
  pipeline->setupproc=setupproc;
  pipeline->stageproc=stageproc;
  pt_create(&(pipeline->slave),_pt_pipeline_slave_code,pipeline,
	    "pt_pipeline_init: create slave");
}

/*************************************************
 * kill the slave, free resources
 */
void pt_pipeline_destroy(pt_pipeline_t *pipeline)
{
  pipeline->terminate=1;
  GATE1(pipeline);
  pt_wait(&(pipeline->slave),NULL,"pt_pipeline_destroy: wait");
  pt_gate_destroy(&(pipeline->gate1));
  pt_gate_destroy(&(pipeline->gate2));
  pipeline->gdata=NULL;
  pipeline->setupproc=NULL;
  pipeline->stageproc=NULL;
}

/*************************************************
 * run the pipeline stage
 */
void pt_pipeline_execute(pt_pipeline_t *pipeline)
{
  GATE1(pipeline);
  SETUP(pipeline);
  GATE2(pipeline);
}

/* EOF pt.c */

