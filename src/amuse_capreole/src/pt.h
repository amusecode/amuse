/***********************************************************************
 * pt.h -- pthreads utility macros
 *
 * Author: Mark Hays <hays@math.arizona.edu>
 */

#ifndef _pt_h_
#define _pt_h_

/* Linux defs:
 *   _REENTRANT to get thread-safe libs
 *   _POSIX_SOURCE to get POSIX semantics
 *   _P is a hack for LinuxThreads -- on my box,
 *      pthread.h includes sched.h. My sched.h
 *      (incorrectly) declares prototypes with
 *      _P instead of __P (which is what everything
 *      else uses... Maybe it's just me.
 */
#ifdef __linux__
#  define _REENTRANT
#  define _POSIX_SOURCE
/* #  define _P __P */
#endif

#include <pthread.h>

#include <string.h>	/* for strerror() */

typedef void *(*pt_startroutine_t)(void *);
typedef void *pt_addr_t;

/*************************************************
 * low level wrappers that die on errors
 */
#define pt_create(t,start,arg,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_create(t, \
			       NULL, \
			       (pt_startroutine_t) (start), \
			       (pt_addr_t) (arg))) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_create_detached(start,arg,msg) \
{ \
    pthread_t t; \
    int errcode; \
    \
    if (errcode=pthread_create(&t, \
			       NULL, \
			       (pt_startroutine_t) (start), \
			       (pt_addr_t) (arg))) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
    if (pthread_detach(t)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_wait(t,exitcode,msg) \
{ \
    pt_addr_t code; \
    int errcode; \
    \
    if (errcode=pthread_join(*(t), \
		      (pt_addr_t) ((exitcode)==NULL ? &code : (exitcode)))) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_exit(status) \
{ \
    pthread_exit(status); \
}

#define pt_mutex_init(m,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_mutex_init(m,NULL)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_mutex_destroy(m,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_mutex_destroy(m)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_cond_init(c,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_cond_init(c,NULL)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_cond_destroy(c,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_cond_destroy(c)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_mutex_lock(m,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_mutex_lock(m)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

/* This one has to do some extra checking so it
 * isn't a macro...
 */
extern int pt_mutex_trylock(pthread_mutex_t *m, char *msg);

#define pt_mutex_unlock(m,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_mutex_unlock(m)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_cond_wait(c,m,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_cond_wait(c,m)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_cond_broadcast(c,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_cond_broadcast(c)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_cond_signal(c,msg) \
{ \
    int errcode; \
    \
    if (errcode=pthread_cond_signal(c)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

/*************************************************
 * N threads simultaneously doing the same thing
 */

typedef struct _pt_arg_t_ {
  int myid;
  int nthreads;
  pthread_t self;
  pt_addr_t data;
} pt_arg_t;

#define pt_myid(th_arg)      ((th_arg)->myid)
#define pt_nthreads(th_arg)  ((th_arg)->nthreads)
#define pt_data(th_arg)      ((th_arg)->data)
#define pt_self(th_arg)      (&((th_arg)->self))
#define pt_thread(th_arg,id) (&(((th_arg)-((th_arg)->myid)+(id))->self))

#define pt_cancel(th_arg,id) \
{ \
    int errcode; \
    \
    if (errcode=pthread_cancel(((th_arg)-((th_arg)->myid)+(id))->self)) { \
        fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
        exit(1); \
    } \
}

#define pt_cancel_all(th_arg) \
{ \
    int myid=(th_arg)->myid,nt=(th_arg)->nthreads,i,errcode; \
    pt_arg_t *base=(th_arg)-myid; \
    \
    for (i=0; i<nt; i++) { \
        if (i==id) continue; \
        if (errcode=pthread_cancel(base[i].self)) { \
            fprintf(stderr,"%s: %s\n",msg,strerror(errcode)); \
            exit(1); \
        } \
    } \
}

extern void _pt_fork(int nthreads,
		     pt_startroutine_t start,
		     pt_addr_t arg,
		     pt_addr_t *exitcodes);

#define pt_fork(nt,start,arg,codes) \
  _pt_fork(nt,(pt_startroutine_t) start, \
	   (pt_addr_t) arg,(pt_addr_t *) codes)

/*************************************************
 * the gate struct (rendezvous point)
 */
typedef struct _pt_gate_t_ {
  int ngate;              
  int nthreads;           
  pthread_mutex_t mutex;  
  pthread_mutex_t block;  
  pthread_cond_t condvar; 
  pthread_cond_t last;    
} pt_gate_t;

extern void pt_gate_init(pt_gate_t *gate,int nthreads);
extern void pt_gate_destroy(pt_gate_t *gate);
extern void pt_gate_sync(pt_gate_t *gate);

/*************************************************
 * the pipeline struct (a single stage)
 */
typedef struct _pt_pipeline_t_ {
  pt_addr_t gdata;
  pt_startroutine_t setupproc;
  pt_startroutine_t stageproc;
  int terminate;
  pthread_t slave;
  pt_gate_t gate1;
  pt_gate_t gate2;
} pt_pipeline_t;

extern void _pt_pipeline_init(pt_pipeline_t *p,
			      pt_addr_t gdata,
			      pt_startroutine_t setup,
			      pt_startroutine_t stage);

#define pt_pipeline_init(p,gdata,setup,stage)  \
  _pt_pipeline_init(p,                         \
		    (pt_addr_t) gdata,         \
		    (pt_startroutine_t) setup, \
		    (pt_startroutine_t) stage)

extern void pt_pipeline_destroy(pt_pipeline_t *p);
extern void pt_pipeline_execute(pt_pipeline_t *p);

/* names of C functions that can be called from
 * FORTRAN without the trailing underscore(s)
 */

#ifdef __linux__
   /* This works on my LinuxThreads + g77-0.5.21 */
#  define F77PIPELINEINIT pipeline_init__
#  define F77PIPELINEDONE pipeline_done__
#  define F77PIPELINEEXEC pipeline_execute__
#else
   /* This seems to work everywhere else */
#  define F77PIPELINEINIT pipeline_init_
#  define F77PIPELINEDONE pipeline_done_
#  define F77PIPELINEEXEC pipeline_execute_
#endif

#endif /* _pt_h_ */

/* EOF pt.h */

