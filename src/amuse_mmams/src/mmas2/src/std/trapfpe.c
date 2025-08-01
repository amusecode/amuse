#ifdef _LINUX_

#include <fenv.h>
void trapfpe(void) {
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

#elif _MACOSX_PPC_

/*
 * fputest.c
 * test for SIGFPE delivery under Darwin PPC
 */


#include <architecture/ppc/fp_regs.h>
#include <mach/mach.h>
#include <pthread.h>

static void *fpu_fpe_enable(void *arg);
#define FE0_MASK (1<<11)
#define FE1_MASK (1<<8) 
/* FE0 FE1 exceptions enabled if either FE0 or FE1 set
 * 0 0 -- floating-point exceptions disabled
 * 0 1 -- floating-point imprecise nonrecoverable
 * 1 0 -- floating-point imprecise recoverable
 * 1 1 -- floating-point precise mode
 */

/* a thread cannot get or set its own MSR bits */
static void * fpu_fpe_enable(void *arg)
{
  thread_t t = *(thread_t *)arg;
  struct ppc_thread_state state;
  unsigned int state_size = PPC_THREAD_STATE_COUNT;
  if (thread_get_state(t, PPC_THREAD_STATE,
		       (natural_t *)&state, &state_size) == KERN_SUCCESS) {
    state.srr1 |= FE1_MASK;
    thread_set_state(t, PPC_THREAD_STATE, (natural_t *)&state, state_size);
  }
  return 0;
}

void trapfpe(void)
{
  static volatile int looping = 0;
  thread_t self = mach_thread_self();
  pthread_t enabler;
  ppc_fp_scr_t r = get_fp_scr();
  /* turn off exception bits to prevent immediate re-fault */
  r.fx = r.fex = r.vx = r.ox = r.ux = r.zx = r.xx = r.vx_snan = r.vx_isi =
    r.vx_idi = r.vx_zdz = r.vx_imz = r.vx_xvc = r.vx_cvi = r.vx_soft = 0;
  /* these only have to be set once, but may as well set anyway */
  r.ve = 1; /* invalid */
  r.oe = 1; /* overflow */
  r.ue = 0; /* underflow */
  r.ze = 1; /* zero divide */
  r.xe = 0; /* inexact */
  if (!looping) {
    looping++;
    set_fp_scr(r);
    /* MSR bit reset after every fault, set it again */
    if (!pthread_create(&enabler, 0, fpu_fpe_enable, &self))
      pthread_join(enabler, 0);
  }
  looping = 0;
}

#else

void trapfpe(void) {};

#endif



