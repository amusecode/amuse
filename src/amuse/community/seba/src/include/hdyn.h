//DUMMY HDYN

#ifndef  STARLAB_HDYN_H
#  define  STARLAB_HDYN_H

#include  "dyn.h"

//-----------------------------------------------------------------------------
//  hdyn  --  a derived class of dynamical particles, with enough information
//             to integrate the equations of motions with a 4th-order Hermite
//             polynomial integrator.
//-----------------------------------------------------------------------------

class  hdyn : public dyn {
    protected:

    real time;

    public:

        hdyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase)
	    : dyn(the_hbpfp, the_sbpfp) {
	  time = 0;
	}

	virtual ~hdyn() {
	}

	inline real get_time()		const	{return time;}
  };
#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/hdyn.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
