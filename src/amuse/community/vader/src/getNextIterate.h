#ifndef _getnextiterate_h_
#define _getnextiterate_h_

#include "vader_common.h"

void getNextIterate(const grid *grd, const bool eos_func, 
		    const wksp *w
#if AA_M > 0 
		    , const unsigned long itCount, 
		    const long nHist, double *residMax
#endif
		    );

#endif
