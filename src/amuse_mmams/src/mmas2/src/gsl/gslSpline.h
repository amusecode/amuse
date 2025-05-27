#ifndef __gslSpline_
#define __gslSpline__

#include "include/stdinc.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

class gslSpline {
protected:
  const gsl_interp_type *spline_type;
  gsl_interp_accel *acc;
  gsl_spline       *spl;
public:
  
  gslSpline(vector<double> &x, vector<double> &y) {
    spline_type = gsl_interp_cspline;
    spl = gsl_spline_alloc(spline_type, x.size());
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(spl, &x[0], &y[0], x.size());
  }
  ~gslSpline() {
    gsl_interp_accel_free(acc);
    gsl_spline_free(spl);
  }
  void init(vector<double> &x, vector<double> &y) {
    gsl_interp_accel_free(acc);
    gsl_spline_free(spl);     
    spl = gsl_spline_alloc(spline_type, x.size());
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(spl, &x[0], &y[0], x.size());
  }

  inline real eval  (double xc)  {return gsl_spline_eval(spl, xc, acc);}
  inline real deriv (double xc)  {return gsl_spline_eval_deriv(spl, xc, acc);}
  inline real deriv2(double xc)  {return gsl_spline_eval_deriv2(spl, xc, acc);}
  inline real integ (double x_lo, double x_up)  {
    return gsl_spline_eval_integ(spl, x_lo, x_up, acc);
  }
};

#endif // __gslSpline__
