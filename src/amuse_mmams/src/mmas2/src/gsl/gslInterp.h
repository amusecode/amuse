#ifndef __gslInterp_
#define __gslInterp__

#include "include/stdinc.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <vector>

class gslInterp {
protected:
  gsl_interp_accel *acc;
  gsl_interp       *interp;
  vector<double>   x, y;
public:
  
  gslInterp(vector<double> &x0, vector<double> &y0) {
    acc   = gsl_interp_accel_alloc();
    for (size_t i = 0; i < x0.size(); i++) {
      x.push_back(x0[i]);
      y.push_back(y0[i]);
    }
    interp = gsl_interp_alloc(gsl_interp_linear, x.size());
    gsl_interp_init(interp, &x[0], &y[0], x.size());
  }
  ~gslInterp() {
    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);     
  }
  void init(vector<double> &x0, vector<double> &y0) {
    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);     
    acc    = gsl_interp_accel_alloc();
    x.clear();
    y.clear();
    for (size_t i = 0; i < x0.size(); i++) {
      x.push_back(x0[i]);
      y.push_back(y0[i]);
    }
    interp = gsl_interp_alloc(gsl_interp_linear, x.size());
    gsl_interp_init(interp, &x[0], &y[0], x.size());
  }

  inline double eval(double xc) {
    return gsl_interp_eval(interp, &x[0], &y[0], xc, acc);
  }
  inline double deriv(double xc) { 
    return gsl_interp_eval_deriv(interp, &x[0], &y[0], xc, acc);
  }
  inline double deriv2(double xc) { 
    return gsl_interp_eval_deriv2(interp, &x[0], &y[0], xc, acc);
  }
  inline double integ(double x_lo, double x_up) {
    return gsl_interp_eval_integ(interp, &x[0], &y[0], x_lo, x_up, acc);
  }
};

#endif // __gslInterp__
