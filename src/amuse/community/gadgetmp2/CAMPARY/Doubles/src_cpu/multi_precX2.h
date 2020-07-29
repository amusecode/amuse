
#ifndef __multi_precX_h__
#define __multi_precX_h__

#include "multi_prec.h"
//#include "gpu_mprec.h"

#define XprecX 2

typedef multi_prec<XprecX> dd_real;

static const dd_real _nan;
static const dd_real _inf;

// 2 * pi
const double myArray_2pi[XprecX] = {0x6.487ed5110b46p0, 0x1.1a62633145c07p-52};
const dd_real _2pi = dd_real(myArray_2pi, XprecX);

// pi
const double myArray_pi[XprecX] = {0x3.243f6a8885a3p0, 0x8.d313198a2e038p-56};
const dd_real _pi = dd_real(myArray_pi, XprecX);

// pi / 2
const double myArray_pi2[XprecX] = {0x1.921fb54442d18p0, 0x4.69898cc51701cp-56};
const dd_real _pi2 = dd_real(myArray_pi2, XprecX);

// pi / 4
const double myArray_pi4[XprecX] = {0xc.90fdaa22168cp-4, 0x2.34c4c6628b80ep-56};
const dd_real _pi4 = dd_real(myArray_pi4, XprecX);

// 3 * pi / 4
const double myArray_3pi4[XprecX] = {0x2.5b2f8fe6643a4p0, 0x6.9e4e5327a2828p-56};
const dd_real _3pi4 = dd_real(myArray_3pi4, XprecX);

// Euler's number
const double myArray_e[XprecX] = {0x2.b7e151628aed2p0, 0xa.6abf7158809dp-56};
const dd_real _e = dd_real(myArray_e, XprecX);

// ln(2)
const double myArray_log2[XprecX] = {0xb.17217f7d1cf78p-4, 0x1.abc9e3b39803fp-56};
const dd_real _log2 = dd_real(myArray_log2, XprecX);

// ln(10)
const double myArray_log10[XprecX] = {0x2.4d763776aaa2cp0, -0xf.a456a4a751f48p-56};
const dd_real _log10 = dd_real(myArray_log10, XprecX);

const double _eps = 0x1.p-104;
const double _min_normalized = 0x1.p-969;

const double myArray_max[XprecX] = {0xf.ffffffffffff8p1020, 0x3.ffffffffffffep968};
const dd_real _max = dd_real(myArray_max, XprecX);
const double myArray_safe_max[XprecX] = {0xf.ffffffffffff8p1020, 0x3.ffffffffffffep968};
const dd_real _safe_max = dd_real(myArray_safe_max, XprecX);
    
const int _ndigits = 31;

#endif
