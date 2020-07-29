
#ifndef __multi_precX_h__
#define __multi_precX_h__

#include "multi_prec.h"
//#include "gpu_mprec.h"

#define XprecX 8

typedef multi_prec<XprecX> dd_real;

static const dd_real _nan;
static const dd_real _inf;

// 2 * pi
const double myArray_2pi[XprecX] = {0x6.487ed5110b46p0, 0x1.1a62633145c07p-52, 
																		-0x1.f1976b7ed8fbcp-108, 0x5.33e63a0105df4p-164, 
																		0x1.31d89cd9128a5p-216, 0x4.3cc71a026ef7cp-276, 
																		0xa.8cd9e69d218d8p-332, 0x1.8158536f92f8ap-384};
const dd_real _2pi = dd_real(myArray_2pi, XprecX);

// pi
const double myArray_pi[XprecX] = {0x3.243f6a8885a3p0, 0x8.d313198a2e038p-56, 
																	 -0xf.8cbb5bf6c7dep-112, 0x2.99f31d0082efap-164, 
																	 0x9.8ec4e6c894528p-220, 0x2.1e638d01377bep-276, 
																	 0x5.466cf34e90c6cp-332, 0xc.0ac29b7c97c5p-388};
const dd_real _pi = dd_real(myArray_pi, XprecX);

// pi / 2
const double myArray_pi2[XprecX] = {0x1.921fb54442d18p0, 0x4.69898cc51701cp-56, 
																		-0x7.c65dadfb63efp-112, 0x1.4cf98e804177dp-164, 
																		0x4.c76273644a294p-220, 0x1.0f31c6809bbdfp-276, 
																		0x2.a33679a748636p-332, 0x6.05614dbe4be28p-388};
const dd_real _pi2 = dd_real(myArray_pi2, XprecX);

// pi / 4
const double myArray_pi4[XprecX] = {0xc.90fdaa22168cp-4, 0x2.34c4c6628b80ep-56, 
																		-0x3.e32ed6fdb1f78p-112, 0xa.67cc74020bbe8p-168, 
																		0x2.63b139b22514ap-220, 0x8.798e3404ddef8p-280, 
																		0x1.519b3cd3a431bp-332, 0x3.02b0a6df25f14p-388};
const dd_real _pi4 = dd_real(myArray_pi4, XprecX);

// 3 * pi / 4
const double myArray_3pi4[XprecX] = {0x2.5b2f8fe6643a4p0, 0x6.9e4e5327a2828p-56, 
																		 0x1.456737b06ea1ap-108, -0x6.0c89aa3f9dcc4p-164, 
																		 -0xd.4ec52e990c22p-224, 0x1.96caa9c0e99cfp-276, 
																		 -0x4.0b2e4985136bp-332, 0x1.90811f49d71d4p-384};
const dd_real _3pi4 = dd_real(myArray_3pi4, XprecX);

// Euler's number
const double myArray_e[XprecX] = {0x2.b7e151628aed2p0, 0xa.6abf7158809dp-56, 
																	-0xb.0c389d18e9f1p-112, 0x3.8b4da56a784dap-164, 
																	-0xf.bae6f3010cdbp-220, -0x1.88c76d93041a1p-272, 
																	0x4.bf8d8d8c31d78p-328, -0x1.c25f937f544eep-380};
const dd_real _e = dd_real(myArray_e, XprecX);

// ln(2)
const double myArray_log2[XprecX] = {0xb.17217f7d1cf78p-4, 0x1.abc9e3b39803fp-56, 
																		 0x2.f6af40f343268p-112, -0xd.6749d275f2e88p-168, 
																		 -0x2.4745505d41848p-220, 0x7.6206debac9854p-276, 
																		 0x1.9552fb4afa1b1p-328, 0xe.d2eae35c1382p-388};
const dd_real _log2 = dd_real(myArray_log2, XprecX);

// ln(10)
const double myArray_log10[XprecX] = {0x2.4d763776aaa2cp0, -0xf.a456a4a751f48p-56, 
																			-0x3.3d75c75c04c18p-108, -0x9.6881bc5f0e788p-164, 
																			0x2.0807c0b5ca58cp-216, -0x3.f4a1395fbe8ccp-272, 
																			-0xe.3cd0ff4e83ca8p-328, 0x2.0b1889061043p-380};
const dd_real _log10 = dd_real(myArray_log10, XprecX);

const double _eps = 0x1.p-416;
const double _min_normalized = 0x1.p-651;

const double myArray_max[XprecX] = {0xf.ffffffffffff8p1020, 0x3.ffffffffffffep968, 
																		0xf.ffffffffffff8p912, 0x3.ffffffffffffep860, 
																		0xf.ffffffffffff8p804, 0x3.ffffffffffffep752, 
																		0xf.ffffffffffff8p696, 0x3.ffffffffffffep644};
const dd_real _max = dd_real(myArray_max, XprecX);
const double myArray_safe_max[XprecX] = {0xf.ffffffffffff8p1020, 0x3.ffffffffffffep968, 
																		     0xf.ffffffffffff8p912, 0x3.ffffffffffffep860, 
																		     0xf.ffffffffffff8p804, 0x3.ffffffffffffep752, 
																		     0xf.ffffffffffff8p696, 0x3.ffffffffffffep644};
const dd_real _safe_max = dd_real(myArray_safe_max, XprecX);
    
const int _ndigits = 124;

#endif
