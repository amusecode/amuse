
#ifndef __multi_precX_h__
#define __multi_precX_h__

#include "multi_prec.h"
//#include "gpu_mprec.h"

#define XprecX 4

typedef multi_prec<XprecX> dd_real;

static const dd_real _nan;
static const dd_real _inf;

// 2 * pi
const double myArray_2pi[XprecX] = {0x6.487ed5110b46p0, 0x1.1a62633145c07p-52, 
																		-0x1.f1976b7ed8fbcp-108, 0x5.33e63a0105df4p-164};
const dd_real _2pi = dd_real(myArray_2pi, XprecX);

// pi
const double myArray_pi[XprecX] = {0x3.243f6a8885a3p0, 0x8.d313198a2e038p-56, 
																	 -0xf.8cbb5bf6c7dep-112, 0x2.99f31d0082efap-164};
const dd_real _pi = dd_real(myArray_pi, XprecX);

// pi / 2
const double myArray_pi2[XprecX] = {0x1.921fb54442d18p0, 0x4.69898cc51701cp-56, 
																		-0x7.c65dadfb63efp-112, 0x1.4cf98e804177dp-164};
const dd_real _pi2 = dd_real(myArray_pi2, XprecX);

// pi / 4
const double myArray_pi4[XprecX] = {0xc.90fdaa22168cp-4, 0x2.34c4c6628b80ep-56, 
																		-0x3.e32ed6fdb1f78p-112, 0xa.67cc74020bbe8p-168};  
const dd_real _pi4 = dd_real(myArray_pi4, XprecX);

// 3 * pi / 4
const double myArray_3pi4[XprecX] = {0x2.5b2f8fe6643a4p0, 0x6.9e4e5327a2828p-56, 
																			0x1.456737b06ea1ap-108, -0x6.0c89aa3f9dcc4p-164};
const dd_real _3pi4 = dd_real(myArray_3pi4, XprecX);

// Euler's number
const double myArray_e[XprecX] = {0x2.b7e151628aed2p0, 0xa.6abf7158809dp-56, 
																	-0xb.0c389d18e9f1p-112, 0x3.8b4da56a784dap-164};
const dd_real _e = dd_real(myArray_e, XprecX);

// ln(2)
const double myArray_log2[XprecX] = {0xb.17217f7d1cf78p-4, 0x1.abc9e3b39803fp-56, 
																		 0x2.f6af40f343268p-112, -0xd.6749d275f2e88p-168};
const dd_real _log2 = dd_real(myArray_log2, XprecX);

// ln(10)
const double myArray_log10[XprecX] = {0x2.4d763776aaa2cp0, -0xf.a456a4a751f48p-56, 
																			-0x3.3d75c75c04c18p-108, -0x9.6881bc5f0e788p-164};
const dd_real _log10 = dd_real(myArray_log10, XprecX);

const double _eps = 0x1.p-208;
const double _min_normalized = 0x1.p-863;

const double myArray_max[XprecX] = {0xf.ffffffffffff8p1020, 0x3.ffffffffffffep968, 
                                    0xf.ffffffffffff8p912, 0x3.ffffffffffffep860};
const dd_real _max = dd_real(myArray_max, XprecX);

const double myArray_safe_max[XprecX] = {0xf.ffffffffffff8p1020, 0x3.ffffffffffffep968, 
                                         0xf.ffffffffffff8p912, 0x3.ffffffffffffep860};
const dd_real _safe_max = dd_real(myArray_safe_max, XprecX);

const int _ndigits = 62;

/*const double myArray_2pi[XprecX] = {6.283185307179586232e+00, 2.449293598294706414e-16,
                                    -5.989539619436679332e-33, 2.224908441726730563e-49};
const dd_real _2pi = dd_real(myArray_2pi, XprecX);
const double myArray_pi[XprecX] = {3.141592653589793116e+00, 1.224646799147353207e-16,
                                   -2.994769809718339666e-33, 1.112454220863365282e-49};
const dd_real _pi = dd_real(myArray_pi, XprecX);
const double myArray_pi2[XprecX] = {1.570796326794896558e+00, 6.123233995736766036e-17,
                                    -1.497384904859169833e-33, 5.562271104316826408e-50};
const dd_real _pi2 = dd_real(myArray_pi2, XprecX);
const double myArray_pi4[XprecX] = {7.853981633974482790e-01, 3.061616997868383018e-17,
                                    -7.486924524295849165e-34, 2.781135552158413204e-50};  
const dd_real _pi4 = dd_real(myArray_pi4, XprecX);
const double myArray_3pi4[XprecX] = {2.356194490192344837e+00, 9.1848509936051484375e-17,
                                     3.9168984647504003225e-33, -2.5867981632704860386e-49};
const dd_real _3pi4 = dd_real(myArray_3pi4, XprecX);
const double myArray_e[XprecX] = {2.718281828459045091e+00, 1.445646891729250158e-16,
                                  -2.127717108038176765e-33, 1.515630159841218954e-49};
const dd_real _e = dd_real(myArray_e, XprecX);

const double myArray_log2[XprecX] = {6.931471805599452862e-01, 2.319046813846299558e-17,
                                     5.707708438416212066e-34, -3.582432210601811423e-50};
const dd_real _log2 = dd_real(myArray_log2, XprecX);
const double myArray_log10[XprecX] = {2.302585092994045901e+00, -2.170756223382249351e-16,
                                      -9.984262454465776570e-33, -4.023357454450206379e-49};
const dd_real _log10 = dd_real(myArray_log10, XprecX);

const double _eps = 1.21543267145725e-63;
const double _min_normalized = 1.6259745436952323e-260;

const double myArray_max[XprecX] = {1.79769313486231570815e+308, 9.97920154767359795037e+291, 
                                    5.53956966280111259858e+275, 3.07507889307840487279e+259};
const dd_real _max = dd_real(myArray_max, XprecX);
const double myArray_safe_max[XprecX] = {1.7976931080746007281e+308, 9.97920154767359795037e+291, 
                                         5.53956966280111259858e+275, 3.07507889307840487279e+259};
const dd_real _safe_max = dd_real(myArray_safe_max, XprecX);

const int _ndigits = 62;*/

#endif
