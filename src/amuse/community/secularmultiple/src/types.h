#include <math.h>
#include <cstdlib>
#include <map>
#include <vector>

/* constants */
/* units (cf. interface.py): 
 * unit_l = units.AU
 * unit_m = units.MSun
 * unit_t = 1.0e6*units.yr
 */
 
#ifndef __FOUND_ROOT
#define ___FOUND_ROOT
#define FOUND_ROOT ((roots_found[i_root] == 1) || (roots_found[i_root] == -1))
#endif 

#ifndef __CONSTANTS
#define __CONSTANTS
#define CONST_G			    (double)	3.94852492465e+13
#define CONST_G_P2          (double)    CONST_G*CONST_G
#define CONST_G_P3          (double)    CONST_G*CONST_G_P2
#define CONST_C_LIGHT		(double)	63239726386.8
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define CONST_MSUN          (double)    1.0
#define CONST_R_SUN         (double)    0.00464913034382
#define CONST_L_SUN         (double)    2.71040410975e+14

#define c_1div2             (double)    1.0/2.0
#define c_1div3             (double)    1.0/3.0
#define c_1div4             (double)    1.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_1div10            (double)    1.0/10.0
#define c_1div15            (double)    1.0/15.0
#define c_1div16            (double)    1.0/16.0
#define c_1div30            (double)    1.0/30.0
#define c_2div3             (double)    2.0/3.0
#define c_3div2             (double)    3.0/2.0
#define c_3div4             (double)    3.0/4.0
#define c_3div5             (double)    3.0/5.0
#define c_3div8             (double)    3.0/8.0
#define c_3div32            (double)    3.0/32.0
#define c_3div1024          (double)    3.0/1024.0
#define c_4div15            (double)    4.0/15.0
#define c_5div2             (double)    5.0/2.0
#define c_5div8             (double)    5.0/8.0
#define c_5div16            (double)    5.0/16.0
#define c_5div64            (double)    5.0/64.0
#define c_7div8             (double)    7.0/8.0
#define c_8div5             (double)    8.0/5.0
#define c_8div7             (double)    8.0/7.0
#define c_9div2             (double)    9.0/2.0
#define c_9div16            (double)    9.0/16.0
#define c_9div32            (double)    9.0/32.0
#define c_11div18           (double)    11.0/18.0
#define c_15div2            (double)    15.0/2.0
#define c_15div4            (double)    15.0/4.0
#define c_15div8            (double)    15.0/8.0
#define c_15div16           (double)    15.0/16.0
#define c_16div5            (double)    16.0/5.0
#define c_25div16           (double)    25.0/16.0
#define c_25div64           (double)    25.0/64.0
#define c_31div2            (double)    31.0/2.0
#define c_32div5            (double)    32.0/5.0
#define c_37div96           (double)    37.0/96.0
#define c_45div8            (double)    45.0/8.0
#define c_64div5            (double)    64.0/5.0
#define c_73div24           (double)    73.0/24.0
#define c_105div4096        (double)    105.0/4096.0
#define c_121div304         (double)    121.0/304.0
#define c_185div16          (double)    185.0/16.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0
#endif


#ifndef __TABLES
#define __TABLES
#define MAX_ORDER (int) 5

#define TABLEWIDTH_A (int) 3
#define TABLELENGTH_A (int) 10
const double A_TABLE[TABLELENGTH_A][TABLEWIDTH_A] =
{{2, 0, -0.500000000000000}, {2, 2, 1.50000000000000}, {3, 
  1, -1.50000000000000}, {3, 3, 2.50000000000000}, {4, 0, 
  0.375000000000000}, {4, 2, -3.75000000000000}, {4, 4, 
  4.37500000000000}, {5, 1, 1.87500000000000}, {5, 
  3, -8.75000000000000}, {5, 5, 7.87500000000000}};
  
#define TABLEWIDTH_B (int) 8
#define TABLELENGTH_B (int) 47
#define HIGHEST_POWER_ECCP2_IN_B_TABLE (int) 3
const double B_TABLE[TABLELENGTH_B][TABLEWIDTH_B] =
{{2, 0, 0, 0, 1.00000000000000, 1.50000000000000, 0, 0}, {2, 1, 1, 
  0, -2.00000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 0, 
  0.500000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 
  2, -0.500000000000000, 0, 0, 0}, {2, 2, 2, 0, 2.50000000000000, 0, 
  0, 0}, {3, 0, 0, 0, 1.00000000000000, 3.00000000000000, 
  0.375000000000000, 0}, {3, 1, 1, 
  0, -2.50000000000000, -1.87500000000000, 0, 0}, {3, 2, 0, 0, 
  0.500000000000000, -0.375000000000000, -0.125000000000000, 0}, {3, 
  2, 0, 2, -0.500000000000000, -0.125000000000000, 0, 0}, {3, 2, 2, 0,
   3.75000000000000, 0.625000000000000, 0, 0}, {3, 3, 1, 
  0, -1.87500000000000, 1.87500000000000, 0, 0}, {3, 3, 1, 2, 
  1.87500000000000, 0, 0, 0}, {3, 3, 3, 0, -4.37500000000000, 0, 0, 
  0}, {4, 0, 0, 0, 1.00000000000000, 5.00000000000000, 
  1.87500000000000, 0}, {4, 1, 1, 
  0, -3.00000000000000, -4.50000000000000, -0.375000000000000, 0}, {4,
   2, 0, 0, 0.500000000000000, -0.125000000000000, -0.375000000000000,
   0}, {4, 2, 0, 2, -0.500000000000000, -0.375000000000000, 0, 0}, {4,
   2, 2, 0, 5.25000000000000, 2.62500000000000, 0, 0}, {4, 3, 1, 
  0, -2.25000000000000, 1.87500000000000, 0.375000000000000, 0}, {4, 
  3, 1, 2, 2.25000000000000, 0.375000000000000, 0, 0}, {4, 3, 3, 
  0, -7.00000000000000, -0.875000000000000, 0, 0}, {4, 4, 0, 0, 
  0.375000000000000, -0.750000000000000, 0.375000000000000, 0}, {4, 4,
   0, 2, -0.750000000000000, 0.750000000000000, 0, 0}, {4, 4, 0, 4, 
  0.375000000000000, 0, 0, 0}, {4, 4, 2, 0, 
  5.25000000000000, -5.25000000000000, 0, 0}, {4, 4, 2, 
  2, -5.25000000000000, 0, 0, 0}, {4, 4, 4, 0, 7.87500000000000, 0, 0,
   0}, {5, 0, 0, 0, 1.00000000000000, 7.50000000000000, 
  5.62500000000000, 0.312500000000000}, {5, 1, 1, 
  0, -3.50000000000000, -8.75000000000000, -2.18750000000000, 0}, {5, 
  2, 0, 0, 0.500000000000000, 
  0.250000000000000, -0.687500000000000, -0.0625000000000000}, {5, 2, 
  0, 2, -0.500000000000000, -0.750000000000000, -0.0625000000000000, 
  0}, {5, 2, 2, 0, 7.00000000000000, 7.00000000000000, 
  0.437500000000000, 0}, {5, 3, 1, 0, -2.62500000000000, 
  1.31250000000000, 1.31250000000000, 0}, {5, 3, 1, 2, 
  2.62500000000000, 1.31250000000000, 0, 0}, {5, 3, 3, 
  0, -10.5000000000000, -3.93750000000000, 0, 0}, {5, 4, 0, 0, 
  0.375000000000000, -0.687500000000000, 0.250000000000000, 
  0.0625000000000000}, {5, 4, 0, 2, -0.750000000000000, 
  0.625000000000000, 0.125000000000000, 0}, {5, 4, 0, 4, 
  0.375000000000000, 0.0625000000000000, 0, 0}, {5, 4, 2, 0, 
  7.00000000000000, -6.12500000000000, -0.875000000000000, 0}, {5, 4, 
  2, 2, -7.00000000000000, -0.875000000000000, 0, 0}, {5, 4, 4, 0, 
  13.1250000000000, 1.31250000000000, 0, 0}, {5, 5, 1, 
  0, -2.18750000000000, 4.37500000000000, -2.18750000000000, 0}, {5, 
  5, 1, 2, 4.37500000000000, -4.37500000000000, 0, 0}, {5, 5, 1, 
  4, -2.18750000000000, 0, 0, 0}, {5, 5, 3, 0, -13.1250000000000, 
  13.1250000000000, 0, 0}, {5, 5, 3, 2, 13.1250000000000, 0, 0, 
  0}, {5, 5, 5, 0, -14.4375000000000, 0, 0, 0}};

#define TABLEWIDTH_D (int) 8
#define TABLELENGTH_D (int) 15
const double D_TABLE[TABLELENGTH_D][TABLEWIDTH_D] = 
{
    {2, 0, 0, 0, 0, 0, 0, 1}, \
    {2, 0, 2, 0, 0, 0, 0, 2}, \
    {2, 0, 2, 0, 2, 0, 0, 3}, \
    {2, 0, 2, 0, 0, 0, 2, 4}, \
    {2, 2, 0, 0, 0, 0, 0, 5}, \
    {2, 2, 0, 2, 0, 0, 0, 6}, \
    {2, 2, 0, 0, 0, 2, 0, 7}, \
    {3, 1, 0, 1, 0, 0, 0, 8}, \
    {3, 1, 2, 0, 1, 1, 1, 9}, \
    {3, 1, 2, 1, 0, 0, 0, 10}, \
    {3, 1, 2, 1, 0, 0, 2, 11}, \
    {3, 1, 2, 1, 2, 0, 0, 12}, \
    {3, 3, 0, 1, 0, 0, 0, 13}, \
    {3, 3, 0, 1, 0, 2, 0, 14}, \
    {3, 3, 0, 3, 0, 0, 0, 15}
};

/* the numbers in each last entry of the D table refer to the functions defined below */

inline double D_TABLE_FUNC1(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(sqrt_ef_p2_minus_one + asec_minus_ef);
}
inline double D_TABLE_FUNC1_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC2(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC2_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC3(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC3_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC4(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC4_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC5(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC5_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC6(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC6_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC7(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC7_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC8(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(c_1div3*one_div_ef_p1*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + ef*asec_minus_ef);
}
inline double D_TABLE_FUNC8_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC9(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div15*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - ef*asec_minus_ef;
}
inline double D_TABLE_FUNC9_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC10(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC10_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC11(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC11_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC12(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_1div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC12_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC13(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC13_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC14(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_3div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC14_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC15(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC15_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}



#endif

/*	ODE solver macros	*/
#ifndef __ODE_MACROS
#define __ODE_MACROS
    #define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
    #define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
    #ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
    #endif
    #ifndef min
        #define min(X,Y) ((X) < (Y) ? (X) : (Y))
    #endif
#endif

/* vector operators */
#ifndef __VECTOR_OPERATORS
#define __VECTOR_OPERATORS
inline void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
inline double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
inline double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}
#endif

/* classes */
#ifndef __Particle
#define __Particle
class Particle
{
    public:
    /* generic properties */
    int index,child1,child2;
    int parent;
    int sibling;
    std::vector<int> parents;
    std::vector<int> connecting_child_in_parents;
    int level,highest_level;
    int is_binary;
    double mass,mass_dot_external,child1_mass,child2_mass,total_system_mass;

    /*******************
    /* body properties *
     * ****************/
    /* general */
    double radius,radius_dot_external,radius_ddot_external;
    double spin_vec_x,spin_vec_y,spin_vec_z;
    double spin_vec_x_dot_external,spin_vec_y_dot_external,spin_vec_z_dot_external;
    int stellar_type;
    
    double position_x,position_y,position_z;
    double velocity_x,velocity_y,velocity_z;

    /* used in ODE solver only */
    double spin_vec[3],dspin_vec_dt[3]; 
    double spin_vec_norm;
    double dmass_dt,dradius_dt;    
    
    void set_ODE_quantities(double delta_time);
    void reset_ODE_quantities();
    
    /*********************
    /* binary properties *
     * ******************/
    /* general */
    double e_vec_x,e_vec_y,e_vec_z;
    double h_vec_x,h_vec_y,h_vec_z;

    double e_vec_x_dot_external,e_vec_y_dot_external,e_vec_z_dot_external;
    double h_vec_x_dot_external,h_vec_y_dot_external,h_vec_z_dot_external;
   
    double true_anomaly;
    
    /* PN terms */
    int include_pairwise_1PN_terms,include_pairwise_25PN_terms;
        
    /* tidal friction */
    int include_tidal_friction_terms,tides_method,include_tidal_bulges_precession_terms,include_rotation_precession_terms;
    double minimum_eccentricity_for_tidal_precession;
    double tides_Q_prime; /* depricated */
    double tides_apsidal_motion_constant, tides_time_lag, tides_gyration_radius;
    double tides_viscous_time_scale;
    int tides_viscous_time_scale_prescription;
    double convective_envelope_mass,convective_envelope_radius,luminosity;

    /* root finding */
    int check_for_secular_breakdown,secular_breakdown_has_occurred;
    int check_for_dynamical_instability,dynamical_instability_has_occurred,dynamical_instability_criterion;
    int dynamical_instability_central_particle;
    double dynamical_instability_K_parameter;
    int check_for_physical_collision_or_orbit_crossing,physical_collision_or_orbit_crossing_has_occurred;
    int check_for_minimum_periapse_distance,minimum_periapse_distance_has_occurred;
    double check_for_minimum_periapse_distance_value;
    int check_for_RLOF_at_pericentre,check_for_RLOF_at_pericentre_use_sepinsky_fit,RLOF_at_pericentre_has_occurred;

    /* used in ODE solver only */
    double e_vec[3],h_vec[3];
    double e_vec_unit[3],h_vec_unit[3];    
    double de_vec_dt[3],dh_vec_dt[3];
    double child1_mass_plus_child2_mass,child1_mass_minus_child2_mass,child1_mass_times_child2_mass;
    double child1_mass_dot_external,child2_mass_dot_external;
    double e,e_p2;
    double j,j_p2,j_p3,j_p4,j_p5; // j=sqrt(1-e^2)
    double h,a;
    
//    Particle(int index, int is_binary, int child1, int child2, double mass, double radius, double spin_vec_x, double spin_vec_y, double spin_vec_z, double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms, int include_tides_terms, double tides_Q_prime, double tides_gyration_radius, int check_for_secular_breakdown, int secular_breakdown_has_occurred, int check_for_dynamical_instability, int dynamical_instability_has_occurred, int dynamical_instability_criterion, int check_for_physical_collision, int physical_collision_has_occurred) : index(index), is_binary(is_binary), child1(child1), child2(child2), mass(mass), radius(radius), spin_vec_x(spin_vec_x), spin_vec_y(spin_vec_y), spin_vec_z(spin_vec_z), e_vec_x(e_vec_x), e_vec_y(e_vec_y), e_vec_z(e_vec_z), h_vec_x(h_vec_x), h_vec_y(h_vec_y), h_vec_z(h_vec_z), include_pairwise_1PN_terms(include_pairwise_1PN_terms), include_pairwise_25PN_terms(include_pairwise_25PN_terms), include_tides_terms(include_tides_terms), tides_Q_prime(tides_Q_prime), tides_gyration_radius(tides_gyration_radius), 

    /* user-specified instantaneous perturbations */
    int sample_orbital_phases_randomly;
    double instantaneous_perturbation_delta_mass;
    double instantaneous_perturbation_delta_position_x,instantaneous_perturbation_delta_position_y,instantaneous_perturbation_delta_position_z;
    double instantaneous_perturbation_delta_velocity_x,instantaneous_perturbation_delta_velocity_y,instantaneous_perturbation_delta_velocity_z;

    Particle(int index, int is_binary) : index(index), is_binary(is_binary)
    {
        /* default values */
        check_for_secular_breakdown = 0;
        check_for_dynamical_instability = 0;
        check_for_physical_collision_or_orbit_crossing = 0;
        check_for_minimum_periapse_distance = 0;
        check_for_RLOF_at_pericentre = 0;

        secular_breakdown_has_occurred = 0;
        dynamical_instability_has_occurred = 0;
        physical_collision_or_orbit_crossing_has_occurred = 0;
        minimum_periapse_distance_has_occurred = 0;
        RLOF_at_pericentre_has_occurred = 0;
        
        
        include_pairwise_1PN_terms = 0;
        include_pairwise_25PN_terms = 0;
        include_tidal_friction_terms = 0;
        include_tidal_bulges_precession_terms = 0;
        include_rotation_precession_terms = 0;

        radius = 1.0e-15; /* this must be set (to nonzero), otherwise ODE solver will have invalid ewt values */
        tides_viscous_time_scale_prescription = 0; /* constant, user-specified t_V */
        minimum_eccentricity_for_tidal_precession = 1.0e-3;
        spin_vec_x_dot_external = spin_vec_y_dot_external = spin_vec_z_dot_external = 0.0;
        mass_dot_external = radius_dot_external = radius_ddot_external = 0.0;
        
        sample_orbital_phases_randomly = 1;
        instantaneous_perturbation_delta_mass = 0.0;
        instantaneous_perturbation_delta_position_x = instantaneous_perturbation_delta_position_y = instantaneous_perturbation_delta_position_z = 0.0;
        instantaneous_perturbation_delta_velocity_x = instantaneous_perturbation_delta_velocity_y = instantaneous_perturbation_delta_velocity_z = 0.0;
        
    }
};
#endif

typedef std::map<int, Particle *> ParticlesMap;
typedef std::map<int, Particle *>::iterator ParticlesMapIterator;




#ifndef __External_Particle
#define __External_Particle
class External_Particle
{
    public:
    /* generic properties */
    int index;
    int mode;
    int path;
    double mass;
    double t_ref,t_passed;
    double eccentricity,periapse_distance;

    double r_vec_x,r_vec_y,r_vec_z;
    
    /* straight line */
    double r0_vec_x,r0_vec_y,r0_vec_z; 
    double rdot_vec_x,rdot_vec_y,rdot_vec_z;
    
    /* hyperbolic orbit */
    double e_hat_vec_x,e_hat_vec_y,e_hat_vec_z;
    double h_hat_vec_x,h_hat_vec_y,h_hat_vec_z;
    
    External_Particle(int index) : index(index) 
    {
        mode = 0;
        path = 0;
    }
    
};
#endif

typedef std::map<int, External_Particle *> External_ParticlesMap;
typedef std::map<int, External_Particle *>::iterator External_ParticlesMapIterator;



/* CVODE UserData */
#ifndef __UserData
#define __UserData
typedef struct {
	ParticlesMap *particlesMap;
    External_ParticlesMap *external_particlesMap;
    double hamiltonian;
    int N_root_finding;
    double start_time;
} *UserData;
#endif
