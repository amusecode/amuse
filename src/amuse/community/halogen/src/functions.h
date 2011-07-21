/* 
** functions.h
**
** Header file for functions.c
*/

DOUBLE Ecosmo(const GI (*));
DOUBLE integral(DOUBLE (*)(DOUBLE, const SI(*)), DOUBLE, DOUBLE, const SI (*));
DOUBLE trapez(DOUBLE (*)(DOUBLE, const SI(*)), DOUBLE, DOUBLE, INT, const SI (*));
DOUBLE integraldf(INT, const GI (*), const SI (*));
DOUBLE trapezdf(INT, INT, const GI (*), const SI (*));
INT locate(INT, const DOUBLE (*), DOUBLE);
DOUBLE lininterpolate(INT, const DOUBLE (*), const DOUBLE (*), DOUBLE);
DOUBLE rand01();
DOUBLE rho(DOUBLE, const SI (*));
DOUBLE drhodr(DOUBLE, const SI (*));
DOUBLE d2rhodr2(DOUBLE, const SI (*));
DOUBLE d2rhodPhi2(DOUBLE, const GI (*), const SI (*));
DOUBLE dlrhodlr(DOUBLE, const SI (*));
DOUBLE eta(DOUBLE, const SI (*));
DOUBLE detadr(DOUBLE, const SI (*));
DOUBLE tau(DOUBLE, const SI (*));
DOUBLE integrandIM(DOUBLE, const SI (*));
DOUBLE integrandIMcutoff(DOUBLE, const SI (*));
DOUBLE integrandMenc(DOUBLE, const SI (*));
DOUBLE integrandPot(DOUBLE, const SI (*));
INT split(INT, DOUBLE, const SI (*));
DOUBLE Menc(DOUBLE, const GI (*));
DOUBLE Pot(DOUBLE, const GI (*));
DOUBLE vescape(DOUBLE, const GI (*));
DOUBLE Tdyn(DOUBLE, const GI (*));
DOUBLE f1(DOUBLE, const SI (*));
DOUBLE f2(DOUBLE, const SI (*));
