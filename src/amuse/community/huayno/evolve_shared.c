/*
 * Reference integrators with single, global shared time step.
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"

static void dkd4_S_m4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void dkd4_S_m5(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void dkd4_S_m6(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk4_S_m4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk4_S_m5(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk4_S_m6(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);

static void dkd6_SS_m11(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void dkd6_SS_m13(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk6_SS_m11(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk6_SS_m13(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);

static void dkd8_SS_m21(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk8_SS_m21(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);

static void dkd10_SS_m35(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void kdk10_SS_m35(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);

static void (*dkd4)(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)=dkd4_S_m6;
static void (*dkd6)(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)=dkd6_SS_m11;
static void (*dkd8)(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)=dkd8_SS_m21;
static void (*dkd10)(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)=dkd10_SS_m35;

void evolve_shared2(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  FLOAT dtsys;
  clevel++;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s,SIGN(dt));
  dtsys=global_timestep(s);
  if(dtsys < fabs(dt))
  {
    evolve_shared2(s,stime, stime+dt/2,dt/2,0);
    evolve_shared2(s,stime+dt/2, etime,dt/2,1);
  }
  else
  {
    deepsteps++;
    simtime+=dt;
    kdk(s,zerosys, stime, etime, dt);
  }
  clevel--;
}

void evolve_shared4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  clevel++;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared4(s,stime, stime+dt/2,dt/2,0);
    evolve_shared4(s,stime+dt/2, etime,dt/2,1);
  } else {
    deepsteps++;
    simtime+=dt;
    dkd4(s, stime, etime, dt);
  }
  clevel--;
}

void evolve_shared6(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  clevel++;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared6(s,stime, stime+dt/2,dt/2,0);
    evolve_shared6(s,stime+dt/2, etime,dt/2,1);
  } else {
    deepsteps++;
    simtime+=dt;
    dkd6(s, stime, etime, dt);
  }
  clevel--;
}

void evolve_shared8(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  clevel++;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared8(s,stime, stime+dt/2,dt/2,0);
    evolve_shared8(s,stime+dt/2, etime,dt/2,1);
  } else {
    deepsteps++;
    simtime+=dt;
    dkd8(s, stime, etime, dt);
  }
  clevel--;
}

void evolve_shared10(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  clevel++;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared10(s,stime, stime+dt/2,dt/2,0);
    evolve_shared10(s,stime+dt/2, etime,dt/2,1);
  } else {
    deepsteps++;
    simtime+=dt;
    dkd10(s, stime, etime, dt);
  }
  clevel--;
}

#define DRIFT(dt) \
  stime += dt; \
  drift(s, stime, dt);

#define KICK(dt) \
  kick(s, s, dt);

#define SPLIT_4TH_S_M6(EVOLVEA,EVOLVEB,dt) \
{ \
  const DOUBLE K1 = 0.0792036964311957L; \
  const DOUBLE K2 = 0.353172906049774L; \
  const DOUBLE K3 = -.0420650803577195L; \
  const DOUBLE K4 = 1 - 2*(K1 + K2 + K3); \
  const DOUBLE D1 = 0.209515106613362L; \
  const DOUBLE D2 = -.143851773179818L; \
  const DOUBLE D3 = 0.5L - D1 - D2; \
  EVOLVEA(K1*dt) \
  EVOLVEB(D1*dt) \
  EVOLVEA(K2*dt) \
  EVOLVEB(D2*dt) \
  EVOLVEA(K3*dt) \
  EVOLVEB(D3*dt) \
  EVOLVEA(K4*dt) \
  EVOLVEB(D3*dt) \
  EVOLVEA(K3*dt) \
  EVOLVEB(D2*dt) \
  EVOLVEA(K2*dt) \
  EVOLVEB(D1*dt) \
  EVOLVEA(K1*dt) \
} 

/* symmetric composition , 4th order, m=5*/
#define SPLIT_4TH_S_M5(EVOLVEA,EVOLVEB,dt) \
{ \
  const DOUBLE K1 = ( (14-sqrt((DOUBLE) 19))/108 ); \
  const DOUBLE K2 = ( (20-7*sqrt((DOUBLE) 19))/108 ); \
  const DOUBLE K3 = ( (((DOUBLE) 1)/2)-(K1+K2) ); \
  const DOUBLE D1 = ( ((DOUBLE) 2)/5 ); \
  const DOUBLE D2 = ( -((DOUBLE) 1)/10 ); \
  const DOUBLE D3 = ( 1-(2*D1+2*D2) ); \
  EVOLVEA(K1*dt) \
  EVOLVEB(D1*dt) \
  EVOLVEA(K2*dt) \
  EVOLVEB(D2*dt) \
  EVOLVEA(K3*dt) \
  EVOLVEB(D3*dt) \
  EVOLVEA(K3*dt) \
  EVOLVEB(D2*dt) \
  EVOLVEA(K2*dt) \
  EVOLVEB(D1*dt) \
  EVOLVEA(K1*dt) \
} 
    
/* symmetric composition , 4th order, m=4*/
#define SPLIT_4TH_S_M4(EVOLVEA,EVOLVEB,dt) \
{ \
  const DOUBLE K1 = ( (642+sqrt( (DOUBLE) 471 ))/3924 ); \
  const DOUBLE K2 = (  121*(12- sqrt( (DOUBLE) 471 ) )/3924 ); \
  const DOUBLE K3 = (  1-2*(K1+K2) ); \
  const DOUBLE D1 = ( ((DOUBLE) 6)/11 ); \
  const DOUBLE D2 = ( ((DOUBLE) 0.5)-D1 ); \
  EVOLVEA(K1*dt) \
  EVOLVEB(D1*dt) \
  EVOLVEA(K2*dt) \
  EVOLVEB(D2*dt) \
  EVOLVEA(K3*dt) \
  EVOLVEB(D2*dt) \
  EVOLVEA(K2*dt) \
  EVOLVEB(D1*dt) \
  EVOLVEA(K1*dt) \
} 

static void dkd4_S_m6(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_4TH_S_M6(DRIFT,KICK,dt)
}

static void kdk4_S_m6(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_4TH_S_M6(KICK,DRIFT,dt)
}

static void dkd4_S_m5(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_4TH_S_M5(DRIFT,KICK,dt)
}

static void kdk4_S_m5(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_4TH_S_M5(KICK,DRIFT,dt)
}

static void dkd4_S_m4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_4TH_S_M4(DRIFT,KICK,dt)
}

static void kdk4_S_m4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_4TH_S_M4(KICK,DRIFT,dt)
}

/* symmetric composition of symmetric maps, 6th order, m=11*/
#define SPLIT_6TH_SS_M11(EVOLVEA,EVOLVEB,dt) \
{ \
const DOUBLE C1 = ((DOUBLE)  0.21375583945878254555518066964857L ); \
const DOUBLE C2 = ((DOUBLE)  0.18329381407425713911385974425217L ); \
const DOUBLE C3 = ((DOUBLE)  0.17692819473098943794898811709929L ); \
const DOUBLE C4 = ((DOUBLE) -0.44329082681170215849622829626258L ); \
const DOUBLE C5 = ((DOUBLE)  0.11728560432865935385403585669136L ); \
const DOUBLE C6 = ((DOUBLE)  0.50405474843802736404832781714239L ); \
EVOLVEA(C1*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA(C1*dt/2) \
}

/* symmetric composition of symmetric maps, 6th order, m=13*/
#define SPLIT_6TH_SS_M13(EVOLVEA,EVOLVEB,dt) \
{ \
const DOUBLE C1 = ((DOUBLE)  0.13861930854051695245808013042625L ); \
const DOUBLE C2 = ((DOUBLE)  0.13346562851074760407046858832209L ); \
const DOUBLE C3 = ((DOUBLE)  0.13070531011449225190542755785015L ); \
const DOUBLE C4 = ((DOUBLE)  0.12961893756907034772505366537091L ); \
const DOUBLE C5 = ((DOUBLE) -0.35000324893920896516170830911323L ); \
const DOUBLE C6 = ((DOUBLE)  0.11805530653002387170273438954049L ); \
const DOUBLE C7 = ((DOUBLE)  0.39907751534871587459988795520665L ); \
EVOLVEA(C1*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C6+C7)*dt/2) \
EVOLVEB(C7*dt) \
EVOLVEA((C6+C7)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA(C1*dt/2) \
}

static void dkd6_SS_m11(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_6TH_SS_M11(DRIFT,KICK,dt)
}

static void kdk6_SS_m11(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_6TH_SS_M11(KICK,DRIFT,dt)
}

static void dkd6_SS_m13(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_6TH_SS_M13(DRIFT,KICK,dt)
}

static void kdk6_SS_m13(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_6TH_SS_M13(KICK,DRIFT,dt)
}

/* symmetric composition of symmetric maps, 8th order, m=21*/
#define SPLIT_8TH_SS_M21(EVOLVEA,EVOLVEB,dt) \
{ \
const DOUBLE C1  = ((DOUBLE)  0.10647728984550031823931967854896L ); \
const DOUBLE C2  = ((DOUBLE)  0.10837408645835726397433410591546L ); \
const DOUBLE C3  = ((DOUBLE)  0.35337821052654342419534541324080L ); \
const DOUBLE C4  = ((DOUBLE) -0.23341414023165082198780281128319L ); \
const DOUBLE C5  = ((DOUBLE) -0.24445266791528841269462171413216L ); \
const DOUBLE C6  = ((DOUBLE)  0.11317848435755633314700952515599L ); \
const DOUBLE C7  = ((DOUBLE)  0.11892905625000350062692972283951L ); \
const DOUBLE C8  = ((DOUBLE)  0.12603912321825988140305670268365L ); \
const DOUBLE C9  = ((DOUBLE)  0.12581718736176041804392391641587L ); \
const DOUBLE C10 = ((DOUBLE)  0.11699135019217642180722881433533L ); \
const DOUBLE C11 = ((DOUBLE) -0.38263596012643665350944670744040L ); \
EVOLVEA(C1*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C6+C7)*dt/2) \
EVOLVEB(C7*dt) \
EVOLVEA((C7+C8)*dt/2) \
EVOLVEB(C8*dt) \
EVOLVEA((C8+C9)*dt/2) \
EVOLVEB(C9*dt) \
EVOLVEA((C9+C10)*dt/2) \
EVOLVEB(C10*dt) \
EVOLVEA((C10+C11)*dt/2) \
EVOLVEB(C11*dt) \
EVOLVEA((C10+C11)*dt/2) \
EVOLVEB(C10*dt) \
EVOLVEA((C9+C10)*dt/2) \
EVOLVEB(C9*dt) \
EVOLVEA((C8+C9)*dt/2) \
EVOLVEB(C8*dt) \
EVOLVEA((C7+C8)*dt/2) \
EVOLVEB(C7*dt) \
EVOLVEA((C6+C7)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA(C1*dt/2) \
}

static void dkd8_SS_m21(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_8TH_SS_M21(DRIFT,KICK,dt)
}

static void kdk8_SS_m21(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_8TH_SS_M21(KICK,DRIFT,dt)
}

/* symmetric composition of symmetric maps, 10th order, m=35*/
#define SPLIT_10TH_SS_M35(EVOLVEA,EVOLVEB,dt) \
{ \
const DOUBLE C1  = ((DOUBLE)  0.078795722521686419263907679337684L ); \
const DOUBLE C2  = ((DOUBLE)  0.31309610341510852776481247192647L ); \
const DOUBLE C3  = ((DOUBLE)  0.027918383235078066109520273275299L ); \
const DOUBLE C4  = ((DOUBLE) -0.22959284159390709415121339679655L ); \
const DOUBLE C5  = ((DOUBLE)  0.13096206107716486317465685927961L ); \
const DOUBLE C6  = ((DOUBLE) -0.26973340565451071434460973222411L ); \
const DOUBLE C7  = ((DOUBLE)  0.074973343155891435666137105641410L ); \
const DOUBLE C8  = ((DOUBLE)  0.11199342399981020488957508073640L ); \
const DOUBLE C9  = ((DOUBLE)  0.36613344954622675119314812353150L ); \
const DOUBLE C10 = ((DOUBLE) -0.39910563013603589787862981058340L ); \
const DOUBLE C11 = ((DOUBLE)  0.10308739852747107731580277001372L ); \
const DOUBLE C12 = ((DOUBLE)  0.41143087395589023782070411897608L ); \
const DOUBLE C13 = ((DOUBLE) -0.0048663605831352617621956593099771L ); \
const DOUBLE C14 = ((DOUBLE) -0.39203335370863990644808193642610L ); \
const DOUBLE C15 = ((DOUBLE)  0.051942502962449647037182904015976L ); \
const DOUBLE C16 = ((DOUBLE)  0.050665090759924496335874344156866L ); \
const DOUBLE C17 = ((DOUBLE)  0.049674370639729879054568800279461L ); \
const DOUBLE C18 = ((DOUBLE)  0.049317735759594537917680008339338L ); \
EVOLVEA(C1*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C6+C7)*dt/2) \
EVOLVEB(C7*dt) \
EVOLVEA((C7+C8)*dt/2) \
EVOLVEB(C8*dt) \
EVOLVEA((C8+C9)*dt/2) \
EVOLVEB(C9*dt) \
EVOLVEA((C9+C10)*dt/2) \
EVOLVEB(C10*dt) \
EVOLVEA((C10+C11)*dt/2) \
EVOLVEB(C11*dt) \
EVOLVEA((C11+C12)*dt/2) \
EVOLVEB(C12*dt) \
EVOLVEA((C12+C13)*dt/2) \
EVOLVEB(C13*dt) \
EVOLVEA((C13+C14)*dt/2) \
EVOLVEB(C14*dt) \
EVOLVEA((C14+C15)*dt/2) \
EVOLVEB(C15*dt) \
EVOLVEA((C15+C16)*dt/2) \
EVOLVEB(C16*dt) \
EVOLVEA((C16+C17)*dt/2) \
EVOLVEB(C17*dt) \
EVOLVEA((C17+C18)*dt/2) \
EVOLVEB(C18*dt) \
EVOLVEA((C17+C18)*dt/2) \
EVOLVEB(C17*dt) \
EVOLVEA((C16+C17)*dt/2) \
EVOLVEB(C16*dt) \
EVOLVEA((C15+C16)*dt/2) \
EVOLVEB(C15*dt) \
EVOLVEA((C14+C15)*dt/2) \
EVOLVEB(C14*dt) \
EVOLVEA((C13+C14)*dt/2) \
EVOLVEB(C13*dt) \
EVOLVEA((C12+C13)*dt/2) \
EVOLVEB(C12*dt) \
EVOLVEA((C11+C12)*dt/2) \
EVOLVEB(C11*dt) \
EVOLVEA((C10+C11)*dt/2) \
EVOLVEB(C10*dt) \
EVOLVEA((C9+C10)*dt/2) \
EVOLVEB(C9*dt) \
EVOLVEA((C8+C9)*dt/2) \
EVOLVEB(C8*dt) \
EVOLVEA((C7+C8)*dt/2) \
EVOLVEB(C7*dt) \
EVOLVEA((C6+C7)*dt/2) \
EVOLVEB(C6*dt) \
EVOLVEA((C5+C6)*dt/2) \
EVOLVEB(C5*dt) \
EVOLVEA((C4+C5)*dt/2) \
EVOLVEB(C4*dt) \
EVOLVEA((C3+C4)*dt/2) \
EVOLVEB(C3*dt) \
EVOLVEA((C2+C3)*dt/2) \
EVOLVEB(C2*dt) \
EVOLVEA((C1+C2)*dt/2) \
EVOLVEB(C1*dt) \
EVOLVEA(C1*dt/2) \
}

static void dkd10_SS_m35(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_10TH_SS_M35(DRIFT,KICK,dt)
}

static void kdk10_SS_m35(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  SPLIT_10TH_SS_M35(KICK,DRIFT,dt)
}

#undef DRIFT
#undef KICK
