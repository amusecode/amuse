#ifdef DEBUG

#define assert(EX) if(.not. ( EX )) call FortranAssert('Assertion',#EX, __FILE__, __LINE__)
#define trace(EX) call FortranTrace(#EX, __FILE__, __LINE__)
#define tracei(EX) call FortranTraceInteger(#EX, EX, __FILE__, __LINE__)
#define tracer(EX) call FortranTraceReal(#EX, EX, __FILE__, __LINE__)
#define cond_trace(cond, EX) if (cond) trace(EX)
#define cond_tracei(cond, EX) if (cond) tracei(EX)
#define cond_tracer(cond, EX) if (cond) tracer(EX)
#define nan_check(EX) call FortranNanCheck(#EX,EX,__FILE__, __LINE__)


! Remove "pure" keyword so that we can debug pure subroutines in debug mode (using trace and assert)
#define pure

#else

#define assert(EX)
#define trace(EX)
#define tracei(EX)
#define tracer(EX)
#define cond_trace(cond, EX)
#define cond_tracei(cond, EX)
#define cond_tracer(cond, EX)
#define nan_check(EX)

#endif

