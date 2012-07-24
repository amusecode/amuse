#ifdef NDEBUG

// the default assert will make empty assert macro's
// when NDEBUG is set
#include <cassert>

#else

#ifndef __LOCALASSERT_H__
#define __LOCALASSERT_H__ 1


#ifdef USE_SYSTEM_ASSERT

#include <cassert>

#else


#if defined __cplusplus
# define __ASSERT_VOID_CAST static_cast<void>
#else
# define __ASSERT_VOID_CAST (void)
#endif

#include <iostream>

class assert_failed
{
    public:
    
    __const char *assertion;
    __const char *file;
    unsigned int line;
    __const char *function;
        
    assert_failed(
        __const char *__assertion,
        __const char *__file,
        unsigned int __line, 
        __const char *__function
    ) : assertion(__assertion), file(__file), line(__line), function(__function)
    {
        
    }
    
    friend std::ostream& operator<<(std::ostream& output, const assert_failed& p)
    {
        output << p.file << ':' << p.line << " assertion failed: '" << p.assertion << "'(" << p.function << ")";
        return output;
    }
};

# define assert(expr)							\
  ((expr)								\
   ? __ASSERT_VOID_CAST (0)						\
   : throw assert_failed (__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION))



/* Version 2.4 and later of GCC define a magical variable `__PRETTY_FUNCTION__'
   which contains the name of the function currently being defined.
   This is broken in G++ before version 2.6.
   C9x has a similar variable called __func__, but prefer the GCC one since
   it demangles C++ function names.  */
# if defined __cplusplus && __GNUC__
#   define __ASSERT_FUNCTION	__PRETTY_FUNCTION__
# else
#  if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#   define __ASSERT_FUNCTION	__func__
#  else
#   define __ASSERT_FUNCTION	((__const char *) 0)
#  endif
# endif

#endif // USE_SYSTEM_ASSERT

#endif //__LOCALASSERT_H__

#endif // NDEBUG
