#ifndef Check_hh
#define Check_hh

#include <string>

/*
  Checks conditions at runtime
  Throws runtime error if condition isn't met
*/
namespace ch_ns
{
    void check(std::string condition,
               std::string file,
               int line);
    
    void check(std::string condition,
               std::string message,
               std::string file,
               int line);
} // namespace ch_ns

// Check only happens in debug
#ifdef NDEBUG
#  define Check(cond)
#  define CheckMsg(cond, desc)
#else
#  define Check(cond)                                           \
    if (!(cond)) ch_ns::check(#cond, __FILE__, __LINE__)
#  define CheckMsg(cond, desc)                                  \
    if (!(cond)) ch_ns::check(#cond, desc, __FILE__, __LINE__)
#endif

// Assert always happens
#define Assert(cond)                                            \
    if (!(cond)) ch_ns::check(#cond, __FILE__, __LINE__)
#define AssertMsg(cond, desc)                                   \
    if (!(cond)) ch_ns::check(#cond, desc, __FILE__, __LINE__)

#endif
