
#ifndef GND_MUTI_MATH_H_
#define GND_MUTI_MATH_H_


#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
#if defined(_INC_MATH) && !defined(_MATH_DEFINES_DEFINED)
#endif

#include <math.h>

#if defined(__linux__) || defined(__MINGW32__)
#define gnd_multi_round(x)	round((x))
#elif defined(_MSC_VER)
#define gnd_multi_round(x)	( (int) ( (x) + 0.5 ) )
#endif


#endif
