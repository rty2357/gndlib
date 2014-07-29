
#ifndef _GND_CROSS_COMPILE_H
#define _GND_CROSS_COMPILE_H

#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
#if defined(_INC_MATH) && !defined(_MATH_DEFINES_DEFINED)
#include <math.h>
#endif

// visual studio secure function
#if defined(_MSC_VER) // visual studio
#pragma warning(disable:4996)
#else
#endif


#endif // _GND_CROSS_COMPILE_H
