#ifndef _config_h
#define _config_h

#include <math.h>
#include <string.h>
#if defined(__cplusplus) || (defined(__STDC_VERSION__) && (__STDC_VERSION__ < 199901L) ) || defined(NOC99COMPLEXH)
typedef double cdouble[2];
#else
#include <complex.h>
typedef double complex cdouble;
#endif
#include <stdlib.h>

#define ALIGNBYTES 16

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifdef PHASERET_SINGLE
#else
#endif

#endif /* _config_h */
