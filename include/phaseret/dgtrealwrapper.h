#ifndef _dgtrealwrapper_h
#define _dgtrealwrapper_h
#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#endif

typedef struct dgtreal_anasyn_plan dgtreal_anasyn_plan;

/** \defgroup dgtrealwrapper Discrete Gabor Transform analysis-synthesis 
 * \addtogroup dgtrealwrapper
 * @{
 */
typedef enum
{
    dgtreal_anasyn_auto,
    dgtreal_anasyn_long,
    dgtreal_anasyn_fb
} dgtreal_anasyn_hint;

/**
 * Note c can be NULL if FFTW_ESTIMATE is used in flags
 */
int
dgtreal_anasyn_init( const double g[], int gl, int L, int W, int a, int M,
                     complex double c[], dgtreal_anasyn_hint hint,
                     ltfat_phaseconvention ptype, unsigned flags,
                     dgtreal_anasyn_plan** p);

int
dgtreal_anasyn_execute_proj(dgtreal_anasyn_plan* p, const complex double cin[],
                            complex double cout[]);

int
dgtreal_anasyn_execute_syn(dgtreal_anasyn_plan* p, const complex double c[],
                           double f[]);

int
dgtreal_anasyn_execute_ana(dgtreal_anasyn_plan* p, const double f[],
                           complex double c[]);

int
dgtreal_anasyn_done(dgtreal_anasyn_plan** p);

/** @} */
#endif
