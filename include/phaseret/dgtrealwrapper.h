#ifndef _dgtrealwrapper_h
#define _dgtrealwrapper_h
#define LTFAT_DOUBLE

#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#endif
#include "ltfat/types.h"

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
 *
 * \returns
 * Status code             | Description
 * ------------------------|--------------------------
 * LTFATERR_SUCCESS        | No error occurred
 * LTFATERR_NULLPOINTER    | \a g or \a p was NULL or \a c was NULL and flags != FFTW_ESTIMATE
 * LTFATERR_BADSIZE        | Signal length L is less or equal to 0.
 * LTFATERR_NOTPOSARG      | At least one of \f W, \f a, \f M, \f gl was less or equal to zero. 
 * LTFATERR_NOTAFRAME      | System does not form a frame
 * LTFATERR_BADTRALEN      | \a L is not divisible by both \a a and \a M. 
 * LTFATERR_INITFAILED     | The FFTW plan creation failed 
 * LTFATERR_NOTSUPPORTED   | This is a non-painless system but its support was not compiled
 * LTFATERR_CANNOTHAPPEN   | \a hint does not have a valid value from \a dgtreal_anasyn_hint or \a ptype is not valid value from \a ltfat_phaseconvention enum
 * LTFATERR_NOMEM          | Signalizes memory allocation error
 */
int
dgtreal_anasyn_init( const LTFAT_REAL g[], int gl, int L, int W, int a, int M,
                     LTFAT_COMPLEX c[], dgtreal_anasyn_hint hint,
                     ltfat_phaseconvention ptype, unsigned flags,
                     dgtreal_anasyn_plan** p);
/**
 * prd
 */
int
dgtreal_anasyn_execute_proj(dgtreal_anasyn_plan* p, const LTFAT_COMPLEX cin[],
                            LTFAT_COMPLEX cout[]);

/**
 * prd
 */
int
dgtreal_anasyn_execute_syn(dgtreal_anasyn_plan* p, const LTFAT_COMPLEX c[],
                           LTFAT_REAL f[]);

/**
 * prd
 */
int
dgtreal_anasyn_execute_ana(dgtreal_anasyn_plan* p, const LTFAT_REAL f[],
                           LTFAT_COMPLEX c[]);

/**
 * prd
 */
int
dgtreal_anasyn_done(dgtreal_anasyn_plan** p);

/** @} */
#endif
