/** \defgroup dgtrealwrapper Discrete Gabor Transform analysis-synthesis
*/

#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "ltfat/types.h"
#include "phaseret/types.h"

#ifndef _phaseret_dgtrealwrapper_h
#define _phaseret_dgtrealwrapper_h
/** \addtogroup dgtrealwrapper
 * @{ */
typedef enum
{
    phaseret_dgtreal_auto,
    phaseret_dgtreal_long,
    phaseret_dgtreal_fb
} phaseret_dgtreal_hint;

typedef struct
{
    ltfat_phaseconvention ptype;
    unsigned fftw_flags;
    phaseret_dgtreal_hint hint;
} phaseret_dgtreal_init_params;


/** Intilialize extra param struct for dgtreal_init
 * 
 */
PHASERET_API int
phaseret_dgtreal_init_params_defaults(phaseret_dgtreal_init_params* params);

/** @} */
#endif

typedef struct PHASERET_NAME(dgtreal_plan) PHASERET_NAME(dgtreal_plan);

/** \addtogroup dgtrealwrapper
 * @{
 *
 */

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
 * LTFATERR_CANNOTHAPPEN   | \a hint does not have a valid value from \a phaseret_dgtreal_hint or \a ptype is not valid value from \a ltfat_phaseconvention enum
 * LTFATERR_NOMEM          | Signalizes memory allocation error
 */
PHASERET_API int
PHASERET_NAME(dgtreal_init)(const LTFAT_REAL g[], int gl, int L, int W, int a, int M,
                            LTFAT_COMPLEX c[], phaseret_dgtreal_init_params* params,
                            PHASERET_NAME(dgtreal_plan)** p);


/** Perform DGTREAL synthesis followed by analysis
 *
 * \note This function CAN work inplace.
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]    p  Transform plan
 * \param[in]  cin  Input coefficients, size M2 x N x W
 * \param[out]   c  Coefficients after projection, size M2 x N x W
 *
 * #### Versions #
 * <tt>
 * phaseret_dgtreal_execute_proj_d(phaseret_dgtreal_plan_d* p,
 *                                        const ltfat_complex_d cin[],
 *                                        ltfat_complex_d c[]);
 *
 * phaseret_dgtreal_execute_proj_s(phaseret_dgtreal_plan_s* p,
 *                                        const ltfat_complex_s cin[],
 *                                        ltfat_complex_s c[]);
 * </tt>
 *
 * \returns
 */
PHASERET_API int
PHASERET_NAME(dgtreal_execute_proj)(PHASERET_NAME(dgtreal_plan)* p,
                                    const LTFAT_COMPLEX cin[], LTFAT_COMPLEX c[]);

/** Perform DGTREAL synthesis
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]    p  Transform plan
 * \param[in]    c  Input coefficients, size M2 x N x W
 * \param[out]   f  Reconstructed signal, size L x W
 *
 * #### Versions #
 * <tt>
 * phaseret_dgtreal_execute_syn_d(phaseret_dgtreal_plan_d* p,
 *                                       const ltfat_complex_d c[],
 *                                       double f[]);
 *
 * phaseret_dgtreal_execute_syn_s(phaseret_dgtreal_plan_s* p,
 *                                       const ltfat_complex_s c[],
 *                                       float f[]);
 * </tt>
 *
 * \returns
 */
PHASERET_API int
PHASERET_NAME(dgtreal_execute_syn)(PHASERET_NAME(dgtreal_plan)* p,
                                   const LTFAT_COMPLEX c[], LTFAT_REAL f[]);

/** Perform DGTREAL analysis
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]    p  Transform plan
 * \param[in]    f  Input signal, size L x W
 * \param[out]   c  Coefficients, size M2 x N x W
 *
 * #### Versions #
 * <tt>
 * phaseret_dgtreal_execute_ana_d(phaseret_dgtreal_plan_d* p,
 *                                       const double f[], ltfat_complex_d c[]);
 *
 * phaseret_dgtreal_execute_ana_s(phaseret_dgtreal_plan_s* p,
 *                                       const float f[], ltfat_complex_s c[]);
 * </tt>
 * \returns
 */
PHASERET_API int
PHASERET_NAME(dgtreal_execute_ana)(PHASERET_NAME(dgtreal_plan)* p,
                                   const LTFAT_REAL f[], LTFAT_COMPLEX c[]);

/** Destroy transform plan
 *
 * \param[in]   p  Transform plan
 *
 * #### Versions #
 * <tt>
 * phaseret_dgtreal_done_d(phaseret_dgtreal_plan_d** p);
 *
 * phaseret_dgtreal_done_s(phaseret_dgtreal_plan_s** p);
 * </tt>
 * \returns
 */
PHASERET_API int
PHASERET_NAME(dgtreal_done)(PHASERET_NAME(dgtreal_plan)** p);

/** @} */



#ifdef __cplusplus
}
#endif
