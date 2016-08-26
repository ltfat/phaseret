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
typedef struct phaseret_dgtreal_params phaseret_dgtreal_params;

/** \addtogroup dgtrealwrapper
 * @{ */
typedef enum
{
    phaseret_dgtreal_auto,
    phaseret_dgtreal_long,
    phaseret_dgtreal_fb
} phaseret_dgtreal_hint;

/** Allocate dgtreal_params structure and initialize to default values
 *
 * \warning The structure must be freed using phaseret_dgtreal_params_free()
 *
 * \returns Allocated struct (or NULL if the memory allocation failed
 * \see phaseret_dgtreal_params_free
 */
PHASERET_API phaseret_dgtreal_params*
phaseret_dgtreal_params_allocdef();

/** Set DGT phase convention
 *
 * \returns
 * Status code          |  Description
 * ---------------------|----------------
 * LTFATERR_SUCESS      |  No error occured
 * LTFATERR_NULLPOINTER |  \a params was NULL
 */
PHASERET_API int
phaseret_dgtreal_params_set_phaseconv(phaseret_dgtreal_params* params, ltfat_phaseconvention ptype);

/** Set FFTW flags
 *
 * \returns
 * Status code          |  Description
 * ---------------------|----------------
 * LTFATERR_SUCESS      |  No error occured
 * LTFATERR_NULLPOINTER |  \a params was NULL
 */
PHASERET_API int
phaseret_dgtreal_params_set_fftwflags(phaseret_dgtreal_params* params, unsigned fftw_flags);

/** Set algorithm hint
 *
 * \returns
 * Status code          |  Description
 * ---------------------|----------------
 * LTFATERR_SUCESS      |  No error occured
 * LTFATERR_NULLPOINTER |  \a params was NULL
 */
PHASERET_API int
phaseret_dgtreal_params_set_hint(phaseret_dgtreal_params* params, phaseret_dgtreal_hint hint);

/** Destroy struct
 *
 * \returns
 * Status code          |  Description
 * ---------------------|----------------
 * LTFATERR_SUCESS      |  No error occured
 * LTFATERR_NULLPOINTER |  \a params was NULL
 */
PHASERET_API int
phaseret_dgtreal_params_free(phaseret_dgtreal_params* params);

/** @} */

// The following function is not part of API
int
phaseret_dgtreal_params_defaults(phaseret_dgtreal_params* params);

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
PHASERET_NAME(dgtreal_init)(const LTFAT_REAL g[], ltfat_int gl, ltfat_int L, ltfat_int W, ltfat_int a, ltfat_int M,
                            LTFAT_COMPLEX c[], phaseret_dgtreal_params* params,
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
