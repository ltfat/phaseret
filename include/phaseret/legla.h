#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#ifndef _phaseret_dgtrealwrapper_h_included
#define _phaseret_dgtrealwrapper_h_included
#include "dgtrealwrapper.h"
#endif
#endif

#include "ltfat/types.h"
#include "phaseret/types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _phaseret_legla_h
#define _phaseret_legla_h
/** \addtogroup legla
 *  @{
 */

typedef enum
{
    EXT_BOTH = (1 << 20), // << DEFAULT
    EXT_UPDOWN = (1 << 21)
} leglaupdate_ext;

typedef enum
{
    ORDER_FWD = (1 << 10), // << DEFAULT
    ORDER_REV = (1 << 11)
} leglaupdate_frameorder;

typedef struct
{
    int height;
    int width;
} phaseret_size;

// typedef struct
// {
//     int y;
//     int x;
// } phaseret_point;

typedef struct
{
    double relthr; ///< Relative threshold for automatic determination of kernel size, default 1e-3
    phaseret_size ksize; ///< Maximum allowed kernel size (default 2*ceil(M/a) -1) or kernel size directly if relthr==0.0
    unsigned leglaflags; ///< LEGLA algorithm flags, default MOD_COEFFICIENTWISE | MOD_MODIFIEDUPDATE
    phaseret_dgtreal_init_params dparams;
} phaseret_legla_init_params;

typedef enum
{
    MOD_STEPWISE = (1 << 0),  // << DEFAULT
    MOD_FRAMEWISE = (1 << 1),
    MOD_COEFFICIENTWISE = (1 << 2),
    MOD_COEFFICIENTWISE_SORTED = (1 << 3),
    MOD_MODIFIEDUPDATE = (1 << 4),
} leglaupdate_mod;

/** Initialize parameter struct for legla_init
 *
 * \param[in,out]  params  Pointer to a structure to be inilialized
 * \returns
 * Status code          |  Description
 * ---------------------|-------------------
 * LTFATERR_SUCCESS     | No error occurred
 * LTFATERR_NULLPOINTER | \a params was NULL
 */
PHASERET_API int
phaseret_legla_init_params_defaults(phaseret_legla_init_params* params);

/** @} */
/* Util */
int phaseret_gcd(int m, int n);
int phaseret_lcm(int m, int n);

#endif

typedef struct PHASERET_NAME(legla_plan) PHASERET_NAME(legla_plan);
typedef struct PHASERET_NAME(leglaupdate_plan) PHASERET_NAME(leglaupdate_plan);
typedef struct PHASERET_NAME(leglaupdate_plan_col) PHASERET_NAME(leglaupdate_plan_col);


/** \addtogroup legla
 *  @{
 */

/** Function prototype for coefficient modification callback
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  \param[in]  userdata   User defined data
 *  \param[in,out]     c   Set of coefficients to by updated, size M2 x N x W
 *  \param[in]         L   Signal length
 *  \param[in]         W   Number of signal channels
 *  \param[in]         a   Time hop factor
 *  \param[in]         M   Number of frequency channels
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_callback_cmod_d(void* userdata, ltfat_complex_d c[],
 *                                ltfatInt L, ltfatInt W, ltfatInt a, ltfatInt M);
 *
 * phaseret_legla_callback_cmod_s(void* userdata, ltfat_complex_s c[],
 *                                ltfatInt L, ltfatInt W, ltfatInt a, ltfatInt M);
 * </tt>
 *  \returns
 *  Status code | Meaning
 *  ------------|-----------------
 *   0          | Signalizes that callback exited without error
 *  <0          | Callback exited with error
 */
typedef int
PHASERET_NAME(legla_callback_cmod)(void* userdata, LTFAT_COMPLEX c[], int L, int W, int a, int M);

/** Function prototype for status callback
 *
 *  The callback is executed at the end of each iteration.
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  \param[in]         p   DGTREAL analysis-synthesis plan
 *  \param[in]  userdata   User defined data
 *  \param[in,out]     c   Set of coefficients at the end of iteration, size M2 x N x W
 *  \param[in]         L   Signal length
 *  \param[in]         W   Number of signal channels
 *  \param[in]         a   Time hop factor
 *  \param[in]         M   Number of frequency channels
 *  \param[in]     alpha   Acceleration parameter
 *  \param[in]      iter   Current iteration
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_callback_status_d(phaseret_dgtreal_plan_d* p,
 *                                 void* userdata, ltfat_complex_d c[],
 *                                 ltfatInt L, ltfatInt W, ltfatInt a,
 *                                 ltfatint M, double* alpha, ltfatInt iter);
 *
 *
 * phaseret_legla_callback_status_s(phaseret_dgtreal_plan_s* p,
 *                                 void* userdata, ltfat_complex_s c[],
 *                                 ltfatInt L, ltfatInt W, ltfatInt a,
 *                                 ltfatint M, double* alpha, ltfatInt iter);
 * </tt>
 *  \returns
 *  Status code | Meaning
 *  ------------|---------------------------------------------------------------------------------------
 *   0          | Signalizes that callback exited without error
 *  >0          | Signalizes that callback exited without error but terminate the algorithm prematurely
 *  <0          | Callback exited with error
 *
 *  \see dgtreal_execute_proj dgtreal_execute_ana dgtreal_execute_syn
 */
typedef int
PHASERET_NAME(legla_callback_status)(PHASERET_NAME(dgtreal_plan)* p, void* userdata, LTFAT_COMPLEX c[],
                                     int L, int W, int a, int M, double* alpha, int iter);

/** Le Roux's Griffin-Lim algorithm
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  \param[in]  cinit   Initial set of coefficients, size M2 x N x W
 *  \param[in]      g   Analysis window, size gl x 1
 *  \param[in]      L   Signal length
 *  \param[in]     gl   Window length
 *  \param[in]      W   Number of signal channels
 *  \param[in]      a   Time hop factor
 *  \param[in]      M   Number of frequency channels
 *  \param[in]   iter   Number of iterations
 *  \param[out]     c   Coefficients with reconstructed phase, size M2 x N x W
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_d(const ltfat_complex_d cinit[], const double g[],
 *                  ltfatInt L, ltfatInt gl, ltfatInt W, ltfatInt a, ltfatInt M,
 *                  ltfatInt iter, ltfat_complex_d c[]);
 *
 * phaseret_legla_s(const ltfat_complex_s cinit[], const float g[],
 *                  ltfatInt L, ltfatInt gl, ltfatInt W, ltfatInt a, ltfatInt M,
 *                  ltfatInt iter, ltfat_complex_s c[]);
 * </tt>
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a cinit or \a g or \a cout was NULL
 *  LTFATERR_NOTPOSARG    | \a iter must be positive
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 */
PHASERET_API int
PHASERET_NAME(legla)(const LTFAT_COMPLEX cinit[], const LTFAT_REAL g[], const int L, const int gl,
                     const int W, const int a, const int M, const int iter, LTFAT_COMPLEX c[]);

/** Le Roux's Griffin-Lim algorithm struct initialization
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  The function can work inplace i.e. \a cinit == \a c
 *
 *  \param[in]   cinit   Initial set of coefficients, size M2 x N x W or NULL
 *  \param[in]       g   Analysis window, size gl x 1
 *  \param[in]       L   Signal length
 *  \param[in]      gl   Window length
 *  \param[in]       W   Number of signal channels
 *  \param[in]       a   Time hop factor
 *  \param[in]       M   Number of frequency channels
 *  \param[in]   alpha   Acceleration parameter
 *  \param[in]  params   Optional parameters
 *  \param[out]      c   Coefficients with reconstructed phase, size M2 x N x W, cannot be NULL
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_init_d(const ltfat_complex_d cinit[], const double g[],
 *                       ltfatInt L, ltfatInt gl, ltfatInt W, ltfatInt a, ltfatInt M,
 *                       double alpha, ltfat_complex_d c[], phaseret_legla_init_params* params,
 *                       phaseret_legla_plan_d** p);
 *
 * phaseret_legla_init_s(const ltfat_complex_s cinit[], const double g[],
 *                       ltfatInt L, ltfatInt gl, ltfatInt W, ltfatInt a, ltfatInt M,
 *                       double alpha, ltfat_complex_s c[], phaseret_legla_init_params* params,
 *                       phaseret_legla_plan_s** p);
 * </tt>
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a g or \a c or \a pout was NULL
 *  LTFATERR_BADARG       | \a alpha must be greater or equal to 0.0
 *  LTFATERR_NOTINRAGE    | \a params->relthr must be in range [0-1]
 *  LTFATERR_BADSIZE      | \a Invalid kernel size: params->ksize
 *  LTFATERR_CANNOTHAPPEN | \a params was not inilialized with phaseret_legla_init_params_defaults
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 */
PHASERET_API int
PHASERET_NAME(legla_init)(const LTFAT_COMPLEX cinit[], const LTFAT_REAL g[], const int L,
                          const int gl, const int W, const int a, const int M,
                          const double alpha, LTFAT_COMPLEX c[],
                          phaseret_legla_init_params* params, PHASERET_NAME(legla_plan)** p);



/** Execute LEGLA plan
 *
 * \param[in]      p   LEGLA plan
 * \param[in]   iter   Number of iterations
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_execute_d(phaseret_legla_plan_d* p, ltfatInt iter);
 *
 * phaseret_legla_execute_s(phaseret_legla_plan_s* p, ltfatInt iter);
 * </tt>
 * \returns
 * Status code          |  Description
 * ---------------------|-------------------
 * LTFATERR_SUCCESS     | No error occurred
 * LTFATERR_NULLPOINTER | \a p or \a p->cinit was NULL.
 * LTFATERR_NOTPOSARG   | \a iter must be positive
 * LTFATERR_BADARG      | \a alpha set in the status callback must be nonnegative
 * LTFATERR_NOMEM       | Memory allocation failed
 * any                  | Status code from any of the callbacks
 */
PHASERET_API int
PHASERET_NAME(legla_execute)(PHASERET_NAME(legla_plan)* p, const int iter);

/** Execute LEGLA plan on new arrays
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]       p   LEGLA plan
 * \param[in]   cinit   Initial coefficient array, size M2 x N x W
 * \param[in]    iter   Number of iterations
 * \param[in]       c   Coefficients with reconstructed phase, size M2 x N x W
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_execute_newarray_d(phaseret_legla_plan_d* p,
 *                                   const ltfat_complex_d cinit[],
 *                                   ltfatInt iter, ltfat_complex_d c[]);
 *
 * phaseret_legla_execute_newarray_s(phaseret_legla_plan_s* p,
 *                                   const ltfat_complex_s cinit[],
 *                                   ltfatInt iter, ltfat_complex_s c[]);
 * </tt>
 * \returns
 * Status code          |  Description
 * ---------------------|-------------------
 * LTFATERR_SUCCESS     | No error occurred
 * LTFATERR_NULLPOINTER | \a p or \a cinit or \a c was NULL.
 * LTFATERR_NOTPOSARG   | \a iter must be positive
 * LTFATERR_BADARG      | \a alpha set in the status callback must be nonnegative
 * LTFATERR_NOMEM       | Memory allocation failed
 * any                  | Status code from any of the callbacks
 */
PHASERET_API int
PHASERET_NAME(legla_execute_newarray)(PHASERET_NAME(legla_plan)* p,
                                      const LTFAT_COMPLEX cinit[],
                                      const int iter, LTFAT_COMPLEX c[]);

/** Delete LEGLA plan
 *
 * \param[in]  p  LEGLA plan
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_done_d(phaseret_legla_plan_d** p);
 *
 * phaseret_legla_done_s(phaseret_legla_plan_s** p);
 * </tt>
 * \returns
 * Status code          |  Description
 * ---------------------|-------------------
 * LTFATERR_SUCCESS     | No error occurred
 * LTFATERR_NULLPOINTER | \a p or \a *p was NULL
 */
PHASERET_API int
PHASERET_NAME(legla_done)(PHASERET_NAME(legla_plan)** p);

/** Register status callback
 *
 *  \param[in]         p   LEGLA plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_set_status_callback_d(phaseret_legla_plan_d* p,
 *                                      phaseret_legla_callback_status_d* callback,
 *                                      void* userdata);
 *
 * phaseret_legla_set_status_callback_s(phaseret_legla_plan_s* p,
 *                                      phaseret_legla_callback_status_s* callback,
 *                                      void* userdata);
 * </tt>
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
PHASERET_API int
PHASERET_NAME(legla_set_status_callback)(PHASERET_NAME(legla_plan)* p,
        PHASERET_NAME(legla_callback_status)* callback, void* userdata);

/** Register coefficient modification callback
 *
 *  \param[in]         p   LEGLA plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *
 * #### Versions #
 * <tt>
 * phaseret_legla_set_cmod_callback_d(phaseret_legla_plan_d* p,
 *                                    phaseret_legla_callback_cmod_d* callback,
 *                                    void* userdata);
 *
 * phaseret_legla_set_cmod_callback_s(phaseret_legla_plan_s* p,
 *                                    phaseret_legla_callback_cmod_s* callback,
 *                                    void* userdata);
 * </tt>
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
PHASERET_API int
PHASERET_NAME(legla_set_cmod_callback)(PHASERET_NAME(legla_plan)* p,
                                       PHASERET_NAME(legla_callback_cmod)* callback,
                                       void* userdata);

/** @}*/

/* Single iteration  */
int
PHASERET_NAME(leglaupdate_init)(const LTFAT_COMPLEX kern[], phaseret_size ksize,
                                int L, int W, int a, int M, int flags,
                                PHASERET_NAME(leglaupdate_plan)** pout);

extern void
PHASERET_NAME(leglaupdate_execute)(PHASERET_NAME(leglaupdate_plan)* plan,
                                   const LTFAT_REAL s[], LTFAT_COMPLEX c[],
                                   LTFAT_COMPLEX cout[]);

void
PHASERET_NAME(leglaupdate_done)(PHASERET_NAME(leglaupdate_plan)** plan);

/* Single col update */
int
PHASERET_NAME(leglaupdate_col_init)(int M, phaseret_size ksize, int flags,
                                    PHASERET_NAME(leglaupdate_plan_col)** pout);

int
PHASERET_NAME(leglaupdate_col_done)( PHASERET_NAME(leglaupdate_plan_col)** p);

void
PHASERET_NAME(leglaupdate_col_execute)(
    PHASERET_NAME(leglaupdate_plan_col)* plan,
    const LTFAT_REAL sCol[],
    const LTFAT_COMPLEX actK[],
    LTFAT_COMPLEX cColFirst[],
    LTFAT_COMPLEX coutrCol[]);

/* Utils */
void
PHASERET_NAME(extendborders)(PHASERET_NAME(leglaupdate_plan_col)* plan,
                             const LTFAT_COMPLEX c[], int N, LTFAT_COMPLEX buf[]);

int
PHASERET_NAME(legla_big2small_kernel)(LTFAT_COMPLEX* bigc, phaseret_size bigsize,
                                      phaseret_size smallsize, LTFAT_COMPLEX* smallc);

int
PHASERET_NAME(legla_findkernelsize)(LTFAT_COMPLEX* bigc, phaseret_size bigsize,
                                    double relthr, phaseret_size* ksize);

/* Modulate kernel */
void
PHASERET_NAME(kernphasefi)(const LTFAT_COMPLEX kern[], phaseret_size ksize,
                           int n, int a, int M, LTFAT_COMPLEX kernmod[]);

/* Format kernel */
void
PHASERET_NAME(formatkernel)(LTFAT_REAL* kernr, LTFAT_REAL* kerni,
                            int kernh, int kernw,
                            int kernwskip, LTFAT_REAL* kernmodr, LTFAT_REAL* kernmodi);



#ifdef __cplusplus
}
#endif
