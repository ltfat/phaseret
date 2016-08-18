#include "config.h"
#include "dgtrealwrapper.h"
#ifndef NOSYSTEMHEADERS
#include <stdlib.h>
#endif


typedef struct legla_plan legla_plan;
typedef struct leglaupdate_plan leglaupdate_plan;
typedef struct leglaupdate_plan_col leglaupdate_plan_col;

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
 *  \returns
 *  Status code | Meaning
 *  ------------|-----------------
 *   0          | Signalizes that callback exited without error
 *  <0          | Callback exited with error
 */
typedef int legla_callback_cmod(void* userdata, complex double c[], int L, int W, int a, int M);

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
 *  \returns
 *  Status code | Meaning
 *  ------------|---------------------------------------------------------------------------------------
 *   0          | Signalizes that callback exited without error
 *  >0          | Signalizes that callback exited without error but terminate the algorithm prematurely
 *  <0          | Callback exited with error
 *
 *  \see dgtreal_anasyn_execute_proj dgtreal_anasyn_execute_ana dgtreal_anasyn_execute_syn
 */
typedef int legla_callback_status(dgtreal_anasyn_plan* p, void* userdata, complex double c[],
                                  int L, int W, int a, int M, double* alpha, int iter);

typedef struct
{
    int height;
    int width;
} phaseret_size;

typedef struct
{
    int y;
    int x;
} phaseret_point;

typedef struct
{
    double relthr; ///< Relative threshold for automatic determination of kernel size, default 1e-3
    phaseret_size ksize; ///< Maximum allowed kernel size (default 2*ceil(M/a) -1) or kernel size directly if relthr==0.0
    dgtreal_anasyn_hint hint; ///< Hint for the alalysis-synthesis functions, default dgtreal_anasyn_auto
    unsigned dgtrealflags; ///< FFTW flags, default FFTW_ESTIMATE
    unsigned leglaflags; ///< LEGLA algorithm flags, default MOD_COEFFICIENTWISE | MOD_MODIFIEDUPDATE
    unsigned private_hash_do_not_use; ///< Checked in legla_init to avoid passing uninitialized struct
} legla_init_params;

typedef enum
{
    MOD_STEPWISE = (1 << 0),  // << DEFAULT
    MOD_FRAMEWISE = (1 << 1),
    MOD_COEFFICIENTWISE = (1 << 2),
    MOD_COEFFICIENTWISE_SORTED = (1 << 3),
    MOD_MODIFIEDUPDATE = (1 << 4),
} leglaupdate_mod;




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
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a cinit or \a g or \a cout was NULL
 *  LTFATERR_NOTPOSARG    | \a iter must be positive
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 */
int
legla(const complex double cinit[], const double g[], const int L, const int gl,
      const int W, const int a, const int M, const int iter, complex double c[]);

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
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a g or \a c or \a pout was NULL
 *  LTFATERR_BADARG       | \a alpha must be greater or equal to 0.0
 *  LTFATERR_NOTINRAGE    | \a params->relthr must be in range [0-1]
 *  LTFATERR_BADSIZE      | \a Invalid kernel size: params->ksize
 *  LTFATERR_CANNOTHAPPEN | \a params was not inilialized with legla_init_params_defaults
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 */
int
legla_init(const complex double cinit[], const double g[], const int L,
           const int gl, const int W, const int a, const int M,
           const double alpha, complex double c[],
           legla_init_params* params, legla_plan** pout);

/** Initialize parameter struct for legla_init
 *
 * \param[in,out]  params  Pointer to a structure to be inilialized
 * \returns
 * Status code          |  Description
 * ---------------------|-------------------
 * LTFATERR_SUCCESS     | No error occurred
 * LTFATERR_NULLPOINTER | \a params was NULL
 */
int
legla_init_params_defaults(legla_init_params* params);


/** Execute LEGLA plan
 *
 * \param[in]      p   LEGLA plan
 * \param[in]   iter   Number of iterations
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
int
legla_execute(legla_plan* p, const int iter);

/** Execute LEGLA plan on new arrays
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]       p   LEGLA plan
 * \param[in]   cinit   Initial coefficient array, size M2 x N x W
 * \param[in]    iter   Number of iterations
 * \param[in]       c   Coefficients with reconstructed phase, size M2 x N x W
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
int
legla_execute_newarray(legla_plan* p, const complex double cinit[], const int iter, complex double c[]);

/** Delete LEGLA plan
 *
 * \param[in]  p  LEGLA plan
 * \returns
 * Status code          |  Description
 * ---------------------|-------------------
 * LTFATERR_SUCCESS     | No error occurred
 * LTFATERR_NULLPOINTER | \a p or \a *p was NULL
 */
int
legla_done(legla_plan** p);

/** Register status callback
 *
 *  \param[in]         p   LEGLA plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
int
legla_set_status_callback(legla_plan* p, legla_callback_status* callback, void* userdata);

/** Register coefficient modification callback
 *
 *  \param[in]         p   LEGLA plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
int
legla_set_cmod_callback(legla_plan* p, legla_callback_cmod* callback, void* userdata);

/** @}*/

/* Single iteration  */
int
leglaupdate_init(const complex double kern[], phaseret_size ksize,
                 int L, int W, int a, int M, int flags, leglaupdate_plan** pout);

extern void
leglaupdate_execute(leglaupdate_plan* plan, const double s[], complex double c[],
                    complex double cout[]);

void
leglaupdate_done(leglaupdate_plan** plan);

/* Single col update */
int
leglaupdate_init_col(int M, phaseret_size ksize,
                     int flags, leglaupdate_plan_col** pout);

void
leglaupdatereal_execute_col(leglaupdate_plan_col* plan,
                            const double sCol[],
                            const complex double actK[],
                            complex double cColFirst[],
                            complex double coutrCol[]);

/* Utils */
void
extendborders(leglaupdate_plan_col* plan, const complex double c[], int N,
              complex double buf[]);

int
legla_big2small_kernel(complex double* bigc, phaseret_size bigsize,
                       phaseret_size smallsize, complex double* smallc);

int
legla_findkernelsize(complex double* bigc, phaseret_size bigsize,
                     double relthr, phaseret_size* ksize);

/* Modulate kernel */
void
kernphasefi(const complex double kern[], phaseret_size ksize,
            int n, int a, int M, complex double kernmod[]);

/* Format kernel */
void
formatkernel(double* kernr, double* kerni,
             int kernh, int kernw,
             int kernwskip, double* kernmodr, double* kernmodi);

/* Util */
int phaseret_gcd(int m, int n);
int phaseret_lcm(int m, int n);
