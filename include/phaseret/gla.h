#include "dgtrealwrapper.h"


typedef struct gla_plan gla_plan;

/** \addtogroup gla
 * @{
 *
 */

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
typedef int gla_callback_status(dgtreal_anasyn_plan* p, void* userdata, complex double c[],
                                int L, int W, int a, int M, double* alpha, int iter);

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
typedef int gla_callback_cmod(void* userdata, complex double c[], int L, int W, int a, int M);

/** Function prototype for signal modification callback
 *
 *  \param[in]  userdata   User defined data
 *  \param[in,out]     f   Time-domain signal to be updated, size L x W
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
typedef int gla_callback_fmod(void* userdata, double f[], int L, int W, int a, int M);


/** Griffin-Lim algorithm
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
 *  \param[out]  cout   Coefficients with reconstructed phase, size M2 x N x W
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a cinit or \a g or \a cout was NULL
 *  LTFATERR_BADSIZE      | Signal length L is less or equal to 0.
 *  LTFATERR_NOTPOSARG    | At least one of \f W, \f a, \f M, \f gl was less or equal to zero. 
 *  LTFATERR_BADTRALEN    | \a L is not divisible by both \a a and \a M. 
 *  LTFATERR_NOTAFRAME    | System does not form a frame
 *  LTFATERR_INITFAILED   | The FFTW plan creation failed 
 *  LTFATERR_NOTSUPPORTED | This is a non-painless system but its support was not compiled
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 */
int
gla(const complex double cinit[], const double g[], const int L, const int gl, const int W,
    const int a, const int M, const int iter, complex double cout[]);

/** Initialize Griffin-Lim algorithm plan
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  \note \a cinit can be NULL if the plan is intended to be used wit the _newarray
 *  execute function. Similarly, \a c can also be NULL only if FFTW_ESTIMATE is passed in flags
 *  (the FFTW planning routine does not touch the array).
 *  On the other hand, the content of c might get overwritten if other FFTW planning flags are used.
 *
 *  \note In-place mode i.e. \a cinit == \a c is allowed.
 *
 *  \param[in]  cinit   Initial set of coefficients, size M2 x N x W or NULL
 *  \param[in]      g   Analysis window, size gl x 1
 *  \param[in]      L   Signal length
 *  \param[in]     gl   Window length
 *  \param[in]      W   Number of signal channels
 *  \param[in]      a   Time hop factor
 *  \param[in]      M   Number of frequency channels
 *  \param[in]  alpha   Acceleration constant
 *  \param[in]      c   Array for holding coefficients with reconstructed phase, size M2 x N x W or NULL if flags == FFTW_ESTIMATE
 *  \param[in]   hint   DGT algorithm hint
 *  \param[in]  flags   FFTW planning flag
 *  \param[out]     p   GLA Plan
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a cinit or \a g or or \a p was NULL or \a c was NULL and flags != FFTW_ESTIMATE
 *  LTFATERR_BADSIZE      | Signal length L is less or equal to 0.
 *  LTFATERR_NOTPOSARG    | At least one of \f W, \f a, \f M, \f gl was less or equal to zero. 
 *  LTFATERR_BADTRALEN    | \a L is not divisible by both \a a and \a M. 
 *  LTFATERR_NOTAFRAME    | System does not form a frame
 *  LTFATERR_INITFAILED   | The FFTW plan creation failed 
 *  LTFATERR_NOTSUPPORTED | This is a non-painless system but its support was not compiled
 *  LTFATERR_BADARG       | \a alpha was set to a negative number 
 *  LTFATERR_CANNOTHAPPEN | \a hint does not have a valid value from \a dgtreal_anasyn_hint or \a ptype is not valid value from \a ltfat_phaseconvention enum
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 */
int
gla_init(const complex double cinit[], const double g[], const int L, const int gl, const int W,
         const int a, const int M, const double alpha,
         complex double c[], dgtreal_anasyn_hint hint, unsigned flags,
         gla_plan** p);

/** Execute Griffin-Lim algorithm plan
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  \param[in]      p   Griffin-lim algorithm plan
 *  \patam[in]   iter   Number of iterations
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p was NULL or the plan was created with \a cinit or \a cout being NULL
 *  LTFATERR_NOTPOSARG    | \a iter was not positive
 *  LTFATERR_BADARG       | \a alpha was set to a negative number in the status callback 
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 *  any                   | Error code from some of the callbacks
 */
int
gla_execute(gla_plan* p, const int iter);

/** Execute Griffin-Lim algorithm plan on a new array
 *
 *  M2 = M/2 + 1, N = L/a
 *
 *  \param[in]      p   Griffin-lim algorithm plan
 *  \param[in]  cinit   Initial set of coefficients, size M2 x N x W
 *  \patam[in]   iter   Number of iterations
 *  \param[in]   cout   Coefficients with reconstructed phase, size M2 x N x W
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a cinit or \a cout was NULL
 *  LTFATERR_NOTPOSARG    | \a iter was not positive
 *  LTFATERR_BADARG       | \a alpha was set to a negative number in the status callback 
 *  LTFATERR_NOMEM        | Memory allocation error occurred
 *  any                   | Error code from some of the callbacks
 */
int
gla_execute_newarray(gla_plan* p, const complex double cinit[], const int iter,
                     complex double cout[]);

/** Destroy Griffin-Lim algorithm plan
 *
 *  \param[in]      p   Griffin-lim algorithm plan
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a *p was NULL
 */
int
gla_done(gla_plan** p);

/** Register status callback
 *
 *  \param[in]         p   Griffin-lim algorithm plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
int
gla_set_status_callback(gla_plan* p, gla_callback_status* callback, void* userdata);

/** Register coefficient modification callback
 *
 *  \param[in]         p   Griffin-lim algorithm plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
int
gla_set_cmod_callback(gla_plan* p, gla_callback_cmod* callback, void* userdata);

/** Register signal modification callback
 *
 *  \param[in]         p   Griffin-lim algorithm plan
 *  \param[in]  callback   Callback function
 *  \param[in]  userdata   User defined data
 *  \returns
 *  Status code           | Description
 *  ----------------------|-----------------------
 *  LTFATERR_SUCCESS      | No error occurred
 *  LTFATERR_NULLPOINTER  | \a p or \a callback was NULL
 */
int
gla_set_fmod_callback(gla_plan* p, gla_callback_fmod* callback, void* userdata);

/** @} */

int
fastupdate(complex double* c, complex double* t, double alpha, int L);
