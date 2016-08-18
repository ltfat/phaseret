#ifndef _pghi_h
#define _pghi_h

// #ifndef NOSYSTEMHEADERS
// #include <ltfat.h>
// #endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pghi_plan pghi_plan;

/** \addtogroup pghi
 * @{
 */

/** Reconstruct complex coefficients from the magnitude using PGHI algorithm
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]        s  Target magnitude, size M2 x N x W
 * \param[in]    gamma  Window specific constant
 * \param[in]        L  Signal length
 * \param[in]        W  Number of channels
 * \param[in]        a  Hop factor
 * \param[in]        M  Number of frequency channels (FFT length)
 * \param[out]       c  Reconstructed coefficients, size M2 x N x W
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Indicates that at least one of the following was NULL: \a c, \a s
 * LTFATERR_NOTPOSARG       | At least one of the followig was not positive: \a L, \a W, \a a
 * LTFATERR_BADARG          | \a gamma was not positive or it was NAN.
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 *
 * \see firwin2gamma
 */
int
pghi(const double s[], double gamma, const int L, const int W,
     const int a, const int M, complex double c[]);

/** Reconstruct complex coefficients from the magnitude using PGHI algorithm 
 * ... using some known coefficients.
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]      cin  Coefficients with initial phase, size M2 x N X W
 * \param[in]     mask  Mask used to select coefficients with known phase, size  M2 x N X W
 * \param[in]    gamma  Window specific constant
 * \param[in]        L  Signal length
 * \param[in]        W  Number of signal channels
 * \param[in]        a  Hop factor
 * \param[in]        M  Number of frequency channels (FFT length)
 * \param[out]       c  Reconstructed coefficients, size M2 x N x W
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Indicates that at least one of the following was NULL: \a cin, \a c, \a mask
 * LTFATERR_NOTPOSARG       | At least one of the followig was not positive: \a L, \a W, \a a
 * LTFATERR_BADARG          | \a gamma was not positive or it was NAN.
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 *
 * \see firwin2gamma
 */
int
pghi_withmask(const complex double cin[], const int mask[],
              double gamma, const int L, const int W,
              const int a, const int M, complex double c[]);

/** Initialize PGHI plan
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]    gamma  Window specific constant
 * \param[in]        L  Signal length
 * \param[in]        W  Number of channels
 * \param[in]        a  Hop factor
 * \param[in]        M  Number of frequency channels (FFT length)
 * \param[in]     tol1  Relative tolerance for the first pass, must be in range [0-1]
 * \param[in]     tol2  Relative tolernace for the second pass, must be in range [0-1] and
 *                      lower or equal to \a tol1. If \a tol2 is NAN or it is equal to
 *                      \a tol1, only the first pass will be done. 
 * \param[out]       p  PGHI plan
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | \a p was NULL.
 * LTFATERR_NOTPOSARG       | At least one of the followig was not positive: \a L, \a W, \a a
 * LTFATERR_BADARG          | \a gamma was not positive
 * LTFATERR_NOTINRANGE      | \a tol1 and \a tol2 were not in range [0-1] or \a tol1 < \a tol2
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 *
 * \see firwin2gamma
 */
int
pghi_init(double gamma, const int L, const int W,
          const int a, const int M, double tol1, double tol2, pghi_plan** p);

/** Execute PGHI plan
 *
 * M2 = M/2 + 1, N = L/a
 *
 * \param[in]      p  PGHI plan
 * \param[in]      s  Target magnitude of coefficients, size M2 x N X W
 * \param[out]     c  Output coefficients with reconstructed phase, size M2 x N X W
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Indicates that at least one of the following was NULL: \a p, \a c, \a s
 */
int
pghi_execute(pghi_plan* p, const double s[], complex double c[]);

/** Execute PGHI plan with respect to mask
 *
 * M2 = M/2 + 1, N = L/a
 *
 * Nonzero values in \a mask represent known coefficient in \cin
 *
 * \param[in]      p  PGHI plan
 * \param[in]    cin  Coefficients with initial phase, size M2 x N X W
 * \param[in]   mask  Mask used to select coefficients with known phase, size  M2 x N X W
 * \param[in] buffer  Work buffer, size M2 x N. Internal heap allocation occurs if it is NULL
 * \param[out]  cout  Output coefficients with reconstructed phase, size M2 x N X W
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Indicates that at least one of the following was NULL: \a p, \a cin, \a mask, \a cout
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
int
pghi_execute_withmask(pghi_plan* p, const complex double cin[], const int mask[],
                      double buffer[], complex double cout[]);

/** Destroy PGHI plan
 *
 * \param[in]   p  PGHI plan
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | \a p was NULL
 */
int
pghi_done(pghi_plan** p);
/** @} */

void
pghimagphase(const double s[], const double phase[], int L, complex double c[]);

void
pghilog(const double in[], int L, double out[]);

void
pghitgrad(const double logs[], double gamma, int a, int M, int N, double tgrad[]);

void
pghifgrad(const double logs[], double gamma, int a, int M, int N, double fgrad[]);


#ifdef __cplusplus
}
#endif
#endif
