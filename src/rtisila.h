/** \addtogroup rtisila
 *  @{
 *
 * @file rtisila.h
 * @author Zdeněk Průša
 * @date 1 Feb 2016
 * @brief RTISI-LA header
 *
 */

#ifndef _rtisila_h
#define _rtisila_h
// complex.h must be before fftw3.h. Then fftw_complex is complex double
#include "config.h"
#include <fftw3.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double* frame; //!< Time domain buffer
    fftw_complex* fftframe; //!< Frequency domain buffer
    fftw_plan fwdp; //!< Real FFT plan
    fftw_plan backp; //!< Real IFFT plan
    const double* g;
    const double* gd;
    const double* specg1;
    const double* specg2;
    const int M;
    const int a;
} rtisilaupdate_plan;

/** Overlay frames to get n-th frame.
 *  \param[in,out] p        RTISILA Update Plan, p.frame contains the overlaid frame
 *  \param[in]     frames  N frames M samples long
 *  \param[in]     g     Analysis window
 *  \param[in]     n        Which frame to overlap
 *  \param[in]     N        Number of frames
 */
void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const double* frames, const double* g, int n, int N);

/** Phase update of a frame.
 * \param[in,out] p      RTISILA Update Plan, p.fftframe contains
 * \param[in]     sframe     Target spectrum magnitude
 * \param[out]    frameupd  Updated frame
 */
void
rtisilaphaseupdate(rtisilaupdate_plan* p, const double* sframe, double* frameupd);

/** Create a RTISILA Update Plan.
 * \param[in,out] p      RTISILA Update Plan, p.fftframe contains
 * \param[in]     sframe     Target spectrum magnitude
 * \param[out]    frameupd  Updated frame
 * \returns RTISILA Update Plan
 */
rtisilaupdate_plan*
rtisilaupdate_init(const double *g, const double* specg1,
                   const double* specg2, const double* gd,
                   int a, int M);

/** Destroy a RTISILA Update Plan.
 * \param[in] p  RTISILA Update Plan
 */
void
rtisilaupdate_done(rtisilaupdate_plan* p);

/** Do maxit iterations of RTISI-LA for a single frame
 *
 * N = lookback + 1 + lookahead
 * M2 = M/2 + 1
 *
 *
 * <em>Note the function can be run inplace i.e. frames and frames2 can
 * point to the same memory location.</em>
 *
 * \param[in,out] p          RTISILA Update Plan
 * \param[in]     frames    N frames M samples long
 * \param[in]     N          Number of frames
 * \param[in]     s    N frames M2 samples long
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of iterations
 * \param[out]    frames2   N output frames M samples long
 */
void
rtisilaupdate_execute(rtisilaupdate_plan* p, const double* frames, int N,
                      const double* s, int lookahead, int maxit, double* frames2);

/** Do maxit iterations of RTISI-LA for a single frame
 *
 * This function just creates a plan, executes it and destroys it.
 *
 * <em>Note the function can be run inplace i.e. frames and frames2 can
 * point to the same memory location.</em>
 *
 * \param[in]     frames    N frames M samples long
 * \param[in]     g       Analysis window
 * \param[in]     specg1     Analysis window used in the first iteration
 *                           for the newest lookahead frame
 * \param[in]     specg2     Analysis window used in the other iterations
 *                           for the newest lookahead frame
 * \param[in]     gd      Synthesis window
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 *                           (possibly zero middlepadded).
 * \param[in]     N          Number of frames N = lookback + 1 + lookahead
 * \param[in]     s    Target magnitude, N frames M samples long
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of iterations
 * \param[out]    frames2   N output frames M samples long
 */
void
rtisilaupdate(const double* frames,
              const double* g, const double* specg1, const double* specg2, const double* gd,
              int a, int M, int N, const double* s, int lookahead, int maxit,
              double* frames2);

/** Do maxit iterations of RTISI-LA for a single frame
 *
 * This function just creates a plan, executes it and destroys it.
 *
 * <em>Note the function can be run inplace i.e. frames and frames2 can
 * point to the same memory location.</em>
 *
 * \param[in]     c          N frames M samples long
 * \param[in]     g       Analysis window
 * \param[in]     specg1     Analysis window used in the first iteration
 *                           for the newest lookahead frame
 * \param[in]     specg2     Analysis window used in the other iterations
 *                           for the newest lookahead frame
 * \param[in]     gd      Synthesis window
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 *                           (possibly zero middlepadded).
 * \param[in]     L          Transform length
 * \param[in]     s    Target magnitude, N frames M samples long
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of iterations
 * \param[out]    frames2   N output frames M samples long
 */
void
rtisila(const double* c,
        const double* g, const double* specg1, const double* specg2, const double* gd,
        int a, int M, int L, const double* s, int lookahead, int maxit, double* c2);



/** ifftshift done in the frequency domain
 *
 * The function modulates the complex spectra such that it effectivelly
 * does ifftshift in the time domain.
 *
 * <em>Note the function can be run inplace i.e. in and out can
 * point to the same memory location.</em>
 *
 * \param[in]     in   Input array
 * \param[in]     L    Length of the arrays
 * \param[out]    out  Output array
 */
void
fftrealifftshift(const double complex* in, int L, double complex* out);


#ifdef __cplusplus
}
#endif

#endif /* _rtisila_h */

/** @}*/
