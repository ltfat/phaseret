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
#include "utils.h"
#include <fftw3.h>


#ifdef __cplusplus
extern "C" {
#endif

// Not storing state between calls
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
 *  \param[in]     frames   N frames M samples long
 *  \param[in]     g        Analysis window
 *  \param[in]     n        Which frame to overlap
 *  \param[in]     N        Number of frames
 */
void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const double* frames, const double* g, int n, int N);

/** Phase update of a frame.
 * \param[in,out] p         RTISILA Update Plan, p.fftframe contains
 * \param[in]     sframe    Target spectrum magnitude
 * \param[out]    frameupd  Updated frame
 */
void
rtisilaphaseupdate(rtisilaupdate_plan* p, const double* sframe, double* frameupd);

/** Create a RTISILA Update Plan.
 * \param[in]     g          Analysis window
 * \param[in]     specg1     Analysis window used in the first iteration
 *                           for the newest lookahead frame
 * \param[in]     specg2     Analysis window used in the other iterations
 *                           for the newest lookahead frame
 * \param[in]     gd         Synthesis window
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 *                           (possibly zero-padded).
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
 * \param[in]     frames     N frames M samples long
 * \param[in]     N          Number of frames
 * \param[in]     s          N frames M2 samples long
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of iterations
 * \param[out]    frames2    N output frames M samples long
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
 * \param[in]     frames     N frames M samples long
 * \param[in]     g          Analysis window
 * \param[in]     specg1     Analysis window used in the first iteration
 *                           for the newest lookahead frame
 * \param[in]     specg2     Analysis window used in the other iterations
 *                           for the newest lookahead frame
 * \param[in]     gd         Synthesis window
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 *                           (possibly zero-padded).
 * \param[in]     N          Number of frames N = lookback + 1 + lookahead
 * \param[in]     s          Target magnitude, N frames M samples long
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of iterations
 * \param[out]    frames2    N output frames M samples long
 */
void
rtisilaupdate(const double* frames,
              const double* g, const double* specg1, const double* specg2, const double* gd,
              int a, int M, int N, const double* s, int lookahead, int maxit,
              double* frames2);

/** Do RTISI-LA for a complete magnitude spectrogram and compensate delay
 *
 * This function just creates a plan, executes it for each col in s and c 
 * and destroys it.
 *
 * \param[in]     s          Magnitude spectrogram  M2 x N array
 * \param[in]     g          Analysis window
 * \param[in]     specg1     Analysis window used in the first iteration
 *                           for the newest lookahead frame
 * \param[in]     specg2     Analysis window used in the other iterations
 *                           for the newest lookahead frame
 * \param[in]     gd         Synthesis window
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 *                           (possibly zero-padded).
 * \param[in]     L          Transform length
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of iterations
 * \param[out]    c          Reconstructed coefficients M2 x N array
 */
void
rtisilaoffline(const double* s,
               const double* g, const double* specg1, const double* specg2, const double* gd,
               int a, int M, int L, int lookahead, int maxit, complex double* c);

/** Plan for rtisila
 *  
 *  Serves for storing state between calls to rtisila_execute.
 */
typedef struct
{
    int maxLookahead;
    int lookahead;
    int lookback;
    int noFrames;
    int maxit;
    double* frames; //!< Buffer for time-domain frames
    double* s; //!< Buffer for target magnitude
    rtisilaupdate_plan* updateplan;
    void** garbageBin;
    int garbageBinSize;
} rtisila_plan;


/** Create a RTISILA Plan.
 * \param[in]     g            Analysis window
 * \param[in]     specg1       Analysis window used in the first iteration
 *                             for the newest lookahead frame
 * \param[in]     specg2       Analysis window used in the other iterations
 *                             for the newest lookahead frame
 * \param[in]     gd           Synthesis window
 * \param[in]     a            Hop size
 * \param[in]     M            FFT length, also length of all the windows
 *                             (possibly zero-padded).
 * \param[in]     lookahead    Number of lookahead frames
 * \param[in]     maxLookahead Maximum number of lookahead frames
 * \param[in]     maxit        Number of iterations. The number of per-frame 
 *                             iterations is (lookahead+1) * maxit. 
 * \returns RTISILA Plan
 */
rtisila_plan*
rtisila_init(const double *g, const double* specg1,
             const double* specg2, const double* gd,
             int a, int M, int lookahead, int maxLookahead, int maxit);

/** Create a RTISILA Plan from a window.
 * \param[in]     g            Analysis window
 * \param[in]     gl           Window length
 * \param[in]     a            Hop size
 * \param[in]     M            FFT length, also length of all the windows
 *                             (possibly zero-padded).
 * \param[in]     lookahead    Number of lookahead frames
 * \param[in]     maxLookahead Maximum number of lookahead frames
 * \param[in]     maxit        Number of iterations. The number of per-frame 
 *                             iterations is (lookahead+1) * maxit. 
 * \returns RTISILA Plan
 */
rtisila_plan*
rtisila_wininit(LTFAT_FIRWIN win, int gl, int a, int M,
                int lookahead, int maxLookahead, int maxit);

/** Destroy a RTISILA Plan.
 * \param[in] p  RTISILA Plan
 */
void
rtisila_done(rtisila_plan* p);


/** Execute RTISILA plan for a single frame
 *  
 *  The function is intedned to be called for consecutive stream of frames
 *  as it reuses some data from the previous frames stored in the plan.
 *
 *  c is lagging behind s by lookahead frames.
 *
 * \param[in]       p   RTISILA plan   
 * \param[in]       s   Target magnitude
 * \param[out]      c   Reconstructed coefficients
 */
void
rtisila_execute(rtisila_plan* p, const double* s, double complex* c);


#ifdef __cplusplus
}
#endif

#endif /* _rtisila_h */

/** @}*/
