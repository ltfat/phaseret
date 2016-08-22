/** \addtogroup rtisila
 *
 * @file rtisila.h
 * @author Zdeněk Průša
 * @date 1 Feb 2016
 * @brief RTISI-LA header
 *
 */

#ifndef _rtisila_h
#define _rtisila_h
#ifndef NOSYSTEMHEADERS
#define LTFAT_DOUBLE
#include "ltfat.h"
#endif
#include "ltfat/types.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct rtisilaupdate_plan rtisilaupdate_plan;

typedef struct rtisila_state rtisila_state;

void
overlaynthframe(const LTFAT_REAL* frames, int gl, int N, int a, int n, LTFAT_REAL* frame);

/** Overlay frames to get n-th frame.
 *  \param[in,out] p        RTISILA Update Plan, p.frame contains the overlaid frame
 *  \param[in]     frames   N frames M samples long
 *  \param[in]     g        Analysis window
 *  \param[in]     n        Which frame to overlap
 *  \param[in]     N        Number of frames
 */
void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const LTFAT_REAL* frames, const LTFAT_REAL* g, int n, int N);

/** Phase update of a frame.
 * \param[in,out] p         RTISILA Update Plan, p.fftframe contains
 * \param[in]     sframe    Target spectrum magnitude
 * \param[out]    frameupd  Updated frame
 */
void
rtisilaphaseupdate(rtisilaupdate_plan* p, const LTFAT_REAL* sframe, LTFAT_REAL* frameupd, LTFAT_COMPLEX* c);

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
int
rtisilaupdate_init(const LTFAT_REAL *g, const LTFAT_REAL* specg1,
                   const LTFAT_REAL* specg2, const LTFAT_REAL* gd,
                   const int gl, int a, int M,
                   rtisilaupdate_plan** p);

/** Destroy a RTISILA Update Plan.
 * \param[in] p  RTISILA Update Plan
 */
int
rtisilaupdate_done(rtisilaupdate_plan** p);

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
rtisilaupdate_execute(rtisilaupdate_plan* p, const LTFAT_REAL* frames, int N,
                      const LTFAT_REAL* s, int lookahead, int maxit, LTFAT_REAL* frames2,
                      LTFAT_COMPLEX* c);

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
rtisilaupdate(const LTFAT_REAL* frames,
              const LTFAT_REAL* g, const LTFAT_REAL* specg1, const LTFAT_REAL* specg2, const LTFAT_REAL* gd,
              const int gl, int a, int M, int N, const LTFAT_REAL* s, int lookahead, int maxit,
              LTFAT_REAL* frames2);

/** \addtogroup rtisila
 *  @{
 *
 */



/** Create a RTISILA state.
 *
 * \param[in]     g            Analysis window
 * \param[in]     gl           Window length
 * \param[in]     W            Number of signal channels
 * \param[in]     a            Hop size
 * \param[in]     M            Number of frequency channels (FFT length)
 * \param[in]     lookahead    (Maximum) number of lookahead frames
 * \param[in]     maxit        Number of iterations. The number of per-frame
 *                             iterations is (lookahead+1) * maxit.
 * \param[out]    p            RTISILA state 
 * \returns 
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | \a p or \a g was NULL
 * LTFATERR_BADARG          | \a lookahead was a negative number
 * LTFATERR_BADSIZE         | \a gl was not positive
 * LTFATERR_NOTPOSARG       | One of the following was not positive: \a W, \a a, \a M, \a maxit
 * LTFATERR_NOTAFRAME       | System is not a frame. 
 * LTFATERR_NOTPAINLESS 	| System is not painless. 
 * LTFATERR_INITFAILED      | FFTW plan creation failed 
 * LTFATERR_NOMEM           | Indentifies that heap allocation failed 
 */
int
rtisila_init(const LTFAT_REAL g[], const int gl, const int W,
             int a, int M, int lookahead, int maxit,
             rtisila_state** p);

/** Create a RTISILA Plan from a window.
 * \param[in]     win          Analysis window
 * \param[in]     gl           Window length
 * \param[in]     W            Number of signal channels
 * \param[in]     a            Hop size
 * \param[in]     M            Number of frequency channels (FFT length)
 * \param[in]     lookahead    (Maximum) number of lookahead frames
 * \param[in]     maxit        Number of iterations. The number of per-frame
 *                             iterations is (lookahead+1) * maxit.
 * \param[out]    p            RTISILA state 
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_CANNOTHAPPEN    | \a win is not a valid value from the \a LTFAT_FIRWIN enum
 * LTFATERR_NULLPOINTER     | \a p or \a g was NULL
 * LTFATERR_BADARG          | \a lookahead was a negative number
 * LTFATERR_BADSIZE         | \a gl was not positive
 * LTFATERR_NOTPOSARG       | One of the following was not positive: \a W, \a a, \a M, \a maxit
 * LTFATERR_NOTAFRAME       | System is not a frame. 
 * LTFATERR_NOTPAINLESS 	| System is not painless. 
 * LTFATERR_INITFAILED      | FFTW plan creation failed 
 * LTFATERR_NOMEM           | Indentifies that heap allocation failed 
 */
int
rtisila_init_win(LTFAT_FIRWIN win, int gl, int W, int a, int M,
                 int lookahead, int maxit, rtisila_state** p);

/** Change number of lookahead frames
 *
 * The number of frames can only be less or equal to the number of lookahead frames
 * specified in the init function.  
 *
 * \note This is not thread safe.
 *
 * \param[in] p          RTISILA Plan
 * \param[in] lookahead  Number of lookahead frame
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | \a p was NULL
 * LTFATERR_BADARG          | \a lookahead was a negative number or greater than max lookahead
 */
int
rtisila_set_lookahead(rtisila_state* p, int lookahead);

/** Execute RTISILA plan for a single time frame
 *
 *  The function is intedned to be called for consecutive stream of frames
 *  as it reuses some data from the previous frames stored in the state.
 *
 *  \a c is lagging behind \a s by \a lookahead frames.
 *
 *  M2=M/2+1
 *
 * \param[in]       p   RTISILA plan
 * \param[in]       s   Target magnitude, size M2 x W
 * \param[out]      c   Reconstructed coefficients, size M2 x W
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | At least one of the following was NULL: \a p, \a s and \a c
 */
int
rtisila_execute(rtisila_state* p, const LTFAT_REAL s[], LTFAT_COMPLEX c[]);

/** Reset buffers of rtisila_state
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | \a p or \a *p was NULL.
 */
int
rtisila_reset(rtisila_state* p);

/** Destroy a RTISILA Plan.
 * \param[in] p  RTISILA Plan
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | \a p or \a *p was NULL.
 */
int
rtisila_done(rtisila_state** p);

/** Do RTISI-LA for a complete magnitude spectrogram and compensate delay
 *
 * This function just creates a plan, executes it for each col in s and c
 * and destroys it.
 *
 * M2 = M/2 + 1, N = L/a
 * The total number of per-frame iterations is: maxit x (lookahead + 1)
 *
 * \param[in]     s          Magnitude spectrogram, size M2 x N x W
 * \param[in]     g          Analysis window, size gl x 1
 * \param[in]     L          Transform length
 * \param[in]     gl         Window length
 * \param[in]     W          Number of signal channels
 * \param[in]     a          Hop size
 * \param[in]     M          Number of frequency channels (FFT length)
 * \param[in]     lookahead  Number of lookahead frames
 * \param[in]     maxit      Number of per-frame iterations
 * \param[out]    c          Reconstructed coefficients M2 x N array
 */
int
rtisilaoffline(const LTFAT_REAL s[], const LTFAT_REAL g[],
               int L, int gl, int W, int a, int M,
               int lookahead, int maxit, LTFAT_COMPLEX c[]);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* _rtisila_h */
