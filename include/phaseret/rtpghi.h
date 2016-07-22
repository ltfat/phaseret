#ifndef _rtpghi_h
#define _rtpghi_h

// #ifndef NOSYSTEMHEADERS
// #include <ltfat.h>
// #endif

#ifdef __cplusplus
extern "C" {
#endif

/** Plan for rtpghi
 * 
 * Serves for storing state between calls to rtpghi_execute.
 *
 */
typedef struct rtpghi_state rtpghi_state;

/** \addtogroup rtpghi
 *  @{
 */

/** Create a RTPGHI Plan.
 * \param[in]     gamma        Window-specific constant Cg*gl^2
 * \param[in]     a            Hop size
 * \param[in]     M            FFT length, also length of all the windows
 *                             (possibly zero-padded).
 * \param[in]     tol          Relative coefficient tolerance.
 * \param[in]     do_causal    Zero delay (1) or 1 frame delay version of the alg. 
 * \returns RTPGHI Plan
 */
int
rtpghi_init(double gamma, int W, int a, int M, double tol, int do_causal, rtpghi_state** pout);

/** Change the version of the algorithm
 *  
 * Either to one-frame-delay version (do_causal==0) or to the 
 * no-delay version (do_causal anything else). 
 *
 * \note This is not thread safe
 *
 * \param[in] p         RTPGHI plan
 * \param[in] do_causal Causal flag
 * \returns Status code
 */
int
rtpghi_set_causal(rtpghi_state* p, int do_causal);

/** Execute RTPGHI plan for a single frame
 *  
 *  The function is intedned to be called for consecutive stream of frames
 *  as it reuses some data from the previous frames stored in the plan.
 *
 *  if do_causal is enebled, c is not lagging, else c is lagging by one 
 *  frame.
 *
 * \param[in]       p   RTPGHI plan   
 * \param[in]       s   Target magnitude
 * \param[out]      c   Reconstructed coefficients
 */
int
rtpghi_execute(rtpghi_state* p, const double s[], complex double c[]);

/** Destroy a RTPGHI Plan.
 * \param[in] p  RTPGHI Plan
 */
int
rtpghi_done(rtpghi_state** p);

/** Do RTPGHI for a complete magnitude spectrogram and compensate delay
 *
 * This function just creates a plan, executes it for each col in s and c 
 * and destroys it.
 *
 * \param[in]     s          Magnitude spectrogram  M2 x N array
 * \param[in]     gamma      Window-specific constant Cg*gl^2
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 * \param[in]     L          Transform length
 *                           (possibly zero-padded).
 * \param[in]     tol        Relative coefficient tolerance.
 * \param[in]     do_causal  Zero delay (1) or 1 frame delay version of the alg.
 * \param[out]    c          Reconstructed coefficients M2 x N array
 */
int
rtpghioffline(const double s[], double gamma, int a, int M, int L, double tol,
              int do_causal, complex double c[]);
/** @}*/

/** Compute phase frequency gradient by differentiation in time
 *
 * \param[in]     logs       Log-magnitude of a 3 x M2 buffer
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 * \param[in]     gamma      Window-specific constant Cg*gl^2
 * \param[in]     do_causal  If true, fgrad is relevant for 3rd buffer col, else it
 *                           is relevant for 2nd buffer. 
 * \param[out]    fgrad      Frequency gradient, array of length M2
 */
void
rtpghifgrad(const double logs[], int a, int M, double gamma,
            int do_causal, double fgrad[]);

/** Compute phase time gradient by differentiation in frequency
 *
 * \param[in]     logs       Log-magnitude, array of length M2
 * \param[in]     a          Hop size
 * \param[in]     M          FFT length, also length of all the windows
 * \param[in]     gamma      Window-specific constant Cg*gl^2
 * \param[out]    tgrad      Time gradient, array of length M2
 */
void
rtpghitgrad(const double logs[], int a, int M, double gamma,
            double tgrad[]);

/** Compute log of input
 * \param[in]   in  Input array of length L
 * \param[in]    L  Length of the arrays
 * \param[out] out  Output array of length L  
 */
void
rtpghilog(const double in[], int L, double out[]);

/** Combine magnitude and phase to a complex array
 * \param[in]        s      Magnitude, array of length L
 * \param[in]    phase      Phase in rad, array of length L
 * \param[in]        L      Length of the arrays
 * \param[out]       c      Output array of length L  
 */
void
rtpghimagphase(const double s[], const double phase[], int L, complex double c[]);

#ifdef __cplusplus
}
#endif


#endif /* _rtpghi_h */

