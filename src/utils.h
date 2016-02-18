/** \defgroup utils Utilities
 * \addtogroup utils
 *  \{
 *
 * \file
 * \author Zdeněk Průša
 * \date 1 Feb 2016
 * \brief UTILS header
 *
 */
#ifndef _utils_h
#define _utils_h

#include "config.h"

/** Just fftshift
 *
 *  *The function can work inplace i.e. in==out*
 *  \param[in]   in   Input array
 *  \param[in]   L    Number of elements in the arrays
 *  \param[out]  out  Output array
 */
void
fftshift(const double* in, int L, double* out);


/** Just ifftshift
 *
 *  *The function can work inplace i.e. in==out*
 *  \param[in]   in   Input array
 *  \param[in]   L    Number of elements in the arrays
 *  \param[out]  out  Output array
 */
void
ifftshift(const double* in, int L, double* out);

/** Just circshift
 *
 *  *The function can work inplace i.e. in==out*
 *  \param[in]   in     Input array
 *  \param[in]   L      Number of elements in the arrays
 *  \param[in]   shift  Shift
 *  \param[out]  out    Output array
 */
void
circshift(const double* in, int L, int shift, double* out);

/** fftshift() in the frequency domain by modulation
 *
 *  The function can work inplace i.e. in==out*
 *  \param[in]   in     Input array
 *  \param[in]   L      Number of elements in the arrays
 *  \param[out]  out    Output array
 */
void
fftrealfftshift(const double complex* in, int L, double complex* out);

/** ifftshift() in the frequency domain by modulation
 *
 *  The function can work inplace i.e. in==out*
 *  \param[in]   in     Input array
 *  \param[in]   L      Number of elements in the arrays
 *  \param[out]  out    Output array
 */
void
fftrealifftshift(const double complex* in, int L, double complex* out);

/** circshift() in the frequency domain by modulation
 *
 *  *The function can work inplace i.e. in==out*
 *  \param[in]   in     Input array
 *  \param[in]   L      Number of elements in the arrays
 *  \param[in]   shift  Shift
 *  \param[out]  out    Output array
 */
void
fftrealcircshift(const complex double* in, int L, double shift,
                 complex double* out);

/** Shifts cols of height x N matrix by one to the left
 *
 *  \param[in,o]   cols     Input/output matrix
 *  \param[in]     height   Height
 *  \param[in]     N        No. of cols
 *  \param[in]     newcol   (optional) Length height vector to be used as the last col.
 *                          If it is NULL, it is set to zeros. 
 */
void
shiftcolsleft(double* cols, int height, int N, const double* newcol);

#endif /* _utils_h */

/** @}*/
