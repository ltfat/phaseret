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

/** Shifts cols of height x N matrix by one to the left
 *
 *  \param[in,o]   cols     Input/output matrix
 *  \param[in]     height   Height
 *  \param[in]     N        No. of cols
 *  \param[in]     newcol   (optional) Length height vector to be used as the last col.
 *                          If it is NULL, it is set to zeros. 
 */
int
shiftcolsleft(double* cols, int height, int N, const double* newcol);

int
force_magnitude(complex double* cin, const double* s, int L, complex double* cout);


#endif /* _utils_h */

/** @}*/
