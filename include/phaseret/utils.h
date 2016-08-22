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

#define LTFAT_DOUBLE
#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#endif
#include "ltfat/types.h"

/** Shifts cols of height x N matrix by one to the left
 *
 *  \param[in,o]   cols     Input/output matrix
 *  \param[in]     height   Height
 *  \param[in]     N        No. of cols
 *  \param[in]     newcol   (optional) Length height vector to be used as the last col.
 *                          If it is NULL, it is set to zeros.
 */
int
shiftcolsleft(LTFAT_REAL* cols, int height, int N, const LTFAT_REAL* newcol);

int
force_magnitude(LTFAT_COMPLEX* cin, const LTFAT_REAL* s, int L, LTFAT_COMPLEX* cout);

void
realimag2absangle(const LTFAT_COMPLEX* cin, const int L, LTFAT_COMPLEX* c);

void
absangle2realimag(const LTFAT_COMPLEX* cin, const int L, LTFAT_COMPLEX* c);


#endif /* _utils_h */

/** @}*/
