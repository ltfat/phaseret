
#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#endif
#include "ltfat/types.h"
#include "phaseret/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Shifts cols of height x N matrix by one to the left
 *
 *  \param[in,o]   cols     Input/output matrix
 *  \param[in]     height   Height
 *  \param[in]     N        No. of cols
 *  \param[in]     newcol   (optional) Length height vector to be used as the last col.
 *                          If it is NULL, it is set to zeros.
 */
int
PHASERET_NAME(shiftcolsleft)(LTFAT_REAL cols[], int height, int N, const LTFAT_REAL newcol[]);

int
PHASERET_NAME(force_magnitude)(LTFAT_COMPLEX cin[], const LTFAT_REAL s[], int L, LTFAT_COMPLEX cout[]);

void
PHASERET_NAME(realimag2absangle)(const LTFAT_COMPLEX cin[], const int L, LTFAT_COMPLEX c[]);

void
PHASERET_NAME(absangle2realimag)(const LTFAT_COMPLEX cin[], const int L, LTFAT_COMPLEX c[]);

void
PHASERET_NAME(absangle2realimag_split2inter)(const LTFAT_REAL s[],
        const LTFAT_REAL phase[], const int L, LTFAT_COMPLEX c[]);

#ifdef __cplusplus
}
#endif


