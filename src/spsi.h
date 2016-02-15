/** \addtogroup spsi
 *  @{
 */
#ifndef _spsi_h
#define _spsi_h
#include "config.h"

void
spsireal_split(const double* cr, const double* ci, int a, int M, int L,
         double* mask, double* coutr, double* couti);

void
spsireal(const double* s, int a, int M, int L, complex double* c, double* initphase);

#endif /* _spsi_h */
/** @} */
