/** \addtogroup pghi
 *  @{
 *
 * @file rtpghi.h
 * @author Zdeněk Průša
 * @date 1 Feb 2016
 * @brief PGHI header
 *
 */

#ifndef _rtpghi_h
#define _rtpghi_h
#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _rtpghiplan rtpghiplan;

rtpghiplan*
rtpghi_init(int M, double tol, int do_causal);

void
rtpghi_done(rtpghiplan* p);

void
rtpghi_execute(rtpghiplan* p, const double* s, complex double* c);

void
rtpghifgrad(const double* logs, int a, int M, double gamma,
            int do_causal, double* fgrad);

void
rtpghitgrad(const double* logs, int a, int M, double gamma,
            double* tgrad);

void
rtpghilog(const double* in, int L, double* out);

void
rtpghimagphase(const double* s, const double* phase, int L, complex double* c);

#ifdef __cplusplus
}
#endif


#endif /* _rtpghi_h */

/** @}*/
