#ifndef NOSYSTEMHEADERS
#include "ltfat.h"
#endif

#ifndef _phaseret_gsrtisilapghi_h
#define _phaseret_gsrtisilapghi_h
// place for non-templated structs, enums, functions etc.
#endif /* _gsrtisila_h */

#include "ltfat/types.h"
#include "phaseret/types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PHASERET_NAME(gsrtisilapghi_state) PHASERET_NAME(gsrtisilapghi_state);

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_init_win)(LTFAT_FIRWIN win, ltfat_int gl, ltfat_int W,
                                      ltfat_int a, ltfat_int M, ltfat_int lookahead,
                                      ltfat_int maxit, double tol, int do_causalrtpghi,
                                      PHASERET_NAME(gsrtisilapghi_state)** pout);

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_init)(const LTFAT_REAL* g, ltfat_int gl, double gamma,
                                  ltfat_int W, ltfat_int a, ltfat_int M, ltfat_int lookahead,
                                  ltfat_int maxit, double tol, int do_causalrtpghi,
                                  PHASERET_NAME(gsrtisilapghi_state)** pout);

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_execute)(PHASERET_NAME(gsrtisilapghi_state)* p,
                                     const LTFAT_REAL s[], LTFAT_COMPLEX c[]);

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_done)(PHASERET_NAME(gsrtisilapghi_state)** p);

PHASERET_API int
PHASERET_NAME(rtisilaoffline)(const LTFAT_REAL s[], const LTFAT_REAL g[],
                              ltfat_int L, ltfat_int gl, ltfat_int W, ltfat_int a, ltfat_int M,
                              ltfat_int lookahead, ltfat_int maxit, LTFAT_COMPLEX c[]);

#ifdef __cplusplus
}
#endif
