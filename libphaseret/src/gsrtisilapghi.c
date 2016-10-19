#include "phaseret/rtisila.h"
#include "phaseret/gsrtisila.h"
#include "phaseret/gsrtisilapghi.h"
#include "phaseret/rtpghi.h"
#include "phaseret/utils.h"
#include "ltfat/macros.h"

struct PHASERET_NAME(gsrtisilapghi_state)
{
    PHASERET_NAME(gsrtisila_state)* gsstate;
    PHASERET_NAME(rtpghi_state)* pghistate;
};

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_init_win)(LTFAT_FIRWIN win, ltfat_int gl,
                                      ltfat_int W, ltfat_int a, ltfat_int M, ltfat_int lookahead,
                                      ltfat_int maxit, double tol, int do_causalrtpghi,
                                      PHASERET_NAME(gsrtisilapghi_state)** pout)
{
    double gamma;
    LTFAT_REAL* g = NULL;
    int status = LTFATERR_SUCCESS;
    int initstatus;
    CHECKMEM(g = LTFAT_NAME_REAL(malloc)(gl));

    // Analysis window
    CHECKSTATUS(LTFAT_NAME(firwin)(win, gl, g), "firwin failed");
    gamma = phaseret_firwin2gamma(win, gl);

    initstatus =
        PHASERET_NAME(gsrtisilapghi_init)(g, gl, gamma, W, a, M, lookahead, maxit, tol,
                                          do_causalrtpghi, pout);

    ltfat_free(g);
    return initstatus;
error:
    if (g) ltfat_free(g);
    return status;
}

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_init)(const LTFAT_REAL* g, ltfat_int gl,
                                  double gamma, ltfat_int W, ltfat_int a, ltfat_int M, ltfat_int lookahead,
                                  ltfat_int maxit, double tol, int do_causalrtpghi,
                                  PHASERET_NAME(gsrtisilapghi_state)** pout)
{
    int status = LTFATERR_SUCCESS;
    int initstatus;
    PHASERET_NAME(gsrtisilapghi_state)* p = NULL;

    CHECKNULL(g); CHECKNULL(pout);
    CHECK(LTFATERR_BADARG, !do_causalrtpghi && lookahead == 0,
          "0 lookahead frames cannot be combined with non-causal RTPGHI");

    if (!do_causalrtpghi) lookahead--;

    CHECKMEM(p = (PHASERET_NAME(gsrtisilapghi_state)*)ltfat_calloc(1, sizeof * p));

    initstatus =
        PHASERET_NAME(rtpghi_init)(gamma, W, a, M, tol, do_causalrtpghi, &p->pghistate);
    CHECKSTATUS(initstatus, "RTPGHI init failed");

    initstatus =
        PHASERET_NAME(gsrtisila_init)(g, gl, W, a, M, lookahead, maxit, &p->gsstate);
    CHECKSTATUS(initstatus, "GSRTISILA init failed");

    PHASERET_NAME(gsrtisila_set_skipinitialization)(p->gsstate, 0);

    return initstatus;
error:
    if (p) PHASERET_NAME(gsrtisilapghi_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_execute)(PHASERET_NAME(gsrtisilapghi_state)* p,
                                     const LTFAT_REAL s[], LTFAT_COMPLEX c[])
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(s); CHECKNULL(c);

    PHASERET_NAME(rtpghi_execute)(p->pghistate, s, c);
    // Here c is also an input
    PHASERET_NAME(gsrtisila_execute)(p->gsstate, s, c);
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(gsrtisilapghi_done)(PHASERET_NAME(gsrtisilapghi_state)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(gsrtisilapghi_state)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    if (pp->gsstate)   PHASERET_NAME(gsrtisila_done)(&pp->gsstate);
    if (pp->pghistate) PHASERET_NAME(rtpghi_done)(&pp->pghistate);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}


/* PHASERET_API int */
/* PHASERET_NAME(gsrtisila_reset)(PHASERET_NAME(gsrtisila_state)* p) */
/* { */
/*     ltfat_int N, W, gl, M2; */
/*     int status = LTFATERR_SUCCESS; */
/*     CHECKNULL(p); */
/*  */
/*     N = p->lookback + 1 + p->maxLookahead; */
/*     W = p->W; */
/*     M2 = p->uplan->M / 2 + 1; */
/*     gl = p->uplan->gl; */
/*  */
/*     memset(p->s, 0, M2 * (1 + p->maxLookahead) * W * sizeof * p->s); */
/*     memset(p->frames, 0, gl * N * W * sizeof * p->frames); */
/*     memset(p->cframes, 0, M2 * N * W * sizeof * p->cframes); */
/* error: */
/*     return status; */
/* } */
/*  */
/* PHASERET_API int */
/* PHASERET_NAME(gsrtisila_set_lookahead)(PHASERET_NAME(gsrtisila_state)* p, */
/*                                        ltfat_int lookahead) */
/* { */
/*     int status = LTFATERR_SUCCESS; */
/*     CHECKNULL(p); */
/*     CHECK(LTFATERR_BADARG, lookahead >= 0 && lookahead <= p->maxLookahead, */
/*           "lookahead can only be in range [0-%d] (passed %d).", p->maxLookahead, */
/*           lookahead); */
/*  */
/*     p->lookahead = lookahead; */
/* error: */
/*     return status; */
/* } */
/*  */
/* PHASERET_API int */
/* PHASERET_NAME(gsrtisila_set_itno)(PHASERET_NAME(gsrtisila_state)* p, */
/*                                   ltfat_int it) */
/* { */
/*     int status = LTFATERR_SUCCESS; */
/*     CHECKNULL(p); */
/*     CHECK(LTFATERR_BADARG, it > 0, "it must be greater than 0."); */
/*  */
/*     p->maxit = it; */
/* error: */
/*     return status; */
/*  */
/* } */
