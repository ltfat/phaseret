#include "phaseret/rtisila.h"
#include "phaseret/gsrtisila.h"
#include "phaseret/utils.h"
#include "ltfat/macros.h"

struct PHASERET_NAME(gsrtisilaupdate_plan)
{
    PHASERET_NAME(rtisilaupdate_plan)* p2;
    const LTFAT_REAL* g;
//    const LTFAT_REAL* gd;
    ltfat_int gl;
    ltfat_int M;
    ltfat_int a;
    ltfat_int gNo;
};

struct PHASERET_NAME(gsrtisila_state)
{
    PHASERET_NAME(gsrtisilaupdate_plan)* uplan;
    ltfat_int maxLookahead;
    ltfat_int lookahead;
    ltfat_int lookback;
    ltfat_int maxit;
    ltfat_int W;
    LTFAT_REAL* frames; //!< Buffer for time-domain frames
    LTFAT_COMPLEX* cframes; //!< Buffer for frequency-domain frames
    LTFAT_REAL* s; //!< Buffer for target magnitude
    void** garbageBin;
    ltfat_int garbageBinSize;
};


PHASERET_API int
PHASERET_NAME(gsrtisilaupdate_init)(const LTFAT_REAL* g, const LTFAT_REAL* gd,
                                    ltfat_int gl, ltfat_int a, ltfat_int M,
                                    ltfat_int gNo,
                                    PHASERET_NAME(gsrtisilaupdate_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(gsrtisilaupdate_plan)* p = NULL;
    ltfat_int M2 = M / 2 + 1;

    CHECKMEM( p = (PHASERET_NAME(gsrtisilaupdate_plan)*)
                  ltfat_calloc(1, sizeof * p));
    p->M = M; p->a = a; p->g = g; p->gl = gl; p->gNo = gNo;

    CHECKSTATUS(
        PHASERET_NAME(rtisilaupdate_init)(NULL, NULL, NULL, gd, gl, a, M, &p->p2),
        "rtisilaupdate_init failed");

    *pout = p;
    return status;
error:
    if (p) PHASERET_NAME(gsrtisilaupdate_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(gsrtisilaupdate_done)(PHASERET_NAME(gsrtisilaupdate_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(gsrtisilaupdate_plan)* pp;
    CHECKNULL(p); CHECKNULL(*p);

    pp = *p;
    if (pp->p2) PHASERET_NAME(rtisilaupdate_done)(&pp->p2);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

PHASERET_API void
PHASERET_NAME(gsrtisilaupdate_execute)(PHASERET_NAME(gsrtisilaupdate_plan)* p,
                                       const LTFAT_REAL* frames, const LTFAT_COMPLEX* cframes, ltfat_int N,
                                       const LTFAT_REAL* s, ltfat_int lookahead, ltfat_int maxit,
                                       LTFAT_REAL* frames2, LTFAT_COMPLEX* cframes2,
                                       LTFAT_COMPLEX* c)
{
    ltfat_int lookback = N - lookahead - 1;
    ltfat_int M = p->M;
    ltfat_int gl = p->gl;
    ltfat_int M2 = M / 2 + 1;

    // If we are not working inplace ...
    if (frames != frames2)
        memcpy(frames2, frames, gl * N  * sizeof * frames);

    if (cframes != cframes2)
        memcpy(cframes2, cframes, M2 * N * sizeof * frames);

    for (ltfat_int it = 0; it < maxit; it++)
    {
        for (ltfat_int nback = lookahead; nback >= 0; nback--)
        {
            ltfat_int indx = lookback + nback;
            ltfat_int nfwd = lookahead - nback;

            PHASERET_NAME(rtisilaoverlaynthframe)(p->p2, frames2,
                                                  p->g + nfwd * gl, indx, N);

            PHASERET_NAME(rtisilaphaseupdate)(p->p2, s + nback * M2,
                                              frames2 +  indx * gl,
                                              cframes2 + indx * M2);
        }
    }

    /* for (ltfat_int ii = lookback; ii < N; ii++) */
    /*     LTFAT_NAME_COMPLEX(fftrealcircshift)(cframes2 + ii * M2, M, -(gl / 2) , */
    /*                                          cframes2 + ii * M2); */

    if (c)
        memcpy(c, cframes2 + lookback * M2, M2 * sizeof * c);
}

PHASERET_API int
PHASERET_NAME(gsrtisila_init)(const LTFAT_REAL* g, ltfat_int gl, ltfat_int W,
                              ltfat_int a, ltfat_int M, ltfat_int lookahead, ltfat_int maxit,
                              PHASERET_NAME(gsrtisila_state)** pout)
{
    int status = LTFATERR_SUCCESS;

    PHASERET_NAME(gsrtisila_state)* p = NULL;
    LTFAT_REAL* wins = NULL;
    LTFAT_REAL* gd = NULL;
    LTFAT_REAL* gana = NULL;
    LTFAT_REAL* gcopy = NULL;

    ltfat_int M2, lookback, winsNo, maxLookahead;
    LTFAT_REAL rellim = 1e-3;

    CHECKNULL(g); CHECKNULL(pout);
    CHECK(LTFATERR_BADSIZE, gl > 0, "gl must be positive (passed %d)", gl);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive (passed %d)", W);
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive (passed %d)", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive (passed %d)", M);
    CHECK(LTFATERR_BADARG, lookahead >= 0,
          "lookahead >=0 failed (passed %d)", lookahead);
    CHECK(LTFATERR_NOTPOSARG, maxit > 0,
          "maxit must be positive (passed %d)", maxit);

    M2 = M / 2 + 1;
    lookback = ceil(((LTFAT_REAL) gl) / a) - 1;
    winsNo = (2 * lookback + 1);
    maxLookahead = lookahead;
    CHECKMEM( p = (PHASERET_NAME(gsrtisila_state)*) ltfat_calloc(1, sizeof * p));

    CHECKMEM( gcopy = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM( gd = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM( wins = LTFAT_NAME_REAL(malloc)(gl * winsNo));
    CHECKMEM( gana = LTFAT_NAME_REAL(malloc)(gl * (lookahead + 1)));

    CHECKSTATUS( LTFAT_NAME(gabdual_painless)(g, gl, a, M, gd),
                 "gabdual failed");

    LTFAT_NAME(fftshift)(g, gl, gcopy);
    LTFAT_NAME(fftshift)(gd, gl, gd);

    for (ltfat_int l = 0; l < gl; l++)
        wins[l] = gcopy[l] * gd[l];

    LTFAT_NAME(periodize_array)(wins, gl, gl * winsNo, wins);

    for (ltfat_int n = 0; n < lookahead + 1; n++)
    {
        LTFAT_REAL* ganachan = gana + (lookahead - n) * gl;
        PHASERET_NAME(overlaynthframe)(wins, gl, winsNo, a, lookback + n,
                                       ganachan);

        for (ltfat_int l = 0; l < gl; l++)
        {
            LTFAT_REAL denom = ganachan[l];
            if (denom == 0) denom = 1;
            if (denom < rellim && denom > 0) denom = rellim;
            if (denom > -rellim && denom < 0) denom = -rellim;
            ganachan[l] = M * gcopy[l] / denom;
        }
    }

    ltfat_free(gcopy); gcopy = NULL;
    ltfat_free(wins); wins = NULL;

    CHECKMEM( p->frames =
                  LTFAT_NAME_REAL(calloc)(gl * (lookback + 1 + maxLookahead) * W));
    CHECKMEM( p->s = LTFAT_NAME_REAL(calloc)( M2 * (1 + maxLookahead) * W));
    CHECKMEM( p->cframes =
                  LTFAT_NAME_COMPLEX(calloc)( M2 * (1 + maxLookahead) * W));

    CHECKSTATUS(
        PHASERET_NAME(gsrtisilaupdate_init)(gana, gd, gl, a, M, lookahead + 1,
                                            &p->uplan),
        "rtisiupdate ini failed");

    p->garbageBinSize = 2;
    CHECKMEM( p->garbageBin = (void**)ltfat_malloc(p->garbageBinSize * sizeof(
                                  void*)));
    p->garbageBin[0] = (void*) gana;
    p->garbageBin[1] = (void*) gd;

    p->lookback = lookback;
    p->lookahead = lookahead;
    p->maxLookahead = maxLookahead;
    p->maxit = maxit;
    p->W = W;

    *pout = p;
    return status;
error:
    if (wins) ltfat_free(wins);
    if (gcopy) ltfat_free(gcopy);
    PHASERET_NAME(gsrtisila_done)(&p);

    *pout = NULL;
    return status;
}

PHASERET_API int
PHASERET_NAME(gsrtisila_done)(PHASERET_NAME(gsrtisila_state)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(gsrtisila_state)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    CHECKSTATUS( PHASERET_NAME(gsrtisilaupdate_done)(&pp->uplan),
                 "gsrtisilaupdate done failed");

    if (pp->s) ltfat_free(pp->s);
    if (pp->frames) ltfat_free(pp->frames);
    if (pp->cframes) ltfat_free(pp->cframes);

    if (pp->garbageBinSize)
    {
        for (ltfat_int ii = 0; ii < pp->garbageBinSize; ii++)
            if (pp->garbageBin[ii])
                ltfat_free(pp->garbageBin[ii]);

        ltfat_free(pp->garbageBin);
    }

    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}
