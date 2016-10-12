#include "ltfat/macros.h"
#include "phaseret/rtpghi.h"
#include "phaseret/utils.h"
#include "float.h"

struct PHASERET_NAME(rtpghi_state)
{
    ltfat_int M;
    ltfat_int a;
    ltfat_int W;
    int* mask;
    int do_causal;
    LTFAT_REAL* slog;
    LTFAT_REAL* s;
    LTFAT_REAL* tgrad; //!< Time gradient buffer
    LTFAT_REAL* fgrad; //!< Frequency gradient buffer
    LTFAT_REAL* phase; //!< Buffer for keeping previously computed frame
    double logtol;
    double tol;
    double gamma;
    LTFAT_REAL* randphase; //!< Precomputed array of random phase
    ltfat_int randphaseLen;
    ltfat_int randphaseId;
    struct LTFAT_NAME(heapinttask)* hit;
};

PHASERET_API int
PHASERET_NAME(rtpghi_set_causal)(PHASERET_NAME(rtpghi_state)* p, int do_causal)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    p->do_causal = do_causal;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghi_init)(double gamma, ltfat_int W, ltfat_int a, ltfat_int M, double tol,
                           int do_causal, PHASERET_NAME(rtpghi_state)** pout)
{
    int status = LTFATERR_SUCCESS;

    ltfat_int M2 = M / 2 + 1;
    PHASERET_NAME(rtpghi_state)* p = NULL;
    CHECKMEM( p = (PHASERET_NAME(rtpghi_state)*) ltfat_calloc(1, sizeof * p));

    p->hit = LTFAT_NAME(heapinttask_init)( M2, 2, 2 * M2 , NULL, 1);
    CHECKMEM( p->slog = LTFAT_NAME_REAL(calloc)(3 * M2 * W));
    CHECKMEM( p->s = LTFAT_NAME_REAL(calloc)(2 * M2 * W));
    CHECKMEM( p->phase = LTFAT_NAME_REAL(calloc)(2 * M2 * W));
    CHECKMEM( p->tgrad = LTFAT_NAME_REAL(calloc)(3 * M2 * W));
    CHECKMEM( p->fgrad = LTFAT_NAME_REAL(calloc)(2 * M2 * W));

    p->randphaseLen = 10 * M2 * W;
    CHECKMEM( p->randphase = LTFAT_NAME_REAL(malloc)(p->randphaseLen));

    for (ltfat_int ii = 0; ii < p->randphaseLen; ii++)
        p->randphase[ii] = 2.0 * M_PI * ((double)rand()) / RAND_MAX;

    CHECKMEM( p->mask = (int*) ltfat_malloc(2 * M2 * W * sizeof * p->mask));

    for (ltfat_int w = 0; w < W; ++w)
    {
        for (ltfat_int ii = 0; ii < M2; ++ii)
        {
            p->mask[ii + w * 2 * M2] = 1;
            p->mask[ii + M2 + w * 2 * M2] = 0;
        }
    }

    p->do_causal = do_causal;
    p->logtol = log(tol + DBL_MIN);
    p->tol = tol;
    p->M = M;
    p->a = a;
    p->gamma = gamma;
    p->W = W;

    *pout = p;
    return status;
error:
    if (p) PHASERET_NAME(rtpghi_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghi_reset)(PHASERET_NAME(rtpghi_state)* p)
{
    int status = LTFATERR_SUCCESS;
    ltfat_int M2, W;
    CHECKNULL(p);
    M2 = p->M / 2 + 1;
    W = p->W;

    p->randphaseId = 0;
    memset(p->slog, 0, 3 * M2 * W * sizeof * p->slog);
    memset(p->s, 0, 2 * M2 * W * sizeof * p->s);
    memset(p->phase, 0, 2 * M2 * W * sizeof * p->phase);
    memset(p->tgrad, 0, 3 * M2 * W * sizeof * p->tgrad);
    memset(p->fgrad, 0, 2 * M2 * W * sizeof * p->tgrad);

error:
    return status;
}




PHASERET_API int
PHASERET_NAME(rtpghi_execute)(PHASERET_NAME(rtpghi_state)* p,
                              const LTFAT_REAL s[], LTFAT_COMPLEX c[])
{
    // n, n-1, n-2 frames
    // s is n-th
    ltfat_int M2 = p->M / 2 + 1;
    ltfat_int W = p->W;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(s); CHECKNULL(c);

    for (ltfat_int w = 0; w < W; ++w)
    {
        LTFAT_REAL* slogCol = p->slog + w * 3 * M2;
        LTFAT_REAL* sCol = p->s + w * 2 * M2;
        LTFAT_REAL* phaseCol = p->phase + w * 2 * M2;
        LTFAT_REAL* tgradCol = p->tgrad + w * 3 * M2;
        LTFAT_REAL* fgradCol = p->fgrad + w * 2 * M2;

        // store log(s)
        PHASERET_NAME(shiftcolsleft)(slogCol, M2, 3, NULL);
        PHASERET_NAME(shiftcolsleft)(sCol, M2, 2, s + w * M2);
        PHASERET_NAME(rtpghilog)(sCol + M2, M2, slogCol + 2 * M2);

        // Store current s

        // Compute and store tgrad for n
        PHASERET_NAME(shiftcolsleft)(tgradCol, M2, 3, NULL);
        PHASERET_NAME(rtpghitgrad)(slogCol + 2 * M2, p->a, p->M, p->gamma,
                                   tgradCol + 2 * M2);

        // Compute fgrad for n or n-1
        PHASERET_NAME(rtpghifgrad)(slogCol, p->a, p->M, p->gamma, p->do_causal,
                                   fgradCol + M2);

        LTFAT_NAME(heapinttask_resetmask)(p->hit,
                                          p->mask, p->do_causal ? slogCol + M2 : slogCol, p->logtol, 1);

        LTFAT_NAME(heapint_execute)(p->hit,
                                    p->do_causal ? tgradCol + M2 : tgradCol, fgradCol, phaseCol);

        PHASERET_NAME(shiftcolsleft)(phaseCol, M2, 2, NULL);

        // Fill in the missing phase from the precomputed random array
        int* donemask = LTFAT_NAME(heapinttask_get_mask)(p->hit);

        for (ltfat_int ii = 0; ii < M2; ii++)
        {
            if (donemask[M2 + ii] <= LTFAT_MASK_UNKNOWN)
            {
                phaseCol[ii] = p->randphase[p->randphaseId++];
                p->randphaseId %= p->randphaseLen;
            }
        }

        // Combine phase with magnitude
        PHASERET_NAME(rtpghimagphase)(p->do_causal ? sCol + M2 : sCol, phaseCol, M2,
                                      c + w * M2);
    }

error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghi_done)(PHASERET_NAME(rtpghi_state)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(rtpghi_state)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    if (pp->hit) LTFAT_NAME(heapinttask_done)(pp->hit);
    if (pp->slog) ltfat_free(pp->slog);
    if (pp->s) ltfat_free(pp->s);
    if (pp->mask) ltfat_free(pp->mask);
    if (pp->phase) ltfat_free(pp->phase);
    if (pp->tgrad) ltfat_free(pp->tgrad);
    if (pp->fgrad) ltfat_free(pp->fgrad);
    if (pp->randphase) ltfat_free(pp->randphase);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghioffline)(const LTFAT_REAL* s, double gamma, ltfat_int L, ltfat_int W,
                             ltfat_int a, ltfat_int M, double tol, int do_causal, LTFAT_COMPLEX* c)
{
    ltfat_int N = L / a;
    ltfat_int M2 = M / 2 + 1;
    PHASERET_NAME(rtpghi_state)* p = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(s); CHECKNULL(c);

    CHECKSTATUS( PHASERET_NAME(rtpghi_init)(gamma, 1, a, M, tol, do_causal, &p),
                 "rtpghi init failed");

    if (do_causal)
    {
        for (ltfat_int w = 0; w < W; w++)
        {
            PHASERET_NAME(rtpghi_reset)(p);

            for (ltfat_int n = 0; n < N; ++n)
            {
                const LTFAT_REAL* sncol = s + n * M2 + w * N * M2;
                LTFAT_COMPLEX* cncol =    c + n * M2 + w * N * M2;
                PHASERET_NAME(rtpghi_execute)(p, sncol, cncol);
            }
        }
    }
    else
    {
        for (ltfat_int w = 0; w < W; w++)
        {
            PHASERET_NAME(rtpghi_reset)(p);
            PHASERET_NAME(rtpghi_execute)(p, s + w * N * M2, c + w * N * M2);

            for (ltfat_int n = 0, nahead = 1; nahead < N; ++n, ++nahead)
            {
                const LTFAT_REAL* sncol = s + nahead * M2 + w * N * M2;
                LTFAT_COMPLEX* cncol =    c + n * M2 + w * N * M2;
                PHASERET_NAME(rtpghi_execute)(p, sncol, cncol);
            }

            PHASERET_NAME(rtpghi_execute)(p, s + w * N * M2, c + (N - 1) * M2 + w * N * M2);
        }
    }

    PHASERET_NAME(rtpghi_done)(&p);
error:
    return status;
}


void
PHASERET_NAME(rtpghifgrad)(const LTFAT_REAL* logs, ltfat_int a, ltfat_int M, double gamma,
                           int do_causal, LTFAT_REAL* fgrad)
{
    ltfat_int M2 = M / 2 + 1;

    const LTFAT_REAL fgradmul = -gamma / (2.0 * a * M);
    const LTFAT_REAL* scol0 = logs;
    const LTFAT_REAL* scol2 = logs + 2 * M2;

    if (do_causal)
    {
        const LTFAT_REAL* scol1 = logs + M2;

        for (ltfat_int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (3.0 * scol2[m] - 4.0 * scol1[m] + scol0[m]);
    }
    else
    {
        for (ltfat_int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (scol2[m] - scol0[m]);
    }
}

void
PHASERET_NAME(rtpghitgrad)(const LTFAT_REAL* logs, ltfat_int a, ltfat_int M, double gamma,
                           LTFAT_REAL* tgrad)
{
    ltfat_int M2 = M / 2 + 1;

    const LTFAT_REAL tgradmul = (a * M) / (gamma * 2.0);
    const LTFAT_REAL tgradplus = 2.0 * M_PI * a / M;

    tgrad[0]      = 0.0;
    tgrad[M2 - 1] = 0.0;

    for (ltfat_int m = 1; m < M2 - 1; m++)
        tgrad[m] = tgradmul * (logs[m + 1] - logs[m - 1]) + tgradplus * m;
}

void
PHASERET_NAME(rtpghilog)(const LTFAT_REAL* in, ltfat_int L, LTFAT_REAL* out)
{
    for (ltfat_int l = 0; l < L; l++)
        out[l] = log(in[l] + DBL_MIN);

}

void
PHASERET_NAME(rtpghimagphase)(const LTFAT_REAL* s, const LTFAT_REAL* phase,
                              ltfat_int L,
                              LTFAT_COMPLEX* c)
{
    for (ltfat_int l = 0; l < L; l++)
        c[l] = s[l] * exp(I * phase[l]);
}
