#include "phaseret/pghi.h"
#include "ltfat/macros.h"
#include <float.h>

struct pghi_plan
{
    double gamma;
    int a;
    int M;
    int W;
    int L;
    double tol1;
    double tol2;
    ltfat_heapinttask_d* hit;
    LTFAT_REAL* tgrad;
    LTFAT_REAL* fgrad;
    /* double* scratch; */
};

int
pghi(const LTFAT_REAL s[], double gamma, const int L, const int W,
     const int a, const int M, LTFAT_COMPLEX c[])
{
    pghi_plan* p = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(s); CHECKNULL(c);
    CHECKSTATUS( pghi_init(gamma, L, W, a, M, 1e-1, 1e-10, &p),
                 "pghi init failed");
    pghi_execute(p, s, c); // This cannot fail

    pghi_done(&p);
error:
    return status;
}

int
pghi_withmask(const LTFAT_COMPLEX cinit[], const int mask[],
              double gamma, const int L, const int W,
              const int a, const int M, LTFAT_COMPLEX c[])
{
    pghi_plan* p = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cinit); CHECKNULL(mask); CHECKNULL(c);

    CHECKSTATUS( pghi_init(gamma, L, W, a, M, 1e-1, 1e-10, &p),
                 "pghi init failed");
    pghi_execute_withmask(p, cinit, mask, NULL, c); // This cannot fail

    pghi_done(&p);
error:
    return status;
}

/* int */
/* pghi_init(double gamma, const int L, const int W, */
/*           const int a, const int M, double tol, pghi_plan** pout) */
/* { */
/*     return pghi_init_twostep(gamma, L, W, a, M, tol, NAN, pout); */
/* } */

int
pghi_init(double gamma, const int L, const int W,
          const int a, const int M, double tol1, double tol2, pghi_plan** pout)
{
    pghi_plan* p = NULL;
    int M2,N;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(pout);
    CHECK(LTFATERR_BADARG, !isnan(gamma) && gamma > 0,
          "gamma cannot be nan and must be positive. (Passed %f).", gamma);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");
    CHECK(LTFATERR_NOTINRANGE, tol1 > 0 && tol1 < 1, "tol1 must be in range ]0,1[");

    if (~isnan(tol2))
    {
        CHECK(LTFATERR_NOTINRANGE, tol2 > 0 && tol2 < 1 && tol2 < tol1,
              "tol2 must be in range ]0,1[ and less or equal to tol1.");
    }

    CHECKMEM( p = (pghi_plan*) calloc(1, sizeof * p));
    p->gamma = gamma; p->a = a; p->M = M; p->W = W; p->L = L; p->tol1 = tol1; 
    p->tol2 = tol2;

    M2 = M / 2 + 1;
    N = L / a;

    CHECKMEM( p->tgrad = LTFAT_NAME_REAL(malloc)(M2 * N));
    CHECKMEM( p->fgrad = LTFAT_NAME_REAL(malloc)(M2 * N));
    /* CHECKMEM( p->scratch = malloc(M2 * N * sizeof * p->scratch)); */
    // Not yet
    p->hit = ltfat_heapinttask_init_d( M2, N, M2 * log((double)M2) , NULL, 1);

    *pout = p;
    return status;
error:
    if (p)
    {
        if (p->tgrad) free(p->tgrad);
        if (p->fgrad) free(p->fgrad);
        /* if (p->scratch) free(p->scratch); */
        if (p->hit) ltfat_heapinttask_done_d(p->hit);
        free(p);
    }
    return status;
}


int
pghi_execute(pghi_plan* p, const LTFAT_REAL s[], LTFAT_COMPLEX c[])
{
    int status = LTFATERR_SUCCESS;
    int M2, W, N;
    CHECKNULL(s); CHECKNULL(c); CHECKNULL(p);

    M2 = p->M / 2 + 1;
    W = p->W;
    N = p->L / p->a;

    for (ltfatInt w = W - 1; w >= 0; --w)
    {
        const LTFAT_REAL* schan = s + w * M2 * N;
        LTFAT_COMPLEX* cchan = c + w * M2 * N;
        const LTFAT_REAL* tgradwchan = p->tgrad + w * M2 * N;
        const LTFAT_REAL* fgradwchan = p->fgrad + w * M2 * N;
        LTFAT_REAL* scratch = ((LTFAT_REAL*)cchan) + M2 * N; // Second half of the output

        pghilog(schan, M2 * N, scratch);
        pghitgrad(scratch, p->gamma, p->a, p->M, N, p->tgrad );
        pghifgrad(scratch, p->gamma, p->a, p->M, N, p->fgrad );

        memset(scratch, 0, M2 * N * sizeof * scratch);

        // Start of without mask
        ltfat_heapinttask_resetmax_d(p->hit, schan, p->tol1);
        ltfat_heapint_execute_d(p->hit, tgradwchan, fgradwchan, scratch);
        int* donemask = ltfat_heapinttask_get_mask_d(p->hit);

        if (!isnan(p->tol2) && p->tol2 < p->tol1)
        {
            // Reuse the just computed mask
            ltfat_heapinttask_resetmask_d(p->hit, donemask, schan, p->tol2, 0);
            ltfat_heapint_execute_d(p->hit, tgradwchan, fgradwchan, scratch);
        }

        // Assign random phase to unused coefficients
        for (int ii = 0; ii < M2 * N; ii++)
            if (donemask[ii] <= LTFAT_MASK_UNKNOWN)
                scratch[ii] = 2.0 * M_PI * ((double)rand()) / RAND_MAX;

        // Combine phase and magnitude
        if (schan != (LTFAT_REAL*) cchan)
        {
            pghimagphase(schan, scratch, M2 * N, cchan);
        }
        else
        {
            // Copy the magnitude first to avoid overwriting it.
            memcpy(p->tgrad, schan, M2 * N * sizeof * schan);
            pghimagphase(p->tgrad, scratch, M2 * N, cchan);
        }
    }
error:
    return status;
}

int
pghi_execute_withmask(pghi_plan* p, const LTFAT_COMPLEX cin[],
                      const int mask[], LTFAT_REAL buffer[], LTFAT_COMPLEX cout[])
{
    LTFAT_REAL* bufferLoc = NULL;
    int freeBufferLoc = 0, M2, W, N;
    LTFAT_REAL* schan;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cin); CHECKNULL(mask); CHECKNULL(cout); CHECKNULL(p);

    M2 = p->M / 2 + 1;
    W = p->W;
    N = p->L / p->a;

    if (buffer)
        bufferLoc = buffer;
    else
    {
        CHECKMEM( bufferLoc = LTFAT_NAME_REAL(malloc)(M2 * N) );
        freeBufferLoc = 1;
    }

    schan = bufferLoc;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_COMPLEX* cinchan = cin + w * M2 * N;
        LTFAT_COMPLEX* coutchan = cout + w * M2 * N;
        const int* maskchan = mask + w * M2 * N;
        const LTFAT_REAL* tgradwchan = p->tgrad + w * M2 * N;
        const LTFAT_REAL* fgradwchan = p->fgrad + w * M2 * N;
        LTFAT_REAL* scratch = ((LTFAT_REAL*)coutchan) + M2 * N; // Second half of the output

        for (int ii = 0; ii < M2 * N; ii++)
            schan[ii] = ltfat_abs(cinchan[ii]);

        pghilog(schan, M2 * N, scratch);
        pghitgrad(scratch, p->gamma, p->a, p->M, N, p->tgrad );
        pghifgrad(scratch, p->gamma, p->a, p->M, N, p->fgrad );

        memset(scratch, 0, M2 * N * sizeof * scratch);

        // Start of without mask
        ltfat_heapinttask_resetmask_d(p->hit, maskchan, schan, p->tol1, 0);
        ltfat_heapint_execute_d(p->hit, tgradwchan, fgradwchan, scratch);
        int* donemask = ltfat_heapinttask_get_mask_d(p->hit);

        if (!isnan(p->tol2))
        {
            // Reuse the just computed mask
            ltfat_heapinttask_resetmask_d(p->hit, donemask, schan, p->tol2, 0);
            ltfat_heapint_execute_d(p->hit, tgradwchan, fgradwchan, scratch);
        }

        // Assign random phase to unused coefficients
        for (int ii = 0; ii < M2 * N; ii++)
            if (donemask[ii] <= LTFAT_MASK_UNKNOWN)
                scratch[ii] = 2.0 * M_PI * ((double)rand()) / RAND_MAX;

        // Combine phase and magnitude
        pghimagphase(schan, scratch, M2 * N, coutchan);

    }
error:
    if (freeBufferLoc) ltfat_free(bufferLoc);
    return status;
}

int
pghi_done(pghi_plan** p)
{
    int status = LTFATERR_SUCCESS;
    pghi_plan* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    ltfat_heapinttask_done_d(pp->hit);
    /* free(pp->scratch); */
    free(pp->fgrad);
    free(pp->tgrad);
    free(pp);
    pp = NULL;
error:
    return status;
}

void
pghimagphase(const LTFAT_REAL s[], const LTFAT_REAL phase[], int L, LTFAT_COMPLEX c[])
{
    for (int l = 0; l < L; l++)
        c[l] = s[l] * exp(I * phase[l]);
}

void
pghilog(const LTFAT_REAL* in, int L, LTFAT_REAL* out)
{
    for (int l = 0; l < L; l++)
        out[l] = log(in[l] + DBL_MIN);

}

void
pghitgrad(const LTFAT_REAL* logs, double gamma, int a, int M, int N, LTFAT_REAL* tgrad)
{
    int M2 = M / 2 + 1;

    const LTFAT_REAL tgradmul = (a * M) / (gamma * 2.0);
    const LTFAT_REAL tgradplus = 2.0 * M_PI * a / M;


    for (int n = 0; n < N; n++)
    {
        LTFAT_REAL* tgradCol = tgrad + n * M2;
        const LTFAT_REAL* logsCol = logs + n * M2;

        tgradCol[0]      = 0.0;
        tgradCol[M2 - 1] = 0.0;

        for (int m = 1; m < M2 - 1; m++)
            tgradCol[m] = tgradmul * (logsCol[m + 1] - logsCol[m - 1]) + tgradplus * m;
    }
}

void
pghifgrad(const LTFAT_REAL* logs, double gamma, int a, int M, int N, LTFAT_REAL* fgrad)
{
    int M2 = M / 2 + 1;

    const double fgradmul = -gamma / (2.0 * a * M);

    for (int n = 1; n < N - 1; n++)
    {
        const LTFAT_REAL* scol0 = logs + (n - 1) * M2;
        const LTFAT_REAL* scol2 = logs + (n + 1) * M2;
        LTFAT_REAL* fgradCol = fgrad + n * M2;

        for (int m = 0; m < M2; ++m)
            fgradCol[m] = fgradmul * (scol2[m] - scol0[m]);
    }

    // Explicit first col
    {
        const LTFAT_REAL* scol0 = logs + (N - 1) * M2;
        const LTFAT_REAL* scol2 = logs + M2;
        LTFAT_REAL* fgradCol = fgrad;

        for (int m = 0; m < M2; ++m)
            fgradCol[m] = fgradmul * (scol2[m] - scol0[m]);
    }

    // Explicit last col
    {
        const LTFAT_REAL* scol0 = logs;
        const LTFAT_REAL* scol2 = logs + (N - 2) * M2;
        LTFAT_REAL* fgradCol = fgrad + (N - 1) * M2;

        for (int m = 0; m < M2; ++m)
            fgradCol[m] = fgradmul * (scol2[m] - scol0[m]);
    }
}
