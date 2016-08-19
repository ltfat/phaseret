#include "ltfat.h"
#include "phaseret/pghi.h"
#include "phaseret/config.h"
#include <math.h>
#include <float.h>
#include "ltfat/errno.h"
#include "ltfat/macros.h"

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
    double* tgrad;
    double* fgrad;
    /* double* scratch; */
};

int
pghi(const double s[], double gamma, const int L, const int W,
     const int a, const int M, complex double c[])
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
pghi_withmask(const complex double cinit[], const int mask[],
              double gamma, const int L, const int W,
              const int a, const int M, complex double c[])
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

    CHECKMEM( p = calloc(1, sizeof * p));

    int M2 = M / 2 + 1;
    int N = L / a;

    *p = (pghi_plan) { .gamma = gamma, .a = a, .M = M, .W = W, .L = L, .tol1 = tol1, .tol2 = tol2 };
    CHECKMEM( p->tgrad = malloc(M2 * N * sizeof * p->tgrad));
    CHECKMEM( p->fgrad = malloc(M2 * N * sizeof * p->fgrad));
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
pghi_execute(pghi_plan* p, const double s[], complex double c[])
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(s); CHECKNULL(c); CHECKNULL(p);

    int M2 = p->M / 2 + 1;
    int W = p->W;
    int N = p->L / p->a;

    for (ltfatInt w = W - 1; w >= 0; --w)
    {
        const double* schan = s + w * M2 * N;
        complex double* cchan = c + w * M2 * N;
        const double* tgradwchan = p->tgrad + w * M2 * N;
        const double* fgradwchan = p->fgrad + w * M2 * N;
        double* scratch = ((double*)cchan) + M2 * N; // Second half of the output

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
        if (schan != (double*) cchan)
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
pghi_execute_withmask(pghi_plan* p, const complex double cin[],
                      const int mask[], double buffer[], complex double cout[])
{
    double* bufferLoc = NULL;
    int freeBufferLoc = 0;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cin); CHECKNULL(mask); CHECKNULL(cout); CHECKNULL(p);

    int M2 = p->M / 2 + 1;
    int W = p->W;
    int N = p->L / p->a;

    if (buffer)
        bufferLoc = buffer;
    else
    {
        CHECKMEM( bufferLoc = ltfat_malloc(M2 * N * sizeof * bufferLoc) );
        freeBufferLoc = 1;
    }

    double* schan = bufferLoc;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const complex double* cinchan = cin + w * M2 * N;
        complex double* coutchan = cout + w * M2 * N;
        const int* maskchan = mask + w * M2 * N;
        const double* tgradwchan = p->tgrad + w * M2 * N;
        const double* fgradwchan = p->fgrad + w * M2 * N;
        double* scratch = ((double*)coutchan) + M2 * N; // Second half of the output

        for (int ii = 0; ii < M2 * N; ii++)
            schan[ii] = cabs(cinchan[ii]);

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
    CHECKNULL(p); CHECKNULL(*p);
    pghi_plan* pp = *p;
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
pghimagphase(const double s[], const double phase[], int L, complex double c[])
{
    for (int l = 0; l < L; l++)
        c[l] = s[l] * cexp(I * phase[l]);
}

void
pghilog(const double* in, int L, double* out)
{
    for (int l = 0; l < L; l++)
        out[l] = log(in[l] + DBL_MIN);

}

void
pghitgrad(const double* logs, double gamma, int a, int M, int N, double* tgrad)
{
    int M2 = M / 2 + 1;

    const double tgradmul = (a * M) / (gamma * 2.0);
    const double tgradplus = 2.0 * M_PI * a / M;


    for (int n = 0; n < N; n++)
    {
        double* tgradCol = tgrad + n * M2;
        const double* logsCol = logs + n * M2;

        tgradCol[0]      = 0.0;
        tgradCol[M2 - 1] = 0.0;

        for (int m = 1; m < M2 - 1; m++)
            tgradCol[m] = tgradmul * (logsCol[m + 1] - logsCol[m - 1]) + tgradplus * m;
    }
}

void
pghifgrad(const double* logs, double gamma, int a, int M, int N, double* fgrad)
{
    int M2 = M / 2 + 1;

    const double fgradmul = -gamma / (2.0 * a * M);

    for (int n = 1; n < N - 1; n++)
    {
        const double* scol0 = logs + (n - 1) * M2;
        const double* scol2 = logs + (n + 1) * M2;
        double* fgradCol = fgrad + n * M2;

        for (int m = 0; m < M2; ++m)
            fgradCol[m] = fgradmul * (scol2[m] - scol0[m]);
    }

    // Explicit first col
    {
        const double* scol0 = logs + (N - 1) * M2;
        const double* scol2 = logs + M2;
        double* fgradCol = fgrad;

        for (int m = 0; m < M2; ++m)
            fgradCol[m] = fgradmul * (scol2[m] - scol0[m]);
    }

    // Explicit last col
    {
        const double* scol0 = logs;
        const double* scol2 = logs + (N - 2) * M2;
        double* fgradCol = fgrad + (N - 1) * M2;

        for (int m = 0; m < M2; ++m)
            fgradCol[m] = fgradmul * (scol2[m] - scol0[m]);
    }
}
