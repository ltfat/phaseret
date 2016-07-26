#include "phaseret/rtisila.h"
#include "phaseret/utils.h"
#include "ltfat/macros.h"
#include <math.h>

struct rtisilaupdate_plan
{
    double* frame; //!< Time domain buffer
    fftw_complex* fftframe; //!< Frequency domain buffer
    fftw_plan fwdp; //!< Real FFT plan
    fftw_plan backp; //!< Real IFFT plan
    const double* g;
    const double* gd;
    const double* specg1;
    const double* specg2;
    int gl;
    int M;
    int a;
};

struct rtisila_state
{
    rtisilaupdate_plan* uplan;
    int maxLookahead;
    int lookahead;
    int lookback;
    int maxit;
    int W;
    double* frames; //!< Buffer for time-domain frames
    double* s; //!< Buffer for target magnitude
    void** garbageBin;
    int garbageBinSize;
};

int
rtisila_set_lookahead(rtisila_state* p, int lookahead)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_BADARG, lookahead >= 0,
          "Zero or more lookahead frames is required. Passed %d.", lookahead);

    if (lookahead < p->lookahead)
        p->lookahead = lookahead;

    return p->lookahead;
error:
    return status;
}


void
overlaynthframe(const double* frames, int gl, int N, int a, int n,
                double* frame)
{
    // Copy n-th frame
    memcpy(frame, frames + gl * n, gl * sizeof * frame);

    // Overlap frames n+1 ... N-1 or until there is overlap
    for (int ii = n + 1; ii < N; ii++)
    {
        int jj = (ii - n) * a;
        int nSamp = gl - jj;

        if (nSamp <= 0) break;

        const double* framestmp = frames + gl * ii;

        for (int kk = 0; kk < nSamp; kk++)
            frame[jj + kk] += framestmp[kk];
    }

    // Overlap frames n-1 ... 0, or until there is no overlap
    for (int ii = n - 1; ii >= 0; ii--)
    {
        int jj = (n - ii) * a;
        int nSamp = gl - jj;

        if (nSamp <= 0) break;

        const double* framestmp = frames + gl * ii;

        for (int kk = 0; kk < nSamp; kk++)
            frame[kk] += framestmp[jj + kk];
    }
}

void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const double* frames,
                       const double* g, int n, int N)
{
    int M = p->M;
    int gl = p->gl;
    int a = p->a;

    overlaynthframe(frames, gl, N, a, n, p->frame);

    /* // Copy n-th frame */
    /* memcpy(p->frame, frames + gl * n, M * sizeof * p->frame); */
    /*  */
    /* // Overlap frames n+1 ... N-1 or until there is overlap */
    /* for (int ii = n + 1; ii < N; ii++) */
    /* { */
    /*     int jj = (ii - n) * a; */
    /*     int nSamp = gl - jj; */
    /*  */
    /*     if (nSamp <= 0) break; */
    /*  */
    /*     const double* framestmp = frames + gl * ii; */
    /*  */
    /*     for (int kk = 0; kk < nSamp; kk++) */
    /*         p->frame[jj + kk] += framestmp[kk]; */
    /* } */
    /*  */
    /* // Overlap frames n-1 ... 0, or until there is no overlap */
    /* for (int ii = n - 1; ii >= 0; ii--) */
    /* { */
    /*     int jj = (n - ii) * a; */
    /*     int nSamp = gl - jj; */
    /*  */
    /*     if (nSamp <= 0) break; */
    /*  */
    /*     const double* framestmp = frames + gl * ii; */
    /*  */
    /*     for (int kk = 0; kk < nSamp; kk++) */
    /*         p->frame[kk] += framestmp[jj + kk]; */
    /* } */

    // Multiply with analysis window
    for (int m = 0; m < gl; m++)
        p->frame[m] *= g[m];

    // Set remaining samples to zeros
    for (int m = gl; m < M; m++)
        p->frame[m] = 0.0;
}

void
rtisilaphaseupdate(rtisilaupdate_plan* p, const double* sframe,
                   double* frameupd, fftw_complex* c)
{
    int M = p->M;
    int gl = p->gl;
    int M2 = M / 2 + 1;

    // FFTREAL
    fftw_execute(p->fwdp);

    // Force the magnitude to target magnitude
    for (int m = 0; m < M2; m++)
        p->fftframe[m] = sframe[m] * cexp(I * carg(p->fftframe[m]));

    // Copy before it gets overwritten
    if (c)
        memcpy(c, p->fftframe, M2 * sizeof * c);

    // IFFTREAL // Overwrites input p->fftframe
    fftw_execute(p->backp);

    // Multiply with the synthesis window
    for (int m = 0; m < gl; m++)
        p->frame[m] *= p->gd[m];

    // Set remaining samples to zeros
    for (int m = gl; m < M; m++)
        p->frame[m] = 0;

    // Copy to the output
    memcpy(frameupd, p->frame, gl * sizeof * frameupd);
}

int
rtisilaupdate_init(const double* g, const double* specg1, const double* specg2,
                   const double* gd, const int gl, int a, int M,
                   rtisilaupdate_plan** pout)
{
    int status = LTFATERR_SUCCESS;
    rtisilaupdate_plan* p = NULL;
    int M2 = M / 2 + 1;

    CHECKMEM( p = calloc(1, sizeof * p));

    *p = (rtisilaupdate_plan)
    {
        .M = M, .a = a, .g = g, .gd = gd, .specg1 = specg1,
         .specg2 = specg2, .gl = gl
    };
    // Real input array for FFTREAL and output array for IFFTREAL
    CHECKMEM( p->frame = fftw_malloc(M * sizeof * p->frame));
    // Complex output array for FFTREAL and input array to IFFTREAL
    CHECKMEM( p->fftframe = fftw_malloc(M2 * sizeof * p->fftframe));

    // FFTREAL plan
    p->fwdp = fftw_plan_dft_r2c_1d(M, p->frame, p->fftframe, FFTW_MEASURE);
    CHECKINIT(p->fwdp, "FFTW plan failed");
    // IFFTREAL plan
    p->backp = fftw_plan_dft_c2r_1d(M, p->fftframe, p->frame, FFTW_MEASURE);
    CHECKINIT(p->backp, "FFTW plan failed");

    *pout = p;
    return status;
error:
    if (p)
    {
        if (p->frame) free(p->frame);
        if (p->fftframe) free(p->fftframe);
        if (p->fwdp) fftw_destroy_plan(p->fwdp);
        if (p->backp) fftw_destroy_plan(p->backp);
        free(p);
    }
    return status;
}

int
rtisilaupdate_done(rtisilaupdate_plan** p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);

    rtisilaupdate_plan* pp = *p;
    fftw_free(pp->frame);
    fftw_free(pp->fftframe);
    fftw_destroy_plan(pp->fwdp);
    fftw_destroy_plan(pp->backp);
    free(pp);
    pp = NULL;
error:
    return status;
}

void
rtisilaupdate_execute(rtisilaupdate_plan* p, const double* frames, int N,
                      const double* s, int lookahead, int maxit, double* frames2,
                      fftw_complex* c)
{
    int lookback = N - lookahead - 1;
    int M = p->M;
    int gl = p->gl;
    int M2 = M / 2 + 1;

    // If we are not working inplace ...
    if (frames != frames2)
        memcpy(frames2, frames, gl * N  * sizeof * frames);

    for (int it = 0; it < maxit; it++)
    {
        for (int nback = lookahead; nback >= 0; nback--)
        {
            int indx = lookback + nback;

            // Newest lookahead frame is treated differently
            if (nback == lookahead)
                if (it == 0)
                    rtisilaoverlaynthframe(p, frames2, p->specg1, indx, N);
                else
                    rtisilaoverlaynthframe(p, frames2, p->specg2, indx, N);
            else
                rtisilaoverlaynthframe(p, frames2, p->g, indx, N);

            if ( nback == 0 && it == (maxit - 1) )
                rtisilaphaseupdate(p, s + nback * M2, frames2 +  indx * gl, c);
            else
                rtisilaphaseupdate(p, s + nback * M2, frames2 +  indx * gl, NULL);
        }
    }

    if (c)
        ltfat_fftrealcircshift_dc(c, M, -(gl / 2) , c);
}

void
rtisilaupdate(const double* frames,
              const double* g, const double* specg1, const double* specg2,
              const double* gd, const int gl,
              int a, int M, int N, const double* s, int lookahead, int maxit,
              double* frames2)
{
    rtisilaupdate_plan* p = NULL;
    rtisilaupdate_init(g, specg1, specg2, gd, gl, a, M, &p);
    rtisilaupdate_execute(p, frames, N, s, lookahead, maxit, frames2, NULL);
    rtisilaupdate_done(&p);
}

void
rtisilaupdatecoef(const double* frames,
                  const double* g, const double* specg1, const double* specg2,
                  const double* gd, const int gl,
                  int a, int M, int N, const double* s, int lookahead, int maxit,
                  double* frames2, double complex* c)
{
    rtisilaupdate_plan* p = NULL;
    rtisilaupdate_init(g, specg1, specg2, gd, gl, a, M, &p);
    rtisilaupdate_execute(p, frames, N, s, lookahead, maxit, frames2, c);
    rtisilaupdate_done(&p);
}

int
rtisila_init(const double* g, const double* gd, const int gl, const int W,
             int a, int M, int lookahead, int maxit, rtisila_state** pout)
{
    int status = LTFATERR_SUCCESS;

    rtisila_state* p = NULL;
    double* wins = NULL;

    CHECKNULL(g); CHECKNULL(gd); CHECKNULL(pout);
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");
    CHECK(LTFATERR_BADARG, lookahead >= 0, "lookahead >=0 failed");
    CHECK(LTFATERR_NOTPOSARG, maxit > 0, "maxit must be positive");

    int M2 = M / 2 + 1;
    int lookback = ceil((double) gl / a) - 1;
    int winsNo = (2 * lookback + 1);
    int maxLookahead = lookahead;
    CHECKMEM( p = calloc(1, sizeof * p));

    CHECKSTATUS(
        rtisilaupdate_init(NULL, NULL, NULL, NULL, gl, a, M, &p->uplan),
        "rtisiupdate ini failed");

    CHECKMEM( p->uplan->g = malloc(gl * sizeof * p->uplan->g));
    CHECKMEM( p->uplan->gd = malloc(gl * sizeof * p->uplan->gd));
    CHECKMEM( p->uplan->specg1 = malloc(gl * sizeof * p->uplan->specg1));
    CHECKMEM( p->uplan->specg2 = malloc(gl * sizeof * p->uplan->specg2));

    ltfat_fftshift_d(g, gl, (double*) p->uplan->g);
    ltfat_fftshift_d(gd, gl, (double*) p->uplan->gd);

    CHECKMEM( wins = malloc(gl * winsNo * sizeof * wins));

    // Spec2
    ltfat_periodize_array_d(p->uplan->gd, gl, gl * winsNo, wins);
    overlaynthframe(wins, gl, winsNo, a, 0, (double*) p->uplan->specg2);

    // Spec1
    memset(wins, 0, gl * sizeof * wins);
    overlaynthframe(wins, gl, winsNo, a, 0, (double*) p->uplan->specg1);

    free(wins); wins = NULL;

    CHECKMEM( p->frames = calloc(gl * (lookback + 1 + maxLookahead) * W,
                                 sizeof * p->frames));
    CHECKMEM( p->s = calloc( M2 * (1 + maxLookahead) * W, sizeof * p->s));

    p->garbageBinSize = 4;
    CHECKMEM( p->garbageBin = malloc(p->garbageBinSize * sizeof(void*)));
    p->garbageBin[0] = (void*) p->uplan->g;
    p->garbageBin[1] = (void*) p->uplan->gd;
    p->garbageBin[2] = (void*) p->uplan->specg1;
    p->garbageBin[3] = (void*) p->uplan->specg2;

    p->lookback = lookback;
    p->lookahead = lookahead;
    p->maxLookahead = maxLookahead;
    p->maxit = maxit;
    p->W = W;

    *pout = p;
    return status;
error:
    if (wins) free(wins);
    if (p)
    {
        if (p->uplan)
        {
            if (p->uplan->g) free((void*)p->uplan->g);
            if (p->uplan->gd) free((void*)p->uplan->gd);
            if (p->uplan->specg1) free((void*)p->uplan->specg1);
            if (p->uplan->specg2) free((void*)p->uplan->specg2);
            rtisilaupdate_done(&p->uplan);
        }
        if (p->s) free(p->s);
        if (p->frames) free(p->frames);
        free(p);
    }
    *pout = NULL;
    return status;
}

int
rtisila_init_win(LTFAT_FIRWIN win, int gl, int W, int a, int M,
                 int lookahead, int maxit,
                 rtisila_state** pout)
{
    double* g = NULL;
    double* gd = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKMEM( g = malloc(gl * sizeof * g));
    CHECKMEM( gd = malloc(gl * sizeof * gd));

    // Analysis window
    CHECKSTATUS( ltfat_firwin_d(win, gl, g), "firwin failed");
    CHECKSTATUS( ltfat_gabdual_painless_d(g, gl, a, M, gd),
                 "gabdual failed");

    int initstatus = rtisila_init(g, gd, gl, W, a, M, lookahead,
                                  maxit, pout);

    free(g); free(gd);
    return initstatus;
error:
    if (g) free(g);
    if (gd) free(gd);
    return status;
}

int
rtisila_done(rtisila_state** p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    rtisila_state* pp = *p;

    CHECKSTATUS( rtisilaupdate_done(&pp->uplan),
                 "rtisilaupdate done failed");

    free(pp->s);
    free(pp->frames);

    if (pp->garbageBinSize)
    {
        for (int ii = 0; ii < pp->garbageBinSize; ii++)
            free(pp->garbageBin[ii]);

        free(pp->garbageBin);
    }

    free(pp);
    pp = NULL;
error:
    return status;
}

int
rtisila_execute(rtisila_state* p, const double* s, double complex* c)
{
    int M = p->uplan->M;
    int gl = p->uplan->gl;
    int M2 = M / 2 + 1;
    int noFrames = p->lookback + 1 + p->lookahead;
    int N = p->lookback + 1 + p->maxLookahead;

    for (int w = 0; w < p->W; w++)
    {
        const double* schan = s + w * M2;
        double complex* cchan = c + w * M2;
        double* frameschan = p->frames + w * N * gl;
        double* sframeschan = p->s  + w * (1 + p->maxLookahead) * M2;
        // Shift frames buffer
        shiftcolsleft( frameschan, gl, noFrames, NULL);

        // Shift scols buffer
        shiftcolsleft( sframeschan, M2, p->lookahead + 1, schan);

        rtisilaupdate_execute(p->uplan, frameschan, noFrames, sframeschan,
                              p->lookahead, p->maxit, frameschan, cchan);
    }

    return LTFATERR_SUCCESS;
}


int
rtisilaoffline(const double* s,
               const double* g, const double* gd,
               const int gl, int a, int M, int L, int lookahead, int maxit,
               complex double* c)
{
    int status = LTFATERR_SUCCESS;
    rtisila_state* p = NULL;

    CHECKNULL(s); CHECKNULL(c);
    int N = L / a;
    int M2 = M / 2 + 1;
    // Just limit lookahead to something sensible
    lookahead = lookahead > N ? N : lookahead;

    CHECKSTATUS(
        rtisila_init(g, gd, gl, 1, a, M, lookahead, maxit, &p),
        "rtisila init failed");

    memcpy(p->s + M2, s, lookahead * M2 * sizeof * p->s);

    for (int n = 0, nahead = lookahead; nahead < N; ++n, ++nahead)
    {
        const double* sncol = s + nahead * M2;
        complex double* cncol = c + n * M2;
        rtisila_execute(p, sncol, cncol);
    }

    for (int n = N - lookahead, nahead = 0; n < N; ++n, ++nahead)
    {
        const double* sncol = s + nahead * M2;
        complex double* cncol = c + n * M2;
        rtisila_execute(p, sncol, cncol);
    }

    rtisila_done(&p);

error:
    return status;
}
