#include "rtisila.h"
#include <math.h>

void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const double* frames,
                       const double* g, int n, int N)
{
    int M = p->M;
    int a = p->a;

    // Copy n-th frame
    memcpy(p->frame, frames + M * n, M * sizeof * p->frame);

    // Overlap frames n+1 ... N-1 or until there is overlap
    for (int ii = n + 1; ii < N; ii++)
    {
        int jj = (ii - n) * a;
        int nSamp = M - jj;

        if (nSamp <= 0) break;

        const double* framestmp = frames + M * ii;

        for (int kk = 0; kk < nSamp; kk++)
            p->frame[jj + kk] += framestmp[kk];
    }

    // Overlap frames n-1 ... 0, or until there is no overlap
    for (int ii = n - 1; ii >= 0; ii--)
    {
        int jj = (n - ii) * a;
        int nSamp = M - jj;

        if (nSamp <= 0) break;

        const double* framestmp = frames + M * ii;

        for (int kk = 0; kk < nSamp; kk++)
            p->frame[kk] += framestmp[jj + kk];
    }

    // Multiply with analysis window
    for (int m = 0; m < M; m++)
        p->frame[m] *= g[m];
}

void
rtisilaphaseupdate(rtisilaupdate_plan* p, const double* sframe,
                   double* frameupd)
{
    int M2 = p->M / 2 + 1;

    // FFTREAL
    fftw_execute(p->fwdp);

    // Force the magnitude to target magnitude
    for (int m = 0; m < M2; m++)
        p->fftframe[m] = sframe[m] * cexp(I * carg(p->fftframe[m]));

    // IFFTREAL
    fftw_execute(p->backp);

    // Multiply with the synthesis window
    for (int m = 0; m < p->M; m++)
        p->frame[m] *= p->gd[m];

    // Copy to the output
    memcpy(frameupd, p->frame, p->M * sizeof * frameupd);
}




rtisilaupdate_plan*
rtisilaupdate_init(const double* g, const double* specg1, const double* specg2,
                   const double* gd, int a, int M)
{
    rtisilaupdate_plan* p = fftw_malloc(sizeof * p);

    int M2 = M / 2 + 1;
    // Real input array for FFTREAL and output array for IFFTREAL
    double*          frame = fftw_malloc(M * sizeof * frame);
    // Complex output array for FFTREAL and input array to IFFTREAL
    // The array also holds the final coefficients of the submit frame
    fftw_complex* fftframe = fftw_malloc(M2 * sizeof * fftframe);

    // FFTREAL plan
    fftw_plan fwdp = fftw_plan_dft_r2c_1d(M, frame, fftframe,
                                          FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    // IFFTREAL plan
    // We want to preserve the complex input as we read from it at the end
    fftw_plan backp = fftw_plan_dft_c2r_1d(M, fftframe, frame,
                                           FFTW_MEASURE | FFTW_PRESERVE_INPUT);

    // This way, we can write to the constant parameters of the structure
    rtisilaupdate_plan pdummy =
    {
        .M = M, .a = a,
        .frame = frame, .fftframe = fftframe,
        .g = g, .specg1 = specg1, .specg2 = specg2, .gd = gd,
        .fwdp = fwdp, .backp = backp
    };

    memcpy(p, &pdummy, sizeof * p);

    return p;
}

void
rtisilaupdate_done(rtisilaupdate_plan* p)
{
    fftw_free(p->frame);
    fftw_free(p->fftframe);
    fftw_destroy_plan(p->fwdp);
    fftw_destroy_plan(p->backp);
    fftw_free(p);
}

void
rtisilaupdate_execute(rtisilaupdate_plan* p, const double* frames, int N,
                      const double* s, int lookahead, int maxit, double* frames2)
{
    int lookback = N - lookahead - 1;
    int M = p->M;
    int M2 = M / 2 + 1;

    // If we are not working inplace ...
    if (frames != frames2)
        memcpy(frames2, frames, M * N  * sizeof * frames);

    for (int it = 0; it < maxit; it++)
    {
        for (int nback = lookahead; nback >= 0; nback--)
        {
            int indx = lookback + nback;

            if (nback == lookahead)
                if (it == 0)
                    rtisilaoverlaynthframe(p, frames2, p->specg1, indx, N);
                else
                    rtisilaoverlaynthframe(p, frames2, p->specg2, indx, N);
            else
                rtisilaoverlaynthframe(p, frames2, p->g, indx, N);

            rtisilaphaseupdate(p, s + nback * M2,
                               frames2 +  indx * M);
        }
    }

    fftrealifftshift(p->fftframe, M, p->fftframe);
}

void
rtisilaupdate(const double* frames,
              const double* g, const double* specg1, const double* specg2,
              const double* gd,
              int a, int M, int N, const double* s, int lookahead, int maxit,
              double* frames2)
{
    rtisilaupdate_plan* p = rtisilaupdate_init(g, specg1, specg2, gd, a, M);
    rtisilaupdate_execute(p, frames, N, s, lookahead, maxit, frames2);
    rtisilaupdate_done(p);
}

void
rtisilaupdatecoef(const double* frames,
                  const double* g, const double* specg1, const double* specg2,
                  const double* gd,
                  int a, int M, int N, const double* s, int lookahead, int maxit,
                  double* frames2, double complex* c)
{
    rtisilaupdate_plan* p = rtisilaupdate_init(g, specg1, specg2, gd, a, M);
    rtisilaupdate_execute(p, frames, N, s, lookahead, maxit, frames2);
    memcpy(c, p->fftframe, (M / 2 + 1)*sizeof * p->fftframe);
    rtisilaupdate_done(p);
}



rtisila_plan*
rtisila_init(const double* g, const double* specg1,
             const double* specg2, const double* gd,
             int a, int M, int lookahead, int maxLookahead, int maxit)
{
    int M2 = M / 2 + 1;
    rtisilaupdate_plan* updplan = rtisilaupdate_init(g, specg1, specg2, gd, a, M);
    int lookback = ceil((double) M / a) - 1;
    double* frames = calloc(M * (lookback + 1 + maxLookahead), sizeof * frames);
    double* s   = calloc(M2 * (1 + maxLookahead), sizeof * s);

    rtisila_plan* p = malloc(sizeof * p);
    rtisila_plan pdummy =
    {
        .updateplan = updplan, .lookahead = lookahead, .lookback = lookback,
        .noFrames = lookback + 1 + maxLookahead, .maxLookahead = maxLookahead,
        .maxit = maxit,
        .s = s, .frames = frames, .garbageBinSize = 0, .garbageBin = NULL

    };
    memcpy(p, &pdummy, sizeof * p);
    return p;
}

rtisila_plan*
rtisila_wininit(LTFAT_FIRWIN win, int gl, int a, int M,
                int lookahead, int maxLookahead, int maxit)
{
    double* g = malloc(gl * sizeof * g);
    double* gd = malloc(gl * sizeof * gd);
    double* specg1 = malloc(gl * sizeof * g);
    double* specg2 = malloc(gl * sizeof * gd);

    firwin_d(win, gl, g);
    gabdual_painless_d(g, gl, a, M, gd);

    // TODO: Fill in specg1 specg2 

    rtisila_plan* p = rtisila_init(g, specg1, specg2, gd, a, M, lookahead,
                                   maxLookahead, maxit);

    p->garbageBinSize = 4;
    p->garbageBin = malloc(p->garbageBinSize * sizeof(void*));
    p->garbageBin[0] = specg1;
    p->garbageBin[1] = specg2;
    p->garbageBin[2] = gd;
    p->garbageBin[3] = g;

    return p;
}

void
rtisila_done(rtisila_plan* p)
{
    free(p->s);
    free(p->frames);
    rtisilaupdate_done(p->updateplan);

    if(p->garbageBinSize)
    {
        for (int ii = 0; ii < p->garbageBinSize; ii++)
            free(p->garbageBin[ii]);
    }

    free(p);
}

void
rtisila_execute(rtisila_plan* p, const double* s, double complex* c)
{
    int M = p->updateplan->M;
    int M2 = M / 2 + 1;
    int noFrames = p->lookback + 1 + p->lookahead;

    // Shift frames buffer
    shiftcolsleft(p->frames, M, noFrames, 0);

    // Shift scols buffer
    shiftcolsleft(p->s, M2, p->lookahead + 1, s);

    rtisilaupdate_execute(p->updateplan, p->frames, noFrames, p->s, p->lookahead,
                          p->maxit, p->frames);

    // Read back the coefficients
    memcpy(c, p->updateplan->fftframe, M2 * sizeof * c);
}


void
rtisilaoffline(const double* s,
               const double* g, const double* specg1, const double* specg2, const double* gd,
               int a, int M, int L, int lookahead, int maxit, complex double* c)
{
    int N = L / a;
    int M2 = M / 2 + 1;
    // Just limit lookahead to something sensible
    lookahead = lookahead > N ? N : lookahead;

    rtisila_plan* p = rtisila_init(g, specg1, specg2, gd, a, M, lookahead,
                                   lookahead, maxit);

    for (int n = 0; n < lookahead; ++n)
    {
        const double* sncol = s + n * M2;
        rtisila_execute(p, sncol, c);
    }

    for (int n = 0, nahead = lookahead; nahead < N; ++n, ++nahead)
    {
        const double* sncol = s + nahead * M2;
        complex double* cncol = c + n * M2;
        rtisila_execute(p, sncol, cncol);
    }

    for (int n = N - lookahead; n < N; ++n)
    {
        complex double* cncol = c + n * M2;
        rtisila_execute(p, NULL, cncol);
    }

    rtisila_done(p);
}
