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
