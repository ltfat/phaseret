#include "rtisila.h"
#include <math.h>

// in might be equal to out
void
fftshift(const double* in, int L, double* out)
{
    const div_t domod = div(L, 2);
    if (domod.rem)
    {
        double first = in[0];
        for (int ii = 0; ii < domod.quot; ii++)
        {
            double tmp = in[ii + 1];
            out[ii] = in[ii + domod.quot + 1];
            out[ii + domod.quot + 1] = tmp;
        }
        out[domod.quot] = first;
    }
    else
    {
        for (int ii = 0; ii < domod.quot; ii++)
        {
            double tmp = in[ii];
            out[ii] = in[ii + domod.quot];
            out[ii + domod.quot] = tmp;
        }
    }
}

// in might be equal to out
void
ifftshift(const double* in, int L, double* out)
{
    const div_t domod = div(L, 2);
    if (domod.rem)
    {
        double middle = in[domod.quot];
        for (int ii = 0; ii < domod.quot; ii++)
        {
            double tmp = in[L - ii - 1];
            out[L - ii - 1] = in[domod.quot - ii - 1];
            out[domod.quot - ii] = tmp;
        }
        out[0] = middle;
    }
    else
    {
        fftshift(in, L, out);
    }
}

void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const double* cframes,
                       const double* gnum, int n, int N)
{
    int M = p->M;
    int a = p->a;

    ifftshift(cframes + M * n, M, p->frame);

    for (int ii = n + 1; ii < N; ii++)
    {
        int jj = (ii - n) * a;
        int nSamp = M - jj;
        if (nSamp <= 0) break;

        ifftshift(cframes + M * ii, M, p->shiftframe);

        for (int kk = 0; kk < nSamp; kk++)
            p->frame[jj + kk] += p->shiftframe[kk];
    }

    for (int ii = n - 1; ii >= 0; ii--)
    {
        int jj = (n - ii) * a;
        int nSamp = M - jj;
        if (nSamp <= 0) break;

        ifftshift(cframes + M * ii, M, p->shiftframe);

        for (int kk = 0; kk < nSamp; kk++)
            p->frame[kk] += p->shiftframe[jj + kk];
    }

    fftshift(p->frame, M, p->frame);

    for (int m = 0; m < M; m++)
        p->frame[m] *= gnum[m];
}

void
rtisilaphaseupdate(rtisilaupdate_plan* p, const double* sframe, double* frameupd)
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
        frameupd[m] = p->frame[m] * p->gd[m];
}




rtisilaupdate_plan*
rtisilaupdate_init(const double *g, const double* specg1, const double* specg2,
                   const double* gd, int a, int M)
{
    rtisilaupdate_plan* p = fftw_malloc(sizeof * p);

    int M2 = M / 2 + 1;
    // Real input array for FFTREAL and output array for IFFTREAL
    double*          frame = fftw_malloc(M * sizeof * frame);
    // Complex output array for FFTREAL and input array to IFFTREAL
    // The array also holds the final coefficients of the submit frame
    fftw_complex* fftframe = fftw_malloc(M2 * sizeof * fftframe);
    // Just a helper array for doing fftshifts
    double*     shiftframe = fftw_malloc(M * sizeof * shiftframe);

    // FFTREAL plan
    fftw_plan fwdp = fftw_plan_dft_r2c_1d(M, frame, fftframe, FFTW_MEASURE);
    // IFFTREAL plan
    // We want to preserve the complex input as we read from it at the end
    fftw_plan backp = fftw_plan_dft_c2r_1d(M, fftframe, frame, FFTW_MEASURE | FFTW_PRESERVE_INPUT);

    // This way, we can write to the constant parameters of the structure
    rtisilaupdate_plan pdummy =
    {
        .M = M, .a = a,
        .frame = frame, .fftframe = fftframe, .shiftframe = shiftframe,
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
    fftw_free(p->shiftframe);
    fftw_destroy_plan(p->fwdp);
    fftw_destroy_plan(p->backp);
    fftw_free(p);
}

void
rtisilaupdate_execute(rtisilaupdate_plan* p, const double* cframes, int N,
                      const double* sframes, int lookahead, int maxit, double* cframes2)
{
    int lookback = N - lookahead - 1;
    int M = p->M;
    int M2 = M / 2 + 1;

    // If we are not working inplace ...
    if (cframes != cframes2)
        memcpy(cframes2, cframes, M * N  * sizeof * cframes);

    for (int it = 0; it < maxit; it++)
    {
        for (int nback = lookahead; nback >= 0; nback--)
        {
            int indx = lookback + nback;

            if (nback == lookahead)
                if (it == 0)
                    rtisilaoverlaynthframe(p, cframes2, p->specg1, indx, N);
                else
                    rtisilaoverlaynthframe(p, cframes2, p->specg2, indx, N);
            else
                rtisilaoverlaynthframe(p, cframes2, p->g, indx, N);

            rtisilaphaseupdate(p, sframes + nback * M2,
                               cframes2 +  indx * M);
        }
    }

}

void
rtisilaupdate(const double* cframes,
              const double* gnum, const double* specg1, const double* specg2, const double* dgnum,
              int a, int M, int N, const double* sframes, int lookahead, int maxit,
              double* cframes2)
{
    rtisilaupdate_plan *p = rtisilaupdate_init(gnum, specg1, specg2, dgnum, a, M);
    rtisilaupdate_execute(p, cframes, N, sframes, lookahead, maxit, cframes2);
    rtisilaupdate_done(p);
}

void
rtisilaupdatecoef(const double* cframes,
                  const double* gnum, const double* specg1, const double* specg2, const double* dgnum,
                  int a, int M, int N, const double* sframes, int lookahead, int maxit,
                  double* cframes2, double complex* c)
{
    rtisilaupdate_plan *p = rtisilaupdate_init(gnum, specg1, specg2, dgnum, a, M);
    rtisilaupdate_execute(p, cframes, N, sframes, lookahead, maxit, cframes2);
    memcpy(c, p->fftframe, (M / 2 + 1)*sizeof * p->fftframe);
    rtisilaupdate_done(p);
}
