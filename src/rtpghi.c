#include "phaseret/rtpghi.h"
#include "phaseret/config.h"
#include "phaseret/utils.h"
#include "float.h"

rtpghi_plan*
rtpghi_init(double gamma, int a, int M, double tol, int do_causal)
{
    int M2 = M / 2 + 1;

    struct ltfat_heapinttask_d* hit = ltfat_heapinttask_init_d( M2, 2, 2 * M2 , NULL, 1);
    double* slog = calloc(3 * M2 , sizeof * slog);
    double* s = calloc(2 * M2 , sizeof * slog);
    double* phase = calloc(2 * M2,  sizeof * phase);
    double* tgrad = calloc(3 * M2,  sizeof * tgrad);
    double* fgrad = calloc(2 * M2,  sizeof * fgrad);

    int randphaseLen = 10 * M2;
    double* randphase = malloc(randphaseLen * sizeof * randphase);

    for (int ii = 0; ii < randphaseLen; ii++)
        randphase[ii] = 2 * M_PI * ((double)rand()) / RAND_MAX;

    int* mask = malloc(2 * M2 * sizeof * mask);

    for (int ii = 0; ii < M2; ++ii)
    {
        mask[ii] = 1;
        mask[ii + M2] = 0;
    }

    // This way, we can write to the constant parameters of the structure
    rtpghi_plan pdummy =
    {
        .hit = hit, .do_causal = do_causal, .logtol = log(tol + DBL_MIN), .tol = tol,
        .slog = slog, .tgrad = tgrad, .fgrad = fgrad, .phase = phase, .M = M, .a = a, .gamma = gamma,
        .randphase = randphase, .randphaseLen = randphaseLen, .randphaseId = 0, .s = s,
        .mask = mask
    };

    rtpghi_plan* p = malloc(sizeof * p);
    memcpy(p, &pdummy, sizeof * p);
    return p;
}


void
rtpghi_execute(rtpghi_plan* p, const double* s, complex double* c)
{
    // n, n-1, n-2 frames
    // s is n-th
    int M2 = p->M / 2 + 1;

    // store log(s)
    shiftcolsleft(p->slog, M2, 3, NULL);
    shiftcolsleft(p->s, M2, 2, s);
    rtpghilog(p->s + M2, M2, p->slog + 2 * M2);

    // Store current s

    // Compute and store tgrad for n
    shiftcolsleft(p->tgrad, M2, 3, NULL);
    rtpghitgrad(p->slog + 2 * M2, p->a, p->M, p->gamma, p->tgrad + 2 * M2);

    // Compute fgrad for n or n-1
    rtpghifgrad(p->slog, p->a, p->M, p->gamma, p->do_causal, p->fgrad + M2);

    ltfat_heapinttask_resetmask_d(p->hit,
                                  p->mask, p->do_causal ? p->slog + M2 : p->slog, p->logtol, 1);

    ltfat_heapint_execute_d(p->hit,
                            p->do_causal ? p->tgrad + M2 : p->tgrad, p->fgrad, p->phase);

    shiftcolsleft(p->phase, M2, 2, NULL);

    // Fill in the missing phase from the precomputed random array
    for (int ii = 0; ii < M2; ii++)
    {
        if (p->hit->donemask[M2 + ii] == 5)
        {
            p->phase[ii] = p->randphase[p->randphaseId++];
            p->randphaseId %= p->randphaseLen;
        }
    }

    // Combine phase with magnitude
    rtpghimagphase(p->do_causal ? p->s + M2 : p->s, p->phase, M2, c);
}


void
rtpghi_done(rtpghi_plan* p)
{
    ltfat_heapinttask_done_d(p->hit);
    free(p->slog);
    free(p->s);
    free(p->mask);
    free(p->phase);
    free(p->tgrad);
    free(p->fgrad);
    free(p->randphase);
    free(p);
}

void
rtpghioffline(const double* s, double gamma, int a, int M, int L, double tol,
              int do_causal, complex double* c)
{
    int N = L / a;
    int M2 = M / 2 + 1;
    rtpghi_plan* p = rtpghi_init(gamma, a, M, tol, do_causal);

    if (do_causal)
    {
        for (int n = 0; n < N; ++n)
        {
            const double* sncol = s + n * M2;
            complex double* cncol = c + n * M2;
            rtpghi_execute(p, sncol, cncol);
        }
    }
    else
    {
        rtpghi_execute(p, s, c);

        for (int n = 0, nahead = 1; nahead < N; ++n, ++nahead)
        {
            const double* sncol = s + nahead * M2;
            complex double* cncol = c + n * M2;
            rtpghi_execute(p, sncol, cncol);
        }

        rtpghi_execute(p, NULL, c + (N - 1) * M2);
    }

    rtpghi_done(p);
}


void
rtpghifgrad(const double* logs, int a, int M, double gamma,
            int do_causal, double* fgrad)
{
    int M2 = M / 2 + 1;

    const double fgradmul = -gamma / (2.0 * a * M);
    const double* scol0 = logs;
    const double* scol2 = logs + 2 * M2;

    if (do_causal)
    {
        const double* scol1 = logs + M2;

        for (int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (3.0 * scol2[m] - 4.0 * scol1[m] + scol0[m]);
    }
    else
    {
        for (int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (scol2[m] - scol0[m]);
    }
}

void
rtpghitgrad(const double* logs, int a, int M, double gamma,
            double* tgrad)
{
    int M2 = M / 2 + 1;

    const double tgradmul = (a * M) / (gamma * 2.0);
    const double tgradplus = 2.0 * M_PI * a / M;

    tgrad[0]      = 0.0;
    tgrad[M2 - 1] = 0.0;

    for (int m = 1; m < M2 - 1; m++)
        tgrad[m] = tgradmul * (logs[m + 1] - logs[m - 1]) + tgradplus * m;
}

void
rtpghilog(const double* in, int L, double* out)
{
    for (int l = 0; l < L; l++)
        out[l] = log(in[l] + DBL_MIN);

}

void
rtpghimagphase(const double* s, const double* phase, int L, complex double* c)
{
    for (int l = 0; l < L; l++)
        c[l] = s[l] * cexp(I * phase[l]);
}
