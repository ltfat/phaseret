#include "rtpghi.h"
#include "float.h"

struct _rtpghiplan
{
    int* mask;
    struct heapinttask_d* hit;
    int do_causal;
    double* phase;
};

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
        c[l] = s[l]*cexp(I*phase[l]);
}
