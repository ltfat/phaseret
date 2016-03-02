#include "spsi.h"
#include "float.h"

void
spsi(const double* s, int a, int M, int N,  double* initphase,
         complex double* c)
{
    int M2 = M / 2 + 1;

    double* tmpphase = initphase;

    if (!initphase)
    {
        tmpphase = malloc( M2 * sizeof * tmpphase);
        memset(tmpphase, 0, M2 * sizeof * tmpphase);
    }

    for (int n = 0; n < N; n++)
    {
        const double* scol = s + n * M2;
        complex double* ccol = c + n * M2;

        spsiupdate(scol, a, M, tmpphase);

        for (int m = 0; m < M2; m++)
            ccol[m] = scol[m] * cexp(I * tmpphase[m]);
    }

    if (!initphase)
        free(tmpphase);
}

void
maskedspsi(const double* s, int a, int M, int N,
               const double* mask, const double* phase, double* initphase, complex double* c)
{
    int M2 = M / 2 + 1;

    double* tmpphase = initphase;

    if (!initphase)
    {
        tmpphase = malloc( M2 * sizeof * tmpphase);
        memset(tmpphase, 0, M2 * sizeof * tmpphase);
    }

    for (int n = 0; n < N; n++)
    {
        const double* scol = s + n * M2;
        complex double* ccol = c + n * M2;
        const double* maskcol = mask + n * M2;
        const double* phasecol = phase + n * M2;

        spsiupdate(scol, a, M, tmpphase);

        /* Overwrite with known phase */
        for (int m = 0; m < M2; m++)
            if (maskcol[m]) tmpphase[m] = phasecol[m];

        for (int m = 0; m < M2; m++)
            ccol[m] = scol[m] * cexp(I * tmpphase[m]);
    }

    if (!initphase)
        free(tmpphase);
}

void
spsiupdate(const double* scol, int a, int M, double* tmpphase)
{
    int M2 = M / 2 + 1;

    for (int m = 1; m < M2 - 1; m++)
    {
        if (scol[m] > scol[m - 1] && scol[m] > scol[m + 1])
        {
            double p; int binup = m, bindown = m;
            double alpha = log(scol[m - 1] + DBL_MIN);
            double beta = log(scol[m] + DBL_MIN);
            double gamma = log(scol[m + 1] + DBL_MIN);
            double denom = alpha - 2.0 * beta + gamma;

            if (denom != 0.0)
                p = 0.5 * (alpha - gamma) / denom;
            else
                p = 0;


            double instf = m + p;
            double peakPhase = tmpphase[m] + 2.0 * M_PI * a * instf / M;
            tmpphase[m] = peakPhase;

            if (p > 0)
            {
                tmpphase[m + 1] = peakPhase;
                binup = m + 2;
                bindown = m - 1;
            }

            if (p < 0)
            {
                tmpphase[m - 1] = peakPhase;
                binup = m + 1;
                bindown = m - 2;
            }

            // Go towards low frequency bins
            int bin = bindown;

            while (bin > 0 && scol[bin] < scol[bin + 1])
            {
                tmpphase[bin] = peakPhase;
                bin--;
            }

            // Go towards high frequency bins
            bin = binup;

            while (bin < M2 - 1 && scol[bin] < scol[bin - 1])
            {
                tmpphase[bin] = peakPhase;
                bin++;
            }
        }
    }


}
