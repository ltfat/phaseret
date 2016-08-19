#include "phaseret/spsi.h"
#include "phaseret/utils.h"
#include "ltfat.h"
#include "ltfat/macros.h"
#include "float.h"

int
spsi(const double* s, int L, int W, int a, int M, double* initphase,
     complex double* c)
{
    int M2 = M / 2 + 1;
    int N = L / a;
    double* tmpphase = initphase;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(s); CHECKNULL(c);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    if (!initphase)
        CHECKMEM( tmpphase = calloc( M2 * W , sizeof * tmpphase) );

    if (s == (double*) c)
    {
        // Inplace, move the abs. values to the second half of the array
        double* chalf = ((double*) c) + W * M2 * N;
        memcpy(chalf, s, W * M2 * N * sizeof * chalf);
        s = chalf;
    }

    for (int w = 0; w < W; w++)
    {
        double* tmpphasecol = tmpphase + w * M2;
        for (int n = 0; n < N; n++)
        {
            const double* scol = s + n * M2 + w * M2 * N;
            complex double* ccol = c + n * M2 + w * M2 * N;

            spsiupdate(scol, 1, a, M, tmpphasecol);

            for (int m = 0; m < M2; m++)
                ccol[m] = scol[m] * cexp(I * tmpphasecol[m]);

        }
    }

error:
    if (!initphase) free(tmpphase);
    return status;
}

int
spsi_withmask(const complex double* cinit, const int* mask, int L, int W, int a,
              int M, double* initphase, complex double* c)
{
    int M2 = M / 2 + 1;
    int N = L / a;
    double* tmpphase = initphase;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(cinit); CHECKNULL(mask); CHECKNULL(c);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    if (!initphase)
        tmpphase = calloc( M2 * W, sizeof * tmpphase);

    for (int w = 0; w < W; w++)
    {
        double* tmpphasecol = tmpphase + w * M2;
        for (int n = 0; n < N; n++)
        {
            complex double* ccol = c + n * M2 + w * M2 * N;
            const complex double* cinitcol = cinit + n * M2 + w * M2 * N;
            const int* maskcol = mask + n * M2 + w * M2 * N;

            realimag2absangle(cinitcol, M2, ccol);
            double* absptr = (double*) ccol;
            double* angleptr = ((double*) ccol) + 1;

            spsiupdate(absptr, 2, a, M, tmpphasecol);

            /* Overwrite with known phase */
            for (int m = 0; m < M2; m++)
                if (maskcol[m]) tmpphasecol[m] = angleptr[2 * m];

            for (int m = 0; m < M2; m++)
                ccol[m] = absptr[2 * m] * cexp(I * tmpphasecol[m]);
        }
    }

    if (!initphase)
        free(tmpphase);

error:
    return status;
}

void
spsiupdate(const double* scol, int stride, int a, int M, double* tmpphase)
{
    int M2 = M / 2 + 1;

    for (int m = 1; m < M2 - 1; m++)
    {
        if (scol[stride * m] > scol[stride * (m - 1)]
            && scol[stride * m] > scol[stride * (m + 1)])
        {
            double p; int binup = m, bindown = m;
            double alpha = log(scol[stride * (m - 1)] + DBL_MIN);
            double beta = log(scol[stride * m] + DBL_MIN);
            double gamma = log(scol[stride * (m + 1)] + DBL_MIN);
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

            while (bin > 0 && scol[stride * bin] < scol[stride * (bin + 1)])
            {
                tmpphase[bin] = peakPhase;
                bin--;
            }

            // Go towards high frequency bins
            bin = binup;

            while (bin < M2 - 1 && scol[stride * bin] < scol[stride * (bin - 1)])
            {
                tmpphase[bin] = peakPhase;
                bin++;
            }
        }
    }


}
