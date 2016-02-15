#include "spsi.h"
#include "float.h"

void
spsireal(const double* s, int a, int M, int N, complex double* c,
         double* initphase)
{
    int M2 = M / 2 + 1;

    double* tmpphase = initphase;

    if (!initphase)
    {
        tmpphase = aligned_alloc(ALIGNBYTES, M2 * sizeof * tmpphase);
        memset(tmpphase, 0, M2 * sizeof * tmpphase);
    }


    for (int n = 0; n < N; n++)
    {
        const double* scol = s + n * M2;
        complex double* ccol = c + n * M2;

        for (int m = 1; m < M2 - 1; m++)
        {
            if (scol[m] > scol[m - 1] && scol[m] > scol[m + 1])
            {
                double p;
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
                // Go towards low frequency bins
                int bin = m - 1;

                while (bin > 0 && scol[bin] < scol[bin + 1])
                {
                    tmpphase[bin] = peakPhase;
                    bin--;
                }

                bin = m + 1;

                while (bin < M2 - 1 && scol[bin] < scol[bin - 1])
                {
                    tmpphase[bin] = peakPhase;
                    bin++;
                }
            }
        }

        for (int m = 0; m < M2; m++)
            ccol[m] = scol[m] * cexp(I * tmpphase[m]);
    }

    if (!initphase)
        free(tmpphase);
}



void
spsireal_split(const double* cr, const double* ci, int a, int M, int N,
               double* mask, double* coutr, double* couti)
{
    int M2 = M / 2 + 1;
    double* cabsPtr;
    double* tmpphase = aligned_alloc(ALIGNBYTES, M2 * sizeof * tmpphase);
    memset(tmpphase, 0, M2 * sizeof * tmpphase);

    if (ci == NULL)
    {
        cabsPtr = (double*)cr;
    }
    else
    {
        cabsPtr = aligned_alloc(ALIGNBYTES, M2 * N * sizeof * cabsPtr);

        for (int w = 0; w < M2 * N; w++)
        {
            cabsPtr[w] = cabs(cr[w] + I * ci[w]);
        }
    }

    for (int n = 0; n < N; n++)
    {
        const double* crcol = cr + n * M2;
        const double* cicol = ci + n * M2;
        double* cabscol = cabsPtr + n * M2;
        double* coutrcol = coutr + n * M2;
        double* couticol = couti + n * M2;

        for (int m = 1; m < M2 - 1; m++)
        {
            if (cabscol[m] > cabscol[m - 1] && cabscol[m] > cabscol[m + 1])
            {
                double p;
                double alpha = cabscol[m - 1];
                double beta = cabscol[m];
                double gamma = cabscol[m + 1];
                double denom = alpha - 2.0 * beta + gamma;
                int bin;

                if (denom != 0.0)
                {
                    p = 0.5 * (alpha - gamma) / denom;
                }
                else
                {
                    p = 0;
                }

                double instf = m + p;
                double peakPhase = tmpphase[m] + 2.0 * M_PI * a * instf / M;
                tmpphase[m] = peakPhase;


                if (p > 0.0)
                {
                    bin = m + 1;
                    tmpphase[bin] = peakPhase;
                    bin = m - 1;

                    while (bin > 0 && cabscol[bin] < cabscol[bin + 1])
                    {
                        tmpphase[bin] = peakPhase;
                        bin--;
                    }

                    bin = m + 2;

                    while (bin < M2 - 1 && cabscol[bin] < cabscol[bin - 1])
                    {
                        tmpphase[bin] = peakPhase;
                        bin++;
                    }
                }
                else if (p < 0.0)
                {
                    bin = m - 1;
                    tmpphase[bin] = peakPhase;
                    bin = m + 1;

                    while (bin < M2 - 1 && cabscol[bin] < cabscol[bin - 1])
                    {
                        tmpphase[bin] = peakPhase;
                        bin++;
                    }

                    bin = m - 2;

                    while (bin > 0 && cabscol[bin] < cabscol[bin + 1])
                    {
                        tmpphase[bin] = peakPhase;
                        bin--;
                    }

                }
            }
        }

        if (mask)
        {
            double* maskcol = mask + n * M2;

            for (int m = 0; m < M2; m++)
            {
                if (maskcol[m])
                {
                    tmpphase[m] = carg(crcol[m] + I * cicol[m]);
                }
            }
        }

        for (int m = 0; m < M2; m++)
        {
            double complex tmpz = cabscol[m] * cexp(I * tmpphase[m]);
            coutrcol[m] = creal(tmpz);
            couticol[m] = cimag(tmpz);
        }

    }

    if (cr != cabsPtr) free(cabsPtr);

    free(tmpphase);
}
