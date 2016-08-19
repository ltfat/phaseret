#include "phaseret/utils.h"
#include <tgmath.h>

int
shiftcolsleft(double* cols, int height, int N, const double* newcol)
{
    for (int n = 0; n < N - 1; n++)
        memcpy(cols + n * height, cols + (n + 1)*height, height * sizeof * cols);

    if (newcol)
        memcpy(cols + (N - 1)*height, newcol, height * sizeof * cols);
    else
        memset(cols + (N - 1)*height, 0, height * sizeof * cols);

    return 0;
}

int
force_magnitude(complex double* cin, const double* s, int L,
                complex double* cout)
{
    for (int ii = 0; ii < L; ii++)
        cout[ii] = s[ii] * cexp(I * carg(cin[ii]));

    return 0;
}


void
realimag2absangle(const complex double* cin, const int L, complex double* c)
{
    double* cplain = (double*) c;

    for (int l = 0; l < L; l++)
    {
        complex double cel = cin[l];
        cplain[2 * l] = cabs(cel);
        cplain[2 * l + 1] = carg(cel);
    }
}

void
absangle2realimag(const complex double* cin, const int L, complex double* c)
{
    double* cinplain = (double*) cin;

    for (int l = 0; l < L; l++)
    {
        double absval = cinplain[2 * l];
        double phaseval = cinplain[2 * l + 1];
        c[l] = absval * cexp(I * phaseval);
    }

}
