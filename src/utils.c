#include "phaseret/utils.h"

void
shiftcolsleft(double* cols, int height, int N, const double* newcol)
{
    for (int n = 0; n < N - 1; n++)
        memcpy(cols + n * height, cols + (n + 1)*height, height * sizeof * cols);

    if (newcol)
        memcpy(cols + (N - 1)*height, newcol, height * sizeof * cols);
    else
        memset(cols + (N - 1)*height, 0, height * sizeof * cols);
}

int
force_magnitude(complex double* cin, const double* s, int L, complex double* cout)
{
    for (int ii = 0; ii < L; ii++)
        cout[ii] = s[ii] * cexp(I * carg(cin[ii]));

}
