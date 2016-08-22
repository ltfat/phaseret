#include "phaseret/utils.h"

int
shiftcolsleft(LTFAT_REAL* cols, int height, int N, const LTFAT_REAL* newcol)
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
force_magnitude(LTFAT_COMPLEX* cin, const LTFAT_REAL* s, int L,
                LTFAT_COMPLEX* cout)
{
    for (int ii = 0; ii < L; ii++)
        cout[ii] = s[ii] * exp(I * ltfat_arg(cin[ii]));

    return 0;
}


void
realimag2absangle(const LTFAT_COMPLEX* cin, const int L, LTFAT_COMPLEX* c)
{
    LTFAT_REAL* cplain = (LTFAT_REAL*) c;

    for (int l = 0; l < L; l++)
    {
        LTFAT_COMPLEX cel = cin[l];
        cplain[2 * l] = ltfat_abs(cel);
        cplain[2 * l + 1] = ltfat_arg(cel);
    }
}

void
absangle2realimag(const LTFAT_COMPLEX* cin, const int L, LTFAT_COMPLEX* c)
{
    LTFAT_REAL* cinplain = (LTFAT_REAL*) cin;

    for (int l = 0; l < L; l++)
    {
        LTFAT_REAL absval = cinplain[2 * l];
        LTFAT_REAL phaseval = cinplain[2 * l + 1];
        c[l] = absval * exp(I * phaseval);
    }

}
