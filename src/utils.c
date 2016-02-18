#include "utils.h"

void
fftrealcircshift(const complex double* in, int L, double shift,
                 complex double* out)
{
    const div_t domod = div(L, 2);

    if (shift != 0)
    {
        double complex phasefact = -I * 2.0 * M_PI * shift / L;

        out[0] = in[0];

        for (int ii = 1; ii < domod.quot + 1; ii++)
            out[ii] = cexp(ii * phasefact) * in[ii];
    }
    else
    {
        if (in != out)
            memcpy(out, in, domod.quot * sizeof * out);
    }
}

void
fftrealfftshift(const double complex* in, int L, double complex* out)
{
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        fftrealcircshift(in, L, domod.quot, out);
    }
    else
    {
        out[0] = in[0];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < domod.quot + 1; ii += 2)
            out[ii] = -1.0 * in[ii];
    }
}

void
fftrealifftshift(const double complex* in, int L, double complex* out)
{
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        fftrealcircshift(in, L, -domod.quot, out);
    else
        fftrealfftshift(in, L, out);
}


void
circshift(const double* in, int L, int shift, double* out)
{
    // Fix shift
    int p = (L - shift) % L;

    if (p < 0) p += L;

    if (in == out)
    {
        int m, count, i, j;

        // Circshift inplace is magic!
        for (m = 0, count = 0; count != L; m++)
        {
            double t = in[m];

            for (i = m, j = m + p; j != m;
                 i = j, j = j + p < L ? j + p : j + p - L, count++)
                out[i] = out[j];

            out[i] = t; count++;
        }
    }
    else
    {
        // Circshit out of place is boring ...
        memcpy(out, in + p, (L - p)*sizeof * out);
        memcpy(out + L - p, in, p * sizeof * out);
    }
}

// in might be equal to out
void
fftshift(const double* in, int L, double* out)
{
    circshift(in, L, (L / 2), out);
}

// in might be equal to out
void
ifftshift(const double* in, int L, double* out)
{
    circshift(in, L, -(L / 2), out);
}

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
