#include "ltfat.h"
#include "phaseret/gla.h"


int
gla_long(double* s, double* g, int L, int W, int a, int M,
         int iter, complex double* c)
{
    int N = L / a;
    int M2 = M / 2 + 1;


    double* gd = malloc(L * sizeof * gd);
    double* f = malloc(L * sizeof * f);

    gabdual_long_d(g, L, 1, a, M, gd);

    ltfat_dgtreal_long_plan_d* dgtPlan = NULL;
    ltfat_idgtreal_long_plan_d* idgtPlan = NULL;

    ltfat_idgtreal_long_init_d(c, gd, L, W, a, M, f, TIMEINV, FFTW_MEASURE,
                               &idgtPlan);
    ltfat_dgtreal_long_init_d(f, g, L, W, a, M, c, TIMEINV, FFTW_MEASURE, &dgtPlan);

    free(gd);

    ltfat_ensurecomplex_array_d(s, N * M2 * W, c);

    for (int ii = 0; ii < iter; ii++)
    {
        ltfat_idgtreal_long_execute(idgtPlan);
        // Callback (userdata,f,L,W,f)
        ltfat_dgtreal_long_execute(dgtPlan);
        force_magnitude(c, s, N * M2 * W, c);
        // Callback (userdata, c, M2, N, c)
    }

    ltfat_dgt_long_done_d(&dgtPlan);
    ltfat_idgt_long_done_d(&idgtPlan);


    free(f);
}

int
fgla_long(double* s, double* g, int L, int W, int a, int M,
          int iter, double alpha, complex double* c)
{
    int N = L / a;
    int M2 = M / 2 + 1;
    double* gd = malloc(L * sizeof * gd);
    double* f = malloc(L * W * sizeof * f);
    complex double* t = malloc(M2 * N  * W * sizeof * t);


    gabdual_long_d(g, L, 1, a, M, gd);

    ltfat_dgtreal_long_plan_d* dgtPlan = NULL;
    ltfat_idgtreal_long_plan_d* idgtPlan = NULL;

    ltfat_idgtreal_long_init_d(c, gd, L, W, a, M, f, TIMEINV, FFTW_MEASURE,
                               &idgtPlan);
    ltfat_dgtreal_long_init_d(f, g, L, W, a, M, c, TIMEINV, FFTW_MEASURE, &dgtPlan);

    free(gd);

    ltfat_ensurecomplex_array_d(s, N * M2 *  W, c);
    memcpy(t, c, (N * M2 *  W) * sizeof * f );

    for (int ii = 0; ii < iter; ii++)
    {
        ltfat_idgtreal_long_execute(idgtPlan);
        ltfat_dgtreal_long_execute(dgtPlan);
        force_magnitude(c, s, N * M2 * W, c);
        fastupdate(c, t, alpha, N * M2 * W );
    }

    ltfat_dgt_long_done_d(&dgtPlan);
    ltfat_idgt_long_done_d(&idgtPlan);


    free(f);
    free(t);
}

int
force_magnitude(complex double* cin, double* s, int L, complex double* cout)
{
    for (int ii = 0; ii < L; ii++)
        cout[ii] = s[ii] * cexp(I * carg(cin[ii]));

}

int
fastupdate(complex double* c, complex double* t, double alpha, int L)
{
    for (int ii = 0; ii < L; ii++)
    {
        complex double cold = c[ii];
        c[ii] = c[ii] + alpha * (c[ii] - t[ii]);
        t[ii] = cold[ii];
    }
}
