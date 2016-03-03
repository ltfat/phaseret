#include "minunit.h"
#include <sndfile.h>
#include <dlfcn.h>
#include <math.h>
#include <phaseret.h>
#include <ltfat.h>

void* lib = NULL;
double* f;
double* frec;
complex double* c;
double* s;
double* s2;
int L;
int a = 128;
int M = 1024;
int M2 = 512 + 1;
int W = 1;
int N;
double* g;
double* gd;


char* test_openwav()
{
    SF_INFO info;
    SNDFILE* sf = sf_open("greasy.wav", SFM_READ, &info);
    mu_assert(sf != NULL, "Failed to load the wav file");

    int Ls = info.frames * info.channels;
    L = (int) (ceil(((double)(Ls)) / M) * M);
    N = L / a;
    debug("L=%d", L);
    debug("N=%d", N);
    g = malloc(L * sizeof * g);
    gd = malloc(L * sizeof * g);
    frec = calloc(L , sizeof * frec);
    c = malloc(M2 * N * sizeof * c);
    s = malloc(M2 * N * sizeof * c);
    s2 = malloc(M2 * N * sizeof * c);


    f = calloc(L, sizeof * f);
    check(f != NULL, "Memory allocation error");

    int num = sf_read_double(sf, f, Ls);

    /* for(int ii=0;ii<L;ii++) */
    /*     debug("%f",f[ii]); */
    /*  */
    sf_close(sf);

error:
    return NULL;
}


char* test_functions()
{
    pgauss_d(L, a * M / L, 0, g);

    // Compute the dual window
    gabdual_long_d(g, L, 1, a, M, gd);

    // Do the transform
    dgtreal_long_d(f, g, L, W, a, M, TIMEINV, c);

    for (int ii = 0; ii < M2 * N; ii++)
        s[ii] = cabs(c[ii]);

    //spsi(s, a, M, 1, NULL , c);
    rtpghioffline(s, a * M, a, M, L, 1e-6, 0, c);

    // Do inverse transform
    idgtreal_long_d(c, gd, L, W, a, M, TIMEINV, frec);

    dgtreal_long_d(frec, g, L, W, a, M, TIMEINV, c);

    for (int ii = 0; ii < M2 * N; ii++)
        s2[ii] = cabs(c[ii]);

    double snorm = 0.0;

    for (ltfatInt l = 0; l < M2 * N; ++l)
    {
        snorm += s[l] * s[l];
    }


    double err = 0.0;

    for (ltfatInt l = 0; l < M2 * N; ++l)
    {
        double dif = fabs( s[l] - s2[l] );
        err += dif * dif;
    }

    printf("Error of the reconstruction is %f\n", 10.0 * log10(err / snorm));

    return NULL;
}

char* test_failures()
{

    return NULL;
}

char* test_free()
{
    LTFAT_SAFEFREEALL(f, g, gd, frec, c, s, s2);

    return NULL;
}

char* all_tests()
{
    mu_suite_start();

    mu_run_test(test_openwav);
    mu_run_test(test_functions);
    mu_run_test(test_failures);
    mu_run_test(test_free);

    return NULL;
}

RUN_TESTS(all_tests);
