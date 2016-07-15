#include "ltfat.h"
#include "ltfat/macros.h"
#include "phaseret/gla.h"
#include "phaseret/utils.h"
#include "dgtrealwrapper_private.h"

struct gla_plan
{
    dgtreal_anasyn_plan* p;
    gla_callback_status* status_callback;
    void* status_callback_userdata;
    gla_callback_cmod* cmod_callback;
    void* cmod_callback_userdata;
    gla_callback_fmod* fmod_callback;
    void* fmod_callback_userdata;
// Used just for fgla
    int do_fast;
    double alpha;
    complex double* t;
};

int
gla(const double s[], const double g[], const int gl, const int L, const int W,
    const int a, const int M, const int iter, complex double c[])
{
    gla_plan* p = NULL;
    gla_init(g, gl, L, W, a, M, 0.99, c, dgtreal_anasyn_auto, FFTW_ESTIMATE, &p);
    gla_execute(p, s, iter);
    gla_done(&p);
}

int
gla_init(const double g[], const int gl, const int L, const int W,
         const int a, const int M, const double alpha,
         complex double c[], dgtreal_anasyn_hint hint, unsigned flags,
         gla_plan** pout)
{
    int N = L / a;
    int M2 = M / 2 + 1;

    gla_plan* p = calloc(1, sizeof * p);
    dgtreal_anasyn_init(g, gl, L, W, a, M, c, hint, flags, &p->p);

    if (alpha > 0.0)
    {
        p->do_fast = 1;
        p->alpha = alpha;
        p->t = malloc(M2 * N * W * sizeof * p->t);
    }

    *pout = p;
}

int
gla_done(gla_plan** p)
{
    gla_plan* pp = *p;
    dgtreal_anasyn_done(&pp->p);
    if (pp->t) free(pp->t);
    free(pp);
    pp = NULL;
}

int
gla_execute_newarray(gla_plan* p, const double s[], const int iter,
                     complex double c[])
{
    // Shallow copy the plan and replace c
    gla_plan p2 = *p;
    dgtreal_anasyn_plan pp2 = *p2.p;
    pp2.c = c;
    p2.p = &pp2;
    gla_execute(&p2, s, iter);
}


int
gla_execute(gla_plan* p, const double s[], const int iter)
{
    dgtreal_anasyn_plan* pp = p->p;
    int M = pp->M;
    int L = pp->L;
    int W = pp->W;
    int a = pp->a;
    int M2 = M / 2 + 1;
    int N = L / a;

    ltfat_ensurecomplex_array_d( s, N * M2 * W, pp->c);

    if (p->do_fast)
        memcpy(p->t, pp->c, (N * M2 * W) * sizeof * p->t );

    for (int ii = 0; ii < iter; ii++)
    {
        // Perform idgtreal
        dgtreal_anasyn_execute_syn(pp, pp->c, pp->f);

        // Optional signal modification
        if (p->fmod_callback)
            p->fmod_callback(p->fmod_callback_userdata, pp->f, L, W, a, M);

        // Perform dgtreal
        dgtreal_anasyn_execute_ana(pp, pp->f, pp->c);

        force_magnitude(pp->c, s, N * M2 * W, pp->c);

        if (p->do_fast)
            fastupdate(pp->c, p->t, p->alpha, N * M2 * W );

        // Optional coefficient modification
        if (p->cmod_callback)
            p->cmod_callback(p->cmod_callback_userdata, pp->c, L, W, a, M);

        // Status callback, optional premature exit
        if (p->status_callback)
        {
            int status = p->status_callback(pp, p->status_callback_userdata,
                                            pp->c, L, W, a, M, &p->alpha, ii);
            if (status > 0)
                break;
        }
    }

}


int
gla_set_status_callback(gla_plan* p,
                        gla_callback_status* callback, void* userdata)
{
    p->status_callback = callback;
    p->status_callback_userdata = userdata;
}

int
gla_set_cmod_callback(gla_plan* p,
                      gla_callback_cmod* callback, void* userdata)
{
    p->cmod_callback = callback;
    p->cmod_callback_userdata = userdata;
}

int
gla_set_fmod_callback(gla_plan* p,
                      gla_callback_fmod* callback, void* userdata)
{
    p->fmod_callback = callback;
    p->fmod_callback_userdata = userdata;
}

int
fastupdate(complex double* c, complex double* t, double alpha, int L)
{
    for (int ii = 0; ii < L; ii++)
    {
        complex double cold = c[ii];
        c[ii] = c[ii] + alpha * (c[ii] - t[ii]);
        t[ii] = cold;
    }
}
