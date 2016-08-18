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
// For storing magnitude
    double* s;
// Storing cinit
    const complex double* cinit;
// Used just for fgla
    int do_fast;
    double alpha;
    complex double* t;
};

int
gla(const complex double cinit[], const double g[],  const int L, const int gl,
    const int W,
    const int a, const int M, const int iter, complex double cout[])
{
    gla_plan* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        gla_init(cinit, g, L, gl, W, a, M, 0.99, cout,
                 dgtreal_anasyn_auto, FFTW_ESTIMATE, &p), "Init failed");
    CHECKSTATUS( gla_execute(p, iter), "Execute failed");

error:
    if (p) gla_done(&p);
    return status;
}

int
gla_init(const complex double cinit[], const double g[], const int L,
         const int gl, const int W, const int a, const int M, const double alpha,
         complex double c[], dgtreal_anasyn_hint hint, unsigned flags,
         gla_plan** pout)
{
    int status = LTFATERR_SUCCESS;
    gla_plan* p = NULL;
    int N = L / a;
    int M2 = M / 2 + 1;

    CHECK(LTFATERR_BADARG, alpha >= 0.0, "alpha cannot be negative");
    CHECKMEM( p = calloc(1, sizeof * p));
    CHECKMEM( p->s = malloc(M2 * N * W * sizeof * p->t));

    if (alpha > 0.0)
    {
        p->do_fast = 1;
        p->alpha = alpha;
        CHECKMEM( p->t = malloc(M2 * N * W * sizeof * p->t));
    }

    CHECKSTATUS(
        dgtreal_anasyn_init(g, gl, L, W, a, M, c, hint, LTFAT_TIMEINV, flags, &p->p),
        "dgtrealwrapper init failed");

    p->cinit = cinit;

    *pout = p;
    return status;
error:
    if (p)
    {
        if (p->t) free(p->t);
        if (p->s) free(p->s);
        free(p);
    }
    *pout = NULL;
    return status;
}

int
gla_done(gla_plan** p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    gla_plan* pp = *p;

    CHECKSTATUS(
        dgtreal_anasyn_done(&pp->p),
        "dgtreal wrapper done failed");

    if (pp->t) free(pp->t);
    free(pp->s);
    free(pp);
    pp = NULL;
error:
    return status;
}

int
gla_execute_newarray(gla_plan* p, const complex double cinit[], const int iter,
                     complex double cout[])
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(cinit); CHECKNULL(cout);
    // Shallow copy the plan and replace c
    gla_plan p2 = *p;
    dgtreal_anasyn_plan pp2 = *p2.p;
    pp2.c = cout;
    p2.p = &pp2;
    p2.cinit = cinit;
    return gla_execute(&p2, iter);
error:
    return status;
}

int
gla_execute(gla_plan* p, const int iter)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, iter > 0,
          "At least one iteration is requred. Passed %d.", iter);

    dgtreal_anasyn_plan* pp = p->p;
    CHECKNULL(pp->c); CHECKNULL(p->cinit);
    int M = pp->M;
    int L = pp->L;
    int W = pp->W;
    int a = pp->a;
    int M2 = M / 2 + 1;
    int N = L / a;

    // Store the magnitude
    for (int ii = 0; ii < N * M2 * W; ii++)
        p->s[ii] = cabs(p->cinit[ii]);

    // Copy to the output array if we are not working inplace
    if (p->cinit != pp->c)
        memcpy(pp->c, p->cinit, (N * M2 * W) * sizeof * pp->c);

    // Inicialize the "acceleration" array
    if (p->do_fast)
        memcpy(p->t, pp->c, (N * M2 * W) * sizeof * p->t );

    for (int ii = 0; ii < iter; ii++)
    {
        // Perform idgtreal
        CHECKSTATUS( dgtreal_anasyn_execute_syn(pp, pp->c, pp->f),
                     "idgtreal failed");

        // Optional signal modification
        if (p->fmod_callback)
            CHECKSTATUS(
                p->fmod_callback(p->fmod_callback_userdata, pp->f, L, W, a, M),
                "fmod callback failed");

        // Perform dgtreal
        CHECKSTATUS( dgtreal_anasyn_execute_ana(pp, pp->f, pp->c),
                     "dgtreal failed");

        force_magnitude(pp->c, p->s, N * M2 * W, pp->c);

        // The acceleration step
        if (p->do_fast)
            fastupdate(pp->c, p->t, p->alpha, N * M2 * W );

        // Optional coefficient modification
        if (p->cmod_callback)
            CHECKSTATUS(
                p->cmod_callback(p->cmod_callback_userdata, pp->c, L, W, a, M),
                "cmod callback failed");

        // Status callback, optional premature exit
        if (p->status_callback)
        {
            int retstatus = p->status_callback(pp, p->status_callback_userdata,
                                               pp->c, L, W, a, M, &p->alpha, ii);
            if (retstatus > 0)
                break;
            else
                CHECKSTATUS(retstatus, "Status callback failed");

            CHECK(LTFATERR_BADARG, p->alpha >= 0.0, "alpha cannot be negative");

            if (p->alpha > 0.0 && !p->do_fast)
            {
                // The plan was not inicialized with acceleration but
                // nonzero alpha was set in the status callback.
                p->do_fast = 1;
                CHECKMEM( p->t = malloc(M2 * N * W * sizeof * p->t));
                memcpy(p->t, pp->c, (N * M2 * W) * sizeof * p->t );
            }
        }
    }

error:
    return status;
}


int
gla_set_status_callback(gla_plan* p,
                        gla_callback_status* callback, void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(callback);

    p->status_callback = callback;
    p->status_callback_userdata = userdata;
error:
    return status;
}

int
gla_set_cmod_callback(gla_plan* p,
                      gla_callback_cmod* callback, void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(callback);

    p->cmod_callback = callback;
    p->cmod_callback_userdata = userdata;
error:
    return status;
}

int
gla_set_fmod_callback(gla_plan* p,
                      gla_callback_fmod* callback, void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(callback);

    p->fmod_callback = callback;
    p->fmod_callback_userdata = userdata;
error:
    return status;
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
    return 0;
}
