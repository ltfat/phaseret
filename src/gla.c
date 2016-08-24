#include "phaseret/gla.h"
#include "phaseret/utils.h"
#include "ltfat/macros.h"
#include "dgtrealwrapper_private.h"

struct PHASERET_NAME(gla_plan)
{
    PHASERET_NAME(dgtreal_plan)* p;
    PHASERET_NAME(gla_callback_status)* status_callback;
    void* status_callback_userdata;
    PHASERET_NAME(gla_callback_cmod)* cmod_callback;
    void* cmod_callback_userdata;
    PHASERET_NAME(gla_callback_fmod)* fmod_callback;
    void* fmod_callback_userdata;
// For storing magnitude
    LTFAT_REAL* s;
// Storing cinit
    const LTFAT_COMPLEX* cinit;
// Used just for fgla
    int do_fast;
    double alpha;
    LTFAT_COMPLEX* t;
};

PHASERET_API int
PHASERET_NAME(gla)(const LTFAT_COMPLEX cinit[], const LTFAT_REAL g[],
                   const int L,
                   const int gl, const int W, const int a, const int M, const int iter,
                   LTFAT_COMPLEX cout[])
{
    PHASERET_NAME(gla_plan)* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        PHASERET_NAME(gla_init)(cinit, g, L, gl, W, a, M, 0.99, cout, NULL, &p),
        "Init failed");

    CHECKSTATUS( PHASERET_NAME(gla_execute)(p, iter), "Execute failed");

error:
    if (p) PHASERET_NAME(gla_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(gla_init)(const LTFAT_COMPLEX cinit[], const LTFAT_REAL g[],
                        const int L, const int gl, const int W, const int a,
                        const int M, const double alpha, LTFAT_COMPLEX c[],
                        phaseret_dgtreal_init_params* params,
                        PHASERET_NAME(gla_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(gla_plan)* p = NULL;
    int N = L / a;
    int M2 = M / 2 + 1;

    CHECK(LTFATERR_BADARG, alpha >= 0.0, "alpha cannot be negative");
    CHECKMEM( p = (PHASERET_NAME(gla_plan)*) ltfat_calloc(1, sizeof * p));
    CHECKMEM( p->s = LTFAT_NAME_REAL(malloc)(M2 * N * W));

    if (alpha > 0.0)
    {
        p->do_fast = 1;
        p->alpha = alpha;
        CHECKMEM( p->t = LTFAT_NAME_COMPLEX(malloc)(M2 * N * W));
    }

    CHECKSTATUS(
        PHASERET_NAME(dgtreal_init)(g, gl, L, W, a, M, c, params, &p->p),
        "dgtrealwrapper init failed");

    p->cinit = cinit;

    *pout = p;
    return status;
error:
    if (p) PHASERET_NAME(gla_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(gla_done)(PHASERET_NAME(gla_plan)** p)
{
    PHASERET_NAME(gla_plan)* pp = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    if (pp->p)
        CHECKSTATUS(
            PHASERET_NAME(dgtreal_done)(&pp->p),
            "dgtreal wrapper done failed");

    ltfat_safefree(pp->t);
    ltfat_safefree(pp->s);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(gla_execute_newarray)(PHASERET_NAME(gla_plan)* p,
                                    const LTFAT_COMPLEX cinit[], const int iter,
                                    LTFAT_COMPLEX cout[])
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(dgtreal_plan) pp2;
    PHASERET_NAME(gla_plan) p2;
    CHECKNULL(p); CHECKNULL(cinit); CHECKNULL(cout);
    // Shallow copy the plan and replace c
    p2 = *p;
    pp2 = *p2.p;
    pp2.c = cout;
    p2.p = &pp2;
    p2.cinit = cinit;
    return PHASERET_NAME(gla_execute)(&p2, iter);
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(gla_execute)(PHASERET_NAME(gla_plan)* p, const int iter)
{
    int status = LTFATERR_SUCCESS;
    int M, L, W, a, M2, N;
    PHASERET_NAME(dgtreal_plan)* pp;
    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, iter > 0,
          "At least one iteration is requred. Passed %d.", iter);

    pp = p->p;
    CHECKNULL(pp->c); CHECKNULL(p->cinit);
    M = pp->M;
    L = pp->L;
    W = pp->W;
    a = pp->a;
    M2 = M / 2 + 1;
    N = L / a;

    // Store the magnitude
    for (int ii = 0; ii < N * M2 * W; ii++)
        p->s[ii] = ltfat_abs(p->cinit[ii]);

    // Copy to the output array if we are not working inplace
    if (p->cinit != pp->c)
        memcpy(pp->c, p->cinit, (N * M2 * W) * sizeof * pp->c);

    // Inicialize the "acceleration" array
    if (p->do_fast)
        memcpy(p->t, pp->c, (N * M2 * W) * sizeof * p->t );

    for (int ii = 0; ii < iter; ii++)
    {
        // Perform idgtreal
        CHECKSTATUS( PHASERET_NAME(dgtreal_execute_syn)(pp, pp->c, pp->f),
                     "idgtreal failed");

        // Optional signal modification
        if (p->fmod_callback)
            CHECKSTATUS(
                p->fmod_callback(p->fmod_callback_userdata, pp->f, L, W, a, M),
                "fmod callback failed");

        // Perform dgtreal
        CHECKSTATUS( PHASERET_NAME(dgtreal_execute_ana)(pp, pp->f, pp->c),
                     "dgtreal failed");

        PHASERET_NAME(force_magnitude)(pp->c, p->s, N * M2 * W, pp->c);

        // The acceleration step
        if (p->do_fast)
            PHASERET_NAME(fastupdate)(pp->c, p->t, p->alpha, N * M2 * W );

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
                CHECKMEM( p->t = LTFAT_NAME_COMPLEX(malloc)(M2 * N * W));
                memcpy(p->t, pp->c, (N * M2 * W) * sizeof * p->t );
            }
        }
    }

error:
    return status;
}


PHASERET_API int
PHASERET_NAME(gla_set_status_callback)(PHASERET_NAME(gla_plan)* p,
                                       PHASERET_NAME(gla_callback_status)* callback,
                                       void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(callback);

    p->status_callback = callback;
    p->status_callback_userdata = userdata;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(gla_set_cmod_callback)(PHASERET_NAME(gla_plan)* p,
                                     PHASERET_NAME(gla_callback_cmod)* callback,
                                     void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(callback);

    p->cmod_callback = callback;
    p->cmod_callback_userdata = userdata;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(gla_set_fmod_callback)(PHASERET_NAME(gla_plan)* p,
                                     PHASERET_NAME(gla_callback_fmod)* callback,
                                     void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(callback);

    p->fmod_callback = callback;
    p->fmod_callback_userdata = userdata;
error:
    return status;
}

int
PHASERET_NAME(fastupdate)(LTFAT_COMPLEX* c, LTFAT_COMPLEX* t, double alpha,
                          int L)
{
    for (int ii = 0; ii < L; ii++)
    {
        LTFAT_COMPLEX cold = c[ii];
        c[ii] = c[ii] + ((LTFAT_REAL)alpha) * (c[ii] - t[ii]);
        t[ii] = cold;
    }
    return 0;
}
