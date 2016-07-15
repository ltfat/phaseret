#include "ltfat.h"
#include "ltfat/macros.h"
#include "phaseret/dgtrealwrapper.h"
#include "dgtrealwrapper_private.h"

int
ltfat_idgtreal_long_execute_d_wrapper(void* plan, const complex double* c,
                                      int L,
                                      int W, double* f)
{
    return ltfat_idgtreal_long_execute_newarray_d(
               (ltfat_idgtreal_long_plan_d*) plan, c, f);
}

int
ltfat_dgtreal_long_execute_d_wrapper(void* plan, const double* f, int L, int W,
                                     complex double* c)
{
    return ltfat_dgtreal_long_execute_newarray_d(
               (ltfat_dgtreal_long_plan_d*) plan, f, c);
}

int
ltfat_idgtreal_fb_execute_d_wrapper(void* plan, const complex double* c, int L,
                                    int W, double* f)
{
    return ltfat_idgtreal_fb_execute_d(
               (ltfat_idgtreal_fb_plan_d*) plan, c, L, W, f);
}

int
ltfat_dgtreal_fb_execute_d_wrapper(void* plan, const double* f, int L, int W,
                                   complex double* c)
{
    return ltfat_dgtreal_fb_execute_d(
               (ltfat_dgtreal_fb_plan_d*) plan, f, L, W, c);
}

int
ltfat_idgtreal_long_done_d_wrapper(void** plan)
{
    return ltfat_idgtreal_long_done_d( (ltfat_idgtreal_long_plan_d**) plan);
}

int
ltfat_dgtreal_long_done_d_wrapper(void** plan)
{
    return ltfat_dgtreal_long_done_d( (ltfat_dgtreal_long_plan_d**) plan);
}

int
ltfat_idgtreal_fb_done_d_wrapper(void** plan)
{
    return ltfat_idgtreal_fb_done_d((ltfat_idgtreal_fb_plan_d**) plan);
}

int
ltfat_dgtreal_fb_done_d_wrapper(void** plan)
{
    return ltfat_dgtreal_fb_done_d((ltfat_dgtreal_fb_plan_d**) plan);
}

int
dgtreal_anasyn_execute_proj(dgtreal_anasyn_plan* p, const complex double cin[],
                            complex double cout[])
{
    p->backtra(p->backtra_userdata, cin, p->L, p->W, p->f);
    p->fwdtra(p->fwdtra_userdata, p->f, p->L, p->W,  cout);
}

int
dgtreal_anasyn_execute_syn(dgtreal_anasyn_plan* p, const complex double c[],
                           double f[])
{
    return p->backtra(p->backtra_userdata, c, p->L, p->W, f);
}

int
dgtreal_anasyn_execute_ana(dgtreal_anasyn_plan* p, const double f[],
                           complex double c[])
{
    return p->fwdtra(p->fwdtra_userdata, f, p->L, p->W,  c);
}

int
dgtreal_anasyn_done(dgtreal_anasyn_plan** p)
{
    dgtreal_anasyn_plan* pp = *p;
    pp->fwddonefunc(&pp->fwdtra_userdata);
    pp->backdonefunc(&pp->backtra_userdata);

    free(pp->f);
    free(pp);
    pp = NULL;
}

int
dgtreal_anasyn_init(const double g[], int gl, int L, int W, int a, int M,
                    complex double c[], dgtreal_anasyn_hint hint,
                    unsigned flags, dgtreal_anasyn_plan** pout)
{
    dgtreal_anasyn_plan* p = NULL;

    int N = L / a;
    int M2 = M / 2 + 1;

    int ispainless = gl <= M;

    double* g2 = NULL;
    int g2l = 0;

    if (ispainless)
    {
        // The length of the dual window is guaranteed to be gl
        g2l = gl;
        g2 = malloc(gl * sizeof * g2);
        ltfat_gabdual_painless_d(g, gl, a, M, g2);
    }
    else
    {
#ifndef NOBLASLAPACK
        g2l = L;
        g2 = malloc(L * sizeof * g2);
        ltfat_fir2long_d(g, gl, L, g2);
        ltfat_gabdual_long_d(g, L, 1, a, M, g2);
#else
        // COMPLAIN !!!
#endif
    }

    p = calloc(1, sizeof * p);
    *p = (dgtreal_anasyn_plan){.M = M, .a = a, .L = L, .W = W, .c = c};

    // Needed for some of the plans
    p->f = malloc(L * W * sizeof * p->f);

    if (dgtreal_anasyn_long == hint)
    {
        // Use _long functions only
        if (ispainless)
        {
            // Make the dual window longer if it is not already
            g2 = realloc(g2, L * sizeof * g2);
            ltfat_fir2long_d(g2, g2l, L, g2);
        }

        ltfat_idgtreal_long_init_d(c, g2, L, W, a, M, p->f, TIMEINV, flags,
                                   (ltfat_idgtreal_long_plan_d**)&p->backtra_userdata);

        p->backtra = &ltfat_idgtreal_long_execute_d_wrapper;
        p->backdonefunc = &ltfat_idgtreal_long_done_d_wrapper;

        // Ensure the original window is long enough
        ltfat_fir2long_d(g, gl, L, g2);
        ltfat_dgtreal_long_init_d(p->f, g2, L, W, a, M, c, TIMEINV, flags,
                                  (ltfat_dgtreal_long_plan_d**)&p->fwdtra_userdata);

        p->fwdtra = &ltfat_dgtreal_long_execute_d_wrapper;
        p->backdonefunc = &ltfat_dgtreal_long_done_d_wrapper;
    }
    else if ( dgtreal_anasyn_fb == hint )
    {
        // Use _fb functions only
        ltfat_idgtreal_fb_init_d( g2, g2l, a, M, TIMEINV, flags,
                                  (ltfat_idgtreal_fb_plan_d**)&p->backtra_userdata);

        p->backtra = &ltfat_idgtreal_fb_execute_d_wrapper;
        p->backdonefunc = &ltfat_idgtreal_fb_done_d_wrapper;

        ltfat_dgtreal_fb_init_d(g, gl, a, M, TIMEINV, flags,
                                (ltfat_dgtreal_fb_plan_d**)&p->fwdtra_userdata);

        p->fwdtra = &ltfat_dgtreal_fb_execute_d_wrapper;
        p->backdonefunc = &ltfat_dgtreal_fb_done_d_wrapper;
    }
    else if ( dgtreal_anasyn_auto == hint )
    {
        // Decide whether to use _fb or _long depending on the window lengths
        if (g2l < L)
        {
            ltfat_idgtreal_fb_init_d( g2, g2l, a, M, TIMEINV, flags,
                                      (ltfat_idgtreal_fb_plan_d**)&p->backtra_userdata);

            p->backtra = &ltfat_idgtreal_fb_execute_d_wrapper;
            p->backdonefunc = &ltfat_idgtreal_fb_done_d_wrapper;
        }
        else
        {
            ltfat_idgtreal_long_init_d(c, g2, L, W, a, M, p->f, TIMEINV, flags,
                                       (ltfat_idgtreal_long_plan_d**)&p->backtra_userdata);
            p->backtra = &ltfat_idgtreal_long_execute_d_wrapper;
            p->backdonefunc = &ltfat_idgtreal_long_done_d_wrapper;
        }

        if (gl < L)
        {
            ltfat_dgtreal_fb_init_d(g, gl, a, M, TIMEINV, flags,
                                    (ltfat_dgtreal_fb_plan_d**)&p->fwdtra_userdata);

            p->fwdtra = &ltfat_dgtreal_fb_execute_d_wrapper;
            p->fwddonefunc = &ltfat_dgtreal_fb_done_d_wrapper;
        }
        else
        {
            ltfat_dgtreal_long_init_d(p->f, g, L, W, a, M, c, TIMEINV, flags,
                                      (ltfat_dgtreal_long_plan_d**)&p->fwdtra_userdata);

            p->fwdtra = &ltfat_dgtreal_long_execute_d_wrapper;
            p->fwddonefunc = &ltfat_dgtreal_long_done_d_wrapper;
        }
    }
    else
    {
        // FAIL
    }

    free(g2);
    *pout = p;
}
