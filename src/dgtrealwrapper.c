#include "ltfat.h"
#include "ltfat/macros.h"
#include "phaseret/dgtrealwrapper.h"
#include "dgtrealwrapper_private.h"

int
ltfat_idgtreal_long_execute_d_wrapper(void* plan, const complex double* c,
                                      int UNUSED(L), int UNUSED(W), double* f)
{
    return ltfat_idgtreal_long_execute_newarray_d(
               (ltfat_idgtreal_long_plan_d*) plan, c, f);
}

int
ltfat_dgtreal_long_execute_d_wrapper(void* plan, const double* f,
                                     int UNUSED(L), int UNUSED(W),
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
    int status = LTFATERR_SUCCESS;
    CHECKSTATUS( p->backtra(p->backtra_userdata, cin, p->L, p->W, p->f),
                 "Back transform failed");
    CHECKSTATUS( p->fwdtra(p->fwdtra_userdata, p->f, p->L, p->W,  cout),
                 "Forward transform failed");
error:
    return status;
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
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    dgtreal_anasyn_plan* pp = *p;
    CHECKSTATUS( pp->fwddonefunc(&pp->fwdtra_userdata),
                 "Forward transform done function failed");
    CHECKSTATUS( pp->backdonefunc(&pp->backtra_userdata),
                 "Back trnasform done function failed");

    free(pp->f);
    free(pp);
    pp = NULL;
error:
    return status;
}

int
dgtreal_anasyn_init(const double g[], int gl, int L, int W, int a, int M,
                    complex double c[], dgtreal_anasyn_hint hint,
                    ltfat_phaseconvention ptype, unsigned flags,
                    dgtreal_anasyn_plan** pout)
{
    dgtreal_anasyn_plan* p = NULL;

    int status = LTFATERR_SUCCESS;

    int ispainless = gl <= M;

    double* g2 = NULL;
    int g2l = 0;

    if (ispainless)
    {
        // The length of the dual window is guaranteed to be gl
        g2l = gl;
        CHECKMEM( g2 = malloc(gl * sizeof * g2));
        CHECKSTATUS( ltfat_gabdual_painless_d(g, gl, a, M, g2),
                     "Gabdual painless call failed");
    }
    else
    {
#ifndef NOBLASLAPACK
        g2l = L;
        CHECKMEM( g2 = malloc(L * sizeof * g2));
        ltfat_fir2long_d(g, gl, L, g2);
        CHECKSTATUS( ltfat_gabdual_long_d(g, L, 1, a, M, g2),
                     "Gabdual long failed");
#else
        CHECK( LTFATERR_BADARG, 0, "Non-painless support was not compiled.");
#endif
    }

    CHECKMEM( p = calloc(1, sizeof * p));
    *p = (dgtreal_anasyn_plan) {.M = M, .a = a, .L = L, .W = W, .c = c};

    // Needed for some of the plans
    CHECKMEM( p->f = malloc(L * W * sizeof * p->f));

    if (dgtreal_anasyn_long == hint)
    {
        // Use _long functions only
        if (ispainless)
        {
            // Make the dual window longer if it is not already
            CHECKMEM( g2 = realloc(g2, L * sizeof * g2));
            ltfat_fir2long_d(g2, g2l, L, g2);
        }

        p->backtra = &ltfat_idgtreal_long_execute_d_wrapper;
        p->backdonefunc = &ltfat_idgtreal_long_done_d_wrapper;

        CHECKSTATUS(
            ltfat_idgtreal_long_init_d(c, g2, L, W, a, M, p->f, ptype, flags,
                                       (ltfat_idgtreal_long_plan_d**)&p->backtra_userdata),
            "idgtreal long init failed"
        );

        p->fwdtra = &ltfat_dgtreal_long_execute_d_wrapper;
        p->fwddonefunc = &ltfat_dgtreal_long_done_d_wrapper;

        // Ensure the original window is long enough
        ltfat_fir2long_d(g, gl, L, g2);

        CHECKSTATUS(
            ltfat_dgtreal_long_init_d(p->f, g2, L, W, a, M, c, ptype, flags,
                                      (ltfat_dgtreal_long_plan_d**)&p->fwdtra_userdata),
            "dgtreal long init failed");

    }
    else if ( dgtreal_anasyn_fb == hint )
    {
        // Use _fb functions only
        p->backtra = &ltfat_idgtreal_fb_execute_d_wrapper;
        p->backdonefunc = &ltfat_idgtreal_fb_done_d_wrapper;

        CHECKSTATUS(
            ltfat_idgtreal_fb_init_d( g2, g2l, a, M, ptype, flags,
                                      (ltfat_idgtreal_fb_plan_d**)&p->backtra_userdata),
            "idgtreal fb init failed");

        p->fwdtra = &ltfat_dgtreal_fb_execute_d_wrapper;
        p->fwddonefunc = &ltfat_dgtreal_fb_done_d_wrapper;

        CHECKSTATUS(
            ltfat_dgtreal_fb_init_d(g, gl, a, M, ptype, flags,
                                    (ltfat_dgtreal_fb_plan_d**)&p->fwdtra_userdata),
            "dgtreal fb init failed");

    }
    else if ( dgtreal_anasyn_auto == hint )
    {
        // Decide whether to use _fb or _long depending on the window lengths
        if (g2l < L)
        {
            p->backtra = &ltfat_idgtreal_fb_execute_d_wrapper;
            p->backdonefunc = &ltfat_idgtreal_fb_done_d_wrapper;

            ltfat_idgtreal_fb_init_d( g2, g2l, a, M, ptype, flags,
                                      (ltfat_idgtreal_fb_plan_d**)&p->backtra_userdata);

        }
        else
        {
            p->backtra = &ltfat_idgtreal_long_execute_d_wrapper;
            p->backdonefunc = &ltfat_idgtreal_long_done_d_wrapper;

            ltfat_idgtreal_long_init_d(c, g2, L, W, a, M, p->f, ptype, flags,
                                       (ltfat_idgtreal_long_plan_d**)&p->backtra_userdata);
        }

        if (gl < L)
        {
            p->fwdtra = &ltfat_dgtreal_fb_execute_d_wrapper;
            p->fwddonefunc = &ltfat_dgtreal_fb_done_d_wrapper;

            ltfat_dgtreal_fb_init_d(g, gl, a, M, ptype, flags,
                                    (ltfat_dgtreal_fb_plan_d**)&p->fwdtra_userdata);

        }
        else
        {
            p->fwdtra = &ltfat_dgtreal_long_execute_d_wrapper;
            p->fwddonefunc = &ltfat_dgtreal_long_done_d_wrapper;

            ltfat_dgtreal_long_init_d(p->f, g, L, W, a, M, c, ptype, flags,
                                      (ltfat_dgtreal_long_plan_d**)&p->fwdtra_userdata);

        }
    }
    else
    {
        CHECKCANTHAPPEN("No such dgtreal_anasyn hint");
    }

    free(g2);
    *pout = p;

    return status;
error:
    if (g2) free(g2);
    if (p)
    {
        if (p->backtra_userdata) p->backdonefunc(&p->backtra_userdata);
        if (p->fwdtra_userdata) p->fwddonefunc(&p->fwdtra_userdata);
        if (p->f) free(p->f);
        free(p);
    }
    return status;
}
