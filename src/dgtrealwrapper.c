#include "phaseret/dgtrealwrapper.h"
#include "ltfat/macros.h"
#include "dgtrealwrapper_private.h"

int
PHASERET_NAME(ltfat_idgtreal_long_execute_wrapper)(void* plan,
        const LTFAT_COMPLEX* c, ltfat_int UNUSED(L), ltfat_int UNUSED(W), LTFAT_REAL* f)
{
    return LTFAT_NAME(idgtreal_long_execute_newarray)(
               (LTFAT_NAME(idgtreal_long_plan)*) plan, c, f);
}

int
PHASERET_NAME(ltfat_dgtreal_long_execute_wrapper)(void* plan,
        const LTFAT_REAL* f, ltfat_int UNUSED(L), ltfat_int UNUSED(W), LTFAT_COMPLEX* c)
{
    return LTFAT_NAME(dgtreal_long_execute_newarray)(
               (LTFAT_NAME(dgtreal_long_plan)*) plan, f, c);
}

int
PHASERET_NAME(ltfat_idgtreal_fb_execute_wrapper)(void* plan,
        const LTFAT_COMPLEX* c, ltfat_int L, ltfat_int W, LTFAT_REAL* f)
{
    return LTFAT_NAME(idgtreal_fb_execute)(
               (LTFAT_NAME(idgtreal_fb_plan)*) plan, c, L, W, f);
}

int
PHASERET_NAME(ltfat_dgtreal_fb_execute_wrapper)(void* plan,
        const LTFAT_REAL* f, ltfat_int L, ltfat_int W, LTFAT_COMPLEX* c)
{
    return LTFAT_NAME(dgtreal_fb_execute)(
               (LTFAT_NAME(dgtreal_fb_plan)*) plan, f, L, W, c);
}

int
PHASERET_NAME(ltfat_idgtreal_long_done_wrapper)(void** plan)
{
    return LTFAT_NAME(idgtreal_long_done)( (LTFAT_NAME(idgtreal_long_plan)**) plan);
}

int
PHASERET_NAME(ltfat_dgtreal_long_done_wrapper)(void** plan)
{
    return LTFAT_NAME(dgtreal_long_done)( (LTFAT_NAME(dgtreal_long_plan)**) plan);
}

int
PHASERET_NAME(ltfat_idgtreal_fb_done_wrapper)(void** plan)
{
    return LTFAT_NAME(idgtreal_fb_done)((LTFAT_NAME(idgtreal_fb_plan)**) plan);
}

int
PHASERET_NAME(ltfat_dgtreal_fb_done_wrapper)(void** plan)
{
    return LTFAT_NAME(dgtreal_fb_done)((LTFAT_NAME(dgtreal_fb_plan)**) plan);
}

PHASERET_API int
PHASERET_NAME(dgtreal_execute_proj)(
    PHASERET_NAME(dgtreal_plan)* p, const LTFAT_COMPLEX cin[],
    LTFAT_COMPLEX cout[])
{
    int status = LTFATERR_SUCCESS;
    CHECKSTATUS( p->backtra(p->backtra_userdata, cin, p->L, p->W, p->f),
                 "Back transform failed");
    CHECKSTATUS( p->fwdtra(p->fwdtra_userdata, p->f, p->L, p->W,  cout),
                 "Forward transform failed");
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(dgtreal_execute_syn)(
    PHASERET_NAME(dgtreal_plan)* p, const LTFAT_COMPLEX c[], LTFAT_REAL f[])
{
    return p->backtra(p->backtra_userdata, c, p->L, p->W, f);
}

PHASERET_API int
PHASERET_NAME(dgtreal_execute_ana)(
    PHASERET_NAME(dgtreal_plan)* p, const LTFAT_REAL f[], LTFAT_COMPLEX c[])
{
    return p->fwdtra(p->fwdtra_userdata, f, p->L, p->W,  c);
}

PHASERET_API int
PHASERET_NAME(dgtreal_done)(PHASERET_NAME(dgtreal_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(dgtreal_plan)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    if (pp->fwdtra_userdata)
    CHECKSTATUS( pp->fwddonefunc(&pp->fwdtra_userdata),
                 "Forward transform done function failed");

    if (pp->backtra_userdata)
    CHECKSTATUS( pp->backdonefunc(&pp->backtra_userdata),
                 "Back trnasform done function failed");

    ltfat_safefree(pp->f);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}


PHASERET_API int
PHASERET_NAME(dgtreal_init)(const LTFAT_REAL g[], ltfat_int gl, ltfat_int L, ltfat_int W,
                            ltfat_int a, ltfat_int M, LTFAT_COMPLEX c[],
                            phaseret_dgtreal_init_params* params, PHASERET_NAME(dgtreal_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(dgtreal_plan)* p = NULL;
    phaseret_dgtreal_init_params paramsLoc;
    int ispainless = gl <= M;
    LTFAT_REAL* g2 = NULL;
    ltfat_int g2l = 0;

    ltfat_int minL = ltfat_lcm(a, M);

    if (params)
        paramsLoc = *params;
    else
        phaseret_dgtreal_init_params_defaults(&paramsLoc);

    CHECK(LTFATERR_BADTRALEN, !(L % minL),
          "L must divisible by lcm(a,M)=%d.", minL);

    if (ispainless)
    {
        // The length of the dual window is guaranteed to be gl
        g2l = gl;
        CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(gl));
        CHECKSTATUS( LTFAT_NAME(gabdual_painless)(g, gl, a, M, g2),
                     "Gabdual painless call failed");
    }
    else
    {
#ifndef NOBLASLAPACK
        g2l = L;
        CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(L));
        LTFAT_NAME(fir2long)(g, gl, L, g2);
        CHECKSTATUS( LTFAT_NAME(gabdual_long)(g, L, a, M, g2),
                     "Gabdual long failed");
#else
        CHECK( LTFATERR_NOTSUPPORTED, 0, "Non-painless support was not compiled.");
#endif
    }

    CHECKMEM( p = (PHASERET_NAME(dgtreal_plan)*) ltfat_calloc(1, sizeof * p));
    p->M = M, p->a = a, p->L = L, p->W = W, p->c = c;

    // Needed for some of the plans
    CHECKMEM( p->f = LTFAT_NAME_REAL(malloc)(L * W));

    if (phaseret_dgtreal_long == paramsLoc.hint)
    {
        // Use _long functions only
        if (ispainless)
        {
            // Make the dual window longer if it is not already
            CHECKMEM( g2 = LTFAT_NAME_REAL(realloc)(g2, g2l, L));
            LTFAT_NAME(fir2long)(g2, g2l, L, g2);
        }

        p->backtra = &PHASERET_NAME(ltfat_idgtreal_long_execute_wrapper);
        p->backdonefunc = &PHASERET_NAME(ltfat_idgtreal_long_done_wrapper);

        CHECKSTATUS(
            LTFAT_NAME(idgtreal_long_init)(c, g2, L, W, a, M, p->f, paramsLoc.ptype,
                                           paramsLoc.fftw_flags,
                                           (LTFAT_NAME(idgtreal_long_plan)**)&p->backtra_userdata),
            "idgtreal long init failed"
        );

        p->fwdtra = &PHASERET_NAME(ltfat_dgtreal_long_execute_wrapper);
        p->fwddonefunc = &PHASERET_NAME(ltfat_dgtreal_long_done_wrapper);

        // Ensure the original window is long enough
        LTFAT_NAME(fir2long)(g, gl, L, g2);

        CHECKSTATUS(
            LTFAT_NAME(dgtreal_long_init)(p->f, g2, L, W, a, M, c, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(dgtreal_long_plan)**)&p->fwdtra_userdata),
            "dgtreal long init failed");

    }
    else if ( phaseret_dgtreal_fb == paramsLoc.hint )
    {
        // Use _fb functions only
        p->backtra = &PHASERET_NAME(ltfat_idgtreal_fb_execute_wrapper);
        p->backdonefunc = &PHASERET_NAME(ltfat_idgtreal_fb_done_wrapper);

        CHECKSTATUS(
            LTFAT_NAME(idgtreal_fb_init)( g2, g2l, a, M, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(idgtreal_fb_plan)**)&p->backtra_userdata),
            "idgtreal fb init failed");

        p->fwdtra = &PHASERET_NAME(ltfat_dgtreal_fb_execute_wrapper);
        p->fwddonefunc = &PHASERET_NAME(ltfat_dgtreal_fb_done_wrapper);

        CHECKSTATUS(
            LTFAT_NAME(dgtreal_fb_init)(g, gl, a, M, paramsLoc.ptype, paramsLoc.fftw_flags,
                                        (LTFAT_NAME(dgtreal_fb_plan)**)&p->fwdtra_userdata),
            "dgtreal fb init failed");

    }
    else if ( phaseret_dgtreal_auto == paramsLoc.hint )
    {
        // Decide whether to use _fb or _long depending on the window lengths
        if (g2l < L)
        {
            p->backtra = &PHASERET_NAME(ltfat_idgtreal_fb_execute_wrapper);
            p->backdonefunc = &PHASERET_NAME(ltfat_idgtreal_fb_done_wrapper);

            LTFAT_NAME(idgtreal_fb_init)( g2, g2l, a, M, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(idgtreal_fb_plan)**)&p->backtra_userdata);

        }
        else
        {
            p->backtra = &PHASERET_NAME(ltfat_idgtreal_long_execute_wrapper);
            p->backdonefunc = &PHASERET_NAME(ltfat_idgtreal_long_done_wrapper);

            LTFAT_NAME(idgtreal_long_init)(c, g2, L, W, a, M, p->f, paramsLoc.ptype,
                                           paramsLoc.fftw_flags,
                                           (LTFAT_NAME(idgtreal_long_plan)**)&p->backtra_userdata);
        }

        if (gl < L)
        {
            p->fwdtra = &PHASERET_NAME(ltfat_dgtreal_fb_execute_wrapper);
            p->fwddonefunc = &PHASERET_NAME(ltfat_dgtreal_fb_done_wrapper);

            LTFAT_NAME(dgtreal_fb_init)(g, gl, a, M, paramsLoc.ptype, paramsLoc.fftw_flags,
                                        (LTFAT_NAME(dgtreal_fb_plan)**)&p->fwdtra_userdata);

        }
        else
        {
            p->fwdtra = &PHASERET_NAME(ltfat_dgtreal_long_execute_wrapper);
            p->fwddonefunc = &PHASERET_NAME(ltfat_dgtreal_long_done_wrapper);

            LTFAT_NAME(dgtreal_long_init)(p->f, g, L, W, a, M, c, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(dgtreal_long_plan)**)&p->fwdtra_userdata);

        }
    }
    else
    {
        CHECKCANTHAPPEN("No such dgtreal hint");
    }

    ltfat_free(g2);
    *pout = p;

    return status;
error:
    ltfat_safefree(g2);
    if (p) PHASERET_NAME(dgtreal_done)(&p);
    return status;
}
