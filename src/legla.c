#include "phaseret/legla.h"
#include "phaseret/gla.h"
#include "phaseret/dgtrealwrapper.h"
#include "dgtrealwrapper_private.h"
#include "ltfat/macros.h"

struct legla_plan
{
    leglaupdate_plan* updateplan;
    dgtreal_anasyn_plan* dgtplan;
    legla_callback_status* status_callback;
    void* status_callback_userdata;
    legla_callback_cmod* cmod_callback;
    void* cmod_callback_userdata;
// Storing magnitude
    LTFAT_REAL* s;
    const LTFAT_COMPLEX* cinit;
// Used just for flegla
    int do_fast;
    double alpha;
    LTFAT_COMPLEX* t;
};

struct leglaupdate_plan
{
    int kNo;
    LTFAT_COMPLEX** k;
    LTFAT_COMPLEX* buf;
    int a;
    int N;
    int W;
    leglaupdate_plan_col* plan_col;
};

struct leglaupdate_plan_col
{
    phaseret_size ksize;
    phaseret_size ksize2;
    int M;
    int flags;
};

#define legla_init_params_hash 968456

int
legla_init_params_defaults(legla_init_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);

    params->relthr = 1e-3;
    params->ksize.width = 0;
    params->ksize.height = 0;
    params->hint = dgtreal_anasyn_auto;
    params->dgtrealflags = FFTW_ESTIMATE;
    params->leglaflags = MOD_COEFFICIENTWISE | MOD_MODIFIEDUPDATE;
    params->private_hash_do_not_use = legla_init_params_hash;
error:
    return status;
}

int
legla(const LTFAT_COMPLEX cinit[], const LTFAT_REAL g[],  const int L,
      const int gl, const int W, const int a, const int M, const int iter,
      LTFAT_COMPLEX cout[])
{
    legla_plan* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        legla_init(cinit, g, gl, L, W, a, M, 0.99, cout, NULL, &p),
        "legla init failed");

    CHECKSTATUS( legla_execute(p, iter), "legla execute failed");

error:
    if (p) legla_done(&p);
    return status;
}


int
legla_init(const LTFAT_COMPLEX cinit[], const LTFAT_REAL g[], const int L,
           const int gl, const int W, const int a, const int M,
           const double alpha, LTFAT_COMPLEX c[],
           legla_init_params* params, legla_plan** pout)
{
    legla_plan* p = NULL;
    legla_init_params pLoc;
    phaseret_size ksize;
    LTFAT_COMPLEX* kernsmall = NULL;
    int M2 = M / 2 + 1;
    int N = L / a;

    int status = LTFATERR_SUCCESS;

    CHECKNULL(c); CHECKNULL(pout);
    CHECK(LTFATERR_BADARG, alpha >= 0, "alpha must be bigger or equal to zero.");

    if (params)
    {
        CHECK(LTFATERR_CANNOTHAPPEN,
              params->private_hash_do_not_use == legla_init_params_hash,
              "params were not initialized with legla_init_params_defaults");
        CHECK(LTFATERR_NOTINRANGE, params->relthr >= 0 && params->relthr <= 1,
              "relthr must be in range [0-1]");

        memcpy(&pLoc, params, sizeof * params);
    }
    else
        legla_init_params_defaults(&pLoc);

    ksize = pLoc.ksize;

    CHECK(LTFATERR_BADSIZE, (ksize.width <= N && ksize.height <= M) ||
          (ksize.width >= 0 && ksize.height >= 0),
          "Bad kernel size {.width=%d,.height=%d}.", ksize.width, ksize.height);

    if (pLoc.relthr == 0.0)
    {
        // Use kernel size directly
        CHECK(LTFATERR_BADSIZE, ksize.width > 0 && ksize.height > 0,
              "Bad kernel size {.width=%d,.height=%d}.", ksize.width, ksize.height);
    }
    else
    {
        CHECK(LTFATERR_BADSIZE, !(ksize.width == 0 && ksize.height > 0) ||
              !(ksize.width > 0 && ksize.height == 0),
              "Kernel size is zero in one direction {.width=%d,.height=%d}.",
              ksize.width, ksize.height);

        if (ksize.width == 0 && ksize.height == 0)
        {
            ksize.width = 2 * ceil(((double)M) / a) - 1;
            ksize.height = 2 * ceil(((double)M) / a) - 1;
        }
    }

    CHECKMEM( p = (legla_plan*) calloc(1, sizeof * p));
    CHECKMEM( p->s = LTFAT_NAME_REAL(malloc)(M2 * N * W));

    CHECKSTATUS(
        dgtreal_anasyn_init(g, gl, L, W, a, M, c, pLoc.hint, LTFAT_FREQINV,
                            pLoc.dgtrealflags, &p->dgtplan), "dgtreal anasyn init failed");

    // Get the "impulse response" and crop it
    memset(c, 0, M2 * N * W * sizeof * c);
    c[0] = 1.0;
    dgtreal_anasyn_execute_proj(p->dgtplan, c, c);
    phaseret_size bigsize;
    bigsize.width = N; bigsize.height = M;

    if ( pLoc.relthr != 0.0)
        legla_findkernelsize(c, bigsize, pLoc.relthr, &ksize);

    DEBUG("Kernel size: {.width=%d,.height=%d}", ksize.width, ksize.height);

    CHECKMEM( kernsmall =
                  LTFAT_NAME_COMPLEX(malloc)(ksize.width * (ksize.height / 2 + 1)));

    legla_big2small_kernel(c, bigsize, ksize, kernsmall);

    CHECKSTATUS(
        leglaupdate_init( kernsmall, ksize, L, W, a, M, pLoc.leglaflags,
                          &p->updateplan),
        "leglaupdate init failed");

    if (alpha > 0.0)
    {
        p->do_fast = 1;
        p->alpha = alpha;
        CHECKMEM( p->t = LTFAT_NAME_COMPLEX(malloc)(M2 * N * W));
    }

    free(kernsmall);
    p->cinit = cinit;
    *pout = p;
    return status;
error:
    if (kernsmall) free(kernsmall);
    if (p)
    {
        if (p->t) free(p->t);
        if (p->s) free(p->s);
        free(p);
        *pout = NULL;
    }
    return status;
}

int
legla_done(legla_plan** p)
{
    int status = LTFATERR_SUCCESS;
    legla_plan* pp; 
    CHECKMEM(p); CHECKMEM(*p);
    pp = *p;

    dgtreal_anasyn_done(&pp->dgtplan);
    leglaupdate_done(&pp->updateplan);
    free(pp->s);
    if (pp->t) free(pp->t);
    free(pp);
    pp = NULL;
error:
    return status;
}

int
legla_execute(legla_plan* p, const int iter)
{
    int status = LTFATERR_SUCCESS;
    int M,L,W,a,M2,N;
    dgtreal_anasyn_plan* pp; 
    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, iter > 0, "At least one iteration is required");
    CHECKNULL(p->cinit);
    pp = p->dgtplan;
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

    if (p->do_fast)
        memcpy(p->t, pp->c, (N * M2 * W) * sizeof * p->t );

    for (int ii = 0; ii < iter; ii++)
    {
        leglaupdate_execute(p->updateplan, p->s, pp->c, pp->c);

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
            int status2 = p->status_callback(pp, p->status_callback_userdata,
                                             pp->c, L, W, a, M, &p->alpha, ii);
            if (status2 > 0)
                break;
            else
                CHECKSTATUS(status2, "status callback failed");

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

int
legla_big2small_kernel(LTFAT_COMPLEX* bigc, phaseret_size bigsize,
                       phaseret_size ksize, LTFAT_COMPLEX* smallc)
{
    div_t wmod = div(ksize.width, 2);
    div_t hmod = div(ksize.height, 2);
    int kernh2 = hmod.quot + 1;

    int M2 = bigsize.height / 2 + 1;

    for (int ii = 0; ii < wmod.quot + wmod.rem; ii++)
    {
        LTFAT_COMPLEX* smallcCol = smallc + kernh2 * ii;
        LTFAT_COMPLEX* bigcCol = bigc + M2 * ii;

        memcpy(smallcCol, bigcCol, kernh2 * sizeof * smallcCol);
    }

    for (int ii = 1; ii < wmod.quot + 1; ii++)
    {
        LTFAT_COMPLEX* smallcCol = smallc + kernh2 * ( ksize.width - ii);
        LTFAT_COMPLEX* bigcCol = bigc + M2 * ( bigsize.width - ii);

        memcpy(smallcCol, bigcCol, kernh2 * sizeof * smallcCol);
    }
    return LTFATERR_SUCCESS;
}

int
legla_findkernelsize(LTFAT_COMPLEX* bigc, phaseret_size bigsize,
                     double relthr, phaseret_size* ksize)

{
    double thr = relthr * ltfat_abs(bigc[0]);
    int realHeight = bigsize.height / 2 + 1;
    div_t wmod = div(ksize->width, 2);
    div_t hmod = div(ksize->height, 2);

    int lastrow = 0;
    for (int n = 0; n < wmod.quot + wmod.rem - 1; n++)
        for (int m = 1; m < hmod.quot + hmod.rem - 1; m++)
            if ( ltfat_abs(bigc[n * realHeight + m]) > thr && m > lastrow )
                lastrow = m;

    ksize->height = 2 * lastrow + 1;
    hmod = div(ksize->height, 2);

    // Kernel is always symmetric in the horizontal direction
    int lastcol = 0;
    for (int m = 0; m <  hmod.quot + hmod.rem - 1; m++)
        for (int n = 1; n < wmod.quot + wmod.rem - 1; n++)
            if ( ltfat_abs(bigc[n * realHeight + m]) > thr && n > lastcol )
                lastcol = n;

    ksize->width = 2 * lastcol + 1;

    return LTFATERR_SUCCESS;
}

int
legla_execute_newarray(legla_plan* p, const LTFAT_COMPLEX cinit[],
                       const int iter, LTFAT_COMPLEX c[])
{
    int status = LTFATERR_SUCCESS;
    legla_plan p2; 
    dgtreal_anasyn_plan pp2; 
    CHECKNULL(p); CHECKNULL(cinit); CHECKNULL(c);

    // Shallow copy the plan and replace c
    p2 = *p;
    pp2 = *p2.dgtplan;
    pp2.c = c;
    p2.dgtplan = &pp2;
    p2.cinit = cinit;
    return legla_execute(&p2, iter);
error:
    return status;
}


int
leglaupdate_init_col( int M, phaseret_size ksize, int flags,
                      leglaupdate_plan_col** pout)
{
    leglaupdate_plan_col* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKMEM( p = (leglaupdate_plan_col*) calloc(1, sizeof * p));
    phaseret_size ksize2;
    ksize2.width = ksize.width / 2 + 1;
    ksize2.height = ksize.height / 2 + 1;
    p->M = M; p->flags = flags; p->ksize = ksize; p->ksize2 = ksize2;

    // Sanitize flags (set defaults)
    if (p->flags & (MOD_FRAMEWISE | MOD_COEFFICIENTWISE))
    {
        p->flags &= ~MOD_STEPWISE; // For safety, clear the default flag.

        // Clear also lower priority flags
        if (p->flags & MOD_COEFFICIENTWISE)
            p->flags &= ~MOD_FRAMEWISE;
    }
    else
        p->flags |= MOD_STEPWISE; // Set the default flag if none of the others is set

    if (p->flags & (ORDER_REV))
        p->flags &= ~ORDER_FWD;
    else
        p->flags |= ORDER_FWD;

    if (p->flags & (EXT_UPDOWN))
        p->flags &= ~EXT_BOTH;
    else
        p->flags |= EXT_BOTH;

    *pout = p;

    return status;
error:
    if (p) free(p);
    return status;
}

int
leglaupdate_init(const LTFAT_COMPLEX kern[], phaseret_size ksize,
                 int L, int W, int a, int M, int flags, leglaupdate_plan** pout)
{
    int N = L / a;
    int M2 = M / 2 + 1;
    int kernh2;
    int status = LTFATERR_SUCCESS;

    leglaupdate_plan* p = NULL;
    LTFAT_COMPLEX* ktmp = NULL;

    CHECKMEM( p = (leglaupdate_plan*) calloc (1, sizeof * p));

    CHECKSTATUS( leglaupdate_init_col( M, ksize, flags, &p->plan_col),
                 "leglaupdate init failed");

    // N
    p->N = (flags & EXT_UPDOWN) ? N - (ksize.width - 1) : N;
    p->a = a;
    p->W = W;

    p->kNo = phaseret_lcm(M, a) / a;

    CHECKMEM( p->k = (LTFAT_COMPLEX**) malloc( p->kNo * sizeof * p->k));

    CHECKMEM( p->buf = 
            LTFAT_NAME_COMPLEX(malloc)( ((M2 + ksize.height - 1) * (p->N + ksize.width - 1))));

    kernh2 = ksize.height / 2 + 1;

    CHECKMEM( ktmp = LTFAT_NAME_COMPLEX(calloc)( ksize.width * kernh2));
    // Involute
    for (int ii = 0; ii < kernh2; ii++)
    {
        const LTFAT_COMPLEX* kRow = kern + ii;
        LTFAT_COMPLEX* kmodRow = ktmp + ii;

        kmodRow[0] = kRow[0];
        for (int jj = 1; jj < ksize.width; jj++)
            kmodRow[jj * kernh2] = conj(kRow[(ksize.width - jj) * kernh2]);
    }

    if (flags & MOD_MODIFIEDUPDATE) ktmp[0] = 0.0 + I * 0.0;

    for (int n = 0; n < p->kNo; n++)
    {
        CHECKMEM( p->k[n] = LTFAT_NAME_COMPLEX(calloc)( ksize.width * kernh2));
        kernphasefi(ktmp, ksize, n, a, M, p->k[n]);
    }

    free(ktmp);

    *pout = p;
    return status;
error:
    if (ktmp) free(ktmp);
    if (p)
    {
        if (p->buf) free(p->buf);
        if (p->k)
        {
            for (int ii = 0; ii < p->kNo; ii++) free(p->k[ii]);
            free(p->k);
        }
        free(p);
    }
    return status;
}

void
leglaupdate_done(leglaupdate_plan** plan)
{

    leglaupdate_plan* pp = *plan;
    for (int n = 0; n < pp->kNo; n++)
        free(pp->k[n]);

    if (pp->k) free(pp->k);

    free(pp->buf);
    free(pp);
    pp = NULL;
}

void
kernphasefi(const LTFAT_COMPLEX kern[], phaseret_size ksize,
            int n, int a, int M, LTFAT_COMPLEX kernmod[])
{
    /* int kernh2 = ksize.height / 2 + 1; */
    /* int kernw2 = ksize.width / 2; */
    div_t wmod = div(ksize.width, 2);
    div_t hmod = div(ksize.height, 2);
    int kernh2 = hmod.quot + 1;
    /*Modulate */
    /* double idx = - ( ksize.height - kernh2); */

    for (int ii = 0; ii < hmod.quot + hmod.rem; ii++)
    {
        const LTFAT_COMPLEX* kRow = kern + ii;
        LTFAT_COMPLEX* kmodRow = kernmod + ii;
        double arg = -2.0 * M_PI * n * a / M * (ii);

        // fftshift
        for (int jj = 0; jj < wmod.quot + wmod.rem - 1; jj++)
            kmodRow[(jj + wmod.quot)*kernh2] = exp(I * arg) * kRow[jj * kernh2];

        for (int jj = 0; jj < wmod.quot; jj++)
            kmodRow[jj * kernh2] = exp(I * arg) * kRow[(jj + wmod.quot + wmod.rem) *
                                   kernh2];
    }
}

int phaseret_lcm(int m, int n)
{
    return m / phaseret_gcd(m, n) * n;
}

int phaseret_gcd(int m, int n)
{
    int tmp;

    while (m)
    {
        tmp = m;
        m = n % m;
        n = tmp;
    }

    return n;
}

extern void
leglaupdate_execute(leglaupdate_plan* plan, const LTFAT_REAL s[],
                    LTFAT_COMPLEX c[], LTFAT_COMPLEX cout[])
{
    leglaupdate_plan_col* p = plan->plan_col;
    int M2 = p->M / 2 + 1;
    int N = plan->N;
    int W = plan->W;
    //int kernh2 = p->ksize2.height;
    //int kernw = p->ksize.width;
    int M2buf = M2 + p->ksize.height - 1;
    int do_onthefly = p->flags & MOD_COEFFICIENTWISE;
    int do_framewise = p->flags & MOD_FRAMEWISE;
    //int do_revorder = p->flags & ORDER_REV;

    LTFAT_COMPLEX** k = plan->k;
    LTFAT_COMPLEX* buf = plan->buf;

    int nfirst;

    for (int w = 0; w < W; w++)
    {
        const LTFAT_REAL* sChan =  s + w * M2 * N;
        LTFAT_COMPLEX* cChan = c + w * M2 * N;
        LTFAT_COMPLEX* coutChan = cout + w * M2 * N;

        extendborders(plan->plan_col, cChan, N, buf);

        /* Outside loop over columns */
        for (nfirst = 0; nfirst < N; nfirst++)
        {
            /* Pick the right kernel */
            LTFAT_COMPLEX* actK = k[nfirst % plan->kNo];
            /* Go to the n-th col in output*/
            LTFAT_COMPLEX* cColFirst = buf + nfirst * M2buf;
            LTFAT_COMPLEX* coutCol = coutChan + nfirst * M2;

            const LTFAT_REAL* sCol = sChan + nfirst * M2;

            leglaupdatereal_execute_col(plan->plan_col, sCol,
                                        actK, cColFirst, coutCol);
        }

        if (!do_onthefly && !do_framewise)
        {
            /* Update the phase only after the projection has been done. */
            for (int n = 0; n < N * M2; n++)
                coutChan[n] = sChan[n] * exp(I * ltfat_arg(coutChan[n]));
        }
    }
}

void
leglaupdatereal_execute_col(leglaupdate_plan_col* plan,
                            const LTFAT_REAL sCol[],
                            const LTFAT_COMPLEX actK[],
                            LTFAT_COMPLEX cColFirst[],
                            LTFAT_COMPLEX coutCol[])
{
    int m, mfirst, mlast;
    int M2 = plan->M / 2 + 1;
    int kernh = plan->ksize.height;
    int kernw = plan->ksize.width;
    int kernh2 = plan->ksize2.height;
    int kernw2 = plan->ksize2.width;
    int M2buf = M2 + kernh - 1;
    int kernhMidId = kernh2 - 1;
    int kernwMidId = kernw2 - 1;

    /* mwSignedIndex Nbuf = N + kernw -1; */
    int do_onthefly = plan->flags & MOD_COEFFICIENTWISE;
    int do_framewise = plan->flags & MOD_FRAMEWISE;

    /* Outside loop over rows */
    for (m = kernh2 - 1, mfirst = 0, mlast = kernh - 1; mfirst < M2;
         m++, mfirst++, mlast++)
    {
        LTFAT_COMPLEX accum = 0.0 + I * 0.0;

        /* inner loop over all cols of the kernel*/
        for (int kn = 0; kn < kernw; kn++)
        {
            /* mexPrintf("kn-loop: %d\n",kn); */
            const LTFAT_COMPLEX* actKCol = actK + kn * kernh2 +  kernh2 - 1;
            LTFAT_COMPLEX* cCol = cColFirst + kn * M2buf;

            /* Inner loop over half of the rows of the kernel excluding the middle row */
            for (int km = 0; km < kernh2 - 1; km++)
            {
                /* Doing the complex conjugated kernel elements simulteneously */
                LTFAT_REAL ar  = ltfat_real(actKCol[-km]);
                LTFAT_REAL ai  = ltfat_imag(actKCol[-km]);
                LTFAT_REAL br  = ltfat_real(cCol[mfirst + km]);
                LTFAT_REAL bi  = ltfat_imag(cCol[mfirst + km]);
                LTFAT_REAL bbr = ltfat_real(cCol[mlast - km]);
                LTFAT_REAL bbi = ltfat_imag(cCol[mlast - km]);
                accum += ar * (br + bbr) - ai * (bi - bbi)
                         + I * ( ar * (bi + bbi) + ai * (br - bbr));
            }

            /* The middle row is real*/
            accum += actKCol[-kernhMidId] * cCol[m];
        }

        coutCol[mfirst] = accum;

        if (do_onthefly)
        {
            /* Update the phase of a coefficient immediatelly */
            coutCol[mfirst] = sCol[mfirst] * exp(I * ltfat_arg(coutCol[mfirst]));
            cColFirst[kernwMidId * M2buf + m] = coutCol[mfirst];
        }

    }

    /* Update the phase of a single column */
    if (do_framewise)
    {
        for (m = kernh2 - 1, mfirst = 0; mfirst < M2; m++, mfirst++)
        {
            coutCol[mfirst] = sCol[mfirst] * exp(I * ltfat_arg(coutCol[mfirst]));
            cColFirst[kernwMidId * M2buf + m] = coutCol[mfirst];
        }
    }
}

void
extendborders(leglaupdate_plan_col* plan, const LTFAT_COMPLEX c[], int N,
              LTFAT_COMPLEX buf[])
{
    int m, n;
    int M2 = plan->M / 2 + 1;
    int kernh = plan->ksize.height;
    int kernw = plan->ksize.width;
    int kernh2 = plan->ksize2.height;
    int kernw2 = plan->ksize2.width;
    int M2buf = M2 + kernh - 1;
    int Nbuf = N + kernw - 1;

    if ( plan->flags & EXT_UPDOWN) Nbuf = N;

    if ( !(plan->flags & EXT_UPDOWN))
    {
        /* Copy input to the center of the buffer */
        for (n = 0; n < N; n++)
        {
            LTFAT_COMPLEX* bufstart = buf + (n + kernw2 - 1) * M2buf + kernh2 - 1;
            const LTFAT_COMPLEX* cstart = c + n * M2;
            memcpy(bufstart, cstart, M2 * sizeof * bufstart);
        }

        /* Periodically extend the left side */
        for (m = 0; m < M2; m++)
        {
            LTFAT_COMPLEX* buftarget = buf + kernh2 - 1 + m;
            LTFAT_COMPLEX* bufsource = buf + (N) * M2buf + kernh2 - 1 + m;

            for (n = 0; n < kernw2 - 1; n++)
                buftarget[n * M2buf] = bufsource[n * M2buf];
        }

        /* Periodically extend the right side*/
        for (m = 0; m < M2; m++)
        {
            LTFAT_COMPLEX* bufsource = buf + (kernw2 - 1) * M2buf + kernh2 - 1 + m;
            LTFAT_COMPLEX* buftarget = buf + (N + kernw2 - 1) * M2buf + kernh2 - 1 + m;

            for (n = 0; n < kernw2 - 1; n++)
                buftarget[n * M2buf] = bufsource[n * M2buf];
        }
    }
    else if (plan->flags & EXT_UPDOWN)
    {
        /* Copy input to the center of the buffer */
        for (n = 0; n < N; n++)
        {
            LTFAT_COMPLEX* bufstart = buf + n * M2buf + kernh2 - 1;
            const LTFAT_COMPLEX* cstart = c + n * M2;
            memcpy(bufstart, cstart, M2 * sizeof * bufstart);
        }
    }

    /* Conjugated odd-symmetric extention of the top border*/
    for (n = 0; n < Nbuf; n++)
    {
        LTFAT_COMPLEX* bufsource = buf + n * M2buf + kernh2;
        LTFAT_COMPLEX* buftarget = bufsource - 2;

        for (m = 0; m < kernh2 - 1; m++)
            buftarget[-m] = conj(bufsource[m]);
    }

    /* Conjugated odd or even symmetry extension of the bottom border.
     * Depending whether M is odd or even
     * */
    for (n = 0; n < Nbuf; n++)
    {
        LTFAT_COMPLEX* buftarget = buf + n * M2buf + kernh2 - 1 + M2;
        LTFAT_COMPLEX* bufsource = buftarget - 2 + plan->M % 2;

        for (m = 0; m < kernh2 - 1; m++)
            buftarget[m] = conj(bufsource[-m]);
    }

}

int
legla_set_status_callback(legla_plan* p, legla_callback_status* callback,
                          void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    p->status_callback = callback;
    p->status_callback_userdata = userdata;
error:
    return status;
}

int
legla_set_cmod_callback(legla_plan* p, legla_callback_cmod* callback,
                        void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    p->cmod_callback = callback;
    p->cmod_callback_userdata = userdata;
error:
    return status;
}
