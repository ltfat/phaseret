#include "ltfat/macros.h"
#include "phaseret/rtpghifb.h"
#include "phaseret/utils.h"
#include "float.h"
#include "rtpghifb_private.h"

#define NORTHFROMW(w,M,N) ((((w) + 1) % (M)) + (w) - (w) % (M))
#define SOUTHFROMW(w,M,N) (((w) - 1 + (M)) % (M) + (w) - (w) % (M))
#define EASTFROMW(w,M,N)  (((w) + (M)) % ((M) * (N)))
#define WESTFROMW(w,M,N)  (((w) - (M) + (M) * (N)) % ((M) * (N)))


PHASERET_API int
PHASERET_NAME(rtpghifb_set_causal)(PHASERET_NAME(rtpghifb_state)* p, int do_causal)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    p->do_causal = do_causal;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifb_set_tol)(PHASERET_NAME(rtpghifb_state)* p, double tol)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_NOTINRANGE, tol > 0 && tol < 1, "tol must be in range ]0,1[");

    p->p->tol = tol;
    p->p->logtol = log(tol + DBL_MIN);
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifb_init)(ltfat_int W, double* fc, ltfat_int a, ltfat_int M,
                           double gamma, double tol, int do_causal,
                           PHASERET_NAME(rtpghifb_state)** pout)
{
    int status = LTFATERR_SUCCESS;


    ltfat_int M2 = M;
    PHASERET_NAME(rtpghifb_state)* p = NULL;

    CHECK(LTFATERR_BADARG, !isnan(gamma) && gamma > 0,
          "gamma cannot be nan and must be positive. (Passed %f).", gamma);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive (passed %d)", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");
    CHECK(LTFATERR_NOTINRANGE, tol > 0 && tol < 1, "tol must be in range ]0,1[");

    CHECKMEM( p = (PHASERET_NAME(rtpghifb_state)*) ltfat_calloc(1, sizeof * p));

    CHECKSTATUS( PHASERET_NAME(rtpghifbupdate_init)( M, W, fc, tol, &p->p));
    CHECKMEM( p->slog =  LTFAT_NAME_REAL(calloc)(3 * M2 * W));
    CHECKMEM( p->tgrad = LTFAT_NAME_REAL(calloc)(3 * M2 * W));
    CHECKMEM( p->s =     LTFAT_NAME_REAL(calloc)(2 * M2 * W));
    CHECKMEM( p->fc = LTFAT_NAME_REAL(calloc)(M2 * W));
    CHECKMEM( p->fgrad = LTFAT_NAME_REAL(calloc)(M2 * W));
    CHECKMEM( p->phase = LTFAT_NAME_REAL(calloc)(M2 * W));

    p->do_causal = do_causal;
    p->M = M;
    p->a = a;
    p->gamma = gamma;
    p->W = W;

    *pout = p;
    return status;
error:
    if (p) PHASERET_NAME(rtpghifb_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifb_reset)(PHASERET_NAME(rtpghifb_state)* p,
                            const LTFAT_REAL** sinit)
{
    int status = LTFATERR_SUCCESS;
    ltfat_int M2, W;
    CHECKNULL(p);
    //M2 = p->M / 2 + 1;
    M2 = p->M;
    W = p->W;

    memset(p->slog, 0,  3 * M2 * W * sizeof * p->slog);
    memset(p->tgrad, 0, 3 * M2 * W * sizeof * p->tgrad);
    memset(p->s, 0,     2 * M2 * W * sizeof * p->s);
    memset(p->fgrad, 0, M2 * W * sizeof * p->tgrad);
    memset(p->phase, 0, M2 * W * sizeof * p->phase);
    memset(p->fc, 0, M2 * W * sizeof * p->fc);

    if (sinit)
        for (ltfat_int w = 0; w < W; w++)
            if (sinit[w])
            {
                memcpy(p->s + 2 * w * M2, sinit[w], 2 * M2 * sizeof * p->s );
                PHASERET_NAME(rtpghifblog)(sinit[w], 2 * M2, p->slog + M2 + w * 3 * M2);
            }

error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifb_execute)(PHASERET_NAME(rtpghifb_state)* p,
                              const LTFAT_REAL s[], LTFAT_COMPLEX c[])
{
    // n, n-1, n-2 frames
    // s is n-th
    //ltfat_int M2 = p->M / 2 + 1;
    ltfat_int M2 = p->M;
    ltfat_int W = p->W;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(s); CHECKNULL(c);

    for (ltfat_int w = 0; w < W; ++w)
    {
        LTFAT_REAL* slogCol = p->slog +   w * 3 * M2;
        LTFAT_REAL* tgradCol = p->tgrad + w * 3 * M2;
        LTFAT_REAL* sCol = p->s +         w * 2 * M2;
        LTFAT_REAL* fgradCol = p->fgrad + w * M2;
        LTFAT_REAL* phaseCol = p->phase + w * M2;

        // store log(s)
        PHASERET_NAME(shiftcolsleft)(slogCol, M2, 3, NULL);
        PHASERET_NAME(shiftcolsleft)(tgradCol, M2, 3, NULL);
        PHASERET_NAME(shiftcolsleft)(sCol, M2, 2, s + w * M2);

        PHASERET_NAME(rtpghifblog)(sCol + M2, M2, slogCol + 2 * M2);

        // Compute and store tgrad for n
        PHASERET_NAME(rtpghifbtgrad)(slogCol + 2 * M2, p->a, p->M, p->gamma,
                                   tgradCol + 2 * M2);

        // Compute fgrad for n or n-1
        PHASERET_NAME(rtpghifbfgrad)(slogCol, p->a, p->M, p->gamma, p->do_causal,
                                   fgradCol);

        PHASERET_NAME(rtpghifbupdate_execute)(p->p,
                                            p->do_causal ? slogCol + M2 : slogCol,
                                                    p->fc,
                                            p->do_causal ? tgradCol + M2 : tgradCol,
                                            fgradCol, phaseCol, phaseCol);

        // Combine phase with magnitude
        PHASERET_NAME(rtpghifbmagphase)(p->do_causal ? sCol + M2 : sCol, phaseCol, M2,
                                      c + w * M2);
    }

error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifb_done)(PHASERET_NAME(rtpghifb_state)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(rtpghifb_state)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    if (pp->p)     PHASERET_NAME(rtpghifbupdate_done)(&pp->p);
    if (pp->slog)  ltfat_free(pp->slog);
    if (pp->s)     ltfat_free(pp->s);
    if (pp->phase) ltfat_free(pp->phase);
    if (pp->tgrad) ltfat_free(pp->tgrad);
    if (pp->fgrad) ltfat_free(pp->fgrad);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifbupdate_init)(ltfat_int M, ltfat_int W, double* fc, double tol,
                                 PHASERET_NAME(rtpghifbupdate_plan)** pout)
{
    int status = LTFATERR_SUCCESS;

    //ltfat_int M2 = M / 2 + 1;
    PHASERET_NAME(rtpghifbupdate_plan)* p = NULL;
    CHECKMEM( p = (PHASERET_NAME(rtpghifbupdate_plan)*) ltfat_calloc(1, sizeof * p));
    //CHECKMEM( p->donemask = (int*) ltfat_calloc(M2, sizeof * p->donemask));
    CHECKMEM( p->donemask = (int*) ltfat_calloc(M, sizeof * p->donemask));
    
    //p->randphaseLen = 10 * M2 * W;
    p->randphaseLen = 10 * M * W;
    CHECKMEM( p->randphase = LTFAT_NAME_REAL(malloc)(p->randphaseLen));

    /* Do this somewhere else */
    for (ltfat_int ii = 0; ii < p->randphaseLen; ii++)
        p->randphase[ii] = (LTFAT_REAL)( 2.0 * M_PI * ((double)rand()) / (RAND_MAX));
      
    //CHECKMEM( p->fc = LTFAT_NAME_REAL(malloc)(M2));
    CHECKMEM( p->fc = LTFAT_NAME_REAL(malloc)(M));
    //for (ltfat_int ii = 0; ii < M2; ii++)
    for (ltfat_int ii = 0; ii < M; ii++)
        p->fc[ii] = fc[ii];
  
    p->logtol = log(tol);
    p->tol = tol;
    p->M = M;
    p->randphaseId = 0;
    //p->h = LTFAT_NAME(heap_init)(2 * M2, NULL);
    p->h = LTFAT_NAME(heap_init)(2 * M, NULL);

    *pout = p;
    return status;
error:
    if (p) PHASERET_NAME(rtpghifbupdate_done)(&p);
    return status;
}

PHASERET_API int
PHASERET_NAME(rtpghifbupdate_execute_withmask)(PHASERET_NAME(rtpghifbupdate_plan)* p,
                                             const LTFAT_REAL slog[],
                                             const LTFAT_REAL fc[],
                                             const LTFAT_REAL tgrad[],
                                             const LTFAT_REAL fgrad[],
                                             const LTFAT_REAL startphase[],
                                             const int mask[], LTFAT_REAL phase[])
{
    //ltfat_int M2 = p->M / 2 + 1;
    ltfat_int M = p->M;
    memcpy(p->donemask, mask, M * sizeof * p->donemask);

    return PHASERET_NAME(rtpghifbupdate_execute_common)(p, slog, fc, tgrad, fgrad, startphase, phase);
}

// slog: M x 2
// tgrad: M x 2
// fgrad: M x 1
// startphase: M x 1
// phase: M x 1
// donemask: M x 1
// heap must be able to hold 2*M values
PHASERET_API int
PHASERET_NAME(rtpghifbupdate_execute)(PHASERET_NAME(rtpghifbupdate_plan)* p,
                                    const LTFAT_REAL slog[],
                                    const LTFAT_REAL fc[],
                                    const LTFAT_REAL tgrad[],
                                    const LTFAT_REAL fgrad[],
                                    const LTFAT_REAL startphase[],
                                    LTFAT_REAL phase[])
{
    //ltfat_int M2 = p->M / 2 + 1;
    ltfat_int M = p->M;
    memset(p->donemask, 0, M * sizeof * p->donemask);
    return PHASERET_NAME(rtpghifbupdate_execute_common)(p, slog, fc, tgrad, fgrad, startphase, phase);
}


PHASERET_API int
PHASERET_NAME(rtpghifbupdate_execute_common)(PHASERET_NAME(rtpghifbupdate_plan)* p,
                                             const LTFAT_REAL slog[],
                                             const LTFAT_REAL fc[],
                                             const LTFAT_REAL tgrad[],
                                             const LTFAT_REAL fgrad[],
                                             const LTFAT_REAL startphase[],
                                             LTFAT_REAL phase[])
{
    LTFAT_NAME(heap)* h = p->h;
    //ltfat_int M2 = p->M / 2 + 1;
    //ltfat_int quickbreak = M2;
    ltfat_int quickbreak = p->M;
    const LTFAT_REAL oneover2 = (LTFAT_REAL) ( 1.0 / 2.0 );

    // We only need to compute M values, so perform quick exit
    // if we have them, but the heap is not yet empty.
    // (deleting from heap involves many operations)
    int* donemask = p->donemask;
    //const LTFAT_REAL* slog2 = slog + M2;
    const LTFAT_REAL* slog2 = slog + p->M;

    // Find max and the absolute thrreshold
    LTFAT_REAL logabstol = slog[0];
    //LTFAT_REAL logabstol = 0;
    //for (ltfat_int m = 1; m < 2 * M2; m++)
    for (ltfat_int m = 1; m < 2 * p->M; m++)
        if (slog[m] > logabstol)
            logabstol = slog[m];

    logabstol += (LTFAT_REAL) p->logtol;
    LTFAT_NAME(heap_reset)(h, slog);

    //for (ltfat_int m = 0; m < M2; m++)
    for (ltfat_int m = 0; m < p->M; m++)
    {
        if ( donemask[m] > 0 )
        {
            // We already know this one
            //LTFAT_NAME(heap_insert)(h, m + M2);
            LTFAT_NAME(heap_insert)(h, m + p->M);
            quickbreak--;
        }
        else
        {
            if ( slog2[m] <= logabstol )
            {
                // This will get a randomly generated phase
                donemask[m] = -1;
                quickbreak--;
            }
            else
            {
                LTFAT_NAME(heap_insert)(h, m);
            }
        }
    }

    ltfat_int w = -1;
    
    while ( (quickbreak > 0) && ( w = LTFAT_NAME(heap_delete)(h) ) >= 0 )
    {
        //if ( w >= M2 )
        if ( w >= p->M )
        {
            // Next frame
            //ltfat_int wprev = w - M2;
            ltfat_int wprev = w - p->M;

            //if ( wprev != M2 - 1 && !donemask[wprev + 1] )
            if ( wprev != p->M - 1 && !donemask[wprev + 1] )
            {
                //frequency step upwards
                LTFAT_REAL step = fc[(wprev+1)] - fc[wprev];
                //printf("f-step up: %f\n", step);

                if (step < 0)
                    step += 2;

                phase[wprev + 1] = phase[wprev] + step*(fgrad[wprev] + fgrad[wprev + 1]) * oneover2;
                donemask[wprev + 1] = 1;

                LTFAT_NAME(heap_insert)(h, w + 1);
                quickbreak--;
            }

            if ( wprev != 0 && !donemask[wprev - 1] )
            {
                //frequency step downwards
                LTFAT_REAL step = fc[(wprev-1)] - fc[wprev];
                //printf("f-step down: %f\n", step);
                if (step > 0)
                    step -= 2;
                
                phase[wprev - 1] = phase[wprev] + step*(fgrad[wprev] + fgrad[wprev - 1]) * oneover2;
                donemask[wprev - 1] = 1;

                LTFAT_NAME(heap_insert)(h, w - 1);
                quickbreak--;
            }
        }
        else
        {
            // Current frame
            if ( !donemask[w] )
            {
                //ltfat_int wnext = w + M2;
                ltfat_int wnext = w + p->M;
                phase[w] = startphase[w] + (tgrad[w] + tgrad[wnext]) * oneover2;
                donemask[w] = 1;

                LTFAT_NAME(heap_insert)(h, wnext);
                quickbreak--;
            }
        }
    }

    // Fill in values below tol
    //for (ltfat_int ii = 0; ii < M2; ii++)
    for (ltfat_int ii = 0; ii < p->M; ii++)
    {
        if (donemask[ii] < 0)
        {
            phase[ii] = p->randphase[p->randphaseId++];
            p->randphaseId %= p->randphaseLen;
        }
    }

    return 0;
}

PHASERET_API int
PHASERET_NAME(rtpghifbupdate_done)(PHASERET_NAME(rtpghifbupdate_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    PHASERET_NAME(rtpghifbupdate_plan)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    if (pp->h)         LTFAT_NAME(heap_done)(pp->h);
    if (pp->donemask)  ltfat_free(pp->donemask);
    if (pp->randphase) ltfat_free(pp->randphase);
    if (pp->fc) ltfat_free(pp->fc);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}



PHASERET_API int
PHASERET_NAME(rtpghifboffline)(const LTFAT_REAL* s, double* fc, ltfat_int L,
                             ltfat_int W, ltfat_int a, ltfat_int M,
                             double gamma, double tol, int do_causal,
                             LTFAT_COMPLEX* c)
{
    ltfat_int N = L / a;
    //ltfat_int M2 = M / 2 + 1;
    ltfat_int M2 = M;
    PHASERET_NAME(rtpghifb_state)* p = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(s); CHECKNULL(c);

    CHECKSTATUS( PHASERET_NAME(rtpghifb_init)(1, fc, a, M, gamma, tol, do_causal, &p));

    if (do_causal)
    {
        for (ltfat_int w = 0; w < W; w++)
        {

            const LTFAT_REAL* schan = s + w * N * M2;
            PHASERET_NAME(rtpghifb_reset)(p, &schan);

            for (ltfat_int n = 0; n < N; ++n)
            {
                const LTFAT_REAL* sncol = schan + n * M2;
                LTFAT_COMPLEX* cncol =    c + n * M2 + w * N * M2;
                PHASERET_NAME(rtpghifb_execute)(p, sncol, cncol);
            }
        }
    }
    else
    {
        for (ltfat_int w = 0; w < W; w++)
        {
            const LTFAT_REAL* schan = s + w * N * M2;
            PHASERET_NAME(rtpghifb_reset)(p, &schan);

            for (ltfat_int n = 0, nahead = 1; nahead < N; ++n, ++nahead)
            {
                const LTFAT_REAL* sncol = schan + nahead * M2;
                LTFAT_COMPLEX* cncol =    c + n * M2 + w * N * M2;
                PHASERET_NAME(rtpghifb_execute)(p, sncol, cncol);
            }

            PHASERET_NAME(rtpghifb_execute)(p, s + w * N * M2, c + (N - 1) * M2 + w * N * M2);
        }
    }

    PHASERET_NAME(rtpghifb_done)(&p);
error:
    return status;
}


void
PHASERET_NAME(rtpghifbfgrad)(const LTFAT_REAL* logs, ltfat_int a, ltfat_int M,
                           double gamma,
                           int do_causal, LTFAT_REAL* fgrad)
{
    //ltfat_int M2 = M / 2 + 1;
    ltfat_int M2 = M;

    const LTFAT_REAL fgradmul = (const LTFAT_REAL)( -gamma / (2.0 * a * M));
    const LTFAT_REAL* scol0 = logs;
    const LTFAT_REAL* scol2 = logs + 2 * M2;

    if (do_causal)
    {
        const LTFAT_REAL* scol1 = logs + M2;

        for (ltfat_int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * ((LTFAT_REAL)(3.0)* scol2[m] - (LTFAT_REAL)(4.0) * scol1[m] + scol0[m]);
    }
    else
    {
        for (ltfat_int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (scol2[m] - scol0[m]);
    }
}

void
PHASERET_NAME(rtpghifbtgrad)(const LTFAT_REAL* logs, ltfat_int a, ltfat_int M,
                           double gamma,
                           LTFAT_REAL* tgrad)
{
    //ltfat_int M2 = M / 2 + 1;
    ltfat_int M2 = M;

    const LTFAT_REAL tgradmul = (const LTFAT_REAL)( (a * M) / (gamma * 2.0));
    const LTFAT_REAL tgradplus = (const LTFAT_REAL)( 2.0 * M_PI * a / ((double)M) );

    tgrad[0]      = 0.0;
    tgrad[M2 - 1] = 0.0;

    for (ltfat_int m = 1; m < M2 - 1; m++)
        tgrad[m] = tgradmul * (logs[m + 1] - logs[m - 1]) + tgradplus * m;
}

void
PHASERET_NAME(rtpghifblog)(const LTFAT_REAL* in, ltfat_int L, LTFAT_REAL* out)
{
#ifdef LTFAT_DOUBLE
    LTFAT_REAL eps = DBL_MIN;
#elif defined(LTFAT_SINGLE)
    LTFAT_REAL eps = FLT_MIN;
#endif

    for (ltfat_int l = 0; l < L; l++)
        out[l] = log(in[l] + eps);
}

void
PHASERET_NAME(rtpghifbmagphase)(const LTFAT_REAL* s, const LTFAT_REAL* phase,
                              ltfat_int L, LTFAT_COMPLEX* c)
{
    for (ltfat_int l = 0; l < L; l++)
        c[l] = s[l] * (cos(phase[l]) + I * sin(phase[l]));
}

