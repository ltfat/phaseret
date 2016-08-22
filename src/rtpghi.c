#include "ltfat/macros.h"
#include "phaseret/rtpghi.h"
#include "phaseret/utils.h"
#include "float.h"

struct rtpghi_state
{
    int M;
    int a;
    int W;
    int* mask;
    int do_causal;
    LTFAT_REAL* slog;
    LTFAT_REAL* s;
    LTFAT_REAL* tgrad; //!< Time gradient buffer
    LTFAT_REAL* fgrad; //!< Frequency gradient buffer
    LTFAT_REAL* phase; //!< Buffer for keeping previously computed frame
    double logtol;
    double tol;
    double gamma;
    LTFAT_REAL* randphase; //!< Precomputed array of random phase
    int randphaseLen;
    int randphaseId;
    struct ltfat_heapinttask_d* hit;
};

int
rtpghi_set_causal(rtpghi_state* p, int do_causal)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    p->do_causal = do_causal;
error:
    return status;
}

int
rtpghi_init(double gamma, int W, int a, int M, double tol, int do_causal,
            rtpghi_state** pout)
{
    int status = LTFATERR_SUCCESS;

    int M2 = M / 2 + 1;
    rtpghi_state* p = NULL;
    CHECKMEM( p = (rtpghi_state*) calloc(1, sizeof * p));

    p->hit = ltfat_heapinttask_init_d( M2, 2, 2 * M2 , NULL, 1);
    CHECKMEM( p->slog = LTFAT_NAME_REAL(calloc)(3 * M2 * W));
    CHECKMEM( p->s = LTFAT_NAME_REAL(calloc)(2 * M2 * W));
    CHECKMEM( p->phase = LTFAT_NAME_REAL(calloc)(2 * M2 * W));
    CHECKMEM( p->tgrad = LTFAT_NAME_REAL(calloc)(3 * M2 * W));
    CHECKMEM( p->fgrad = LTFAT_NAME_REAL(calloc)(2 * M2 * W));

    p->randphaseLen = 10 * M2 * W;
    CHECKMEM( p->randphase = LTFAT_NAME_REAL(malloc)(p->randphaseLen));

    for (int ii = 0; ii < p->randphaseLen; ii++)
        p->randphase[ii] = 2.0 * M_PI * ((double)rand()) / RAND_MAX;

    CHECKMEM( p->mask = (int*) malloc(2 * M2 * W * sizeof * p->mask));

    for (int w = 0; w < W; ++w)
    {
        for (int ii = 0; ii < M2; ++ii)
        {
            p->mask[ii + w * 2 * M2] = 1;
            p->mask[ii + M2 + w * 2 * M2] = 0;
        }
    }

    p->do_causal = do_causal;
    p->logtol = log(tol + DBL_MIN);
    p->tol = tol;
    p->M = M;
    p->a = a;
    p->gamma = gamma;
    p->W = W;

    *pout = p;
    return status;
error:
    if (p)
    {
        if (p->slog) free(p->slog);
        if (p->s) free(p->s);
        if (p->tgrad) free(p->tgrad);
        if (p->fgrad) free(p->fgrad);
        if (p->phase) free(p->phase);
        if (p->randphase) free(p->randphase);
        if (p->hit) ltfat_heapinttask_done_d(p->hit);
        free(p);
    }
    return status;
}


int
rtpghi_execute(rtpghi_state* p, const LTFAT_REAL s[], LTFAT_COMPLEX c[])
{
    // n, n-1, n-2 frames
    // s is n-th
    int M2 = p->M / 2 + 1;
    int W = p->W;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(s); CHECKNULL(c);

    for (int w = 0; w < W; ++w)
    {
        LTFAT_REAL* slogCol = p->slog + w * 3 * M2;
        LTFAT_REAL* sCol = p->s + w * 2 * M2;
        LTFAT_REAL* phaseCol = p->phase + w * 2 * M2;
        LTFAT_REAL* tgradCol = p->tgrad + w * 3 * M2;
        LTFAT_REAL* fgradCol = p->fgrad + w * 2 * M2;

        // store log(s)
        shiftcolsleft(slogCol, M2, 3, NULL);
        shiftcolsleft(sCol, M2, 2, s + w * M2);
        rtpghilog(sCol + M2, M2, slogCol + 2 * M2);

        // Store current s

        // Compute and store tgrad for n
        shiftcolsleft(tgradCol, M2, 3, NULL);
        rtpghitgrad(slogCol + 2 * M2, p->a, p->M, p->gamma, tgradCol + 2 * M2);

        // Compute fgrad for n or n-1
        rtpghifgrad(slogCol, p->a, p->M, p->gamma, p->do_causal, fgradCol + M2);

        ltfat_heapinttask_resetmask_d(p->hit,
                                      p->mask, p->do_causal ? slogCol + M2 : slogCol, p->logtol, 1);

        ltfat_heapint_execute_d(p->hit,
                                p->do_causal ? tgradCol + M2 : tgradCol, fgradCol, phaseCol);

        shiftcolsleft(phaseCol, M2, 2, NULL);

        // Fill in the missing phase from the precomputed random array
        int* donemask = ltfat_heapinttask_get_mask_d(p->hit);

        for (int ii = 0; ii < M2; ii++)
        {
            if (donemask[M2 + ii] <= LTFAT_MASK_UNKNOWN)
            {
                phaseCol[ii] = p->randphase[p->randphaseId++];
                p->randphaseId %= p->randphaseLen;
            }
        }

        // Combine phase with magnitude
        rtpghimagphase(p->do_causal ? sCol + M2 : sCol, phaseCol, M2, c + w * M2);
    }

error:
    return status;
}

int
rtpghi_done(rtpghi_state** p)
{
    int status = LTFATERR_SUCCESS;
    rtpghi_state* pp; 
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    ltfat_heapinttask_done_d(pp->hit);
    free(pp->slog);
    free(pp->s);
    free(pp->mask);
    free(pp->phase);
    free(pp->tgrad);
    free(pp->fgrad);
    free(pp->randphase);
    free(pp);
    pp = NULL;
error:
    return status;
}

int
rtpghioffline(const LTFAT_REAL* s, double gamma, int L, int W, int a, int M,
              double tol, int do_causal, LTFAT_COMPLEX* c)
{
    int N = L / a;
    int M2 = M / 2 + 1;
    rtpghi_state* p = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(s); CHECKNULL(c);

    CHECKSTATUS( rtpghi_init(gamma, W, a, M, tol, do_causal, &p),
                 "rtpghi init failed");

    if (do_causal)
    {
        for (int n = 0; n < N; ++n)
        {
            const LTFAT_REAL* sncol = s + n * M2;
            LTFAT_COMPLEX* cncol = c + n * M2;
            rtpghi_execute(p, sncol, cncol);
        }
    }
    else
    {
        rtpghi_execute(p, s, c);

        for (int n = 0, nahead = 1; nahead < N; ++n, ++nahead)
        {
            const LTFAT_REAL* sncol = s + nahead * M2;
            LTFAT_COMPLEX* cncol = c + n * M2;
            rtpghi_execute(p, sncol, cncol);
        }

        rtpghi_execute(p, NULL, c + (N - 1) * M2);
    }

    rtpghi_done(&p);
error:
    return status;
}


void
rtpghifgrad(const LTFAT_REAL* logs, int a, int M, double gamma,
            int do_causal, LTFAT_REAL* fgrad)
{
    int M2 = M / 2 + 1;

    const LTFAT_REAL fgradmul = -gamma / (2.0 * a * M);
    const LTFAT_REAL* scol0 = logs;
    const LTFAT_REAL* scol2 = logs + 2 * M2;

    if (do_causal)
    {
        const LTFAT_REAL* scol1 = logs + M2;

        for (int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (3.0 * scol2[m] - 4.0 * scol1[m] + scol0[m]);
    }
    else
    {
        for (int m = 0; m < M2; ++m)
            fgrad[m] = fgradmul * (scol2[m] - scol0[m]);
    }
}

void
rtpghitgrad(const LTFAT_REAL* logs, int a, int M, double gamma,
            LTFAT_REAL* tgrad)
{
    int M2 = M / 2 + 1;

    const LTFAT_REAL tgradmul = (a * M) / (gamma * 2.0);
    const LTFAT_REAL tgradplus = 2.0 * M_PI * a / M;

    tgrad[0]      = 0.0;
    tgrad[M2 - 1] = 0.0;

    for (int m = 1; m < M2 - 1; m++)
        tgrad[m] = tgradmul * (logs[m + 1] - logs[m - 1]) + tgradplus * m;
}

void
rtpghilog(const LTFAT_REAL* in, int L, LTFAT_REAL* out)
{
    for (int l = 0; l < L; l++)
        out[l] = log(in[l] + DBL_MIN);

}

void
rtpghimagphase(const LTFAT_REAL* s, const LTFAT_REAL* phase, int L, LTFAT_COMPLEX* c)
{
    for (int l = 0; l < L; l++)
        c[l] = s[l] * exp(I * phase[l]);
}


double
firwin2gamma(LTFAT_FIRWIN win, int gl)
{
    if( gl <= 0)
    {
        return NAN;
    }
    double gamma;

    switch (win)
    {
    case LTFAT_HANN:
    case LTFAT_HANNING:
    case LTFAT_NUTTALL10:
        gamma = 0.25645; break;
    case LTFAT_SQRTHANN:
    case LTFAT_COSINE:
    case LTFAT_SINE:
        gamma = 0.41532; break;
    case LTFAT_HAMMING:
        gamma = 0.29794; break;
    case LTFAT_NUTTALL01:
        gamma = 0.29610; break;
    case LTFAT_SQUARE:
    case LTFAT_RECT:
        gamma = 0.85732; break;
    case LTFAT_TRIA:
    case LTFAT_TRIANGULAR:
    case LTFAT_BARTLETT:
        gamma = 0.27561; break;
    case LTFAT_SQRTTRIA:
        gamma = 0.48068; break;
    case LTFAT_BLACKMAN:
        gamma = 0.17954; break;
    case LTFAT_BLACKMAN2:
        gamma = 0.18465; break;
    case LTFAT_NUTTALL:
    case LTFAT_NUTTALL12:
        gamma = 0.12807; break;
    case LTFAT_OGG:
    case LTFAT_ITERSINE:
        gamma = 0.35744; break;
    case LTFAT_NUTTALL20:
        gamma = 0.14315; break;
    case LTFAT_NUTTALL11:
        gamma = 0.17001; break;
    case LTFAT_NUTTALL02:
        gamma = 0.18284; break;
    case LTFAT_NUTTALL30:
        gamma = 0.09895; break;
    case LTFAT_NUTTALL21:
        gamma = 0.11636; break;
    case LTFAT_NUTTALL03:
        gamma = 0.13369; break;
    default:
        return NAN;
    };

    gamma *= gl * gl;

    return gamma;

}
