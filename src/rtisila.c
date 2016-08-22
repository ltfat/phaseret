#include "phaseret/rtisila.h"
#include "phaseret/utils.h"
#include "ltfat/macros.h"

struct rtisilaupdate_plan
{
    LTFAT_REAL* frame; //!< Time domain buffer
    LTFAT_COMPLEX* fftframe; //!< Frequency domain buffer
    fftw_plan fwdp; //!< Real FFT plan
    fftw_plan backp; //!< Real IFFT plan
    const LTFAT_REAL* g;
    const LTFAT_REAL* gd;
    const LTFAT_REAL* specg1;
    const LTFAT_REAL* specg2;
    int gl;
    int M;
    int a;
};

struct rtisila_state
{
    rtisilaupdate_plan* uplan;
    int maxLookahead;
    int lookahead;
    int lookback;
    int maxit;
    int W;
    LTFAT_REAL* frames; //!< Buffer for time-domain frames
    LTFAT_REAL* s; //!< Buffer for target magnitude
    void** garbageBin;
    int garbageBinSize;
};

int
rtisila_set_lookahead(rtisila_state* p, int lookahead)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_BADARG, lookahead >= 0 && lookahead <= p->lookahead,
          "lookahead can only be in range [0-%d] (passed %d).",
          p->lookahead, lookahead);

    p->lookahead = lookahead;
error:
    return status;
}


void
overlaynthframe(const LTFAT_REAL* frames, int gl, int N, int a, int n,
                LTFAT_REAL* frame)
{
    // Copy n-th frame
    memcpy(frame, frames + gl * n, gl * sizeof * frame);

    // Overlap frames n+1 ... N-1 or until there is overlap
    for (int ii = n + 1; ii < N; ii++)
    {
        int jj = (ii - n) * a;
        int nSamp = gl - jj;

        if (nSamp <= 0) break;

        const LTFAT_REAL* framestmp = frames + gl * ii;

        for (int kk = 0; kk < nSamp; kk++)
            frame[jj + kk] += framestmp[kk];
    }

    // Overlap frames n-1 ... 0, or until there is no overlap
    for (int ii = n - 1; ii >= 0; ii--)
    {
        int jj = (n - ii) * a;
        int nSamp = gl - jj;

        if (nSamp <= 0) break;

        const LTFAT_REAL* framestmp = frames + gl * ii;

        for (int kk = 0; kk < nSamp; kk++)
            frame[kk] += framestmp[jj + kk];
    }
}

void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const LTFAT_REAL* frames,
                       const LTFAT_REAL* g, int n, int N)
{
    int M = p->M;
    int gl = p->gl;
    int a = p->a;

    overlaynthframe(frames, gl, N, a, n, p->frame);

    /* // Copy n-th frame */
    /* memcpy(p->frame, frames + gl * n, M * sizeof * p->frame); */
    /*  */
    /* // Overlap frames n+1 ... N-1 or until there is overlap */
    /* for (int ii = n + 1; ii < N; ii++) */
    /* { */
    /*     int jj = (ii - n) * a; */
    /*     int nSamp = gl - jj; */
    /*  */
    /*     if (nSamp <= 0) break; */
    /*  */
    /*     const double* framestmp = frames + gl * ii; */
    /*  */
    /*     for (int kk = 0; kk < nSamp; kk++) */
    /*         p->frame[jj + kk] += framestmp[kk]; */
    /* } */
    /*  */
    /* // Overlap frames n-1 ... 0, or until there is no overlap */
    /* for (int ii = n - 1; ii >= 0; ii--) */
    /* { */
    /*     int jj = (n - ii) * a; */
    /*     int nSamp = gl - jj; */
    /*  */
    /*     if (nSamp <= 0) break; */
    /*  */
    /*     const double* framestmp = frames + gl * ii; */
    /*  */
    /*     for (int kk = 0; kk < nSamp; kk++) */
    /*         p->frame[kk] += framestmp[jj + kk]; */
    /* } */

    // Multiply with analysis window
    for (int m = 0; m < gl; m++)
        p->frame[m] *= g[m];

    // Set remaining samples to zeros
    for (int m = gl; m < M; m++)
        p->frame[m] = 0.0;
}

void
rtisilaphaseupdate(rtisilaupdate_plan* p, const LTFAT_REAL* sframe,
                   LTFAT_REAL* frameupd, LTFAT_COMPLEX* c)
{
    int M = p->M;
    int gl = p->gl;
    int M2 = M / 2 + 1;

    // FFTREAL
    fftw_execute(p->fwdp);

    // Force the magnitude to target magnitude
    for (int m = 0; m < M2; m++)
        p->fftframe[m] = sframe[m] * exp(I * ltfat_arg(p->fftframe[m]));

    // Copy before it gets overwritten
    if (c)
        memcpy(c, p->fftframe, M2 * sizeof * c);

    // IFFTREAL // Overwrites input p->fftframe
    fftw_execute(p->backp);

    // Multiply with the synthesis window
    for (int m = 0; m < gl; m++)
        p->frame[m] *= p->gd[m];

    // Set remaining samples to zeros
    for (int m = gl; m < M; m++)
        p->frame[m] = 0;

    // Copy to the output
    memcpy(frameupd, p->frame, gl * sizeof * frameupd);
}

int
rtisilaupdate_init(const LTFAT_REAL* g, const LTFAT_REAL* specg1, const LTFAT_REAL* specg2,
                   const LTFAT_REAL* gd, const int gl, int a, int M,
                   rtisilaupdate_plan** pout)
{
    int status = LTFATERR_SUCCESS;
    rtisilaupdate_plan* p = NULL;
    int M2 = M / 2 + 1;

    CHECKMEM( p = (rtisilaupdate_plan*) calloc(1, sizeof * p));
    p->M = M; p->a = a; p->g = g; p->gd = gd; p->specg1 = specg1; 
    p->specg2 = specg2; p->gl = gl;

    // Real input array for FFTREAL and output array for IFFTREAL
    CHECKMEM( p->frame = LTFAT_NAME_REAL(malloc)(M));
    // Complex output array for FFTREAL and input array to IFFTREAL
    CHECKMEM( p->fftframe = LTFAT_NAME_COMPLEX(malloc)(M2));

    // FFTREAL plan
    p->fwdp = fftw_plan_dft_r2c_1d(M, p->frame,
                                   (LTFAT_FFTW(complex)*) p->fftframe, FFTW_MEASURE);
    CHECKINIT(p->fwdp, "FFTW plan failed");
    // IFFTREAL plan
    p->backp = fftw_plan_dft_c2r_1d(M,(LTFAT_FFTW(complex)*) p->fftframe,
                                    p->frame, FFTW_MEASURE);
    CHECKINIT(p->backp, "FFTW plan failed");

    *pout = p;
    return status;
error:
    if (p)
    {
        if (p->frame) free(p->frame);
        if (p->fftframe) free(p->fftframe);
        if (p->fwdp) fftw_destroy_plan(p->fwdp);
        if (p->backp) fftw_destroy_plan(p->backp);
        free(p);
    }
    return status;
}

int
rtisilaupdate_done(rtisilaupdate_plan** p)
{
    int status = LTFATERR_SUCCESS;
    rtisilaupdate_plan* pp; 
    CHECKNULL(p); CHECKNULL(*p);

    pp = *p;
    fftw_free(pp->frame);
    fftw_free(pp->fftframe);
    fftw_destroy_plan(pp->fwdp);
    fftw_destroy_plan(pp->backp);
    free(pp);
    pp = NULL;
error:
    return status;
}

void
rtisilaupdate_execute(rtisilaupdate_plan* p, const LTFAT_REAL* frames, int N,
                      const LTFAT_REAL* s, int lookahead, int maxit, LTFAT_REAL* frames2,
                      LTFAT_COMPLEX* c)
{
    int lookback = N - lookahead - 1;
    int M = p->M;
    int gl = p->gl;
    int M2 = M / 2 + 1;

    // If we are not working inplace ...
    if (frames != frames2)
        memcpy(frames2, frames, gl * N  * sizeof * frames);

    for (int it = 0; it < maxit; it++)
    {
        for (int nback = lookahead; nback >= 0; nback--)
        {
            int indx = lookback + nback;

            // Newest lookahead frame is treated differently
            if (nback == lookahead)
                if (it == 0)
                    rtisilaoverlaynthframe(p, frames2, p->specg1, indx, N);
                else
                    rtisilaoverlaynthframe(p, frames2, p->specg2, indx, N);
            else
                rtisilaoverlaynthframe(p, frames2, p->g, indx, N);

            if ( nback == 0 && it == (maxit - 1) )
                rtisilaphaseupdate(p, s + nback * M2, frames2 +  indx * gl, c);
            else
                rtisilaphaseupdate(p, s + nback * M2, frames2 +  indx * gl, NULL);
        }
    }

    if (c)
        ltfat_fftrealcircshift_dc(c, M, -(gl / 2) , c);
}

void
rtisilaupdate(const LTFAT_REAL* frames,
              const LTFAT_REAL* g, const LTFAT_REAL* specg1, const LTFAT_REAL* specg2,
              const LTFAT_REAL* gd, const int gl,
              int a, int M, int N, const LTFAT_REAL* s, int lookahead, int maxit,
              LTFAT_REAL* frames2)
{
    rtisilaupdate_plan* p = NULL;
    rtisilaupdate_init(g, specg1, specg2, gd, gl, a, M, &p);
    rtisilaupdate_execute(p, frames, N, s, lookahead, maxit, frames2, NULL);
    rtisilaupdate_done(&p);
}

void
rtisilaupdatecoef(const LTFAT_REAL* frames,
                  const LTFAT_REAL* g, const LTFAT_REAL* specg1, const LTFAT_REAL* specg2,
                  const LTFAT_REAL* gd, const int gl,
                  int a, int M, int N, const LTFAT_REAL* s, int lookahead, int maxit,
                  LTFAT_REAL* frames2, LTFAT_COMPLEX* c)
{
    rtisilaupdate_plan* p = NULL;
    rtisilaupdate_init(g, specg1, specg2, gd, gl, a, M, &p);
    rtisilaupdate_execute(p, frames, N, s, lookahead, maxit, frames2, c);
    rtisilaupdate_done(&p);
}

int
rtisila_reset(rtisila_state* p)
{
    int N, W, gl, M2;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    N = p->lookback + 1 + p->maxLookahead;
    W = p->W;
    M2 = p->uplan->M / 2 + 1;
    gl = p->uplan->gl;

    memset(p->s, 0, M2 * (1 + p->maxLookahead) * W * sizeof * p->s);
    memset(p->frames, 0, gl * N * W * sizeof * p->frames);
error:
    return status;
}

int
rtisila_init(const LTFAT_REAL* g, const int gl, const int W,
             int a, int M, int lookahead, int maxit, rtisila_state** pout)
{
    int status = LTFATERR_SUCCESS;

    rtisila_state* p = NULL;
    LTFAT_REAL* wins = NULL;

    int M2, lookback, winsNo, maxLookahead;
    CHECKNULL(g); CHECKNULL(pout);
    CHECK(LTFATERR_BADSIZE, gl > 0, "gl must be positive (passed %d)", gl);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive (passed %d)", W);
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive (passed %d)", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive (passed %d)", M);
    CHECK(LTFATERR_BADARG, lookahead >= 0,
          "lookahead >=0 failed (passed %d)", lookahead);
    CHECK(LTFATERR_NOTPOSARG, maxit > 0,
          "maxit must be positive (passed %d)", maxit);

    M2 = M / 2 + 1;
    lookback = ceil(((LTFAT_REAL) gl) / a) - 1;
    winsNo = (2 * lookback + 1);
    maxLookahead = lookahead;
    CHECKMEM( p = (rtisila_state*) calloc(1, sizeof * p));

    CHECKSTATUS(
        rtisilaupdate_init(NULL, NULL, NULL, NULL, gl, a, M, &p->uplan),
        "rtisiupdate ini failed");

    CHECKMEM( p->uplan->g = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM( p->uplan->gd = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM( p->uplan->specg1 = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM( p->uplan->specg2 = LTFAT_NAME_REAL(malloc)(gl));

    CHECKSTATUS( ltfat_gabdual_painless_d(g, gl, a, M, (LTFAT_REAL*) p->uplan->gd),
                 "gabdual failed");

    ltfat_fftshift_d(g, gl, (LTFAT_REAL*) p->uplan->g);
    ltfat_fftshift_d(p->uplan->gd, gl, (LTFAT_REAL*) p->uplan->gd);

    CHECKMEM( wins = LTFAT_NAME_REAL(malloc)(gl));

    // Spec2
    ltfat_periodize_array_d(p->uplan->gd, gl, gl * winsNo, wins);
    overlaynthframe(wins, gl, winsNo, a, 0, (LTFAT_REAL*) p->uplan->specg2);

    // Spec1
    memset(wins, 0, gl * sizeof * wins);
    overlaynthframe(wins, gl, winsNo, a, 0, (LTFAT_REAL*) p->uplan->specg1);

    free(wins); wins = NULL;

    CHECKMEM( p->frames = 
            LTFAT_NAME_REAL(calloc)(gl * (lookback + 1 + maxLookahead) * W));
    CHECKMEM( p->s = LTFAT_NAME_REAL(calloc)( M2 * (1 + maxLookahead) * W));

    p->garbageBinSize = 4;
    CHECKMEM( p->garbageBin = (void**)malloc(p->garbageBinSize * sizeof(void*)));
    p->garbageBin[0] = (void*) p->uplan->g;
    p->garbageBin[1] = (void*) p->uplan->gd;
    p->garbageBin[2] = (void*) p->uplan->specg1;
    p->garbageBin[3] = (void*) p->uplan->specg2;

    p->lookback = lookback;
    p->lookahead = lookahead;
    p->maxLookahead = maxLookahead;
    p->maxit = maxit;
    p->W = W;

    *pout = p;
    return status;
error:
    if (wins) free(wins);
    if (p)
    {
        if (p->uplan)
        {
            if (p->uplan->g) free((void*)p->uplan->g);
            if (p->uplan->gd) free((void*)p->uplan->gd);
            if (p->uplan->specg1) free((void*)p->uplan->specg1);
            if (p->uplan->specg2) free((void*)p->uplan->specg2);
            rtisilaupdate_done(&p->uplan);
        }
        if (p->s) free(p->s);
        if (p->frames) free(p->frames);
        free(p);
    }
    *pout = NULL;
    return status;
}

int
rtisila_init_win(LTFAT_FIRWIN win, int gl, int W, int a, int M,
                 int lookahead, int maxit, rtisila_state** pout)
{
    LTFAT_REAL* g = NULL;
    int status = LTFATERR_SUCCESS;
    int initstatus;
    CHECKMEM( g = LTFAT_NAME_REAL(malloc)(gl));

    // Analysis window
    CHECKSTATUS( ltfat_firwin_d(win, gl, g), "firwin failed");

    initstatus = rtisila_init(g, gl, W, a, M, lookahead, maxit, pout);

    free(g);
    return initstatus;
error:
    if (g) free(g);
    return status;
}

int
rtisila_done(rtisila_state** p)
{
    int status = LTFATERR_SUCCESS;
    rtisila_state* pp; 
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    CHECKSTATUS( rtisilaupdate_done(&pp->uplan),
                 "rtisilaupdate done failed");

    free(pp->s);
    free(pp->frames);

    if (pp->garbageBinSize)
    {
        for (int ii = 0; ii < pp->garbageBinSize; ii++)
            free(pp->garbageBin[ii]);

        free(pp->garbageBin);
    }

    free(pp);
    pp = NULL;
error:
    return status;
}

int
rtisila_execute(rtisila_state* p, const LTFAT_REAL* s, LTFAT_COMPLEX* c)
{
    int M, gl, M2, noFrames, N;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(s); CHECKNULL(c);

    M = p->uplan->M;
    gl = p->uplan->gl;
    M2 = M / 2 + 1;
    noFrames = p->lookback + 1 + p->lookahead;
    N = p->lookback + 1 + p->maxLookahead;

    for (int w = 0; w < p->W; w++)
    {
        const LTFAT_REAL* schan = s + w * M2;
        LTFAT_COMPLEX* cchan = c + w * M2;
        LTFAT_REAL* frameschan = p->frames + w * N * gl;
        LTFAT_REAL* sframeschan = p->s  + w * (1 + p->maxLookahead) * M2;
        // Shift frames buffer
        shiftcolsleft( frameschan, gl, noFrames, NULL);

        // Shift scols buffer
        shiftcolsleft( sframeschan, M2, p->lookahead + 1, schan);

        rtisilaupdate_execute(p->uplan, frameschan, noFrames, sframeschan,
                              p->lookahead, p->maxit, frameschan, cchan);
    }

error:
    return status;
}


int
rtisilaoffline(const LTFAT_REAL s[], const LTFAT_REAL g[],
               int L, int gl, int W, int a, int M,
               int lookahead, int maxit, LTFAT_COMPLEX c[])
{
    int status = LTFATERR_SUCCESS;
    rtisila_state* p = NULL;
    int N = L / a;
    int M2 = M / 2 + 1;

    CHECKNULL(s); CHECKNULL(g); CHECKNULL(c);
    // Just limit lookahead to something sensible
    lookahead = lookahead > N ? N : lookahead;

    CHECKSTATUS(
        rtisila_init(g, gl, 1, a, M, lookahead, maxit, &p),
        "rtisila init failed");

    for (int w = 0; w < W; w++)
    {
        rtisila_reset(p);

        memcpy(p->s + M2 + w * N * M2, s + w * N * M2, lookahead * M2 * sizeof * p->s);

        for (int n = 0, nahead = lookahead; nahead < N; ++n, ++nahead)
        {
            const LTFAT_REAL* sncol = s + nahead * M2 + w * N * M2;
            LTFAT_COMPLEX* cncol = c + n * M2 + w * N * M2;
            rtisila_execute(p, sncol, cncol);
        }

        for (int n = N - lookahead, nahead = 0; n < N; ++n, ++nahead)
        {
            const LTFAT_REAL* sncol = s + nahead * M2 + w * N * M2;
            LTFAT_COMPLEX* cncol = c + n * M2 + w * N * M2;
            rtisila_execute(p, sncol, cncol);
        }
    }
error:
    if (p) rtisila_done(&p);
    return status;
}
