#include "mex_helper.h"
#define LTFAT_DOUBLE
#include "phaseret/legla.h"

/*
 Calling convention                   0         1 2 3           4
 cout = comp_leglaupdaterealsinglecol(c,kernsmall,s,M,do_onthefly)
 *
 * kernsmall size is kernh x kernw
 * c size is M2 x kernw
 * s size is M2 x 1
 * cout size M2 x 1
 * */



void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nlhs); UNUSED(nrhs);
    double* cr, *ci, *coutr, *couti, *s, *kernr, *kerni;
    mwSignedIndex M, M2, N, kernh2, kernh, kernw;
    M = (mwSignedIndex) mxGetScalar(prhs[3]);
    M2 = M / 2 + 1;
    N  = mxGetN(prhs[0]);
    cr = mxGetPr(prhs[0]);
    ci = mxGetPi(prhs[0]);

    int isreal = !mxIsComplex(prhs[0]);

    /* This is real only input. Create an empty array so that the code does not
     * explode */
    if (isreal)
        ci = mxCalloc(M2 * N, sizeof * ci);

    int do_onthefly = (int)mxGetScalar(prhs[4]);

    kernr = mxGetPr(prhs[1]);
    kerni = mxGetPi(prhs[1]);
    kernh2 = mxGetM(prhs[1]);
    kernh = 2 * kernh2 - 1;
    kernw = mxGetN(prhs[1]);
    int iskernreal = !mxIsComplex(prhs[1]);

    if (iskernreal)
        kerni = mxCalloc(kernw * kernh2, sizeof * kerni);

    s = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(M2, 1, mxCOMPLEX);
    coutr = mxGetPr(plhs[0]);
    couti = mxGetPi(plhs[0]);

    leglaupdate_mod modflag = do_onthefly ? MOD_COEFFICIENTWISE : MOD_FRAMEWISE;
    phaseret_leglaupdate_plan_col_d* plan = NULL;

    phaseret_leglaupdate_col_init_d( M, (phaseret_size) {.width = kernw, .height = kernh},
    EXT_UPDOWN | modflag , &plan);

    complex double* kern = mxMalloc(kernw * kernh2 * sizeof * kern);
    complex double* c = mxMalloc(M2 * N * sizeof * c);
    split2complex(cr, ci, M2 * N, c);
    split2complex(kernr, kerni, kernh2 * kernw, kern);

    complex double* cout = mxMalloc(M2 * sizeof * cout);
    complex double* buf = mxMalloc(((M2 + kernh - 1) * N) * sizeof * buf);
    phaseret_extendborders_d(plan, c, N, buf);

    phaseret_leglaupdate_col_execute_d(plan, s, kern, buf, cout);

    complex2split(cout, M2, coutr, couti);

    mxFree(kern);
    mxFree(buf);
    mxFree(c);
    mxFree(cout);
    phaseret_leglaupdate_col_done_d(&plan);

    if (isreal) mxFree(ci);
    if (iskernreal) mxFree(kerni);
}
