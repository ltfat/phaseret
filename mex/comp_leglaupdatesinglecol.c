#include "mex_helper.h"
#include "legla.h"

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
    mwSignedIndex M, M2, N, kernh, kernw;
    M = (mwSignedIndex) mxGetScalar(prhs[3]);
    M2 = M / 2 + 1;  N  = mxGetN(prhs[0]);
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
    kernh = mxGetM(prhs[1]);
    kernw = mxGetN(prhs[1]);
    int iskernreal = !mxIsComplex(prhs[1]);

    if (iskernreal)
        kerni = mxCalloc(kernw * kernh, sizeof * kerni);

    s = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(M2, 1, mxCOMPLEX);
    coutr = mxGetPr(plhs[0]);
    couti = mxGetPi(plhs[0]);

    leglaupdate_mod modflag = do_onthefly ? MOD_COEFFICIENTWISE : MOD_FRAMEWISE;
    leglaupdate_plan_col* plan = NULL;

    leglaupdate_init_col( M, (phaseret_size) {.width = kernw, .height = kernh},
    EXT_UPDOWN | modflag , &plan);

    complex double* kern = mxMalloc(kernw * kernh * sizeof * kern);
    complex double* c = mxMalloc(M2 * N * sizeof * c);
    split2complex(cr, ci, M2 * N, c);
    split2complex(kernr, kerni, M2 * N, kern);

    complex double* cout = mxMalloc(M2 * N * sizeof * cout);
    complex double* buf = mxMalloc(((M2 + kernh - 1) * N) * sizeof * buf);
    extendborders(plan, c, N, buf);

    leglaupdatereal_execute_col(plan, s, kern, buf, cout);

    complex2split(cout, M2 * N, coutr, couti);

    mxFree(kern);
    mxFree(buf);
    mxFree(c);
    mxFree(cout);

    if (isreal) mxFree(ci);
    if (iskernreal) mxFree(kerni);
}
