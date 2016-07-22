/** \addtogroup legla
 *  @{
 */
#include "mex_helper.h"
#include "legla.h"

/** MEX interface for LeRoux's version of GLA
*
* Matlab calling convention
*
*     c = comp_leglaupdatereal(c,kernsmall,s,a,M,flags.do_onthefly)
*
*/

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nlhs);
    double* cr, *ci, *coutr, *couti, *s, *kernr, *kerni;
    mwSignedIndex M, a, M2, N, kernh2, kernw;
    int do_onthefly = 0;
    a = (mwSignedIndex) mxGetScalar(prhs[3]);
    M = (mwSignedIndex) mxGetScalar(prhs[4]);
    M2 = M / 2 + 1;
    N  = mxGetN(prhs[0]);
    cr = mxGetPr(prhs[0]);
    ci = mxGetPi(prhs[0]);
    int L = N * a;

    if (nrhs > 5)
        do_onthefly = (int) mxGetScalar(prhs[5]);

    int isreal = !mxIsComplex(prhs[0]);

    /* This is real only input. Create an empty array so that the code does not
     * explode */
    if (isreal)
        ci = mxCalloc(M2 * N, sizeof * ci);

    kernr = mxGetPr(prhs[1]);
    kerni = mxGetPi(prhs[1]);
    kernh2 = mxGetM(prhs[1]);
    kernw = mxGetN(prhs[1]);

    int iskernreal = !mxIsComplex(prhs[1]);
    if (iskernreal)
        kerni = mxCalloc(kernh2 * kernw, sizeof * kerni);

    s = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(M2, N, mxCOMPLEX);
    coutr = mxGetPr(plhs[0]);
    couti = mxGetPi(plhs[0]);

    complex double* kern = mxMalloc(kernh2 * kernw * sizeof * kern);
    complex double* c = mxMalloc(M2 * N * sizeof * c);
    complex double* cout = mxMalloc(M2 * N * sizeof * cout);

    split2complex(cr, ci, M2 * N, c);
    split2complex(kernr, kerni, kernh2 * kernw, kern);

    leglaupdate_plan* plan = NULL;
    leglaupdate_init(kern, (phaseret_size) {.width = kernw, .height = 2 * kernh2 - 1},
    L, 1, a, M, do_onthefly ? MOD_COEFFICIENTWISE : 0, &plan);

    leglaupdate_execute(plan, s, c, cout);

    leglaupdate_done(&plan);

    complex2split(cout, M2 * N, coutr, couti);

    mxFree(cout);
    mxFree(c);
    mxFree(kern);
    if (isreal) mxFree(ci);
    if (iskernreal) mxFree(kerni);
}
// @} */
