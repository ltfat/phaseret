#ifdef __STDC_UTF_16__
#include <uchar.h>
#endif
#include "mex.h"
#include "legla.h"

/*
   Calling convention
                            0         1 2 3 4                 5
   c = comp_leglaupdatereal(c,kernsmall,s,a,M,flags.do_onthefly)
 * */

    void
mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

    double *cr,*ci,*coutr,*couti,*s,*kernr,*kerni;
    mwSignedIndex M,a,M2,N,kernh,kernw;
    int do_onthefly = 0;
    a = (mwSignedIndex) mxGetScalar(prhs[3]);
    M = (mwSignedIndex) mxGetScalar(prhs[4]);
    M2 = M/2 + 1;
    N  = mxGetN(prhs[0]);
    cr = mxGetPr(prhs[0]);
    ci = mxGetPi(prhs[0]);

    if(nrhs>5)
    {
        do_onthefly = (int) mxGetScalar(prhs[5]);
    }

    int isreal = !mxIsComplex(prhs[0]);

    if(isreal)
    {
        /* This is real only input. Create an empty array so that the code does not
         * explode */
        ci = mxCalloc(M2*N,sizeof*ci);
    }

    kernr = mxGetPr(prhs[1]);
    kerni = mxGetPi(prhs[1]);
    kernh = mxGetM(prhs[1]);
    kernw = mxGetN(prhs[1]);
    s = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(M2,N,mxCOMPLEX);
    coutr = mxGetPr(plhs[0]);
    couti = mxGetPi(plhs[0]);

    leglaupdate_plan* plan = leglaupdate_init(s,a,M,N,kernr,kerni,kernh,kernw,
                                              do_onthefly?MOD_COEFFICIENTWISE:0);

    leglaupdatereal_execute(plan,cr,ci,coutr,couti);

    leglaupdate_done(plan);


    if(isreal) mxFree(ci);
}
