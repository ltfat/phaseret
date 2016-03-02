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
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    UNUSED(nlhs);UNUSED(nrhs);
    double *cr,*ci,*coutr,*couti,*s,*kernr,*kerni;
    mwSignedIndex M,M2,N,kernh,kernw;
    M = (mwSignedIndex) mxGetScalar(prhs[3]);
    M2 = M/2 + 1;  N  = mxGetN(prhs[0]);
    cr = mxGetPr(prhs[0]);
    ci = mxGetPi(prhs[0]);

    int isreal = !mxIsComplex(prhs[0]);

    if(isreal)
    {
        /* This is real only input. Create an empty array so that the code does not
         * explode */
        ci = mxCalloc(M2*N,sizeof*ci);
    }

    int do_onthefly = (int)mxGetScalar(prhs[4]);

    kernr = mxGetPr(prhs[1]);
    kerni = mxGetPi(prhs[1]);
    kernh = mxGetM(prhs[1]);
    kernw = mxGetN(prhs[1]);
    int iskernreal = !mxIsComplex(prhs[1]);

    if(iskernreal)
    {
        /* This is real only input. Create an empty array so that the code does not
         * explode */
        kerni = mxCalloc(kernw*kernh,sizeof*kerni);
    }

    s = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(M2,1,mxCOMPLEX);
    coutr = mxGetPr(plhs[0]);
    couti = mxGetPi(plhs[0]);

    leglaupdate_mod modflag = do_onthefly?MOD_COEFFICIENTWISE:MOD_FRAMEWISE;
    leglaupdate_plan_col plan = leglaupdate_init_col( M, kernh, kernw, EXT_UPDOWN | modflag );

    double* kr = mxMalloc(plan.kernw*plan.kernwskip);
    double* ki = mxMalloc(plan.kernw*plan.kernwskip);
    double* bufr = mxMalloc(((M2+kernh-1)*N)*sizeof*bufr);
    double* bufi = mxMalloc(((M2+kernh-1)*N)*sizeof*bufi);

    extendborders(&plan,cr,ci,N,bufr,bufi);

    formatkernel(kernr,kerni,kernh, kernw, plan.kernwskip, kr, ki);

    leglaupdatereal_execute_col(&plan,bufr,bufi,kr,ki,s,coutr,couti);

    mxFree(kr);
    mxFree(ki);
    mxFree(bufr);
    mxFree(bufi);

    if(isreal) mxFree(ci);
}
