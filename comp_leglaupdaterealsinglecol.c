#ifdef __STDC_UTF_16__
#include <uchar.h>
#endif
#include "mex.h"
#include "leglalib.h"

/*
 Calling convention                   0         1 2 3 
 cout = comp_leglaupdaterealsinglecol(c,kernsmall,s,M)
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
    double *cr,*ci,*coutr,*couti,*s,*kernr,*kerni;
    mwSignedIndex M,M2,N,kernh,kernw;
    M = (mwSignedIndex) mxGetScalar(prhs[3]);
    M2 = M/2 + 1;
    N  = mxGetN(prhs[0]);
    cr = mxGetPr(prhs[0]);
    ci = mxGetPi(prhs[0]);

    if(ci == NULL)
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

    plhs[0] = mxCreateDoubleMatrix(M2,1,mxCOMPLEX);
    coutr = mxGetPr(plhs[0]);
    couti = mxGetPi(plhs[0]);

    leglaupdate_plan_col plan = leglaupdate_init_col(s, M, kernh, kernw, EXT_UPDOWN | MOD_FRAMEWISE );

    double* kr = aligned_alloc(ALIGNBYTES,plan.kernw*plan.kernwskip);
    double* ki = aligned_alloc(ALIGNBYTES,plan.kernw*plan.kernwskip);
    double* bufr = aligned_alloc(ALIGNBYTES,((M2+kernh-1)*N)*sizeof*bufr);
    double* bufi = aligned_alloc(ALIGNBYTES,((M2+kernh-1)*N)*sizeof*bufi);

    extendborders(&plan,cr,ci,N,bufr,bufi);

    formatkernel(kernr,kerni,kernh, kernw, plan.kernwskip, kr, ki);

    leglaupdatereal_execute_col(&plan,0,bufr,bufi,kr,ki,coutr,couti);

    free(kr);
    free(ki);
    free(bufr);
    free(bufi);

    if(ci) mxFree(ci);
}
