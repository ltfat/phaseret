#ifdef __STDC_UTF_16__
#include <uchar.h>
#endif
#include "mex.h"
#include "rtisila.h"

/*
   Calling convention
                                       0    1       2      3     4 5 6       7         8     9
   cframes2 = comp_rtisilaupdate(cframes,gnum,scpecg1,specg2,dgnum,a,M,sframes,lookahead,maxit)
 * */

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    const mxArray* mxcframes = prhs[0];

    int N = mxGetDimensions(mxcframes)[1];
    int W = 1;

    if (mxGetNumberOfDimensions(mxcframes) > 2)
        W = mxGetDimensions(mxcframes)[2];

    const double* cframes = mxGetData(mxcframes);
    const double* gnum = mxGetData(prhs[1]);
    const double* specg1 = mxGetData(prhs[2]);
    const double* specg2 = mxGetData(prhs[3]);
    const double* dgnum = mxGetData(prhs[4]);
    int a = (int) mxGetScalar(prhs[5]);
    int M = (int) mxGetScalar(prhs[6]);
    const double* sframes = mxGetData(prhs[7]);
    int lookahead = (int) mxGetScalar(prhs[8]);
    int maxit = (int) mxGetScalar(prhs[9]);
    int M2 = M / 2 + 1;

    plhs[0] = mxDuplicateArray(mxcframes);
    double* cframes2 = mxGetData(plhs[0]);

    double* cr = NULL;
    double* ci = NULL;

    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(M2, W, mxCOMPLEX);
        cr = mxGetData(plhs[1]);
        ci = mxGetImagData(plhs[1]);
    }

    rtisilaupdate_plan* p = rtisilaupdate_init(gnum, specg1, specg2, dgnum, a, M);

    for (int w = 0; w < W; w++)
    {
        rtisilaupdate_execute(p, cframes + w * M * N, N,
                              sframes + w * M2 * N,
                              lookahead, maxit, cframes2 + w * M * N);

        if (nlhs > 1)
        {
            double complex* cc = (double complex*) p->fftframe;
            double* crChan = cr + w * M2;
            double* ciChan = ci + w * M2;

            for (int m = 0; m < M2; m++)
            {
                crChan[m] = creal(cc[m]);
                ciChan[m] = cimag(cc[m]);
            }
        }
    }


    rtisilaupdate_done(p);

}
