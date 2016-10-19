/** \addtogroup mex
 *  \{
 *  \file
 *  \brief MEX file interface for rtisilaupdate()
 *
 *  **This is a computaional routine. The input arguments are not checked for correctness!!**
 *
 *  Matlab calling convention:
 *  --------------------------
 *
 *      [frames2,cframes2,c] = comp_rtisilaupdate(frames,cframes,g,gd,a,M,s,lookahead,maxit)
 *
 *  Input arg.   | Description
 *  ------------ | -------------------------------------------------------------
 *  frames       | M x N real matrix, frame buffer, N = lookback + 1 + lookahead
 *  cframes      | M2 x N real matrix, frame buffer, N = lookback + 1 + lookahead
 *  g            | M el. vector, analysis window
 *  gd           | M el. vector, synthesis window
 *  a            | Hop factor
 *  M            | Number of channels (FFT length)
 *  s            | M2 x lookahead+1 matrix
 *
 *  Output arg.  | Description
 *  -----------  | -------------------------------------------------------------
 *  frames2      | frames with updated last lookahead + 1 cols.
 *  cframes2     | cframes with updated last lookahead + 1 cols.
 *  c            | Complex coefficients of the submit frame.
 *
 *  \author Zdeněk Průša
 *  \date 17.02.2016
 *  \}
 */

#include "mex_helper.h"
#define LTFAT_DOUBLE
#include "phaseret/gsrtisila.h"

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nrhs);
    const mxArray* mxcframes = prhs[0];
    const mxArray* mxcoefbuf = prhs[1];

    int N = mxGetDimensions(mxcframes)[1];
    int W = 1;

    if (mxGetNumberOfDimensions(mxcframes) > 2)
        W = mxGetDimensions(mxcframes)[2];

    const double* cframes = mxGetData(mxcframes);
    const double* coefbufr = mxGetData(mxcoefbuf);
    const double* coefbufi = mxGetImagData(mxcoefbuf);
    const double* gnums = mxGetData(prhs[2]);
    const double* dgnum = mxGetData(prhs[3]);
    int a = (int) mxGetScalar(prhs[4]);
    int M = (int) mxGetScalar(prhs[5]);
    const double* sframes = mxGetData(prhs[6]);
    int lookahead = (int) mxGetScalar(prhs[7]);
    int maxit = (int) mxGetScalar(prhs[8]);
    int M2 = M / 2 + 1;
    int gl = mxGetM(prhs[2]);

    plhs[0] = mxDuplicateArray(mxcframes);
    double* cframes2 = mxGetData(plhs[0]);

    double* cr = NULL;
    double* ci = NULL;

    double* cbufr = NULL;
    double* cbufi = NULL;

    complex double* coefbuf2 = mxMalloc(M2 * N * sizeof * coefbuf2);
    split2complex(coefbufr,coefbufi, M2*N, coefbuf2);

    complex double* cout = mxMalloc(M2 * sizeof * cout);

    if (nlhs > 1)
    {
        plhs[1] = mxDuplicateArray(mxcoefbuf);
        cbufr = mxGetData(plhs[1]);
        cbufi = mxGetImagData(plhs[1]);
    }

    if (nlhs > 2)
    {
        plhs[2] = mxCreateDoubleMatrix(M2, W, mxCOMPLEX);
        cr = mxGetData(plhs[2]);
        ci = mxGetImagData(plhs[2]);
    }

    phaseret_gsrtisilaupdate_plan_d* p = NULL;
    phaseret_gsrtisilaupdate_init_d(gnums, dgnum, gl, a, M, (lookahead + 1), 0, &p);

    for (int w = 0; w < W; w++)
    {
        phaseret_gsrtisilaupdate_execute_d(p, cframes + w * M * N,
                                           coefbuf2 + w * M2 * N,
                                           N,
                                           sframes + w * M2 * N,
                                           lookahead, maxit,
                                           cframes2 + w * M * N,
                                           coefbuf2 + w * M2 * N,
                                           cout);

        if (nlhs > 1)
            complex2split(coefbuf2, N*M2, cbufr + w * N * M2, cbufi + w * N * M2);

        if (nlhs > 2)
            complex2split(cout, M2,  cr + w * M2, ci + w * M2);
    }



    mxFree(coefbuf2);
    mxFree(cout);

    phaseret_gsrtisilaupdate_done_d(&p);
}
