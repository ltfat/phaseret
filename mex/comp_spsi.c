/** \addtogroup mex
 *  \{
 *  \file
 *  \brief MEX file interface for spsi()
 *
 *  **This is a computaional routine. The input arguments are not checked for correctness!!**
 *
 *  Matlab calling convention:
 *  --------------------------
 *
 *      [c,endphase] = comp_spsi(s,a,M[,startphase])
 *
 *  Input arg.   | Description
 *  ------------ | -------------------------------------------------------------
 *  s            | M2 x N real matrix, tarhet magnitude
 *  a            | Hop factor
 *  M            | Number of channels (FFT length)
 *  startphase   | (optional) Length M2 vector, phase of -1 col
 *
 *  Output arg.  | Description
 *  -----------  | -------------------------------------------------------------
 *  c            | M2 x N complex coefficients
 *  endphase     | length M2 vector phase of N-1 col
 *
 *  \author Zdeněk Průša
 *  \date 17.02.2016
 *  \}
 */


#include "mex_helper.h"
#include "spsi.h"

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nrhs);
    double* s, *cr, *ci, *initphase = (void*) 0;
    mwSignedIndex N, M, a, M2;
    a = (mwSignedIndex) mxGetScalar(prhs[1]);
    M = (mwSignedIndex) mxGetScalar(prhs[2]);
    M2 = M / 2 + 1;
    N  = mxGetN(prhs[0]);
    s = mxGetData(prhs[0]);

    if (nrhs > 3)
        initphase = mxGetData(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix(M2, N, mxCOMPLEX);
    cr = mxGetData(plhs[0]);
    ci = mxGetImagData(plhs[0]);

    double complex* cc = mxMalloc( M2 * N * sizeof * cc);

    spsi(s, a, M, N, initphase, cc);

    complex2split(cc, M2 * N, cr, ci);

    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(M2, 1, mxREAL);
        double* endphase = mxGetData(plhs[1]);
        memcpy(endphase, initphase, M2 * sizeof * initphase);
    }

    mxFree(cc);
}
