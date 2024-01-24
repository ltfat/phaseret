/** \addtogroup mex
 *  \{
 *  \file
 *  \brief MEX file interface for rtpghifb()
 *
 *  **This is a computational routine. The input arguments are not checked for correctness!!**
 *
 *  Matlab calling convention:
 *  --------------------------
 *
 *      [phase] = comp_rtpghifbupdate(slog,fc, tgrad, fgrad, prevphase, tol, M)
 *
 *  Input arg.   | Description
 *  ------------ | -------------------------------------------------------------
 *  slog         | M x 2 real matrix, tarhet magnitude
 *  fc           | M x 1 real
 *  tgrad        | M x 2 real
 *  fgrad        | M x 1 real
 *  prevphase    | M x 1 real
 *  tol          | scalar 
 *  M            | Number of channels (FFT length)
 *
 *  Output arg.  | Description
 *  -----------  | -------------------------------------------------------------
 *  phase        | M x 1 real
 *
 *  \author Clara Hollomey
 *  \date 12.05.2023
 *  \}
 */

#include "mex_helper.h"
#define LTFAT_DOUBLE
#include "phaseret/rtpghifb.h"
#include "../libltfat/modules/libphaseret/src/rtpghifb_private.h"

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nrhs); UNUSED(nlhs);
    double* slog = mxGetData(prhs[0]);
    double* fc = mxGetData(prhs[1]);
    double* tgrad = mxGetData(prhs[2]);
    double* fgrad = mxGetData(prhs[3]);
    double* prevphase = mxGetData(prhs[4]);
    double tol = mxGetScalar(prhs[5]);
    mwSignedIndex M = (mwSignedIndex) mxGetScalar(prhs[6]);
    

    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
    double* phase = mxGetData(plhs[0]);
    

     phaseret_rtpghifbupdate_plan_d* plan = NULL;
     phaseret_rtpghifbupdate_init_d(M,1, fc, tol,&plan);
     phaseret_rtpghifbupdate_execute_d(plan, slog, fc, tgrad, fgrad, prevphase, phase);
     phaseret_rtpghifbupdate_done_d(&plan);
}

