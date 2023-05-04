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
 *  slog         | M2 x 2 real matrix, tarhet magnitude
 *  tgrad        | M2 x 2 real
 *  fgrad        | M2 x 1 real
 *  prevphase    | M2 x 1 real
 *  tol          | scalar 
 *  M            | Number of channels (FFT length)
 *
 *  Output arg.  | Description
 *  -----------  | -------------------------------------------------------------
 *  phase        | M2 x 1 real
 *
 *  \author Zdeněk Průša
 *  \date 17.02.2016
 *  \}
 */

#include "mex_helper.h"
#define LTFAT_DOUBLE
#include "phaseret/rtpghifb.h"
#include "../libltfat/modules/libphaseret/src/rtpghifb_private.h"
#include "ltfat_mex_template_helper.h"
#include "mex.h"

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nrhs); UNUSED(nlhs);
    const mxArray* mxs  = prhs[0];
    double* slog = mxGetData(prhs[0]);
    double* fc = mxGetData(prhs[1]);
    double* tgrad = mxGetData(prhs[2]);
    double* fgrad = mxGetData(prhs[3]);
    double* prevphase = mxGetData(prhs[4]);
    double tol = mxGetScalar(prhs[5]);
    mwSignedIndex M = (mwSignedIndex) mxGetScalar(prhs[6]);
    mwSignedIndex M2 = M / 2 + 1;
    mwSize N = mxGetN(mxs);

    plhs[0] = mxCreateDoubleMatrix(M2, 1, mxREAL);
    double* phase = mxGetData(plhs[0]);
    

     phaseret_rtpghifbupdate_plan_d* plan = NULL;
     phaseret_rtpghifbupdate_init_d(N,M,1, fc, tol,&plan);
     phaseret_rtpghifbupdate_execute_d(plan, slog, fc, tgrad, fgrad, prevphase, phase);
     phaseret_rtpghifbupdate_done_d(&plan);
}

