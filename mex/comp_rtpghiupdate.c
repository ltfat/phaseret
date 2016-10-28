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
#include "phaseret/rtpghi.h"

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nrhs); UNUSED(nlhs);
    double* slog = mxGetData(prhs[0]);
    double* tgrad = mxGetData(prhs[1]);
    double* fgrad = mxGetData(prhs[2]);
    double* prevphase = mxGetData(prhs[3]);
    double tol = mxGetScalar(prhs[4]);
    mwSignedIndex M = (mwSignedIndex) mxGetScalar(prhs[5]);
    mwSignedIndex M2 = M / 2 + 1;

    plhs[0] = mxCreateDoubleMatrix(M2, 1, mxREAL);
    double* phase = mxGetData(plhs[0]);

     phaseret_rtpghiupdate_plan_d* plan = NULL;
     phaseret_rtpghiupdate_init_d(M,1,tol,&plan);
     phaseret_rtpghiupdate_execute_d(plan, slog, tgrad, fgrad, prevphase, phase);
     phaseret_rtpghiupdate_done_d(&plan);
}
