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
    double* a = mxGetData(prhs[1]);
    double* fc = mxGetData(prhs[2]);
    double* tgrad = mxGetData(prhs[3]);
    double* fgrad = mxGetData(prhs[4]);
    double* prevphase = mxGetData(prhs[5]);
    double* neighs = mxGetData(prhs[6]);
    double* posInfo = mxGetData(prhs[7]);
    double tol = mxGetScalar(prhs[8]);
    mwSignedIndex M = (mwSignedIndex) mxGetScalar(prhs[9]);
    

    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
    double* phase = mxGetData(plhs[0]);
    

     phaseret_rtpghifbupdate_plan_d* plan = NULL;
     
     /*
     for (int ii = 0; ii < M; ii++)
     {
        printf("a: %f\n ", a[ii]);
        printf("neighs: %f\n", neighs[ii]);
        printf("pos: %f\n", posInfo[ii]);
     }
      */
     
     phaseret_rtpghifbupdate_init_d(M,1, fc, tol,&plan);
     phaseret_rtpghifbupdate_execute_d(plan, slog, fc, tgrad, fgrad, prevphase, phase);
     phaseret_rtpghifbupdate_done_d(&plan);
}

