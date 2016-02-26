#include "mex_helper.h"

void
mexFunction(int nlhs, mxArray* plhs[],
            int nrhs, const mxArray* prhs[])
{
    UNUSED(nlhs);UNUSED(nrhs);UNUSED(prhs);
    plhs[0] = mxCreateString("mexok");
}
