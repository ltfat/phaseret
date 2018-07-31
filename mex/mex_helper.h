#ifndef _mex_helper_h
#define _mex_helper_h

#include <complex.h>
// #ifdef __STDC_UTF_16__
// #include <uchar.h>
// #endif
#include "mex.h"

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

void
complex2split(const complex double* in, int L, double* outr, double* outi)
{
    if(outi)
    {
        for (int ii = 0; ii < L; ++ii)
        {
            outr[ii] = creal(in[ii]);
            outi[ii] = cimag(in[ii]);
        }
    }
    else
    {   
        for (int ii = 0; ii < L; ++ii)
        {
            outr[ii] = creal(in[ii]);
        }  
    }

}

void
split2complex(const double* inr, const double* ini, int L, complex double* out)
{
    
    if(ini)
    {
        for (int ii = 0; ii < L; ++ii)
        {
            out[ii] = inr[ii] + I*ini[ii];
        }
    }
    else
    {
        for (int ii = 0; ii < L; ++ii)
        {
            double* out2 = (double*) out;
            out2[2*ii] = inr[ii];
            out2[2*ii+1] = 0.0;
        }       
    }

}

#endif
