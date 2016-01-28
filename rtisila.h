#ifndef _rtisila_h
#define _rtisila_h
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>

#define EXTERN

typedef struct
{
    double* frame;
    double* shiftframe;
    fftw_complex* fftframe;
    fftw_plan fwdp;
    fftw_plan backp;
    const double* g;
    const double* gd;
    const double* specg1;
    const double* specg2;
    const int M;
    const int a;
} rtisilaupdate_plan;

// in might be equal to out
void
fftshift(const double* in, int L, double* out);

// in might be equal to out
void
ifftshift(const double* in, int L, double* out);

void
rtisilaoverlaynthframe(rtisilaupdate_plan* p, const double* cframes, const double* gnum, int n, int N);

void
rtisilaphaseupdate(rtisilaupdate_plan* p, const double* sframe, double* frameupd);



rtisilaupdate_plan*
rtisilaupdate_init(const double *g, const double* specg1, const double* specg2, const double* gd,
                   int a, int M);

void
rtisilaupdate_done(rtisilaupdate_plan* p);

void
rtisilaupdate_execute(rtisilaupdate_plan* p, const double* cframes, int N,
                      const double* sframes, int lookahead, int maxit, double* cframes2);

void
rtisilaupdate(const double* cframes,
              const double* gnum, const double* specg1, const double* specg2, const double* dgnum,
              int a, int M, int N, const double* sframes, int lookahead, int maxit,
              double* cframes2);
#endif
