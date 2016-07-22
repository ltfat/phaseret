#ifndef _rtpghi_h
#define _rtpghi_h

// #ifndef NOSYSTEMHEADERS
// #include <ltfat.h>
// #endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pghi_plan pghi_plan;

/** \addtogroup pghi
 * @{
 */
int
pghi(const double s[], double gamma, const int L, const int W,
     const int a, const int M, complex double c[]);

int
pghi_init(double gamma, const int L, const int W,
          const int a, const int M, double tol, pghi_plan** pout);

int
pghi_init_twostep(double gamma, const int L, const int W,
                  const int a, const int M, double tol1, double tol2, pghi_plan** pout);

int
pghi_execute(pghi_plan* p, const double s[], complex double c[]);

int
pghi_done(pghi_plan** p);
/** @} */

void
pghimagphase(const double s[], const double phase[], int L, complex double c[]);

void
pghilog(const double in[], int L, double out[]);

void
pghitgrad(const double logs[], double gamma, int a, int M, int N, double tgrad[]);

void
pghifgrad(const double logs[], double gamma, int a, int M, int N, double fgrad[]);


#ifdef __cplusplus
}
#endif
#endif
