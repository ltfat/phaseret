#include <stdlib.h>
#include <stdalign.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif


#if ! ( defined(mex_h) || defined(MEX_H) )
typedef int mwSignedIndex;
#endif


#define ALIGNBYTES 32 /* 4 doubles */

typedef enum
{
    MOD_STEPWISE = (1 << 0),  // << DEFAULT
    MOD_FRAMEWISE = (1 << 1),
    MOD_COEFFICIENTWISE = (1 << 2)
} leglaupdate_mod;

typedef enum
{
    ORDER_FWD = (1 << 10), // << DEFAULT
    ORDER_REV = (1 << 11)
} leglaupdate_frameorder;

typedef enum
{
    EXT_BOTH = (1 << 20), // << DEFAULT
    EXT_UPDOWN = (1 << 21)
} leglaupdate_ext;

typedef struct
{
    mwSignedIndex kernh;
    mwSignedIndex kernw;
    mwSignedIndex kernwskip;
    mwSignedIndex kernh2;
    mwSignedIndex kernw2;
    mwSignedIndex M;
    int flags;
} leglaupdate_plan_col;

typedef struct
{
    double** kr;
    double** ki;
    double* bufr;
    double* bufi;
    double*  s;
    mwSignedIndex a;
    mwSignedIndex N;
    mwSignedIndex kernelNo;
    leglaupdate_plan_col plan_col;
} leglaupdate_plan;



/* Planning */
leglaupdate_plan*
leglaupdate_init(double* s,mwSignedIndex a, mwSignedIndex M,
        mwSignedIndex N, double* kernr, double* kerni,
        mwSignedIndex kernh, mwSignedIndex kernw,
        int flags);

leglaupdate_plan_col
leglaupdate_init_col(mwSignedIndex M,
        mwSignedIndex kernh, mwSignedIndex kernw,
        int flags);

/* Executing */
void
leglaupdatereal_execute(leglaupdate_plan* plan,  double* cr, double* ci, double* coutr, double* couti);

void
leglaupdatereal_execute_col(leglaupdate_plan_col* plan,
        double* crColFirst, double* ciColFirst,
        double* actKr, double* actKi,
        double* sCol,
        double* coutrCol, double* coutiCol);

/* Cleanup */
void
leglaupdate_done(leglaupdate_plan* plan);

/*  */
void
extendborders(leglaupdate_plan_col* plan, const double* cr, const double* ci, mwSignedIndex N, double* bufr, double* bufi);

/* Modulate kernel */
void
kernphasefi(double* kernr, double* kerni, mwSignedIndex kernh, mwSignedIndex kernw,
        mwSignedIndex kernwskip,
        mwSignedIndex n, mwSignedIndex a, mwSignedIndex M,
        double* kernmodr, double* kernmodi);

/* Format kernel */
void
formatkernel(double* kernr, double* kerni,
        mwSignedIndex kernh, mwSignedIndex kernw,
        mwSignedIndex kernwskip, double* kernmodr, double* kernmodi);

/* Util */
mwSignedIndex gcd(mwSignedIndex m, mwSignedIndex n);
mwSignedIndex lcm(mwSignedIndex m, mwSignedIndex n);

