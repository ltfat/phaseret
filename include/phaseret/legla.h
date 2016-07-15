/** \addtogroup legla
 *  @{
 */
#include "config.h"

// #if ! ( defined(mex_h) || defined(MEX_H) )
// typedef int mwSignedIndex;
// #endif

typedef struct legla_plan legla_plan;

typedef struct
{
    int height;
    int width;
} phaseret_size;

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
    int kernh;
    int kernw;
    int kernwskip;
    int kernh2;
    int kernw2;
    int M;
    int flags;
} leglaupdate_plan_col;

typedef struct
{
    double** kr;
    double** ki;
    double* bufr;
    double* bufi;
    double*  s;
    int a;
    int N;
    int kernelNo;
    leglaupdate_plan_col plan_col;
} leglaupdate_plan;

/* Planning */
leglaupdate_plan*
leglaupdate_init( int a, int M,
                  int N, double* kernr, double* kerni,
                  int kernh, int kernw,
                  int flags);

leglaupdate_plan_col
leglaupdate_init_col(int M,
                     int kernh, int kernw,
                     int flags);

/* Executing */
extern void
leglaupdatereal_execute(leglaupdate_plan* plan, const double* s, double* cr,
                        double* ci, double* coutr, double* couti);

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
extendborders(leglaupdate_plan_col* plan, const double* cr, const double* ci, int N, double* bufr, double* bufi);

/* Modulate kernel */
void
kernphasefi(double* kernr, double* kerni, int kernh, int kernw,
            int kernwskip,
            int n, int a, int M,
            double* kernmodr, double* kernmodi);

/* Format kernel */
void
formatkernel(double* kernr, double* kerni,
             int kernh, int kernw,
             int kernwskip, double* kernmodr, double* kernmodi);

/* Util */
int phaseret_gcd(int m, int n);
int phaseret_lcm(int m, int n);

/** @}*/
