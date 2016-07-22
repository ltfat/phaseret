#include "config.h"
#include "dgtrealwrapper.h"


typedef struct legla_plan legla_plan;
typedef struct leglaupdate_plan leglaupdate_plan;
typedef struct leglaupdate_plan_col leglaupdate_plan_col;

/** \addtogroup legla
 *  @{
 */
typedef int legla_callback_cmod(void* userdata, complex double* c, int L, int W, int a, int M);
typedef int legla_callback_status(dgtreal_anasyn_plan* p, void* userdata, complex double* c,
                                  int L, int W, int a, int M, double* alpha, int iter);

typedef struct
{
    int height;
    int width;
} phaseret_size;

typedef enum
{
    MOD_STEPWISE = (1 << 0),  // << DEFAULT
    MOD_FRAMEWISE = (1 << 1),
    MOD_COEFFICIENTWISE = (1 << 2),
    MOD_MODIFIEDUPDATE = (1 << 3)
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

/* Top level */
int
legla(const complex double cinit[], const double g[], const int gl, const int L,
      const int W, const int a, const int M, const int iter, complex double cout[]);

/* High level level */
int
legla_init(const double g[], const int gl, const int L, const int W,
           const int a, const int M, phaseret_size ksize, const double alpha,
           complex double c[], dgtreal_anasyn_hint hint, unsigned dgtrealflags,
           unsigned leglaflags, legla_plan** pout);

int
legla_execute(legla_plan* p, const complex double cinit[], const int iter);


int
legla_execute_newarray(legla_plan* p, const complex double cinit[], const int iter, complex double cout[]);

int
legla_done(legla_plan** p);

int
legla_set_status_callback(legla_plan* p, legla_callback_status* callback, void* userdata);

int
legla_set_cmod_callback(legla_plan* p, legla_callback_cmod* callback, void* userdata);

/** @}*/

/* Single iteration  */
int
leglaupdate_init(const complex double kern[], phaseret_size ksize,
                 int L, int W, int a, int M, int flags, leglaupdate_plan** pout);

extern void
leglaupdate_execute(leglaupdate_plan* plan, const double s[], complex double c[],
                    complex double cout[]);

void
leglaupdate_done(leglaupdate_plan** plan);

/* Single col update */
int
leglaupdate_init_col(int M, phaseret_size ksize,
                     int flags, leglaupdate_plan_col** pout);

void
leglaupdatereal_execute_col(leglaupdate_plan_col* plan,
                            const double sCol[],
                            const complex double actK[],
                            complex double cColFirst[],
                            complex double coutrCol[]);

/* Utils */
void
extendborders(leglaupdate_plan_col* plan, const complex double c[], int N,
              complex double buf[]);

int
legla_big2small_kernel(complex double* bigc, phaseret_size bigsize,
                       phaseret_size smallsize, complex double* smallc);

/* Modulate kernel */
void
kernphasefi(const complex double kern[], phaseret_size ksize,
            int n, int a, int M, complex double kernmod[]);

/* Format kernel */
void
formatkernel(double* kernr, double* kerni,
             int kernh, int kernw,
             int kernwskip, double* kernmodr, double* kernmodi);

/* Util */
int phaseret_gcd(int m, int n);
int phaseret_lcm(int m, int n);
