#include "dgtrealwrapper.h"

/** \addtogroup gla
 * @{
 */

typedef struct gla_plan gla_plan;

typedef int gla_callback_cmod(void* userdata, complex double* c, int L, int W, int a, int M);

typedef int gla_callback_fmod(void* userdata, double* f, int L, int W, int a, int M);

typedef int gla_callback_status(dgtreal_anasyn_plan* p, void* userdata, complex double* c,
                                int L, int W, int a, int M, double* alpha, int iter);

int
gla(const double s[], const double g[], const int gl, const int L, const int W,
    const int a, const int M, const int iter, complex double c[]);

int
gla_init(const double g[], const int gl, const int L, const int W,
         const int a, const int M, const double alpha,
         complex double c[], dgtreal_anasyn_hint hint, unsigned flags,
         gla_plan** pout);

int
gla_execute(gla_plan* p, const double s[], const int iter);

int
gla_execute_newarray(gla_plan* p, const double s[], const int iter,
                     complex double c[]);

int
gla_done(gla_plan** p);

int
gla(const double s[], const double g[], const int gl, const int L, const int W,
        const int a, const int M, const int iter, complex double c[]);

int
gla_set_status_callback(gla_plan* p, gla_callback_status* callback, void* userdata);

int
gla_set_cmod_callback(gla_plan* p, gla_callback_cmod* callback, void* userdata);

int
gla_set_fmod_callback(gla_plan* p, gla_callback_fmod* callback, void* userdata);

/** @} */

int
fastupdate(complex double* c, complex double* t, double alpha, int L);
