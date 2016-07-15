typedef int complextorealtransform(void* userdata, const complex double* c, int L, int W, double* f);
typedef int realtocomplextransform(void* userdata, const double* f, int L, int W, complex double* c);
typedef int donefunc(void** pla);

struct dgtreal_anasyn_plan
{
    int L;
    int W;
    int a;
    int M;
    double* f;
    complex double* c;
    complextorealtransform* backtra;
    void* backtra_userdata;
    donefunc* backdonefunc;
    realtocomplextransform* fwdtra;
    void* fwdtra_userdata;
    donefunc* fwddonefunc;
};


int
ltfat_idgtreal_long_execute_d_wrapper(void* plan, const complex double* c, int L,
                                      int W, double* f);

int
ltfat_dgtreal_long_execute_d_wrapper(void* plan, const double* f, int L, int W,
                                     complex double* c);

int
ltfat_idgtreal_fb_execute_d_wrapper(void* plan, const complex double* c, int L,
                                    int W, double* f);

int
ltfat_dgtreal_fb_execute_d_wrapper(void* plan, const double* f, int L, int W,
                                   complex double* c);

int
ltfat_idgtreal_long_done_d_wrapper(void** plan);

int
ltfat_dgtreal_long_done_d_wrapper(void** plan);

int
ltfat_idgtreal_fb_done_d_wrapper(void** plan);

int
ltfat_dgtreal_fb_done_d_wrapper(void** plan);
