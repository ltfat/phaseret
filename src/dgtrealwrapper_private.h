typedef int complextorealtransform(void* userdata, const LTFAT_COMPLEX* c, int L, int W, LTFAT_REAL* f);
typedef int realtocomplextransform(void* userdata, const LTFAT_REAL* f, int L, int W, LTFAT_COMPLEX* c);
typedef int donefunc(void** pla);

struct dgtreal_anasyn_plan
{
    int L;
    int W;
    int a;
    int M;
    LTFAT_REAL* f;
    LTFAT_COMPLEX* c;
    complextorealtransform* backtra;
    void* backtra_userdata;
    donefunc* backdonefunc;
    realtocomplextransform* fwdtra;
    void* fwdtra_userdata;
    donefunc* fwddonefunc;
};


int
ltfat_idgtreal_long_execute_d_wrapper(void* plan, const LTFAT_COMPLEX* c, int L,
                                      int W, LTFAT_REAL* f);

int
ltfat_dgtreal_long_execute_d_wrapper(void* plan, const LTFAT_REAL* f, int L, int W,
                                     LTFAT_COMPLEX* c);

int
ltfat_idgtreal_fb_execute_d_wrapper(void* plan, const LTFAT_COMPLEX* c, int L,
                                    int W, LTFAT_REAL* f);

int
ltfat_dgtreal_fb_execute_d_wrapper(void* plan, const LTFAT_REAL* f, int L, int W,
                                   LTFAT_COMPLEX* c);

int
ltfat_idgtreal_long_done_d_wrapper(void** plan);

int
ltfat_dgtreal_long_done_d_wrapper(void** plan);

int
ltfat_idgtreal_fb_done_d_wrapper(void** plan);

int
ltfat_dgtreal_fb_done_d_wrapper(void** plan);
