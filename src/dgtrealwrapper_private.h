typedef int PHASERET_NAME(complextorealtransform)(void* userdata, const LTFAT_COMPLEX* c, int L, int W, LTFAT_REAL* f);
typedef int PHASERET_NAME(realtocomplextransform)(void* userdata, const LTFAT_REAL* f, int L, int W, LTFAT_COMPLEX* c);
typedef int PHASERET_NAME(donefunc)(void** pla);

struct PHASERET_NAME(dgtreal_plan)
{
    int L;
    int W;
    int a;
    int M;
    LTFAT_REAL* f;
    LTFAT_COMPLEX* c;
    PHASERET_NAME(complextorealtransform)* backtra;
    void* backtra_userdata;
    PHASERET_NAME(donefunc)* backdonefunc;
    PHASERET_NAME(realtocomplextransform)* fwdtra;
    void* fwdtra_userdata;
    PHASERET_NAME(donefunc)* fwddonefunc;
};

#ifdef __cplusplus
extern "C" {
#endif

int
PHASERET_NAME(ltfat_idgtreal_long_execute_wrapper)(void* plan, const LTFAT_COMPLEX* c,
        int L, int W, LTFAT_REAL* f);

int
PHASERET_NAME(ltfat_dgtreal_long_execute_wrapper)(void* plan, const LTFAT_REAL* f,
        int L, int W, LTFAT_COMPLEX* c);

int
PHASERET_NAME(ltfat_idgtreal_fb_execute_wrapper)(void* plan, const LTFAT_COMPLEX* c, int L,
        int W, LTFAT_REAL* f);

int
PHASERET_NAME(ltfat_dgtreal_fb_execute_wrapper)(void* plan, const LTFAT_REAL* f, int L, int W,
        LTFAT_COMPLEX* c);

int
PHASERET_NAME(ltfat_idgtreal_long_done_wrapper)(void** plan);

int
PHASERET_NAME(ltfat_dgtreal_long_done_wrapper)(void** plan);

int
PHASERET_NAME(ltfat_idgtreal_fb_done_wrapper)(void** plan);

int
PHASERET_NAME(ltfat_dgtreal_fb_done_wrapper)(void** plan);

#ifdef __cplusplus
}
#endif

