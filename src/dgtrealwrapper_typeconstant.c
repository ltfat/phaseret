#include "phaseret/dgtrealwrapper.h"
#include "ltfat/macros.h"


PHASERET_API int
phaseret_dgtreal_init_params_defaults(phaseret_dgtreal_init_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->ptype = LTFAT_FREQINV;
    params->fftw_flags = FFTW_ESTIMATE;
    params->hint = phaseret_dgtreal_auto;
error:
    return status;
}
