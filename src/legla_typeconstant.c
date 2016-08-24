#include "phaseret/legla.h"
#include "ltfat/macros.h"

PHASERET_API int
phaseret_legla_init_params_defaults(phaseret_legla_init_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);

    params->relthr = 1e-3;
    params->ksize.width = 0;
    params->ksize.height = 0;
    params->leglaflags = MOD_COEFFICIENTWISE | MOD_MODIFIEDUPDATE;
    phaseret_dgtreal_init_params_defaults(&params->dparams);
error:
    return status;
}
