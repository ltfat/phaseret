#include "phaseret/dgtrealwrapper.h"
#include "ltfat/macros.h"
#include "dgtrealwrapper_private.h"

int
phaseret_dgtreal_params_defaults(phaseret_dgtreal_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->ptype = LTFAT_FREQINV;
    params->fftw_flags = FFTW_ESTIMATE;
    params->hint = phaseret_dgtreal_auto;
error:
    return status;
}

PHASERET_API phaseret_dgtreal_params*
phaseret_dgtreal_params_allocdef()
{
    phaseret_dgtreal_params* params =
        (phaseret_dgtreal_params*) ltfat_calloc(1, sizeof * params);

    phaseret_dgtreal_params_defaults(params);

    return params;
}

PHASERET_API int
phaseret_dgtreal_params_set_phaseconv(phaseret_dgtreal_params* params,
                                      ltfat_phaseconvention ptype)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->ptype = ptype;
error:
    return status;
}

PHASERET_API int
phaseret_dgtreal_params_set_fftwflags(phaseret_dgtreal_params* params,
                                      unsigned fftw_flags)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->fftw_flags = fftw_flags;
error:
    return status;

}

PHASERET_API int
phaseret_dgtreal_params_set_hint(phaseret_dgtreal_params* params,
                                 phaseret_dgtreal_hint hint)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->hint = hint;
error:
    return status;
}

PHASERET_API int
phaseret_dgtreal_params_free(phaseret_dgtreal_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    ltfat_free(params);
error:
    return status;
}
