#ifndef _PHASERET_RTPGHIFB_PRIVATE_H
#define _PHASERET_RTPGHIFB_PRIVATE_H


struct PHASERET_NAME(rtpghifb_state)
{
    PHASERET_NAME(rtpghifbupdate_plan)* p;
    ltfat_int M;
    ltfat_int N;
    ltfat_int a;
    ltfat_int W;
    int do_causal;
    LTFAT_REAL* slog;
    LTFAT_REAL* fc;
    LTFAT_REAL* s;
    LTFAT_REAL* tgrad; //!< Time gradient buffer
    LTFAT_REAL* fgrad; //!< Frequency gradient buffer
    LTFAT_REAL* phase; //!< Buffer for keeping previously computed frame
    double gamma;
};

struct PHASERET_NAME(rtpghifbupdate_plan)
{
    LTFAT_NAME(heap)* h;
    int* donemask;
    double logtol;
    double tol;
    ltfat_int M;
    ltfat_int N;
    LTFAT_REAL* randphase; //!< Precomputed array of random phase
    ltfat_int randphaseLen;
    ltfat_int randphaseId;
    LTFAT_REAL* fc;
};

#endif