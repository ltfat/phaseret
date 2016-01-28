#include "legla.h"

leglaupdate_plan_col
leglaupdate_init_col( mwSignedIndex M,
                      mwSignedIndex kernh, mwSignedIndex kernw,
                      int flags)
{
    leglaupdate_plan_col retval = (leglaupdate_plan_col)
    {
        .kernh = kernh, .kernw = kernw,
         .M = M, .flags = flags,
          .kernh2 = kernh / 2 + 1, .kernw2 = kernw / 2 + 1
    };

    retval.kernwskip = ((retval.kernh2 * sizeof(double)) / ALIGNBYTES + 1) *
                       ALIGNBYTES;

    // Sanitize flags (set defaults)
    if (retval.flags & (MOD_FRAMEWISE | MOD_COEFFICIENTWISE))
    {
        retval.flags &= ~MOD_STEPWISE; // For safety, clear the default flag.

        // Clear also lower priority flags
        if (retval.flags & MOD_COEFFICIENTWISE)
            retval.flags &= ~MOD_FRAMEWISE;
    }
    else
    {
        retval.flags |=
            MOD_STEPWISE; // Set the default flag if none of the others is set
    }

    if (retval.flags & (ORDER_REV))
        retval.flags &= ~ORDER_FWD;
    else
        retval.flags |= ORDER_FWD;

    if (retval.flags & (EXT_UPDOWN))
        retval.flags &= ~EXT_BOTH;
    else
        retval.flags |= EXT_BOTH;

    return retval;
}


leglaupdate_plan*
leglaupdate_init(double* s, mwSignedIndex a, mwSignedIndex M,
                 mwSignedIndex N, double* kernr, double* kerni,
                 mwSignedIndex kernh, mwSignedIndex kernw,
                 int flags)
{
    leglaupdate_plan* plan = aligned_alloc(ALIGNBYTES, sizeof(leglaupdate_plan));
    plan->plan_col = leglaupdate_init_col( M, kernh, kernw, flags);
    plan->s = s;

    // N
    plan->N = (plan->plan_col.flags & EXT_UPDOWN) ? N - (kernw - 1) : N;

    plan->kernelNo = lcm(M, a) / a;

    int M2 = M / 2 + 1;


    plan->kr = aligned_alloc(ALIGNBYTES, plan->kernelNo * sizeof(double*));
    plan->ki = aligned_alloc(ALIGNBYTES, plan->kernelNo * sizeof(double*));

    plan->bufr = aligned_alloc(ALIGNBYTES,
                               ((M2 + kernh - 1) * (plan->N + kernw - 1)) * sizeof * plan->bufr);
    plan->bufi = aligned_alloc(ALIGNBYTES,
                               ((M2 + kernh - 1) * (plan->N + kernw - 1)) * sizeof * plan->bufi);

    for (int n = 0; n < plan->kernelNo; n++)
    {
        plan->kr[n] = aligned_alloc(ALIGNBYTES,
                                    plan->plan_col.kernw * plan->plan_col.kernwskip);
        plan->ki[n] = aligned_alloc(ALIGNBYTES,
                                    plan->plan_col.kernw * plan->plan_col.kernwskip);
        kernphasefi(kernr, kerni, kernh, kernw, plan->plan_col.kernwskip, n, a, M,
                    plan->kr[n], plan->ki[n]);
    }

    return plan;
}


void
extendborders(leglaupdate_plan_col* plan, const double* cr, const double* ci,
              mwSignedIndex N, double* bufr, double* bufi)
{
    mwSignedIndex m, n;
    mwSignedIndex M2 = plan->M / 2 + 1;
    mwSignedIndex kernh = plan->kernh;
    mwSignedIndex kernw = plan->kernw;
    mwSignedIndex kernw2 = plan->kernw2;
    mwSignedIndex kernh2 = plan->kernh2;
    /* mwSignedIndex kernskipd = plan->kernwskip/(sizeof*cr); */
    mwSignedIndex M2buf = M2 + kernh - 1;
    mwSignedIndex Nbuf = N + kernw - 1;

    if ( plan->flags & EXT_UPDOWN) Nbuf = N;

    if ( !(plan->flags & EXT_UPDOWN))
    {
        /* Copy input to the center of the buffer */
        for (n = 0; n < N; n++)
        {
            double* bufrstart = bufr + (n + kernw2 - 1) * M2buf + kernh2 - 1;
            double* bufistart = bufi + (n + kernw2 - 1) * M2buf + kernh2 - 1;
            const double* crstart = cr + n * M2;
            const double* cistart = ci + n * M2;
            memcpy(bufrstart, crstart, M2 * sizeof * bufrstart);
            memcpy(bufistart, cistart, M2 * sizeof * bufistart);
        }


        /* Periodically extend the left side */
        for (m = 0; m < M2; m++)
        {
            double* bufrtarget = bufr + kernh2 - 1 + m;
            double* bufitarget = bufi + kernh2 - 1 + m;
            double* bufrsource = bufr + (N) * M2buf + kernh2 - 1 + m;
            double* bufisource = bufi + (N) * M2buf + kernh2 - 1 + m;

            for (n = 0; n < kernw2 - 1; n++)
            {
                *bufrtarget = *bufrsource;
                *bufitarget = *bufisource;
                bufrtarget += M2buf;
                bufitarget += M2buf;
                bufrsource += M2buf;
                bufisource += M2buf;
            }
        }

        /* Periodically extend the right side*/
        for (m = 0; m < M2; m++)
        {
            double* bufrsource = bufr + (kernw2 - 1) * M2buf + kernh2 - 1 + m;
            double* bufisource = bufi + (kernw2 - 1) * M2buf + kernh2 - 1 + m;

            double* bufrtarget = bufr + (N + kernw2 - 1) * M2buf + kernh2 - 1 + m;
            double* bufitarget = bufi + (N + kernw2 - 1) * M2buf + kernh2 - 1 + m;

            for (n = 0; n < kernw2 - 1; n++)
            {
                *bufrtarget = *bufrsource;
                *bufitarget = *bufisource;
                bufrtarget += M2buf;
                bufitarget += M2buf;
                bufrsource += M2buf;
                bufisource += M2buf;
            }
        }
    }
    else if (plan->flags & EXT_UPDOWN)
    {
        /* Copy input to the center of the buffer */
        for (n = 0; n < N; n++)
        {
            double* bufrstart = bufr + n * M2buf + kernh2 - 1;
            double* bufistart = bufi + n * M2buf + kernh2 - 1;
            const double* crstart = cr + n * M2;
            const double* cistart = ci + n * M2;
            memcpy(bufrstart, crstart, M2 * sizeof * bufrstart);
            memcpy(bufistart, cistart, M2 * sizeof * bufistart);
        }
    }

    /* Conjugated odd-symmetric extention of the top border*/
    for (n = 0; n < Nbuf; n++)
    {
        double* bufrsource = bufr + n * M2buf + kernh2;
        double* bufisource = bufi + n * M2buf + kernh2;

        double* bufrtarget = bufr + n * M2buf + kernh2 - 2;
        double* bufitarget = bufi + n * M2buf + kernh2 - 2;

        for (m = 0; m < kernh2 - 1; m++)
        {
            *bufrtarget-- = *bufrsource++;
            *bufitarget-- = -*bufisource++;
        }
    }

    /* Conjugated odd or even symmetry extension of the bottom border.
     * Depending whether M is odd or even
     * */
    for (n = 0; n < Nbuf; n++)
    {
        double* bufrsource = bufr + n * M2buf + kernh2 - 1 + M2 - 2 + plan->M % 2;
        double* bufisource = bufi + n * M2buf + kernh2 - 1 + M2 - 2 + plan->M % 2;

        double* bufrtarget = bufr + n * M2buf + kernh2 - 1 + M2;
        double* bufitarget = bufi + n * M2buf + kernh2 - 1 + M2;

        for (m = 0; m < kernh2 - 1; m++)
        {
            *bufrtarget++ = *bufrsource--;
            *bufitarget++ = -*bufisource--;
        }
    }

}


void
leglaupdatereal_execute_col(leglaupdate_plan_col* plan,
                            double* crColFirst, double* ciColFirst,
                            double* actKr, double* actKi,
                            double* sCol,
                            double* coutrCol, double* coutiCol)
{
    mwSignedIndex m, mfirst, mlast;
    mwSignedIndex M2 = plan->M / 2 + 1;
    mwSignedIndex kernh = plan->kernh;
    mwSignedIndex kernw = plan->kernw;
    mwSignedIndex kernw2 = plan->kernw2;
    mwSignedIndex kernh2 = plan->kernh2;
    mwSignedIndex kernskipd = plan->kernwskip / (sizeof * crColFirst);
    mwSignedIndex M2buf = M2 + kernh - 1;
    mwSignedIndex kernhMidId = kernh2 - 1;
    mwSignedIndex kernwMidId = kernw2 - 1;
    /* double* bufr = plan->bufr; */
    /* double* bufi = plan->bufi; */

    /* mwSignedIndex Nbuf = N + kernw -1; */
    int do_onthefly = plan->flags & MOD_COEFFICIENTWISE;
    int do_framewise = plan->flags & MOD_FRAMEWISE;

    /* Outside loop over rows */
    for (m = kernh2 - 1, mfirst = 0, mlast = kernh - 1; mfirst < M2;
         m++, mfirst++, mlast++)
    {
        double accumr = 0.0, accumi = 0.0;

        /* mexPrintf("m-loop: %d,%d,%d\n",mfirst,m,mlast); */
        /* inner loop over all cols of the kernel*/
        for (mwSignedIndex kn = 0; kn < kernw; kn++)
        {
            /* mexPrintf("kn-loop: %d\n",kn); */
            double* actKrCol = actKr + kn * kernskipd;
            double* actKiCol = actKi + kn * kernskipd;

            double* crCol = crColFirst + kn * M2buf;
            double* ciCol = ciColFirst + kn * M2buf;

            /* Inner loop over half of the rows of the kernel excluding the middle row */
            /* Doing the complex conjugated kernel elements simulteneously */
            for (mwSignedIndex km = 0; km < kernh2 - 1; km++)
            {
                double ar  = actKrCol[km];
                double ai  = actKiCol[km];
                double br  = crCol[mfirst + km];
                double bi  = ciCol[mfirst + km];
                double bbr = crCol[mlast - km];
                double bbi = ciCol[mlast - km];
                accumr += ar * (br + bbr) - ai * (bi - bbi);
                accumi += ar * (bi + bbi) + ai * (br - bbr);
            }

            /* The middle row is real*/
            accumr += actKrCol[kernhMidId] * crCol[m];
            accumi += actKrCol[kernhMidId] * ciCol[m];
        }

        coutrCol[mfirst] = accumr;
        coutiCol[mfirst] = accumi;

        /* Update the phase of a coefficient immediatelly */
        if (do_onthefly)
        {
            double arg = atan2(coutiCol[mfirst], coutrCol[mfirst]);
            coutrCol[mfirst] = sCol[mfirst] * cos(arg);
            coutiCol[mfirst] = sCol[mfirst] * sin(arg);
            crColFirst[kernwMidId * M2buf + m] = coutrCol[mfirst];
            ciColFirst[kernwMidId * M2buf + m] = coutiCol[mfirst];
        }

    }

    /* Update the phase of a single column */
    if (do_framewise)
    {
        for (m = kernh2 - 1, mfirst = 0; mfirst < M2; m++, mfirst++)
        {
            double arg = atan2(coutiCol[mfirst], coutrCol[mfirst]);
            coutrCol[mfirst] = sCol[mfirst] * cos(arg);
            coutiCol[mfirst] = sCol[mfirst] * sin(arg);
            crColFirst[kernwMidId * M2buf + m] = coutrCol[mfirst];
            ciColFirst[kernwMidId * M2buf + m] = coutiCol[mfirst];
        }
    }

}

void
leglaupdatereal_execute(leglaupdate_plan* plan,  double* cr, double* ci,
                        double* coutr, double* couti)
{
    mwSignedIndex M2 = plan->plan_col.M / 2 + 1;
    mwSignedIndex N = plan->N;
    mwSignedIndex kernh = plan->plan_col.kernh;
    mwSignedIndex kernw2 = plan->plan_col.kernw2;
    mwSignedIndex M2buf = M2 + kernh - 1;
    /* mwSignedIndex Nbuf = N + kernw -1; */
    int do_onthefly = plan->plan_col.flags & MOD_COEFFICIENTWISE;
    int do_framewise = plan->plan_col.flags & MOD_FRAMEWISE;
    int do_revorder = plan->plan_col.flags & ORDER_REV;


    double* s = plan->s;
    /* Clear the output */
    memset(coutr, 0, M2 * N * sizeof * coutr);
    memset(couti, 0, M2 * N * sizeof * couti);

    double** kr = plan->kr;
    double** ki = plan->ki;

    double* bufr = plan->bufr;
    double* bufi = plan->bufi;

    /* mexPrintf(" M2 %d, N %d, kernh %d, kernw %d, kernskipd %d \n",M2,N,kernh,kernw,kernskipd); */

    mwSignedIndex n, nfirst;

    /* Copy input to the buffer and do explicit periodic extensions at the left
     * and right borders and conjugated symmetric extensions at the top and bottom borders. */
    extendborders(&plan->plan_col, cr, ci, N, bufr, bufi);

    /* Outside loop over columns */
    for (n = kernw2 - 1, nfirst = 0; nfirst < N; n++, nfirst++)
    {

        /* mexPrintf("n-loop: %d,%d\n",nfirst,n); */
        /* Pick the right kernel */
        double* actKr = kr[nfirst % plan->kernelNo];
        double* actKi = ki[nfirst % plan->kernelNo];
        /* Go to the n-th col in output*/
        double* coutrCol = coutr + nfirst * M2;
        double* coutiCol = couti + nfirst * M2;

        double* crColFirst = bufr + nfirst * M2buf;
        double* ciColFirst = bufi + nfirst * M2buf;

        double* sCol = s + nfirst * M2;

        leglaupdatereal_execute_col(&plan->plan_col, crColFirst, ciColFirst,
                                    actKr, actKi, sCol, coutrCol, coutiCol);

    }

    if (!do_onthefly && !do_framewise)
    {
        /* Update the phase only after the projection has been done. */
        for (mwSignedIndex n = 0; n < N * M2; n++)
        {
            double arg = atan2(couti[n], coutr[n]);
            coutr[n] = s[n] * cos(arg);
            couti[n] = s[n] * sin(arg);
        }
    }

}

void
leglaupdate_done(leglaupdate_plan* plan)
{
    for (int n = 0; n < plan->kernelNo; n++)
    {
        free(plan->kr[n]);
        free(plan->ki[n]);
    }

    if (plan->kr) free(plan->kr);

    if (plan->ki) free(plan->ki);

    free(plan->bufr);
    free(plan->bufi);
    free(plan);
}

/* kernmod is kernh/2+1 X kernw */
void
kernphasefi(double* kernr, double* kerni, mwSignedIndex kernh,
            mwSignedIndex kernw,
            mwSignedIndex kernwskip,
            mwSignedIndex n, mwSignedIndex a, mwSignedIndex M,
            double* kernmodr, double* kernmodi)
{
    formatkernel(kernr, kerni, kernh, kernw, kernwskip, kernmodr, kernmodi);

    if (n != 0)
    {
        mwSignedIndex kernh2 = kernh / 2 + 1;
        mwSignedIndex kernwskipd = kernwskip / (sizeof * kernr);
        /*Modulate */
        double idx = - ( kernh - kernh2);

        for (mwSignedIndex ii = 0; ii < kernh2; ii++)
        {
            double* kernmodrRow = kernmodr + ii;
            double* kernmodiRow = kernmodi + ii;
            double arg = 2.0 * M_PI * n * a / M * (idx + ii);
            double mulr = cos(arg);
            double muli = sin(arg);

            for (mwSignedIndex jj = 0; jj < kernw; jj++)
            {
                double aRe = kernmodrRow[jj * kernwskipd];
                double aIm = kernmodiRow[jj * kernwskipd];
                //double complex tmp = CMPLX(mulr,muli)*CMPLX(kernmodrRow[jj*kernwskipd],kernmodiRow[jj*kernwskipd]);
                kernmodrRow[jj * kernwskipd] = aRe * mulr - aIm * muli;
                kernmodiRow[jj * kernwskipd] = aRe * muli + aIm * mulr;
            }
        }

    }
}

void
formatkernel(double* kernr, double* kerni,
             mwSignedIndex kernh, mwSignedIndex kernw,
             mwSignedIndex kernwskip, double* kernmodr, double* kernmodi)
{
    mwSignedIndex kernh2 = kernh / 2 + 1;
    mwSignedIndex kernwskipd = kernwskip / (sizeof * kernr);

    memset(kernmodr, 0, kernwskip * kernw);
    memset(kernmodr, 0, kernwskip * kernw);

    for (mwSignedIndex ii = 0; ii < kernw; ii++)
    {
        memcpy(kernmodr + ii * kernwskipd, kernr + ii * kernh,
               kernh2 * sizeof * kernmodr);
        memcpy(kernmodi + ii * kernwskipd, kerni + ii * kernh,
               kernh2 * sizeof * kernmodi);
    }

}


mwSignedIndex lcm(mwSignedIndex m, mwSignedIndex n)
{
    return m / gcd(m, n) * n;
}

mwSignedIndex gcd(mwSignedIndex m, mwSignedIndex n)
{
    mwSignedIndex tmp;

    while (m)
    {
        tmp = m;
        m = n % m;
        n = tmp;
    }

    return n;
}
