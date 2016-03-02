/** \addtogroup spsi
 *  \{
 *  \file
 *  \author Zdeněk Průša
 *  \date 1 Feb 2016
 *  \brief SPSI header

 *  
 */
#ifndef _spsi_h
#define _spsi_h
#include "config.h"

/** SPSI algorithm implementation
 *  \param[in] s  M x N arrays, target magnitude
 *  \param[in] a  Hop factor
 *  \param[in] M  Number of frequency channels
 *  \param[in] N  Number of time frames
 *  \param[in,out] initphase Length M vector, [in] phase of -1 frame,
 *                           [out] phase of [N-1] frame.
 *                           Can be NULL, but then a memory allocation will occur.
 *  \param[out] c M x N array of coefficients
 */
void
spsi(const double* s, int a, int M, int N, double* initphase, complex double* c);

/** Masked SPSI algorithm implementation
 *
 *  Works as spsi() except arg(c[ii])=phase[ii] if mask[ii] evaluates to true.
 *
 *  \param[in] s  M x N arrays, target magnitude
 *  \param[in] a  Hop factor
 *  \param[in] M  Number of frequency channels
 *  \param[in] N  Number of time frames
 *  \param[in] mask M x N array, mask of coefficients with known phase
 *  \param[in] phase M x N array, known phase
 *  \param[in,out] initphase Length M vector, [in] phase of -1 frame,
 *                           [out] phase of [N-1] frame.
 *                           Can be NULL, but then a memory allocation will occur.
 *  \param[out] c M x N array of coefficients
 */
void
maskedspsi(const double* s, int a, int M, int N,
               const double* mask, const double* phase, double* initphase, complex double* c);


void
spsiupdate(const double* scol, int a, int M, double* tmpphase);

#endif /* _spsi_h */
/** \} */
