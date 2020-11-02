/*
 * Calls the blas routine to compute a norm of a vector.
 */

#include <math.h>
#include "declarations.h"

double norm2(n,x)
     int n;
     double *x;
{
  double nrm;
  int incx=1;

#ifdef ARMA_BLAS_UNDERSCORE
#ifdef ARMA_BLAS_CAPITALS
  nrm=DNRM2(&n,x,&incx);
#else
  nrm=dnrm2(&n,x,&incx);
#endif
#else
#ifdef ARMA_BLAS_CAPITALS
  nrm=DNRM2_(&n,x,&incx);
#else
  nrm=dnrm2_(&n,x,&incx);
#endif
#endif
  
  return(nrm);
}

double norm1(n,x)
     int n;
     double *x;
{
  double nrm;
  int incx=1;

#ifdef ARMA_BLAS_UNDERSCORE
#ifdef ARMA_BLAS_CAPITALS
  nrm=DASUM(&n,x,&incx);
#else
  nrm=dasum(&n,x,&incx);
#endif
#else
#ifdef ARMA_BLAS_CAPITALS
  nrm=DASUM_(&n,x,&incx);
#else
  nrm=dasum_(&n,x,&incx);
#endif
#endif
  
  return(nrm);
}

double norminf(n,x)
     int n;
     double *x;
{
  int i;
  double nrm;
  int incx=1;

#ifdef ARMA_BLAS_UNDERSCORE
#ifdef ARMA_BLAS_CAPITALS
  i=IDAMAX(&n,x,&incx);
  nrm=fabs(x[i-1]);
#else
  i=idamax(&n,x,&incx);
  nrm=fabs(x[i-1]);
#endif
#else
#ifdef ARMA_BLAS_CAPITALS
  i=IDAMAX_(&n,x,&incx);
  nrm=fabs(x[i-1]);
#else
  i=idamax_(&n,x,&incx);
  nrm=fabs(x[i-1]);
#endif
#endif
  
  return(nrm);
}



