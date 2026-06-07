#include        <stdlib.h>
#include        <math.h>
#include        "Eigen.h"

#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

/*
*       Routines adapted from Numerical Recipes in C
*       (Press, et al. Cambridge University Press, 1988)
*       Modifications by David A. Clausi
*       Department of Systems Design Engineering
*       University of Waterloo
*       Waterloo, Ontario
*       March, 1993
*/


//* Eigenvalue Decomposition Routine *//
//* A,V are [1..n][1..n], d is [1..n] *//
void eigen( double** A, int n, double* d, double** V )
{

  jacobi(A,n,d,V);

}

/* jacobi */
/*
* Computes all eigenvalues and eigenvalues of a real symmetric matrix
*  A[1..n][1..n]. On output, elements of a above the diagonal are destroyed.  d[1..n]
*  returns the eigenvalues of A.  V[1..n][1..n] is a matrix whose columns contain,
*  on output, the normalized eigenvectors of A. 
*
* pp. 364-366, Numerical Recipes in C
*
*/
void jacobi( double** A, int n, double* d, double** V )
 {

  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;
  double* b = (double*) malloc( (unsigned) n * sizeof(double) ); b--;
  double* z = (double*) malloc( (unsigned) n * sizeof(double) ); z--;

  for ( ip = 1; ip <= n; ip++ )  
  {
	  for ( iq = 1; iq <= n; iq++ ) 
		  V[ip][iq] = 0.0;
	  V[ip][ip] = 1.0;
  }

  for ( ip = 1; ip <= n; ip++ ) 
  {
    b[ip] = d[ip] = A[ip][ip];
    z[ip] = 0.0;
  }

  for (i=1;i<=50;i++)  {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) { /* sum off-diag elements */
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(A[ip][iq]);
    }
    if (sm == 0.0)  {
      free(++z);
      free(++b);
      return;
    }
    if (i<4)  tresh=0.2*sm/(n*n);
    else tresh=0.0;
    for (ip=1;ip<=n-1;ip++)  {
      for (iq=ip+1;iq<=n;iq++)  {
        g=100.0*fabs(A[ip][iq]);
        if (i>4 && fabs(d[ip])+g == fabs(d[ip])
	  && fabs(d[iq])+g == fabs(d[iq]))
	  A[ip][iq]=0.0;
        else if (fabs(A[ip][iq])>tresh)  {
	  h=d[iq]-d[ip];
	  if (fabs(h)+g == fabs(h))
	    t=(A[ip][iq])/h;
          else  {
	    theta=0.5*h/(A[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t= -t;
          }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*A[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  A[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {  ROTATE(A,j,ip,j,iq) }
	  for (j=ip+1;j<=iq-1;j++) {  ROTATE(A,ip,j,j,iq) }
	  for (j=iq+1;j<=n;j++) {  ROTATE(A,ip,j,iq,j) }
	  for (j=1;j<=n;j++) {  ROTATE(V,j,ip,j,iq) }
        }
      }
    }

    for(ip=1;ip<=n;ip++)  {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }

  free(++z);
  free(++b);
}

