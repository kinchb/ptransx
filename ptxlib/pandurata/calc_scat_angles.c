#include "panhead.h"

void calc_scat_angles(double deg, double *mu, double *psi)
{
  double cth,lambda, za, zb, zq, xlo,xmid,xhi,ylo,ymid,yhi;
  int i;

  lambda = (double)rand()/(RAND_MAX);
  zq = 1.-2.*lambda;
  za = -2.*zq+sqrt(4.*zq*zq+1.);
  zb = -2.*zq-sqrt(4.*zq*zq+1.);
  cth = pow(za,1./3.)-pow(-zb,1./3.);
  lambda = (double)rand()/(RAND_MAX);
  za = (cth*cth-1.)/(cth*cth+1.);
  xlo = 0;
  xhi = 2*PI;
  ylo = 0-lambda;
  yhi = 1-lambda;
  for (i=0;i<=15;i++) {
    xmid = 0.5*(xlo+xhi);
    ymid = xmid/(2*PI)+deg/(4.*PI)*za*sin(2.*xmid)-lambda;
    //printf("%d %12.5e %12.5e %12.5e\n",i,xlo,xmid,xhi);
    //printf("%d %12.5e %12.5e %12.5e\n",i,ylo,ymid,yhi);
    if (ymid <= 0) {
      xlo = xmid;
      ylo = ymid;
    }
    if (ymid > 0) {
      xhi = xmid;
      yhi = ymid;
    }
  }
  *psi = xmid;
  *mu = cth;
}

