//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
//#include "bhdisk.h"

double calc_mag(double x[])
{
  double mag;
  mag = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  return mag;
}

void cross(double a[],double b[],double c[])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

void normalize(double x[])
{
  double xmag;
  xmag = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  if (xmag > 0) {
    x[0]=x[0]/xmag;
    x[1]=x[1]/xmag;
    x[2]=x[2]/xmag;
  }
}

double dot(double a[],double b[])
{
  double dotp;
  dotp = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return dotp;
}

void cart_to_spher(double x[], double r[])
{
  r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  if (r[0] > 0) {
    r[1]=acos(x[2]/r[0]);
  } else r[1]=0;
  r[2]=atan2(x[1],x[0]);
}

void spher_to_cart(double r[], double x[])
{
  x[0]=r[0]*sin(r[1])*cos(r[2]);
  x[1]=r[0]*sin(r[1])*sin(r[2]);
  x[2]=r[0]*cos(r[1]);
}

void cartv_to_spherv(double x[], double v[], double p[])
{
  double r,ph_hat[3];
  double r_hat[3]={0,0,1},th_hat[3]={1,0,0};
  int j;
  r = calc_mag(x);
  if (r != 0) {
    for (j=0;j<=2;j++) r_hat[j]=x[j]/r;
  }
  r = sqrt(r_hat[0]*r_hat[0]+r_hat[1]*r_hat[1]);
  if (r != 0) {
    th_hat[0] = r_hat[0]*r_hat[2]/r;
    th_hat[1] = r_hat[1]*r_hat[2]/r;
    th_hat[2] = -r;
  } else {
    r = sqrt(v[0]*v[0]+v[1]*v[1]);
    if (r !=0) {
      th_hat[0] = v[0]/r;
      th_hat[1] = v[1]/r;
    }
    th_hat[2] = 0;
  }
  cross(r_hat, th_hat, ph_hat);
  p[0]=dot(v,r_hat);
  p[1]=dot(v,th_hat);
  p[2]=dot(v,ph_hat);
}

void spherv_to_cartv(double r[], double p[], double v[])
{
  v[0]=-r[0]*sin(r[1])*sin(r[2])*p[2]+
    (r[0]*cos(r[1])*p[1]+sin(r[1])*p[0])*cos(r[2]);
  v[1]=r[0]*sin(r[1])*cos(r[2])*p[2]+
    (r[0]*cos(r[1])*p[1]+sin(r[1])*p[0])*sin(r[2]);
  v[2]=-r[0]*sin(r[1])*p[1]+cos(r[1])*p[0];
}

