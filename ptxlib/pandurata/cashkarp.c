//#include "bhdisk.h"
//#include "panhead.h"

void cashkarp(double dt, double y0[], double yn[], double del[])
{
  int i;
  double k1[8],k2[8],k3[8],k4[8],k5[8],k6[8],yns[8];
  double a2=0.2,a3=0.3,a4=0.6,a5=1,a6=0.875;
  double b21=0.2,b31=0.075,b32=0.225,b41=0.3,b42=-0.9,b43=1.2,
    b51=-11./54.,b52=2.5,b53=-70./27.,b54=35./27.,
    b61=1631./55296.,b62=175./512.,b63=575./13824.,b64=44275./110592.,
    b65=253./4096.;
  double c1=37./378.0,c2=0.0,c3=250./621.,c4=125./594.,c5=0.0,c6=512./1771.;
  double cs1=2825./27648.,cs2=0,cs3=18575./48384.,cs4=13525./55296.,
    cs5=277./14336.,cs6=0.25;
  accel(y0,k1);
  for (i=0;i<=7;i++) {
    k1[i]=k1[i]*dt;
    yn[i]=y0[i]+b21*k1[i];
  }
  accel(yn,k2); 
  for (i=0;i<=7;i++) {
    k2[i]=k2[i]*dt;
    yn[i]=y0[i]+b31*k1[i]+b32*k2[i];
  } 
  accel(yn,k3); 
  for (i=0;i<=7;i++) {
    k3[i]=k3[i]*dt;
    yn[i]=y0[i]+b41*k1[i]+b42*k2[i]+b43*k3[i];
  } 
  accel(yn,k4); 
  for (i=0;i<=7;i++) {
    k4[i]=k4[i]*dt;
    yn[i]=y0[i]+b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i];
  } 
  accel(yn,k5); 
  for (i=0;i<=7;i++) {
    k5[i]=k5[i]*dt;
    yn[i]=y0[i]+b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i];
  } 
  accel(yn,k6); 

  for (i=0;i<=7;i++) { 
    k6[i]=dt*k6[i];
    yn[i]=y0[i]+c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]; 
    yns[i]=y0[i]+cs1*k1[i]+cs2*k2[i]+cs3*k3[i]+cs4*k4[i]+cs5*k5[i]+cs6*k6[i]; 
    del[i]=yn[i]-yns[i];
  } 
}

