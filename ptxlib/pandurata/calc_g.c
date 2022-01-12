//#include "bhdisk.h"
#include "panhead.h"

void calc_g(double g_dn[4][4], double g_up[4][4], double y[])
{
  double r,r2,a2,Sig,Del,cth,sth2,cth2;
  int i,j;
  r = y[1];
  r2 = r*r;
  a2 = aa*aa;
  cth = cos(y[2]);
  cth2 = cth*cth;
  sth2 = 1.-cth2;
  //sth2 = sin(y[2])*sin(y[2]);
  //cth2 = cos(y[2])*cos(y[2]);
  Sig = r2+a2*cth2;
  Del = r2-2*M*r+a2;
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      g_dn[i][j] = 0;
      g_up[i][j] = 0;
    }
  }
  g_dn[0][0] = -(1-2*M*r/Sig);
  g_dn[1][1] = Sig/Del;
  g_dn[2][2] = Sig; 
  g_dn[3][3] = (r2+a2+2*M*a2*r*sth2/Sig)*sth2;
  g_dn[0][3] = -2*aa*M*r*sth2/Sig;
  g_dn[3][0] = g_dn[0][3];

  g_up[0][0] = -((r2+a2)*(r2+a2)-a2*Del*sth2)/(Sig*Del);
  g_up[1][1] = Del/Sig;
  g_up[2][2] = 1/Sig;
  if (sth2 != 0) {
    g_up[3][3] = (Del-a2*sth2)/(Sig*Del*sth2);
  } else g_up[3][3] = 0;
  g_up[0][3] = -2*aa*M*r/(Sig*Del);
  g_up[3][0] = g_up[0][3];

}
