#include <math.h>

double dot_g4(double g_ab[4][4],double a[],double b[])
{
  double dotp;
  dotp = g_ab[0][0]*a[0]*b[0]+g_ab[0][3]*a[0]*b[3]+
    g_ab[1][1]*a[1]*b[1]+g_ab[2][2]*a[2]*b[2]+
    g_ab[3][0]*a[3]*b[0]+g_ab[3][3]*a[3]*b[3];
  return dotp;
}

//given a time-like vector a and metric g_ab, returns a normalized vector 
//a with g_ab.a.a=-1
void norml_tl(double g_ab[4][4], double a[])
{
  double nrm;
  int i;
  nrm = dot_g4(g_ab,a,a);
  nrm = sqrt(fabs(nrm));
  for (i=0;i<=3;i++) a[i]=a[i]/nrm;
}

//given a space-like vector a and metric g_ab, returns a normalized vector 
//a with g_ab.a.a=+1
void norml_sl(double g_ab[4][4], double a[])
{
  double nrm;
  int i;
  nrm = dot_g4(g_ab,a,a);
  nrm = sqrt(fabs(nrm));
  for (i=0;i<=3;i++) a[i]=a[i]/nrm;
}

void boost(double beta, double n_[], double p_[])
{
  double gamma,Lambda[4][4],p_p[4];
  int i,j;
  gamma = 1./sqrt(1.-beta*beta);
  Lambda[0][0]=gamma;
  for (i=1;i<=3;i++) {
    Lambda[0][i] = -beta*gamma*n_[i-1];
    Lambda[i][0] = Lambda[0][i];
    for (j=1;j<=3;j++) {
      Lambda[i][j] = (gamma-1.)*n_[i-1]*n_[j-1];
    }
    Lambda[i][i] = Lambda[i][i]+1.;
  }
  for (i=0;i<=3;i++) {
    p_p[i]=0;
    for (j=0;j<=3;j++) p_p[i]=p_p[i]+Lambda[i][j]*p_[j];
  }
  for (i=0;i<=3;i++) p_[i]=p_p[i];
}

