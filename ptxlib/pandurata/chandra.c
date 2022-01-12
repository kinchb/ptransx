#include <math.h>
#define eps 1.38e-16

void chandra_limit(double mu, double *deg, double *darken)
{
  double dmu,mu_[21],subtot,subtot2,wlo,whi;
  double I_mu[21] = {0.41441,0.47490,0.52397,0.57001,0.61439,0.65770,
		      0.70029,0.74234,0.78398,0.82530,0.86637,
		      0.90722,0.94789,0.98842,1.02882,1.06911,
		      1.10931,1.14943,1.18947,1.22945,1.26938},
    deg_mu[21] = {0.11713,0.08979,0.07448,0.06311,0.05410,0.04667,
		  0.04041,0.03502,0.03033,0.02619,0.02252,
		  0.01923,0.01627,0.01358,0.011123,0.008880,
		  0.006818,0.004919,0.003155,0.001522,0.0};
  int j,mlo,mhi,Nmu=21;
  
  dmu = 1.0/(Nmu-1.);
  for (j=0;j<Nmu;j++) mu_[j]=j*dmu;
  mlo = (int)(mu/dmu);
  mhi = mlo+1;
  if (mhi < Nmu) wlo = (mu_[mhi]-mu)/dmu;
  if (mhi == Nmu) {
    mlo = Nmu-2;
    mhi = Nmu-1;
    wlo = 0;
  }
  whi = 1.-wlo;
  *deg = wlo*deg_mu[mlo]+whi*deg_mu[mhi];
  *darken = wlo*I_mu[mlo]+whi*I_mu[mhi];
  //printf("%d %d %g %g %g\n",mlo,mhi,mu,*deg,*darken);

  return;
}

void diffuse_reflection(double mu0, double ph0, double mu, double ph,
			double *I_l, double *I_r, double *U_)
{
  double dmu,mu_[21],wlo,whi;
  double Cpsi,Cphi,Cchi,Czeta,CuH1,CH2,
    Cpsi0,Cphi0,Cchi0,Czeta0,CuH10,CH20,dph0,cos2dph0,sin2dph0,mu2,mu02,
    tmpI_l,tmpI_r,tmpU_;
  double S1[3][3], S2[3][3], S3[3][3], S4[3][3], QQ[3][3], SS[3][3], S_[3][3];
  double ch_psi[21]={0.000000, 0.0407500, .0914400, .151050, 0.219160,
		     0.295570, 0.380140, 0.472760, 0.573380, 0.681950,
		     0.798430, 0.922770, 1.05497,  1.19501,  1.34286,
		     1.498520, 1.661980, 1.83321,  2.01223,  2.19902, 2.39357},
         ch_phi[21]={1.00000, 1.12988,  1.20976,  1.26850,  1.31108,
		     1.33973, 1.35569,  1.35971,  1.35228,  1.33374,
		     1.30437, 1.26432,  1.21375,  1.15279,  1.08153,
		     1.00003, 0.908360, 0.806550, 0.694680, 0.572760, 0.440830},
	 ch_chi[21]={1.00000, 1.10352, 1.18638, 1.26329, 1.33687,
		     1.40821, 1.47801, 1.54664, 1.61435, 1.68132,
		     1.74772, 1.81362, 1.87911, 1.94425, 2.00907,
		     2.07365, 2.13799, 2.20213, 2.26609, 2.32990, 2.39356},
	 ch_zeta[21]={0.00000,  0.01824, 0.03764, 0.0578000, 0.078520,
		      0.099690, 0.12121, 0.14303, 0.165100,  0.187380,
		      0.209840, 0.23247, 0.25523, 0.278120,  0.301120,
		      0.324210, 0.34739, 0.37065, 0.393980,  0.417380, 0.44083},
	 ch_uH1[21]={1.00000, 1.07167, 1.11602, 1.14837, 1.17155,
		     1.18685, 1.19487, 1.19599, 1.19030, 1.17774,
		     1.15816, 1.13118, 1.09624, 1.05256, 0.99899,
		     0.93381, 0.85435, 0.75611, 0.63033, 0.45471, 0.00000},
	 ch_H2[21]={1.00000, 1.04967, 1.08621, 1.11762, 1.14552,
		    1.17075, 1.19383, 1.21508, 1.23476, 1.25308,
		    1.27019, 1.28624, 1.30132, 1.31554, 1.32895,
		    1.34166, 1.35371, 1.36515, 1.37601, 1.38638, 1.39625};

  int i,j,mlo,mhi,Nmu=21;
  dmu = 1.0/(Nmu-1.);
  for (j=0;j<Nmu;j++) mu_[j]=j*dmu;
  dph0 = ph0-ph;
  cos2dph0 = cos(2.*dph0);
  sin2dph0 = sin(2.*dph0);
  mu2 = mu*mu;
  mu02 = mu0*mu0;
  mlo = (int)(mu/dmu);
  mhi = mlo+1;
  if (mhi < Nmu) wlo = (mu_[mhi]-mu)/dmu;
  if (mhi == Nmu) {
    mlo = Nmu-2;
    mhi = Nmu-1;
    wlo = 0;
  }
  whi = 1.-wlo;
  Cpsi = wlo*ch_psi[mlo]+whi*ch_psi[mhi];
  Cphi = wlo*ch_phi[mlo]+whi*ch_phi[mhi];
  Cchi = wlo*ch_chi[mlo]+whi*ch_chi[mhi];
  Czeta = wlo*ch_zeta[mlo]+whi*ch_zeta[mhi];
  CuH1 = wlo*ch_uH1[mlo]+whi*ch_uH1[mhi];
  CH2 = wlo*ch_H2[mlo]+whi*ch_H2[mhi];

  mlo = (int)(mu0/dmu);
  mhi = mlo+1;
  if (mhi < Nmu) wlo = (mu_[mhi]-mu0)/dmu;
  if (mhi == Nmu) {
    mlo = Nmu-2;
    mhi = Nmu-1;
    wlo = 0;
  }
  whi = 1.-wlo;
  Cpsi0 = wlo*ch_psi[mlo]+whi*ch_psi[mhi];
  Cphi0 = wlo*ch_phi[mlo]+whi*ch_phi[mhi];
  Cchi0 = wlo*ch_chi[mlo]+whi*ch_chi[mhi];
  Czeta0 = wlo*ch_zeta[mlo]+whi*ch_zeta[mhi];
  CuH10 = wlo*ch_uH1[mlo]+whi*ch_uH1[mhi];
  CH20 = wlo*ch_H2[mlo]+whi*ch_H2[mhi];

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      S1[i][j]=0;
      S2[i][j]=0;
      S3[i][j]=0;
      S4[i][j]=0;
      QQ[i][j]=0;
      SS[i][j]=0;
      S_[i][j]=0;
    }
  }
  QQ[0][0]=1;
  QQ[1][1]=1;
  QQ[2][2]=2;

  S1[0][0]=Cpsi;
  S1[0][1]=sqrt(2.)*Cphi;
  S1[1][0]=Cchi;
  S1[1][1]=sqrt(2.)*Czeta;
  S2[0][0]=Cpsi0;
  S2[0][1]=Cchi0;
  S2[1][0]=sqrt(2.)*Cphi0;
  S2[1][1]=sqrt(2.)*Czeta0;
  S3[0][0]=-4.*mu*mu0*cos(dph0);
  S3[0][2]=2.*mu*sin(dph0);
  S3[2][0]=2.*mu0*sin(dph0);
  S3[2][2]=cos(dph0);
  S4[0][0]=mu2*mu02*cos2dph0;
  S4[0][1]=-mu2*cos2dph0;
  S4[0][2]=-mu2*mu0*sin2dph0;
  S4[1][0]=-mu02*cos2dph0;
  S4[1][1]=cos2dph0;
  S4[1][2]=mu0*sin2dph0;
  S4[2][0]=-mu*mu02*sin2dph0;
  S4[2][1]=mu*sin2dph0;
  S4[2][2]=-mu*mu0*cos2dph0;
  m_mult3(S1,S2,S_);
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      SS[i][j]=3./8.*mu/(mu+mu0+eps)*(S_[i][j]+S3[i][j]*CuH1*CuH10+S4[i][j]*CH2*CH20);
    }
  }
  m_mult3(QQ,SS,S_);
  tmpI_l=S_[0][0]*(*I_l)+S_[0][1]*(*I_r)+S_[0][2]*(*U_);
  tmpI_r=S_[1][0]*(*I_l)+S_[1][1]*(*I_r)+S_[1][2]*(*U_);
  tmpU_=S_[2][0]*(*I_l)+S_[2][1]*(*I_r)+S_[2][2]*(*U_);
  *I_l=tmpI_l;
  *I_r=tmpI_r;
  *U_=tmpU_;
  return;
}

