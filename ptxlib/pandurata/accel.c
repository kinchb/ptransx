//#include "bhdisk.h"
#include "panhead.h"

void accel(double y[], double dy[])
{
  double r,r2,a2,M2,Sig,sinth,sth2,costh,cth2,omg,obr,alp,Del,p_mod,
    domg_dr,domg_dt,dalp_dr,dalp_dt,dobr_dr,dobr_dt,al2,al2Sig,
    o_Sig,o_Sig2,Sig2,SigDel,sthcth,z1,o_p_mod,o_obr;
  r = y[1];
  r2 = r*r;
  a2 = aa*aa;
  M2 = 2.*M;
  sinth = sin(y[2]);
  sth2 = sinth*sinth;
  costh = cos(y[2]);
  cth2 = 1.-sth2;
  sthcth = sinth*costh;
  Sig = r2+a2*cth2;
  Del = r2-M2*r+a2;
  Sig2 = Sig*Sig;
  o_Sig = 1./Sig;
  o_Sig2 = o_Sig*o_Sig;
  SigDel = Sig*Del;
  obr = (r2+a2+M2*r*a2*sth2*o_Sig)*sth2;
  o_obr = 1./obr;
  z1 = 1./(SigDel+M2*r*(a2+r2));
  omg = M2*aa*r*z1;// /(SigDel+2*M*r*(a2+r2));
  al2 = SigDel*z1;// /(SigDel+2*M*r*(a2+r2));
  alp = sqrt(al2);
  al2Sig = al2*o_Sig;
  p_mod = y[4]+omg*y[7];
  o_p_mod = 1./p_mod;

  //domg_dr = -omg*r/(SigDel+2*M*r*(a2+r2))*(3.*r2+a2*(1+cth2)-a2*a2*cth2/r2);
  //domg_dt = -omg*r/(SigDel+2*M*r*(a2+r2))*(2.*M*a2-a2*r-a2*a2/r)*sinth*costh;
  domg_dr = -omg*r*z1;// /(SigDel+2*M*r*(a2+r2));
  domg_dt = domg_dr*(M2*a2-a2*r-a2*a2/r)*sthcth;
  domg_dr = domg_dr*(3.*r2+a2*(1+cth2)-a2*a2*cth2/r2);

  dalp_dr = -alp*al2Sig*M/Del*((a2*a2-r2*r2)/Del-2.*r2*a2*sth2*o_Sig);
  dalp_dt = -alp*al2Sig*M2*a2*r*sthcth*(a2+r2)/(SigDel);

  dy[1] = -y[5]*o_p_mod*al2Sig*Del;
  dy[2] = -y[6]*o_p_mod*al2Sig;
  if (sth2 != 0) {
    dobr_dr = -(sth2*(2.*r+M2*a2*sth2*(a2*cth2-r2)*o_Sig2))*(o_obr*o_obr);
    dobr_dt = -(4.*M*a2*r*sth2*sthcth*(r2+a2+Sig)*o_Sig2 +
		2.*sthcth*(r2+a2))*(o_obr*o_obr);
    dy[3] = (omg*y[4]+(omg*omg-al2*o_obr)*y[7])*o_p_mod;
 
    dy[5] = -domg_dr*y[7]+0.5*al2*o_p_mod
      *(2*o_Sig*(r-M-r*Del*o_Sig)*y[5]*y[5]
	-2*r*o_Sig2*y[6]*y[6]+dobr_dr*y[7]*y[7]) 
      +p_mod/alp*dalp_dr;
    dy[6] = -domg_dt*y[7]+0.5*al2*o_p_mod
      *(2.*a2*sthcth*o_Sig2*(Del*y[5]*y[5]+y[6]*y[6])
	+dobr_dt*y[7]*y[7])
      +p_mod/alp*dalp_dt;
  } else {
    dy[3] = omg;
    dy[5] = 0.5*al2/y[4]*(2*o_Sig*(r-M-r*Del*o_Sig)*y[5]*y[5]
			  -2*r*o_Sig2*y[6]*y[6])
      +y[4]/alp*dalp_dr;
    dy[6] = 0;
  }
  dy[0] = 0;
  dy[4] = 0;
  dy[7] = 0;
}
