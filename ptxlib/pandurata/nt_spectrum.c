#include "panhead.h"
#define indexre(a,b) ((Nr+1)*(b)+(a))

void nt_spectrum(double Risco, double rr[], double drr[], 
		 double nu[], double Inur[], double Ts_r[], double qnur[])
{
  double Mbh,Mdot,Mdot_17,Z1,Z2,Re,dr,a2,r2,r3,u_star,du,
    subtot,Lisco,Eisco,nom_eff,rhos,Ts,Teff,ps,Lum,Ledd,x,kap_es,kap_ff,f_hard,
    sdens_r[Nr+1],v_r[Nr+1],h_r[Nr+1],g_r[Nr+1],
    alpha_r[Nr+1],rhoc_r[Nr+1],Tc_r[Nr+1],presc_r[Nr+1],rhos_r[Nr+1],
    r_cm[Nr+1],Q_int[Nr+1],flux_r[Nr+1],Hatm_r[Nr+1],
    ntA[Nr+1],ntB[Nr+1],ntC[Nr+1],ntD[Nr+1],ntE[Nr+1],ntF[Nr+1],ntG[Nr+1],
    ntJ[Nr+1],ntsig[Nr+1],ntL[Nr+1],ntQ[Nr+1],subu[Nr+1],rhB[Nr+1],rhC[Nr+1],
    akR1[Nr+1],akR2[Nr+1],akQ[Nr+1];
  //  double rr[Nr+1];
  long i,j,k,ir_in;

  Mbh = Mstar*3.*1.99e33;
  kap_es = 0.4;
  f_hard = 1.8;
  a2 = aa*aa;
  
  //Z1 = 1.+pow((1.-aa*aa),1./3.)*(pow((1.+aa),1./3.)+pow((1.-aa),1./3.));
  //Z2 = sqrt(3.*aa*aa+Z1*Z1);
  //Risco = (3.+Z2-sqrt((3.-Z1)*(3+Z1+2.*Z2)));
  Re = 1.001*(1.+sqrt(1.-aa*aa));
  //  printf("%g %12.5e\n",aa,Risco);
  Lisco = 2./sqrt(3.)*(3.*sqrt(Risco)-2.*aa)/sqrt(Risco);
  Eisco = (Risco*Risco-2*M*Risco+aa*sqrt(M*Risco))/
    (Risco*sqrt(Risco*Risco-3*M*Risco+2*aa*sqrt(M*Risco)));
  nom_eff = del_eta+(1.-Eisco);
  //printf("nominal eff: %g %g\n",Eisco,nom_eff);
  Mdot = 1.38e38*Mstar*3.*L_Edd/(nom_eff*9e20);
  //For emission inside ISCO, scale Mdot so that total luminosity is roughly the same as NT
  if (em_model >= 3) Mdot=Mdot*(0.2+nom_eff);
  //Mdot = 8e17;
  Mdot_17 = Mdot/1e17;
  printf("Mdot: %g %g\n",Mdot,Mdot_17);

  j=-1;
  for (i=0;i<=Nr;i++) {
    if (rr[i]<Risco) j=i;
  }
  ir_in = j+1;
  for (i=ir_in;i<=Nr;i++) {
    //    rr[i] = Risco+((double)i+0.5)*dr;
    r_cm[i] = rr[i]*Mstar*4.4e5;
    r2 = rr[i]*rr[i];
    r3 = rr[i]*r2;
    ntA[i] = 1.+a2/r2+2.*a2/r3;
    ntB[i] = 1.+aa*pow(rr[i],-1.5);
    ntC[i] = 1.-3./rr[i]+2.*aa*pow(rr[i],-1.5);
    ntD[i] = 1.-2./rr[i]+a2/r2;
    ntE[i] = 1.+4.*a2/r2-4.*a2/r3+3.*a2*a2/(r2*r2);
    ntF[i] = 1.-2.*aa*pow(rr[i],-1.5)+a2/r2;
    ntL[i] = ntF[i]/sqrt(ntC[i])-Lisco/sqrt(rr[i]);
    ntsig[i] = 0.75*sqrt(Mbh*Gn*pow(r_cm[i],-3.))*ntD[i]/ntC[i];
    u_star = 1./rr[i];
    du = u_star/Nr;
    subtot = 0.0;
    for (j=0;j<=Nr;j++) {
      subu[j] = (double)j*du;
      subtot += (1.-2.*aa*pow(subu[j],1.5)+a2*subu[i]*subu[i])/
	((1.+aa*pow(subu[j],1.5))*(1.-3.*subu[j]+2.*aa*pow(subu[j],1.5)));
    }
    ntJ[i] = exp(1.5*subtot*du);
    subtot = 0.;
    for (j=ir_in;j<=i;j++) {
      Q_int[j]=ntL[j]*ntF[j]/(ntB[j]*ntC[j]*ntJ[j]*pow(rr[j],1.5));
      //subtot += Q_int[j];
      subtot += Q_int[j]*drr[j];
    }
    //subtot = subtot-0.5*(Q_int[0]+Q_int[i]);
    subtot = subtot-0.5*(Q_int[ir_in]*drr[ir_in]+Q_int[i]*drr[i]);

    //ADD NON-ZERO TORQUE AT ISCO, AS IN AGOL+KROLIK(1999)
    //ntQ[i] = ntL[i]-3.*ntJ[i]/(2.*sqrt(rr[i]))*dr*subtot;
    ntQ[i] = ntL[i]-3.*ntJ[i]/(2.*sqrt(rr[i]))*subtot;
    akR1[i] = ntQ[i]/(ntB[i]*sqrt(ntC[i]));
    akR2[i] = del_eta*Risco*sqrt(Risco)*sqrt(ntC[ir_in])/ntC[i]/sqrt(rr[i]);
    akQ[i] = ntB[i]*sqrt(ntC[i])*(akR1[i]+akR2[i]);
    ntQ[i] = akQ[i];

    rhB[i] = 1.-3./rr[i]+2.*aa*pow(rr[i],-1.5);
    rhC[i] = 1.-4.*aa*pow(rr[i],-1.5)+3.*a2/r2;
    g_r[i] = Gn*Mbh*rhB[i]*rhC[i]*pow(r_cm[i],-3.);
    flux_r[i] = 6.e25*(Mdot_17/(Mstar*Mstar))*ntQ[i]/
      (pow(rr[i],3.)*ntB[i]*sqrt(ntC[i]));
    sdens_r[i] = 20.*(Mstar/Mdot_17/nt_alpha)*pow(rr[i],1.5)*pow(ntB[i],3.)*
      sqrt(ntC[i])*ntE[i]/(ntA[i]*ntA[i]*ntQ[i]);
    h_r[i] = 1.e5*Mdot_17*ntA[i]*ntA[i]*sqrt(ntC[i])*ntQ[i]/
      (pow(ntB[i],3.)*ntD[i]*ntE[i]);
    rhoc_r[i] = sdens_r[i]/(2.*h_r[i]);
    Tc_r[i] = 4.e7*pow(nt_alpha*Mstar,-0.25)*sqrt(ntB[i])*pow(ntE[i],0.25)/
      (pow(rr[i],3./8.)*sqrt(ntA[i]));
    Teff = pow(flux_r[i]/sigma,0.25);
    Ts_r[i] = Teff*f_hard;

    subtot = 0;
    for (j=0;j<=Ne;j++) {
      x = 1.e3*nu[j]/(kB_ev*Ts_r[i]);
      if (rr[i] <= Risco) {
	Inur[indexre(i,j)]=0;
	qnur[indexre(i,j)]=0;
      }
      if (rr[i] > Risco) {
	Inur[indexre(i,j)]=pow(nu[j],3.)/(exp(x)-1.)/pow(f_hard,4.);
	//Inur[indexre(i,j)]=sqrt(rhoc_r[i])*pow(Ts_r[i],5./4.)*
	//pow(x,1.5)*exp(-x/2.)/sqrt(exp(x)-1.);
	//Inur[indexre(i,j)]=1e10*pow(Hatm_r[i],-1./3.)*pow(Ts_r[i],11./6.)*
	//x*x*exp(-x)*pow((1.-exp(-x)),-2./3.);
	kap_ff = 1.5e25*rhoc_r[i]*pow(Ts_r[i],-3.5)*pow(x,-3.)*(1.-exp(-x));
	qnur[indexre(i,j)]=kap_es/(kap_es+kap_ff);
      }
      //printf("%g %g %g\n",rr[i],nu[j],Inur[indexre(i,j)]);
      if (j > 0) {
	subtot=subtot+(nu[j]-nu[j-1])*Inur[indexre(i,j)];
      }
    }
    //printf("%11.4e %11.4e %11.4e\n",rhoc_r[i],Ts_r[i],flux_r[i]/subtot);
  }
  /*
  subtot=0.;
  for (i=0;i<=Nr;i++) {
    subtot+=4.*PI*flux_r[i]*rr[i]*dr*pow(4.4e5*Mstar,2.)/sqrt(ntD[i])/sqrt(ntA[i]/ntD[i]);
  }
  Lum = subtot;
  Ledd = 3.6e38*Mstar;
  printf("Luminosity = %6.3e Ledd\n",Lum/Ledd);
  printf("efficiency = %6.3e\n", Lum/(cc*cc*Mdot));
  */
  return;
}

