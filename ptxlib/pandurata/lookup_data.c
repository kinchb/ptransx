#include "panhead.h"

void lookup_data(double part_x0[], 
		 double rr[], double tt[], double pp[], double g_dn[4][4],
		 double rho_ijk[], double T_ijk[], double bb_ijk[], 
		 double ut_ijk[], double ur_ijk[], double uz_ijk[], 
		 double up_ijk[], double weights[], double *rho0, 
		 double *T_e0, double *Bmag, double part_v[])
{
  double r,t,f,dth,dph,w_r[2],w_z[2],w_p[2],wght,a2,Sig,Delta,omg,alpha,omgb2,
    za,zb,zc,nrm,det,g_up[4][4];
  int i,j,k,i_r[2],i_z[2],i_p[2];
  int ingrid;

  int r_ndx, t_ndx, p_ndx;

  ingrid = 0;
  r = part_x0[1];
  t = part_x0[2];
  t = acos(cos(t));
  f = part_x0[3];
  while (f < 0.0) f+=2.*PI;
  f = fmod(f,PI/2.);
  dph = pp[1]-pp[0];
  dth = tt[1]-tt[0];
  if ((r < rr[Nr])&&(r > rr[0])&&(t > tt[0])&&(t < tt[Nth])) {
    ingrid = 1;
  }

  if (ingrid == 0) {
    *rho0 = 0;
    *T_e0 = 0;
    *Bmag = 0;
    a2=aa*aa;
    Sig = r*r+a2*cos(t)*cos(t);      //rho^2 in some texts
    Delta = r*r-2*M*r+a2;
    alpha = sqrt(Sig*Delta/(Sig*Delta+2*M*r*(a2+r*r)));
    omg = 2.*M*r*aa/(Sig*Delta+2.*M*r*(a2+r*r));
    part_v[0]=1./alpha;
    part_v[1]=0;
    part_v[2]=0;
    part_v[3]=omg/alpha;
  }
  if (ingrid == 1) {
    i_r[0]=(int)(log(r/rr[0])/log(rr[1]/rr[0]));
    i_r[1]=i_r[0]+1;
    w_r[0]=(rr[i_r[1]]-r)/(rr[i_r[1]]-rr[i_r[0]]);
    w_r[1]=1.-w_r[0];
    i_z[0]=(int)((t-tt[0])/dth);
    i_z[1]=i_z[0]+1;
    w_z[0]=(tt[i_z[1]]-t)/dth;
    w_z[1]=1.-w_z[0];
    i_p[0]=(int)((f+pp[0])/dph);
    //    printf("%d %d %d %d %d %d\n",i_r[0],i_r[1],i_z[0],i_z[1],i_p[0],i_p[1]);
    if (i_p[0] == 0) {
      i_p[0] = Nph;
      i_p[1] = 0;
      w_p[0] = (pp[0]-f)/dph;
    } else if (i_p[0] == Nph+1) {
      i_p[0] = Nph;
      i_p[1] = 0;
      w_p[0] = (pp[Nph]+dph-f)/dph;
    } else {
      i_p[0] = i_p[0]-1;
      i_p[1] = i_p[0]+1;
      w_p[0] = (pp[i_p[1]]-f)/dph;
    }
    w_p[1] = 1.-w_p[0];
    if (i_p[0] > Nph) {
      printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",r,t,f,w_r[0],w_z[0],w_p[0]);
      printf("%d %d %d %d %d %d\n",i_r[0],i_r[1],i_z[0],i_z[1],i_p[0],i_p[1]);
    }
    //PERFORM BILINEAR INTERPOLATION ON GRID
    *rho0 = 0.;
    *T_e0 = 0.;
    *Bmag = 0.;
    part_v[0] = 0.;
    part_v[1] = 0.;
    part_v[2] = 0.;
    part_v[3] = 0.;
    for (i=0;i<=1;i++) {
        for (j=0;j<=1;j++) {
	        for (k=0;k<=1;k++) {
                wght = w_r[i]*w_z[j]*w_p[k];
                *rho0+=rho_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
                *T_e0+=T_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
                *Bmag+=bb_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
                part_v[0]+=ut_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
                part_v[1]+=ur_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
                part_v[2]+=uz_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
                part_v[3]+=up_ijk[indexijk(i_r[i],i_z[j],i_p[k])]*wght;
	        }
        }
    }

    //NORMALIZE THE 4-VELOCITY
    za = g_dn[0][0];
    zb = 2.*g_dn[0][3]*part_v[3];
    zc = g_dn[1][1]*part_v[1]*part_v[1]+g_dn[2][2]*part_v[2]*part_v[2] +
      g_dn[3][3]*part_v[3]*part_v[3]+1.;
    det = zb*zb-4.*za*zc;
    if (det >= 0) {
      part_v[0]=(-zb-sqrt(det))/(2.*za);
    }
    //IF NOT NORMALIZABLE, SET TO ZAMO
    if ((det <= 0)||(isnan(det) != 0)) {
      *rho0 = 0;
      *T_e0 = 0;
      *Bmag = 0;
      a2=aa*aa;
      Sig = r*r+a2*cos(t)*cos(t);      //rho^2 in some texts
      Delta = r*r-2*M*r+a2;
      alpha = sqrt(Sig*Delta/(Sig*Delta+2*M*r*(a2+r*r)));
      omg = 2.*M*r*aa/(Sig*Delta+2.*M*r*(a2+r*r));
      omgb2 = (Sig*Delta+2.*M*r*(a2+r*r))/Sig*sin(t)*sin(t);
      part_v[0]=1./alpha;
      part_v[1]=0;
      part_v[2]=0;
      part_v[3]=omg/alpha;
    }
    //if (isnan(part_v[3]) != 0) printf("%g %g\n",alpha,det);
    nrm = g_dn[0][0]*part_v[0]*part_v[0]+g_dn[1][1]*part_v[1]*part_v[1]+
      g_dn[2][2]*part_v[2]*part_v[2]+g_dn[3][3]*part_v[3]*part_v[3]+
      2.*g_dn[0][3]*part_v[0]*part_v[3];
    //    calc_g(g_dn,g_up,part_x0);
    if (nrm > -0.99) {
      printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	     part_v[0],part_v[1],part_v[2],part_v[3],w_r[1],w_z[1],w_p[1]);
      printf("%12.5e %12.5e %g %g\n",nrm,det,1./alpha,omg/alpha);
    }
    for (i=0;i<=1;i++) {
      weights[i]=i_r[i];
      weights[i+2]=i_z[i];
      weights[i+4]=i_p[i];
      weights[i+6]=w_r[i];
      weights[i+8]=w_z[i];
      weights[i+10]=w_p[i];
    }

  r_ndx = (int)(weights[0] + 0.1);
  t_ndx = (int)(weights[2] + 0.1);

  if ((r - rr[r_ndx]) > (rr[r_ndx+1] - r)) {
      r_ndx++;
  }
  if ((t - tt[t_ndx]) > (tt[t_ndx+1] - t)) {
      t_ndx++;
  }
  if (f < pp[0]) {
      p_ndx = 0;
  } else if (f >= pp[Nph]) {
      p_ndx = Nph;
  } else {
      p_ndx = (int)((f - pp[0])/dph);
      if ((f - pp[p_ndx]) > (pp[p_ndx+1] - f)) {
          p_ndx++;
      }
  }

  *T_e0 = T_ijk[indexijk(r_ndx,t_ndx,p_ndx)];

  weights[0] = r_ndx;
  weights[2] = t_ndx;
  weights[4] = p_ndx;

  }
  return;
}

