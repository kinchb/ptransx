#include <math.h>
#define TINY 1.0e-20

void ludcmp_js(double **a, int n, int *indx, double *d, int call, int ir, double radius)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=(double *)malloc((n+1)*sizeof(double));
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) {
		  printf("Singular matrix in routine ludcmp: %d, %d, %lf\n", call, ir, radius);

		  printf("\n");
		  printf("%12.5e %12.5e %12.5e %12.5e\n",
			 a[1][1],a[1][2],a[1][3],a[1][4]);
		  printf("%12.5e %12.5e %12.5e %12.5e\n",
			 a[2][1],a[2][2],a[2][3],a[2][4]);
		  printf("%12.5e %12.5e %12.5e %12.5e\n",
			 a[3][1],a[3][2],a[3][3],a[3][4]);
		  printf("%12.5e %12.5e %12.5e %12.5e\n",
			 a[4][1],a[4][2],a[4][3],a[4][4]);

		}
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free(vv);
}
#undef TINY
