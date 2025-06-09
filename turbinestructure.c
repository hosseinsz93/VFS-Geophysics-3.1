/**
 * Turbine blade structure functions
 * The structure equations are solved using real dimensions, although the flow may be solved without dimension 
 * The equations derived in Kallesoe's paper are solved (Wind Energy 2007; 10:209-230)
 */

#include "variables.h"
#include <algorithm>
using namespace std;

int Nsmooth=1;

/*
 * Direct solver from Recipe using LU decomposition
 */

#define NR_END 1
#define FREE_ARG char*

PetscErrorCode resetMatrixCoeff(Beam *blade);

PetscErrorCode WriteForces_Term(char *file,int NumOfElments, Beam *blade, vector<double> &VeF)
{
  FILE *f;
  char filen[80];
  
  sprintf(filen, file);
  f = fopen(filen, "w");

  for (int i=1; i<NumOfElments; i++)
  {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le \n",blade->r[i],VeF.at(i));
  }
  fclose(f);

  return(0);

}

PetscErrorCode Write_Displacement(char *file,int NumOfElments, int iblade, Beam *blade, IBMNodes *ibm)
{
  FILE *f;
  char filen[80];
  
  sprintf(filen, file);
  f = fopen(filen, "w");
  int ibi=0;

  for (int i=1; i<NumOfElments; i++)
  {
    double u = ibm[ibi].bladestructure[iblade].u[i] ;
    double v = ibm[ibi].bladestructure[iblade].v[i] ;
    double dudt = ibm[ibi].bladestructure[iblade].dudt[i] ;
    double dvdt = ibm[ibi].bladestructure[iblade].dvdt[i] ;
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le\n",blade->r[i],u,v,dudt,dvdt);
  }
  fclose(f);

  return(0);

}

PetscErrorCode Write_Amatrix(char *file,int NumOfElments, int iblade, Beam *blade, IBMNodes *ibm)
{
  FILE *f;
  char filen[80];
  
  sprintf(filen, file);
  f = fopen(filen, "w");
  int ibi=0;


	for (int i=0;i<4*NumOfElments;i++)
	{	
		for (int j=0;j<4*NumOfElments;j++)
		{

         PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n",blade->A_LU[i][j]);
		}
         }	

  fclose(f);
  return(0);

}

PetscErrorCode Write_Bmatrix(char *file,int NumOfElments, int iblade, Beam *blade, IBMNodes *ibm)
{
  FILE *f;
  char filen[80];
  
  sprintf(filen, file);
  f = fopen(filen, "w");
  int ibi=0;


	for (int i=0;i<4*NumOfElments;i++)
	{	

         PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n",blade->b[i]);
         }	

  fclose(f);
  return(0);

}

PetscErrorCode Init_TurbineStructure(IBMNodes *ibm)
{

  for (int ibi=0; ibi<NumberOfTurbines; ibi++) {
    for (int iblade=0; iblade<ibm[ibi].num_blade; iblade++)  {
      int NElmts = ibm[ibi].bladestructure[iblade].NumberofElements;
      for (int i=0; i<NElmts; i++) {
        ibm[ibi].bladestructure[iblade].dudt[i] = 0.;
        ibm[ibi].bladestructure[iblade].dvdt[i] = 0.;
        ibm[ibi].bladestructure[iblade].dudt_o[i] = 0.;
        ibm[ibi].bladestructure[iblade].dvdt_o[i] = 0.;

        ibm[ibi].bladestructure[iblade].coefu_dudr[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefu_d2udr2[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefu_d3udr3[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefu_d4udr4[i] = 0.;

        ibm[ibi].bladestructure[iblade].coefu_dvdr[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefu_d2vdr2[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefu_d3vdr3[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefu_d4vdr4[i] = 0.;

        ibm[ibi].bladestructure[iblade].coefv_dudr[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefv_d2udr2[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefv_d3udr3[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefv_d4udr4[i] = 0.;

        ibm[ibi].bladestructure[iblade].coefv_dvdr[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefv_d2vdr2[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefv_d3vdr3[i] = 0.;
        ibm[ibi].bladestructure[iblade].coefv_d4vdr4[i] = 0.;
      }
    }
  }
  return(0);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *fvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_fvector(double *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

/**
 * LU decomposition
 */
#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	//vv=fvector(1,n);
	vv=new double[n];
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != (n-1)) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	delete[] vv;
}
#undef TINY

/**
 * Direct solver: forward substitution and backward substitution
 */

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum=0.;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii-1;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];


	}
}


/** 
 * Second order finite difference stencils 
 */

/** 
 * Center difference for first order derivative
 */
double dfdr(double fm1, double fp1, double h) 
{
	if (fabs((0.5*fp1-0.5*fm1)/(fm1+1.e-6))<1.e-6) 
		return 0;
	else
		return (0.5*fp1-0.5*fm1)/h;
}

/** 
 * Center difference for second order derivative
 */
double d2fdr2(double fm1, double f, double fp1, double h) 
{
	if (fabs((fp1-2.0*f+fm1)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (fp1-2.0*f+fm1)/pow(h,2);
}

/** 
 * Center difference for third order derivative
 */
double d3fdr3(double fm2, double fm1, double f, double fp1, double fp2, double h) 
{
	if (fabs((-0.5*fm2+1.0*fm1-1.0*fp1+0.5*fp2)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-0.5*fm2+1.0*fm1-1.0*fp1+0.5*fp2)/pow(h,3);
}

/** 
 * Center difference for fourth order derivative
 */
double d4fdr4(double fm2, double fm1, double f, double fp1, double fp2, double h) 
{
	if (fabs((fm2-4.0*fm1+6.0*f-4.0*fp1+fp2)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (fm2-4.0*fm1+6.0*f-4.0*fp1+fp2)/pow(h,4);
}

/** 
 * Approximate first order derivative at j using values at j, j+1 and j+2
 */
double dfdr_forward(double f, double fp1, double fp2, double h)
{	
	if (fabs((-1.5*f+2.0*fp1-0.5*fp2)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-1.5*f+2.0*fp1-0.5*fp2)/h;
}

/** 
 * Approximate first order derivative at j using values at j-2, j-1 and j
 */
double dfdr_backward(double fm2, double fm1, double f, double h)
{
	if (fabs((0.5*fm2-2.0*fm1+1.5*f)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (0.5*fm2-2.0*fm1+1.5*f)/h; 
}

/** 
 * Approximate second order derivative at j using values at j, j+1, j+2 and j+3
 */
double d2fdr2_forward(double f, double fp1, double fp2, double fp3, double h)
{	
	if (fabs((2.0*f-5.0*fp1+4.0*fp2-fp3)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (2.0*f-5.0*fp1+4.0*fp2-fp3)/pow(h,2);
}

/** 
 * approximate second order derivative at j using values at j-3, j-2, j-1 and j
 */
double d2fdr2_backward(double fm3, double fm2, double fm1, double f, double h)
{
	if (fabs((-fm3+4.0*fm2-5.0*fm1+2.0*f)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-fm3+4.0*fm2-5.0*fm1+2.0*f)/pow(h,2); 
}

/** 
 * approximate third order derivative at j using values at j-1, j, j+1, j+2 and j+3
 */
double d3fdr3_forward1(double fm1, double f, double fp1, double fp2, double fp3, double h)
{
	if (fabs((-1.5*fm1+5.0*f-6.0*fp1+3.0*fp2-0.5*fp3)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-1.5*fm1+5.0*f-6.0*fp1+3.0*fp2-0.5*fp3)/pow(h,3); 
}

/** 
 * approximate third order derivative at j using values at j, j+1, j+2, j+3 and j+4
 */
double d3fdr3_forward2(double f, double fp1, double fp2, double fp3, double fp4, double h)
{
	if (fabs((-2.5*f+9.0*fp1-12.0*fp2+7.0*fp3-1.5*fp4)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-2.5*f+9.0*fp1-12.0*fp2+7.0*fp3-1.5*fp4)/pow(h,3); 
}

/** 
 * approximate third order derivative at j using values at j-3, j-2, j-1, j and j+1
 */
double d3fdr3_backward1(double fm3, double fm2, double fm1, double f, double fp1, double h)
{
	if (fabs((0.5*fm3-3.0*fm2+6.0*fm1-5.0*f+1.5*fp1)/(f+1.e-6))<1.e-6) 
		return 0;
	else
		return (0.5*fm3-3.0*fm2+6.0*fm1-5.0*f+1.5*fp1)/pow(h,3); 
}

/** 
 * approximate third order derivative at j using values at j-4, j-3, j-2, j-1 and j
 */
double d3fdr3_backward2(double fm4, double fm3, double fm2, double fm1, double f, double h)
{
	if (fabs((1.5*fm4-7.0*fm3+12.0*fm2-9.0*fm1+2.5*f)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (1.5*fm4-7.0*fm3+12.0*fm2-9.0*fm1+2.5*f)/pow(h,3); 
}

/** 
 * approximate fourth order derivative at j using values at j-1, j, j+1, j+2, j+3 and j+4
 */
double d4fdr4_forward1(double fm1, double f, double fp1, double fp2, double fp3, double fp4, double h)
{
	if (fabs((2.0*fm1-9.0*f+16.0*fp1-14.0*fp2+6.0*fp3-fp4)/(f+1.e-6))<1.e-6)
		return 0;
	else 
		return (2.0*fm1-9.0*f+16.0*fp1-14.0*fp2+6.0*fp3-fp4)/pow(h,4); 
}

/** 
 * approximate fourth order derivative at j using values at j, j+1, j+2, j+3, j+4 and j+5
 */
double d4fdr4_forward2(double f, double fp1, double fp2, double fp3, double fp4, double fp5, double h)
{
	if (fabs((3.0*f-14.0*fp1+26.0*fp2-24.0*fp3+11.0*fp4-2.0*fp5)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (3.0*f-14.0*fp1+26.0*fp2-24.0*fp3+11.0*fp4-2.0*fp5)/pow(h,4); 
}

/** 
 * approximate fourth order derivative at j using values at j-4, j-3, j-2, j-1, j and j+1
 */
double d4fdr4_backward1(double fm4, double fm3, double fm2, double fm1, double f, double fp1, double h)
{
	if (fabs((-1.0*fm4+6.0*fm3-14.0*fm2+16.0*fm1-9.0*f+2.0*fp1)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-1.0*fm4+6.0*fm3-14.0*fm2+16.0*fm1-9.0*f+2.0*fp1)/pow(h,4); 
}

/** 
 * approximate fourth order derivative at j using values at j-4, j-3, j-2, j-1, j and j+1
 */
double d4fdr4_backward2(double fm5, double fm4, double fm3, double fm2, double fm1, double f, double h)
{
	if (fabs((-2.0*fm5+11.0*fm4-24.0*fm3+26.0*fm2-14.0*fm1+3.0*f)/(f+1.e-6))<1.e-6)
		return 0;
	else
		return (-2.0*fm5+11.0*fm4-24.0*fm3+26.0*fm2-14.0*fm1+3.0*f)/pow(h,4); 
}


/**
 * Use least squares to interpolate from one profile to the other profile 

 * A = 
 
 * [ a00, a01, a02]
 * [ a10, a11, a12]
 * [ a20, a21, a22]
 
 * B = 
 
 * f0
 * f1
 * f2

 * A*X=B

 * X =
 
 * (a01*a12*f2 - a02*a11*f2 - a01*a22*f1 + a02*a21*f1 + a11*a22*f0 - a12*a21*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20)
 * -(a00*a12*f2 - a02*a10*f2 - a00*a22*f1 + a02*a20*f1 + a10*a22*f0 - a12*a20*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20)
 * (a00*a11*f2 - a01*a10*f2 - a00*a21*f1 + a01*a20*f1 + a10*a21*f0 - a11*a20*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20)

 * for 5-point, second-order least squares

 * a00=10; a01=a10=2*\sum_i x_i; a02=a20=a11=2*\sum_i x_i^2; a12=a21=2*\sum_i x_i^3; a22=2*\sum_i x_i^4

 * f0=2*\sum_i F_i; f1=2*\sum_i F_i x_i; f2=2*\sum_i F_i x_i^2

 */ 


PetscErrorCode interpolation_leastsquares(vector<double> r0, vector<double> g0, vector<double> r1, vector<double> g1, int N0, int N1)
{
	int jclose, j0, j1, jstart, jend;
	double a00, a01, a02, a10, a11, a12, a20, a21, a22, f0, f1, f2;
	double c0, c1, c2;

	for (j1=0;j1<N1;j1++) {

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d A \n", j1);
		double rr=r1.at(j1);

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d B \n", j1);
		double dd=1.e6;
		for (j0=0; j0<N0; j0++) {
			double dis = pow(rr-r0.at(j0),2);
			if (dis<dd) 
			{
				dd=dis; jclose=j0;
			}
		}	

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d C jclose=%d\n", j1, jclose);
//		if (jclose==0 || jclose==1) {
	//		jstart=0; jend=5;
	//	} else if (jclose==N0-1 || jclose==N0-2 || jclose==N0-3) {
	//		jstart=N0-5; jend=N0;
	//	} else {
			jstart=jclose-2; jend=jclose+3;
	//	}

		if (jstart<0) 
		{
			jstart = 0; jend = jstart+5;
		}

		if (jend>N0-1) 
		{
			jend = N0; jstart =jend-5;
		}

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d jstart=%d jend=%d \n", j1, jstart, jend);
		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d D \n", j1);
		a00=10;
		a01=0; a02=0; a12=0; a22=0; 
		f0=0; f1=0; f2=0;
		for (j0=jstart; j0<jend; j0++) {
			a01+=2*r0.at(j0);	
			a02+=2*pow(r0.at(j0),2);	
			a12+=2*pow(r0.at(j0),3);	
			a22+=2*pow(r0.at(j0),4);	
			
			f0+=2*g0.at(j0);
			f1+=2*g0.at(j0)*r0.at(j0);
			f2+=2*g0.at(j0)*r0.at(j0)*r0.at(j0);
		}				

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d E \n", j1);
		a10=a01; a20=a02; a11=a02; a21=a12;
	
		c0=(a01*a12*f2 - a02*a11*f2 - a01*a22*f1 + a02*a21*f1 + a11*a22*f0 - a12*a21*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
		c1=-(a00*a12*f2 - a02*a10*f2 - a00*a22*f1 + a02*a20*f1 + a10*a22*f0 - a12*a20*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
  		c2=(a00*a11*f2 - a01*a10*f2 - a00*a21*f1 + a01*a20*f1 + a10*a21*f0 - a11*a20*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);

		g1.at(j1)=c0+c1*rr+c2*pow(rr,2);
	
	}

	return(0);
}


PetscErrorCode _interpolation_leastsquares(vector<double> r0, vector<double> g0, vector<double> r1, vector<double> &g1, int N0, int N1)
{
	int jclose, j0, j1, jstart, jend;
	double a00, a01, a02, a10, a11, a12, a20, a21, a22, f0, f1, f2;
	double c0, c1, c2;

	for (j1=0;j1<N1;j1++) {

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d A \n", j1);
		double rr=r1.at(j1);

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d B \n", j1);
		double dd=1.e6;
		for (j0=0; j0<N0; j0++) {
			double dis = pow(rr-r0.at(j0),2);
			if (dis<dd) 
			{
				dd=dis; jclose=j0;
			}
		}	

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d C jclose=%d\n", j1, jclose);
//		if (jclose==0 || jclose==1) {
	//		jstart=0; jend=5;
	//	} else if (jclose==N0-1 || jclose==N0-2 || jclose==N0-3) {
	//		jstart=N0-5; jend=N0;
	//	} else {
			jstart=jclose-2; jend=jclose+3;
	//	}

		if (jstart<0) 
		{
			jstart = 0; jend = jstart+5;
		}

		if (jend>N0-1) 
		{
			jend = N0; jstart =jend-5;
		}

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d jstart=%d jend=%d \n", j1, jstart, jend);
		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d D \n", j1);
		a00=10;
		a01=0; a02=0; a12=0; a22=0; 
		f0=0; f1=0; f2=0;
		for (j0=jstart; j0<jend; j0++) {
			a01+=2*r0.at(j0);	
			a02+=2*pow(r0.at(j0),2);	
			a12+=2*pow(r0.at(j0),3);	
			a22+=2*pow(r0.at(j0),4);	
			
			f0+=2*g0.at(j0);
			f1+=2*g0.at(j0)*r0.at(j0);
			f2+=2*g0.at(j0)*r0.at(j0)*r0.at(j0);
		}				

		//PetscPrintf(PETSC_COMM_WORLD, "j1=%d E \n", j1);
		a10=a01; a20=a02; a11=a02; a21=a12;
	
		c0=(a01*a12*f2 - a02*a11*f2 - a01*a22*f1 + a02*a21*f1 + a11*a22*f0 - a12*a21*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
		c1=-(a00*a12*f2 - a02*a10*f2 - a00*a22*f1 + a02*a20*f1 + a10*a22*f0 - a12*a20*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
  		c2=(a00*a11*f2 - a01*a10*f2 - a00*a21*f1 + a01*a20*f1 + a10*a21*f0 - a11*a20*f0)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);

		g1.at(j1)=c0+c1*rr+c2*pow(rr,2);

//		cout << "test: j1: " << j1 << "\t" << g1.at(j1) << endl;
	
	}

	return(0);
}

/*
 * Smooth profiles
 * */
extern  double dfunc_4h(double r);
PetscErrorCode smoothprofile(vector<double> &f, int N_smooth)
{
	int i,j,k;
	int NElmts = f.size();
	vector<double> f0;
	f0.resize(NElmts);

	for (i=0;i<N_smooth;i++)
	{
		for (j=0;j<NElmts;j++)
			f0.at(j) = f.at(j);

		for (j=0;j<NElmts;j++)
		{
			int k1=max(j-2,0);
			int k2=min(j+3,NElmts);
			double sum_f=0.0;
			double sum_d=0.0;
			for (k=k1; k<k2; k++) 
			{			
					double weight=(double)(k-j);
					sum_f+=f0[k]*dfunc_4h(weight);
					sum_d+=dfunc_4h(weight);
			}
			f.at(j)=sum_f/(sum_d+1.e-19);
		}
	}
	
	return(0);
}

/**
 * Coordinate transformation 
 * tower-top and yaw position are assumed fixed 
 * (X, Y, Z) inertial frame, Y downwind, Z upward 
 * (xhat, yhat, zhat) rotates with the hub, zhat is aligned with the pitch axis of the blade, yhat is aligned with the Y axis, the angle between (xhat, yhat, zhat) and (X, Y, Z) is \Phi 
 * (x, y, z) describe the position of the blade, rotates /beta (the pitch angle) around the zhat axis 
 * (eta, xi) elastic principle, rotate twistpre+theta along elastic axis ea
 * elastic axis is computed as (u+lpi,v,w) 
 * lpi: distance from pitch axis to elastic axis 
 */

/**
 * Transformation from XYZ to xyzhat
 */ 
PetscErrorCode xyzhat_fromXYZ(double phi, double fX, double fY, double fZ, double *fxhat, double *fyhat, double *fzhat)
{
	*fxhat = cos(phi)*fX - sin(phi)*fZ;
	*fyhat = fY;
	*fzhat = sin(phi)*fX + cos(phi)*fZ;

	return(0);
}

/**
 * Transformation from xyzhat to XYZ
 */ 
PetscErrorCode XYZ_fromxyzhat(double phi, double fxhat, double fyhat, double fzhat, double *fX, double *fY, double *fZ)
{
	*fX =  cos(phi)*fxhat + sin(phi)*fzhat;
	*fY =  fyhat;
	*fZ = -sin(phi)*fxhat + cos(phi)*fzhat;
	return(0);
}

/**
 * Transformation from xyzhat to xyz
 */ 
PetscErrorCode xyz_fromxyzhat(double beta, double fxhat, double fyhat, double fzhat, double *fx, double *fy, double *fz)
{
	*fx =  cos(beta)*fxhat + sin(beta)*fyhat;
	*fy = -sin(beta)*fxhat + cos(beta)*fyhat;
	*fz =  fzhat;

	return(0);
}

/**
 * Transformation from xyz to xyzhat
 */ 
PetscErrorCode xyzhat_fromxyz(double beta, double fx, double fy, double fz, double *fxhat, double *fyhat, double *fzhat)
{
	*fxhat =  cos(beta)*fx - sin(beta)*fy;
	*fyhat =  sin(beta)*fx + cos(beta)*fy;
	*fzhat =  fz;
	return(0);
}

/**
 *  Compute spatial derivatives of blade property 
 *  Only the derivatives of the internal points are computed 
 */
PetscErrorCode computespatialderivatives0_bladestructure(IBMNodes *ibm)
{
	int i, j;
	double fm5, fm4, fm3, fm2, fm1, f, fp1, fp2, fp3, fp4, fp5;

	int ibi, iblade;

	for (ibi=0;ibi<NumberOfTurbines;ibi++) 
	
	for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {
//		cout << "\t **** "<< ibi << " " << iblade << endl;
		double h=ibm[ibi].bladestructure[iblade].dh;
		
		int NElmts = ibm[ibi].bladestructure[iblade].NumberofElements;	

		/* first derivatives */
		for (j=0;j<NElmts;j++) {

			if (j==0) 
			{
				f=ibm[ibi].bladestructure[iblade].mass[j];
				fp1=ibm[ibi].bladestructure[iblade].mass[j+1];
				fp2=ibm[ibi].bladestructure[iblade].mass[j+2];
				ibm[ibi].bladestructure[iblade].dmassdr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].E[j];
				fp1=ibm[ibi].bladestructure[iblade].E[j+1];
				fp2=ibm[ibi].bladestructure[iblade].E[j+2];
				ibm[ibi].bladestructure[iblade].dEdr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].G[j];
				fp1=ibm[ibi].bladestructure[iblade].G[j+1];
				fp2=ibm[ibi].bladestructure[iblade].G[j+2];
				ibm[ibi].bladestructure[iblade].dGdr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].Iksi[j];
				fp1=ibm[ibi].bladestructure[iblade].Iksi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Iksi[j+2];
				ibm[ibi].bladestructure[iblade].dIksidr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].Ieta[j];
				fp1=ibm[ibi].bladestructure[iblade].Ieta[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Ieta[j+2];
				ibm[ibi].bladestructure[iblade].dIetadr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].Icg[j];
				fp1=ibm[ibi].bladestructure[iblade].Icg[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Icg[j+2];
				ibm[ibi].bladestructure[iblade].dIcgdr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].twistpre[j];
				fp1=ibm[ibi].bladestructure[iblade].twistpre[j+1];
				fp2=ibm[ibi].bladestructure[iblade].twistpre[j+2];
				ibm[ibi].bladestructure[iblade].dtwistpredr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].lcg[j];
				fp1=ibm[ibi].bladestructure[iblade].lcg[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lcg[j+2];
				ibm[ibi].bladestructure[iblade].dlcgdr[j]=dfdr_forward(f, fp1, fp2, h);
        	
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				ibm[ibi].bladestructure[iblade].dlpidr[j]=dfdr_forward(f, fp1, fp2, h);
         	
				f=ibm[ibi].bladestructure[iblade].Ietaetaksi[j];
				fp1=ibm[ibi].bladestructure[iblade].Ietaetaksi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Ietaetaksi[j+2];
				ibm[ibi].bladestructure[iblade].dIetaetaksidr[j]=dfdr_forward(f, fp1, fp2, h);
          	
				f=ibm[ibi].bladestructure[iblade].Ietaksiksi[j];
				fp1=ibm[ibi].bladestructure[iblade].Ietaksiksi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Ietaksiksi[j+2];
				ibm[ibi].bladestructure[iblade].dIetaksiksidr[j]=dfdr_forward(f, fp1, fp2, h);
        		
				f=ibm[ibi].bladestructure[iblade].J[j];
				fp1=ibm[ibi].bladestructure[iblade].J[j+1];
				fp2=ibm[ibi].bladestructure[iblade].J[j+2];
				ibm[ibi].bladestructure[iblade].dJdr[j]=dfdr_forward(f, fp1, fp2, h);
			}	
			else if (j==NElmts-1)
			{
				f=ibm[ibi].bladestructure[iblade].mass[j];
				fm1=ibm[ibi].bladestructure[iblade].mass[j-1];
				fm2=ibm[ibi].bladestructure[iblade].mass[j-2];
				ibm[ibi].bladestructure[iblade].dmassdr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].E[j];
				fm1=ibm[ibi].bladestructure[iblade].E[j-1];
				fm2=ibm[ibi].bladestructure[iblade].E[j-2];
				ibm[ibi].bladestructure[iblade].dEdr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].G[j];
				fm1=ibm[ibi].bladestructure[iblade].G[j-1];
				fm2=ibm[ibi].bladestructure[iblade].G[j-2];
				ibm[ibi].bladestructure[iblade].dGdr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].Iksi[j];
				fm1=ibm[ibi].bladestructure[iblade].Iksi[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Iksi[j-2];
				ibm[ibi].bladestructure[iblade].dIksidr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].Ieta[j];
				fm1=ibm[ibi].bladestructure[iblade].Ieta[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Ieta[j-2];
				ibm[ibi].bladestructure[iblade].dIetadr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].Icg[j];
				fm1=ibm[ibi].bladestructure[iblade].Icg[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Icg[j-2];
				ibm[ibi].bladestructure[iblade].dIcgdr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].twistpre[j];
				fm1=ibm[ibi].bladestructure[iblade].twistpre[j-1];
				fm2=ibm[ibi].bladestructure[iblade].twistpre[j-2];
				ibm[ibi].bladestructure[iblade].dtwistpredr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].lcg[j];
				fm1=ibm[ibi].bladestructure[iblade].lcg[j-1];
				fm2=ibm[ibi].bladestructure[iblade].lcg[j-2];
				ibm[ibi].bladestructure[iblade].dlcgdr[j]=dfdr_backward(fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				ibm[ibi].bladestructure[iblade].dlpidr[j]=dfdr_backward(fm2, fm1, f, h);
         	
				f=ibm[ibi].bladestructure[iblade].Ietaetaksi[j];
				fm1=ibm[ibi].bladestructure[iblade].Ietaetaksi[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Ietaetaksi[j-2];
				ibm[ibi].bladestructure[iblade].dIetaetaksidr[j]=dfdr_backward(fm2, fm1, f, h);
          	
				f=ibm[ibi].bladestructure[iblade].Ietaksiksi[j];
				fm1=ibm[ibi].bladestructure[iblade].Ietaksiksi[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Ietaksiksi[j-2];
				ibm[ibi].bladestructure[iblade].dIetaksiksidr[j]=dfdr_backward(fm2, fm1, f, h);
        		
				f=ibm[ibi].bladestructure[iblade].J[j];
				fm1=ibm[ibi].bladestructure[iblade].J[j-1];
				fm2=ibm[ibi].bladestructure[iblade].J[j-2];
				ibm[ibi].bladestructure[iblade].dJdr[j]=dfdr_backward(fm2, fm1, f, h);

			} else
			{
				fm1=ibm[ibi].bladestructure[iblade].mass[j-1];
				fp1=ibm[ibi].bladestructure[iblade].mass[j+1];
				ibm[ibi].bladestructure[iblade].dmassdr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].E[j-1];
				fp1=ibm[ibi].bladestructure[iblade].E[j+1];
				ibm[ibi].bladestructure[iblade].dEdr[j]=dfdr(fm1, fp1, h);
        
//				cout << "E\t " << j << " " << fm1 << " " << fp1 << " " << fp1-fm1 << endl; 	
				fm1=ibm[ibi].bladestructure[iblade].G[j-1];
				fp1=ibm[ibi].bladestructure[iblade].G[j+1];
				ibm[ibi].bladestructure[iblade].dGdr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].Iksi[j-1];
				fp1=ibm[ibi].bladestructure[iblade].Iksi[j+1];
				ibm[ibi].bladestructure[iblade].dIksidr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].Ieta[j-1];
				fp1=ibm[ibi].bladestructure[iblade].Ieta[j+1];
				ibm[ibi].bladestructure[iblade].dIetadr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].Icg[j-1];
				fp1=ibm[ibi].bladestructure[iblade].Icg[j+1];
				ibm[ibi].bladestructure[iblade].dIcgdr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].twistpre[j-1];
				fp1=ibm[ibi].bladestructure[iblade].twistpre[j+1];
				ibm[ibi].bladestructure[iblade].dtwistpredr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].lcg[j-1];
				fp1=ibm[ibi].bladestructure[iblade].lcg[j+1];
				ibm[ibi].bladestructure[iblade].dlcgdr[j]=dfdr(fm1, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				ibm[ibi].bladestructure[iblade].dlpidr[j]=dfdr(fm1, fp1, h);
         	
				fm1=ibm[ibi].bladestructure[iblade].Ietaetaksi[j-1];
				fp1=ibm[ibi].bladestructure[iblade].Ietaetaksi[j+1];
				ibm[ibi].bladestructure[iblade].dIetaetaksidr[j]=dfdr(fm1, fp1, h);
          	
				fm1=ibm[ibi].bladestructure[iblade].Ietaksiksi[j-1];
				fp1=ibm[ibi].bladestructure[iblade].Ietaksiksi[j+1];
				ibm[ibi].bladestructure[iblade].dIetaksiksidr[j]=dfdr(fm1, fp1, h);
        		
				fm1=ibm[ibi].bladestructure[iblade].J[j-1];
				fp1=ibm[ibi].bladestructure[iblade].J[j+1];
				ibm[ibi].bladestructure[iblade].dJdr[j]=dfdr(fm1, fp1, h);
			}
		}
        	
		/* second derivatives */
		for (j=0;j<NElmts;j++) {
        	
			if (j==0) 
			{
				f=ibm[ibi].bladestructure[iblade].mass[j];
				fp1=ibm[ibi].bladestructure[iblade].mass[j+1];
				fp2=ibm[ibi].bladestructure[iblade].mass[j+2];
				fp3=ibm[ibi].bladestructure[iblade].mass[j+3];
				ibm[ibi].bladestructure[iblade].d2massdr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].E[j];
				fp1=ibm[ibi].bladestructure[iblade].E[j+1];
				fp2=ibm[ibi].bladestructure[iblade].E[j+2];
				fp3=ibm[ibi].bladestructure[iblade].E[j+3];
				ibm[ibi].bladestructure[iblade].d2Edr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].G[j];
				fp1=ibm[ibi].bladestructure[iblade].G[j+1];
				fp2=ibm[ibi].bladestructure[iblade].G[j+2];
				fp3=ibm[ibi].bladestructure[iblade].G[j+3];
				ibm[ibi].bladestructure[iblade].d2Gdr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].Iksi[j];
				fp1=ibm[ibi].bladestructure[iblade].Iksi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Iksi[j+2];
				fp3=ibm[ibi].bladestructure[iblade].Iksi[j+3];
				ibm[ibi].bladestructure[iblade].d2Iksidr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].Ieta[j];
				fp1=ibm[ibi].bladestructure[iblade].Ieta[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Ieta[j+2];
				fp3=ibm[ibi].bladestructure[iblade].Ieta[j+3];
				ibm[ibi].bladestructure[iblade].d2Ietadr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);

				f=ibm[ibi].bladestructure[iblade].Icg[j];
				fp1=ibm[ibi].bladestructure[iblade].Icg[j+1];
				fp2=ibm[ibi].bladestructure[iblade].Icg[j+2];
				fp3=ibm[ibi].bladestructure[iblade].Icg[j+3];
				ibm[ibi].bladestructure[iblade].d2Icgdr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].twistpre[j];
				fp1=ibm[ibi].bladestructure[iblade].twistpre[j+1];
				fp2=ibm[ibi].bladestructure[iblade].twistpre[j+2];
				fp3=ibm[ibi].bladestructure[iblade].twistpre[j+3];
				ibm[ibi].bladestructure[iblade].d2twistpredr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].lcg[j];
				fp1=ibm[ibi].bladestructure[iblade].lcg[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lcg[j+2];
				fp3=ibm[ibi].bladestructure[iblade].lcg[j+3];
				ibm[ibi].bladestructure[iblade].d2lcgdr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
        	
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				fp3=ibm[ibi].bladestructure[iblade].lpi[j+3];
				ibm[ibi].bladestructure[iblade].d2lpidr2[j]=d2fdr2_forward(f, fp1, fp2, fp3, h);
			}
			else if (j==NElmts-1)
			{
				f=ibm[ibi].bladestructure[iblade].mass[j];
				fm1=ibm[ibi].bladestructure[iblade].mass[j-1];
				fm2=ibm[ibi].bladestructure[iblade].mass[j-2];
				fm3=ibm[ibi].bladestructure[iblade].mass[j-3];
				ibm[ibi].bladestructure[iblade].d2massdr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].E[j];
				fm1=ibm[ibi].bladestructure[iblade].E[j-1];
				fm2=ibm[ibi].bladestructure[iblade].E[j-2];
				fm3=ibm[ibi].bladestructure[iblade].E[j-3];
				ibm[ibi].bladestructure[iblade].d2Edr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].G[j];
				fm1=ibm[ibi].bladestructure[iblade].G[j-1];
				fm2=ibm[ibi].bladestructure[iblade].G[j-2];
				fm3=ibm[ibi].bladestructure[iblade].G[j-3];
				ibm[ibi].bladestructure[iblade].d2Gdr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].Iksi[j];
				fm1=ibm[ibi].bladestructure[iblade].Iksi[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Iksi[j-2];
				fm3=ibm[ibi].bladestructure[iblade].Iksi[j-3];
				ibm[ibi].bladestructure[iblade].d2Iksidr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].Ieta[j];
				fm1=ibm[ibi].bladestructure[iblade].Ieta[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Ieta[j-2];
				fm3=ibm[ibi].bladestructure[iblade].Ieta[j-3];
				ibm[ibi].bladestructure[iblade].d2Ietadr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);

				f=ibm[ibi].bladestructure[iblade].Icg[j];
				fm1=ibm[ibi].bladestructure[iblade].Icg[j-1];
				fm2=ibm[ibi].bladestructure[iblade].Icg[j-2];
				fm3=ibm[ibi].bladestructure[iblade].Icg[j-3];
				ibm[ibi].bladestructure[iblade].d2Icgdr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].twistpre[j];
				fm1=ibm[ibi].bladestructure[iblade].twistpre[j-1];
				fm2=ibm[ibi].bladestructure[iblade].twistpre[j-2];
				fm3=ibm[ibi].bladestructure[iblade].twistpre[j-3];
				ibm[ibi].bladestructure[iblade].d2twistpredr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].lcg[j];
				fm1=ibm[ibi].bladestructure[iblade].lcg[j-1];
				fm2=ibm[ibi].bladestructure[iblade].lcg[j-2];
				fm3=ibm[ibi].bladestructure[iblade].lcg[j-3];
				ibm[ibi].bladestructure[iblade].d2lcgdr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
        	
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm3=ibm[ibi].bladestructure[iblade].lpi[j-3];
				ibm[ibi].bladestructure[iblade].d2lpidr2[j]=d2fdr2_backward(fm3, fm2, fm1, f, h);
			}
			else 
			{	
				fm1=ibm[ibi].bladestructure[iblade].mass[j-1];
				f=ibm[ibi].bladestructure[iblade].mass[j];
				fp1=ibm[ibi].bladestructure[iblade].mass[j+1];
				ibm[ibi].bladestructure[iblade].d2massdr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].E[j-1];
				f=ibm[ibi].bladestructure[iblade].E[j];
				fp1=ibm[ibi].bladestructure[iblade].E[j+1];
				ibm[ibi].bladestructure[iblade].d2Edr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].G[j-1];
				f=ibm[ibi].bladestructure[iblade].G[j];
				fp1=ibm[ibi].bladestructure[iblade].G[j+1];
				ibm[ibi].bladestructure[iblade].d2Gdr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].Iksi[j-1];
				f=ibm[ibi].bladestructure[iblade].Iksi[j];
				fp1=ibm[ibi].bladestructure[iblade].Iksi[j+1];
				ibm[ibi].bladestructure[iblade].d2Iksidr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].Ieta[j-1];
				f=ibm[ibi].bladestructure[iblade].Ieta[j];
				fp1=ibm[ibi].bladestructure[iblade].Ieta[j+1];
				ibm[ibi].bladestructure[iblade].d2Ietadr2[j]=d2fdr2(fm1, f, fp1, h);

				fm1=ibm[ibi].bladestructure[iblade].Icg[j-1];
				f=ibm[ibi].bladestructure[iblade].Icg[j];
				fp1=ibm[ibi].bladestructure[iblade].Icg[j+1];
				ibm[ibi].bladestructure[iblade].d2Icgdr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].twistpre[j-1];
				f=ibm[ibi].bladestructure[iblade].twistpre[j];
				fp1=ibm[ibi].bladestructure[iblade].twistpre[j+1];
				ibm[ibi].bladestructure[iblade].d2twistpredr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].lcg[j-1];
				f=ibm[ibi].bladestructure[iblade].lcg[j];
				fp1=ibm[ibi].bladestructure[iblade].lcg[j+1];
				ibm[ibi].bladestructure[iblade].d2lcgdr2[j]=d2fdr2(fm1, f, fp1, h);
        	
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				ibm[ibi].bladestructure[iblade].d2lpidr2[j]=d2fdr2(fm1, f, fp1, h);
        		}	
		}


	
	        /* third derivatives */
		for (j=0;j<NElmts;j++) 
		{
        		if (j==0) 
			{
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				fp3=ibm[ibi].bladestructure[iblade].lpi[j+3];
				fp4=ibm[ibi].bladestructure[iblade].lpi[j+4];
				ibm[ibi].bladestructure[iblade].d3lpidr3[j]=d3fdr3_forward2(f, fp1, fp2, fp3, fp4, h);

			}
			else if (j==1)
			{
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				fp3=ibm[ibi].bladestructure[iblade].lpi[j+3];
				ibm[ibi].bladestructure[iblade].d3lpidr3[j]=d3fdr3_forward1(fm1, f, fp1, fp2, fp3, h);
			}
			else if (j==NElmts-1)
			{
				fm4=ibm[ibi].bladestructure[iblade].lpi[j-4];
				fm3=ibm[ibi].bladestructure[iblade].lpi[j-3];
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				ibm[ibi].bladestructure[iblade].d3lpidr3[j]=d3fdr3_backward2(fm4, fm3, fm2, fm1, f, h);
			}
			else if (j==NElmts-2)
			{
				fm3=ibm[ibi].bladestructure[iblade].lpi[j-3];
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				ibm[ibi].bladestructure[iblade].d3lpidr3[j]=d3fdr3_backward1(fm3, fm2, fm1, f, fp1, h);
			}
			else
			{
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
	        		fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				ibm[ibi].bladestructure[iblade].d3lpidr3[j]=d3fdr3(fm2, fm1, f, fp1, fp2, h); 
			}
		}
                
		/* fourth derivatives */
		for (j=0;j<NElmts;j++) 
		{
			if (j==0)
			{
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
		        	fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				fp3=ibm[ibi].bladestructure[iblade].lpi[j+3];
				fp4=ibm[ibi].bladestructure[iblade].lpi[j+4];
				fp5=ibm[ibi].bladestructure[iblade].lpi[j+5];
				ibm[ibi].bladestructure[iblade].d4lpidr4[j]=d4fdr4_forward2(f, fp1, fp2, fp3, fp4, fp5, h);
			}
			else if (j==1)
			{
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
		        	fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				fp3=ibm[ibi].bladestructure[iblade].lpi[j+3];
				fp4=ibm[ibi].bladestructure[iblade].lpi[j+4];
				ibm[ibi].bladestructure[iblade].d4lpidr4[j]=d4fdr4_forward1(fm1, f, fp1, fp2, fp3, fp4, h);
			}
			else if (j==NElmts-1)
			{
		        	fm5=ibm[ibi].bladestructure[iblade].lpi[j-5];
		        	fm4=ibm[ibi].bladestructure[iblade].lpi[j-4];
		        	fm3=ibm[ibi].bladestructure[iblade].lpi[j-3];
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
		        	ibm[ibi].bladestructure[iblade].d4lpidr4[j]=d4fdr4_backward2(fm5, fm4, fm3, fm2, fm1, f, h);
			}
			else if (j==NElmts-2)
			{
		        	fm4=ibm[ibi].bladestructure[iblade].lpi[j-4];
		        	fm3=ibm[ibi].bladestructure[iblade].lpi[j-3];
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
				f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
		        	ibm[ibi].bladestructure[iblade].d4lpidr4[j]=d4fdr4_backward1(fm4, fm3, fm2, fm1, f, fp1, h);
			}
			else 
			{
				fm2=ibm[ibi].bladestructure[iblade].lpi[j-2];
				fm1=ibm[ibi].bladestructure[iblade].lpi[j-1];
		        	f=ibm[ibi].bladestructure[iblade].lpi[j];
				fp1=ibm[ibi].bladestructure[iblade].lpi[j+1];
				fp2=ibm[ibi].bladestructure[iblade].lpi[j+2];
				ibm[ibi].bladestructure[iblade].d4lpidr4[j]=d4fdr4(fm2, fm1, f, fp1, fp2, h); 
			}
		}
	}        	

	int rank;	
  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (!rank) {
		FILE *fd;
    	    	char filen[80]; 
        
		for (ibi=0;ibi<NumberOfTurbines;ibi++) {
		for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {
    	    	        sprintf(filen,"%s/bladeproperty_1stderivative%3.3d_%3.3d.out" , path, ibi, iblade);
    	    	        fd = fopen(filen, "w"); 
    	    		PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = r dmassdr dEdr dGdr dIksidr dIetadr dIcgdr dtwistpredr dlcgdr dlpidr dIetaetaksidr dIetaksiksidr dJdr \n");
    	    		
			for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le %le %le %le %le %le %le %le %le %le %le %le \n", ibm[ibi].bladestructure[iblade].r[i], 
                                	ibm[ibi].bladestructure[iblade].dmassdr[i], ibm[ibi].bladestructure[iblade].dEdr[i], ibm[ibi].bladestructure[iblade].dGdr[i], 
					ibm[ibi].bladestructure[iblade].dIksidr[i], ibm[ibi].bladestructure[iblade].dIetadr[i], ibm[ibi].bladestructure[iblade].dIcgdr[i], 
					ibm[ibi].bladestructure[iblade].dtwistpredr[i], ibm[ibi].bladestructure[iblade].dlcgdr[i], ibm[ibi].bladestructure[iblade].dlpidr[i] , 
					ibm[ibi].bladestructure[iblade].dIetaetaksidr[i], ibm[ibi].bladestructure[iblade].dIetaksiksidr[i], ibm[ibi].bladestructure[iblade].dJdr[i] );
			}
			fclose(fd);
        
    	    	        sprintf(filen,"%s/bladeproperty_2ndderivative%3.3d_%3.3d.out" , path, ibi, iblade);
    	    	        fd = fopen(filen, "w"); 
    	    		PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = r d2massdr2 d2Edr2 d2Gdr2 d2Iksidr2 d2Ietadr2 d2Icgdr2 d2twistpredr2 d2lcgdr2 d2lpidr2 \n");
    	    		
			for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le %le %le %le %le %le %le %le %le \n", ibm[ibi].bladestructure[iblade].r[i], 
                                	ibm[ibi].bladestructure[iblade].d2massdr2[i], ibm[ibi].bladestructure[iblade].d2Edr2[i], ibm[ibi].bladestructure[iblade].d2Gdr2[i], 
					ibm[ibi].bladestructure[iblade].d2Iksidr2[i], ibm[ibi].bladestructure[iblade].d2Ietadr2[i], ibm[ibi].bladestructure[iblade].d2Icgdr2[i], 
					ibm[ibi].bladestructure[iblade].d2twistpredr2[i], ibm[ibi].bladestructure[iblade].d2lcgdr2[i], ibm[ibi].bladestructure[iblade].d2lpidr2[i]);
			}
			fclose(fd);
        
    	    	        sprintf(filen,"%s/bladeproperty_3rdderivative%3.3d_%3.3d.out" , path, ibi, iblade);
    	    	        fd = fopen(filen, "w"); 
    	    		PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = r d3lpidr3 \n");
    	    		
			for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le \n", ibm[ibi].bladestructure[iblade].r[i], ibm[ibi].bladestructure[iblade].d3lpidr3[i]);
			}
			fclose(fd);
        
    	    	        sprintf(filen,"%s/bladeproperty_4thderivative%3.3d_%3.3d.out" , path, ibi, iblade);
    	    	        fd = fopen(filen, "w"); 
    	    		PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = r d4lpidr4 \n");
    	    		
			for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le \n", ibm[ibi].bladestructure[iblade].r[i], ibm[ibi].bladestructure[iblade].d4lpidr4[i]);
			}
			fclose(fd);
		}
		}
	}
	return(0);
}       	



/** 
 * Read blade property, use least squares to compute the blade material properties at node locations, allocate variables 
 */
PetscErrorCode init_bladestructure(IBMNodes *ibm)
{

	int	rank;
	int	i,j,ibi;

	int 	iblade;

  	char string[80];

    	PetscPrintf(PETSC_COMM_WORLD, "READ blade structure property! \n");

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	for (ibi=0;ibi<NumberOfTurbines;ibi++) {
	for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {

  		if(!rank) { // root processor read in the data
    			FILE *fd;
    			char filen[80]; 

    		        sprintf(filen,"%s/bladeproperty%3.3d_%3.3d.inp" , path, ibi, iblade);
    		        fd = fopen(filen, "r"); 
		        if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open %s", filen), exit(0);
                 	
		        if (fd) {
                 	
    				PetscPrintf(PETSC_COMM_SELF, "iturbine=%d iblade=%d \n", ibi, iblade);
	
				// read number of elements 
      		        	fgets(string, 80, fd);
	//			std::cout << string << std::endl;
      		 		fscanf(fd, "%i \n",&(ibm[ibi].bladestructure[iblade].NumberofElements));

				// read bctype
      			        fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			        fscanf(fd, "%i \n",&(ibm[ibi].bladestructure[iblade].bctype));
                                
				// read u boundary values 
      			 	fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			        fscanf(fd, "%le %le %le %le \n",&(ibm[ibi].bladestructure[iblade].u_bc0), &(ibm[ibi].bladestructure[iblade].dudr_bc0), &(ibm[ibi].bladestructure[iblade].d2udr2_bc0), &(ibm[ibi].bladestructure[iblade].d3udr3_bc0));
                         	
				// read u boundary values 
      			        fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			 	fscanf(fd, "%le %le %le %le \n",&(ibm[ibi].bladestructure[iblade].u_bc1), &(ibm[ibi].bladestructure[iblade].dudr_bc1), &(ibm[ibi].bladestructure[iblade].d2udr2_bc1), &(ibm[ibi].bladestructure[iblade].d3udr3_bc1));
                                
				// read v boundary values 
      			        fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			        fscanf(fd, "%le %le %le %le \n",&(ibm[ibi].bladestructure[iblade].v_bc0), &(ibm[ibi].bladestructure[iblade].dvdr_bc0), &(ibm[ibi].bladestructure[iblade].d2vdr2_bc0), &(ibm[ibi].bladestructure[iblade].d3vdr3_bc0));
                                
				// read v boundary values 
      			 	fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			        fscanf(fd, "%le %le %le %le \n",&(ibm[ibi].bladestructure[iblade].v_bc1), &(ibm[ibi].bladestructure[iblade].dvdr_bc1), &(ibm[ibi].bladestructure[iblade].d2vdr2_bc1), &(ibm[ibi].bladestructure[iblade].d3vdr3_bc1));
                         	
				// read theta boundary values 
      			        fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			 	fscanf(fd, "%le %le \n",&(ibm[ibi].bladestructure[iblade].theta_bc0), &(ibm[ibi].bladestructure[iblade].dthetadr_bc0));
                                
				// read theta boundary values 
      			        fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			        fscanf(fd, "%le %le \n",&(ibm[ibi].bladestructure[iblade].theta_bc1), &(ibm[ibi].bladestructure[iblade].dthetadr_bc1));
                                
                                
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].NumberofElements), 1, MPI_INT, 0, PETSC_COMM_WORLD); 
                                
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].bctype), 1, MPI_INT, 0, PETSC_COMM_WORLD); 
                                
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].u_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dudr_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2udr2_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3udr3_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                         	
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].u_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dudr_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2udr2_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3udr3_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                                
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].v_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dvdr_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2vdr2_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3vdr3_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                                
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].v_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dvdr_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2vdr2_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3vdr3_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                         	
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].theta_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dthetadr_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                                
			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].theta_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dthetadr_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                                

      			        fgets(string, 80, fd);
//				std::cout << string << std::endl;
      			        fscanf(fd, "%i \n ",&(ibm[ibi].bladestructure[iblade].NumberofInp));
                                

      			 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].NumberofInp), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                               
    				PetscPrintf(PETSC_COMM_SELF, "\t %d \n", ibm[ibi].bladestructure[iblade].NumberofElements);
    				PetscPrintf(PETSC_COMM_SELF, "\t %d \n", ibm[ibi].bladestructure[iblade].bctype);
    				PetscPrintf(PETSC_COMM_SELF, "\t %f %f %f %f\n", ibm[ibi].bladestructure[iblade].u_bc0, ibm[ibi].bladestructure[iblade].dudr_bc0, ibm[ibi].bladestructure[iblade].d2udr2_bc0, ibm[ibi].bladestructure[iblade].d3udr3_bc0);
    				PetscPrintf(PETSC_COMM_SELF, "\t %f %f %f %f\n", ibm[ibi].bladestructure[iblade].u_bc1, ibm[ibi].bladestructure[iblade].dudr_bc1, ibm[ibi].bladestructure[iblade].d2udr2_bc1, ibm[ibi].bladestructure[iblade].d3udr3_bc1);
    				PetscPrintf(PETSC_COMM_SELF, "\t %f %f %f %f\n", ibm[ibi].bladestructure[iblade].v_bc0, ibm[ibi].bladestructure[iblade].dvdr_bc0, ibm[ibi].bladestructure[iblade].d2vdr2_bc0, ibm[ibi].bladestructure[iblade].d3vdr3_bc0);
    				PetscPrintf(PETSC_COMM_SELF, "\t %f %f %f %f\n", ibm[ibi].bladestructure[iblade].v_bc1, ibm[ibi].bladestructure[iblade].dvdr_bc1, ibm[ibi].bladestructure[iblade].d2vdr2_bc1, ibm[ibi].bladestructure[iblade].d3vdr3_bc1);
    				PetscPrintf(PETSC_COMM_SELF, "\t %f %f \n", ibm[ibi].bladestructure[iblade].theta_bc0, ibm[ibi].bladestructure[iblade].dthetadr_bc0);
    				PetscPrintf(PETSC_COMM_SELF, "\t %f %f \n", ibm[ibi].bladestructure[iblade].theta_bc1, ibm[ibi].bladestructure[iblade].dthetadr_bc1);
    				PetscPrintf(PETSC_COMM_SELF, "\t %d \n", ibm[ibi].bladestructure[iblade].NumberofInp);


                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].r_inp));  
                                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].mass_inp));  
                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].E_inp));  
                                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].G_inp));  
                                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Iksi_inp));  
                                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ieta_inp));  
                                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Icg_inp));  
                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].twistpre_inp));  
                                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].lcg_inp));  
                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].lpi_inp));  
                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ietaetaksi_inp));  
                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ietaksiksi_inp));  
                         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].J_inp));  
                                
    				PetscPrintf(PETSC_COMM_SELF, "\t Read Input blade structure \n");

			 	for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofInp; i++) {
					ibm[ibi].bladestructure[iblade].E_inp[i] = 0.0;
					ibm[ibi].bladestructure[iblade].G_inp[i] = 0.0;
				}
			 	for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofInp; i++) {
					double r_inp, mass_inp, E_inp, G_inp, Iksi_inp, Ieta_inp, Icg_inp, twistpre_inp, lcg_inp, lpi_inp, Ietaetaksi_inp, Ietaksiksi_inp, J_inp;
			        	fscanf(fd, "%le %le %le %le %le %le %le %le %le %le %le %le %le \n", &r_inp, &mass_inp, &E_inp, &G_inp, &Iksi_inp, &Ieta_inp, &Icg_inp, &twistpre_inp, &lcg_inp, &lpi_inp, &Ietaetaksi_inp, &Ietaksiksi_inp, &J_inp );

					twistpre_inp *= -(acos(-1.0)/180.0);

					ibm[ibi].bladestructure[iblade].r_inp[i] = r_inp;
					ibm[ibi].bladestructure[iblade].mass_inp[i] = mass_inp; 
					ibm[ibi].bladestructure[iblade].E_inp[i] = E_inp;
					ibm[ibi].bladestructure[iblade].G_inp[i] = G_inp;
					ibm[ibi].bladestructure[iblade].Iksi_inp[i] = Iksi_inp; 
					ibm[ibi].bladestructure[iblade].Ieta_inp[i] = Ieta_inp;
					ibm[ibi].bladestructure[iblade].Icg_inp[i] = Icg_inp;
					ibm[ibi].bladestructure[iblade].twistpre_inp[i] = twistpre_inp; 
					ibm[ibi].bladestructure[iblade].lcg_inp[i] = lcg_inp;
					ibm[ibi].bladestructure[iblade].lpi_inp[i] = lpi_inp;
					ibm[ibi].bladestructure[iblade].Ietaetaksi_inp[i] = Ietaetaksi_inp;
					ibm[ibi].bladestructure[iblade].Ietaksiksi_inp[i] = Ietaksiksi_inp;
					ibm[ibi].bladestructure[iblade].J_inp[i] = J_inp;
		

//					printf("%le %le %le %le %le %le %le %le %le %le %le %le %le \n", (ibm[ibi].bladestructure[iblade].r_inp[i]), (ibm[ibi].bladestructure[iblade].mass_inp[i]), (ibm[ibi].bladestructure[iblade].E_inp[i]), (ibm[ibi].bladestructure[iblade].G_inp[i]), (ibm[ibi].bladestructure[iblade].Iksi_inp[i]), (ibm[ibi].bladestructure[iblade].Ieta_inp[i]), (ibm[ibi].bladestructure[iblade].Icg_inp[i]), (ibm[ibi].bladestructure[iblade].twistpre_inp[i]), (ibm[ibi].bladestructure[iblade].lcg_inp[i]), (ibm[ibi].bladestructure[iblade].lpi_inp[i]) , (ibm[ibi].bladestructure[iblade].Ietaetaksi_inp[i]), (ibm[ibi].bladestructure[iblade].Ietaksiksi_inp[i]), (ibm[ibi].bladestructure[iblade].J_inp[i]) );

			        } 
 
    				PetscPrintf(PETSC_COMM_SELF, "\t Broadcast Input blade structure \n");

			        MPI_Bcast(ibm[ibi].bladestructure[iblade].r_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(ibm[ibi].bladestructure[iblade].mass_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].E_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(ibm[ibi].bladestructure[iblade].G_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].Iksi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			 	MPI_Bcast(ibm[ibi].bladestructure[iblade].Ieta_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].Icg_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].twistpre_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].lcg_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].lpi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].Ietaetaksi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].Ietaksiksi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			        MPI_Bcast(ibm[ibi].bladestructure[iblade].J_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
//				MPI_Barrier(PETSC_COMM_WORLD);
			}               
		        PetscPrintf(PETSC_COMM_SELF, "Close %s \n", filen);
			fclose(fd);
		}
		else 
		{

		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].NumberofElements), 1, MPI_INT, 0, PETSC_COMM_WORLD); 
                        
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].bctype), 1, MPI_INT, 0, PETSC_COMM_WORLD); 
                 	
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].u_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dudr_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2udr2_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3udr3_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 

		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].u_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dudr_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2udr2_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3udr3_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                 	
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].v_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dvdr_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2vdr2_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3vdr3_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                        
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].v_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dvdr_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d2vdr2_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].d3vdr3_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                        
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].theta_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dthetadr_bc0), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                 	
		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].theta_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(&(ibm[ibi].bladestructure[iblade].dthetadr_bc1), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                 	
      		        MPI_Bcast(&(ibm[ibi].bladestructure[iblade].NumberofInp), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                 	
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].r_inp));  
                 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].mass_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].E_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].G_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Iksi_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ieta_inp));  
                 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Icg_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].twistpre_inp));  
                 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].lcg_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].lpi_inp));  
                        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ietaetaksi_inp));  
                       	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ietaksiksi_inp));  
                       	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofInp)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].J_inp));  
                   	
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].r_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].mass_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].E_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].G_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(ibm[ibi].bladestructure[iblade].Iksi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].Ieta_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(ibm[ibi].bladestructure[iblade].Icg_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].twistpre_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		 	MPI_Bcast(ibm[ibi].bladestructure[iblade].lcg_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].lpi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
			MPI_Bcast(ibm[ibi].bladestructure[iblade].Ietaetaksi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].Ietaksiksi_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].J_inp, ibm[ibi].bladestructure[iblade].NumberofInp, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                
	//		MPI_Barrier(PETSC_COMM_WORLD);
		}
	}
	}

	MPI_Barrier(PETSC_COMM_WORLD);


	if (!rank) {

		FILE *fd;
    	    	char filen[80]; 
        
		for (ibi=0;ibi<NumberOfTurbines;ibi++) {
		for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {
    	    	        sprintf(filen,"%s/bladeproperty_input%3.3d_%3.3d.out" , path, ibi, iblade);
    	    	        fd = fopen(filen, "w"); 
    	    		PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = r mass E G Iksi Ieta Icg twistpre lcg lpi Ietaetaksi Ietaksiksi J \n");
    	    		
			for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofInp; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le %le %le %le %le %le %le %le %le %le %le %le \n", ibm[ibi].bladestructure[iblade].r_inp[i], ibm[ibi].bladestructure[iblade].mass_inp[i], ibm[ibi].bladestructure[iblade].E_inp[i], ibm[ibi].bladestructure[iblade].G_inp[i], ibm[ibi].bladestructure[iblade].Iksi_inp[i], ibm[ibi].bladestructure[iblade].Ieta_inp[i], ibm[ibi].bladestructure[iblade].Icg_inp[i], ibm[ibi].bladestructure[iblade].twistpre_inp[i], ibm[ibi].bladestructure[iblade].lcg_inp[i], ibm[ibi].bladestructure[iblade].lpi_inp[i], ibm[ibi].bladestructure[iblade].Ietaetaksi_inp[i], ibm[ibi].bladestructure[iblade].Ietaksiksi_inp[i], ibm[ibi].bladestructure[iblade].J_inp[i] );
			}
			fclose(fd);
        
		}
		}
	}

	MPI_Barrier(PETSC_COMM_WORLD);

    	PetscPrintf(PETSC_COMM_WORLD, "Allocate variables! \n");

	/* allocate variables */
	for (ibi=0;ibi<NumberOfTurbines;ibi++) {
	for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {

		int NElmts = ibm[ibi].bladestructure[iblade].NumberofElements;
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].r));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].mass));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].E));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].G));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Iksi));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ieta));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Icg));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].twistpre));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].lcg));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].lpi));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ietaetaksi));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].Ietaksiksi));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].J));  

	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dmassdr));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dEdr));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dGdr));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dIksidr));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dIetadr));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dIcgdr));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dtwistpredr));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dlcgdr));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dlpidr));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dIetaetaksidr));  
        	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dIetaksiksidr));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dJdr));  

                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2massdr2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2Edr2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2Gdr2));  
 	       	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2Iksidr2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2Ietadr2));  
 	       	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2Icgdr2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2twistpredr2));  
 	       	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2lcgdr2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d2lpidr2));  

                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d3lpidr3));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].d4lpidr4));  
 	
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].u));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].v));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].theta));  
 	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].u_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].v_o));  
 	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].theta_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].u_om1));  
 	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].v_om1));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].theta_om1));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].u_om2));  
 	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].v_om2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].theta_om2));  
 
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dudt));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dvdt));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dthetadt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dudt_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dvdt_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dthetadt_o));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dudt_om1));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dvdt_om1));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dthetadt_om1));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dudt_om2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dvdt_om2));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].dthetadt_om2));  

	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].fu));  
         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].fv));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].ftheta));  
         	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].M));  
                
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].fu_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].fv_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].ftheta_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].M_o));  
                
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsu));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsv));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhstheta));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdudt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdvdt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdthetadt));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsu_o));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsv_o));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhstheta_o));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdudt_o));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdvdt_o));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdthetadt_o));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsu_om1));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsv_om1));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhstheta_om1));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdudt_om1));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdvdt_om1));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdthetadt_om1));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsu_om2));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsv_om2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhstheta_om2));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdudt_om2));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdvdt_om2));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].rhsdthetadt_om2));  
       
 	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].err_u));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].err_v));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].err_theta));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].err_dudt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].err_dvdt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].err_dthetadt));  

  	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].ucg));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].vcg));  
 	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].ucg_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].vcg_o));  
                PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].ucg_om1));  
 	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].vcg_om1));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].resu));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].resv));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].restheta));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].resdudt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].resdvdt));  
	        PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].resdthetadt));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_dudr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d2udr2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d3udr3));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d4udr4));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_dvdr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d2vdr2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d3vdr3));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d4vdr4));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_dthetadr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefu_d2thetadr2));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_dudr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d2udr2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d3udr3));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d4udr4));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_dvdr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d2vdr2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d3vdr3));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d4vdr4));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_dthetadr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefv_d2thetadr2));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_dudr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d2udr2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d3udr3));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d4udr4));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_dvdr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d2vdr2));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d3vdr3));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d4vdr4));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_dthetadr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coeftheta_d2thetadr2));  

	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefdudt_ddt_dudr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefdudt_ddt_dvdr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefdvdt_ddt_dudr));  
	 	PetscMalloc((ibm[ibi].bladestructure[iblade].NumberofElements)*sizeof(PetscReal), &(ibm[ibi].bladestructure[iblade].coefdvdt_ddt_dvdr));  

		if (torsion_turbinestructure) 
		{
			ibm[ibi].bladestructure[iblade].indx = new int[6*NElmts];
			ibm[ibi].bladestructure[iblade].b = new double[6*NElmts];
			ibm[ibi].bladestructure[iblade].bb = new double[6*NElmts];
			ibm[ibi].bladestructure[iblade].A = new double*[6*NElmts];
			ibm[ibi].bladestructure[iblade].A_LU = new double*[6*NElmts];
			for (i=0; i<6*NElmts; i++)
			{
				ibm[ibi].bladestructure[iblade].A[i] = new double[6*NElmts];
				ibm[ibi].bladestructure[iblade].A_LU[i] = new double[6*NElmts];
			}
		} 
		else
		{
			ibm[ibi].bladestructure[iblade].indx = new int[4*NElmts];
			ibm[ibi].bladestructure[iblade].b = new double[4*NElmts];
			ibm[ibi].bladestructure[iblade].bb = new double[4*NElmts];
			ibm[ibi].bladestructure[iblade].A = new double*[4*NElmts];
			ibm[ibi].bladestructure[iblade].A_LU = new double*[4*NElmts];
			for (i=0; i<4*NElmts; i++)
			{
				ibm[ibi].bladestructure[iblade].A[i] = new double[4*NElmts];
				ibm[ibi].bladestructure[iblade].A_LU[i] = new double[4*NElmts];
			}
		}

		/* initialize */
	    	PetscPrintf(PETSC_COMM_WORLD, "Initialize variables! \n");
		for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
		{
			ibm[ibi].bladestructure[iblade].u[i]=0.0;
			ibm[ibi].bladestructure[iblade].u_o[i]=0.0;
			ibm[ibi].bladestructure[iblade].u_om1[i]=0.0;
			ibm[ibi].bladestructure[iblade].u_om2[i]=0.0;

			ibm[ibi].bladestructure[iblade].v[i]=0.0;
			ibm[ibi].bladestructure[iblade].v_o[i]=0.0;
			ibm[ibi].bladestructure[iblade].v_om1[i]=0.0;
			ibm[ibi].bladestructure[iblade].v_om2[i]=0.0;

			ibm[ibi].bladestructure[iblade].theta[i]=0.0;
			ibm[ibi].bladestructure[iblade].theta_o[i]=0.0;
			ibm[ibi].bladestructure[iblade].theta_om1[i]=0.0;
			ibm[ibi].bladestructure[iblade].theta_om2[i]=0.0;

			ibm[ibi].bladestructure[iblade].fu[i]=0.0;
			ibm[ibi].bladestructure[iblade].fv[i]=0.0;
			ibm[ibi].bladestructure[iblade].M[i]=0.0;

			ibm[ibi].bladestructure[iblade].err_u[i]=0.0;
			ibm[ibi].bladestructure[iblade].err_v[i]=0.0;
			ibm[ibi].bladestructure[iblade].err_theta[i]=0.0;

			ibm[ibi].bladestructure[iblade].err_dudt[i]=0.0;
			ibm[ibi].bladestructure[iblade].err_dvdt[i]=0.0;
			ibm[ibi].bladestructure[iblade].err_dthetadt[i]=0.0;

			ibm[ibi].bladestructure[iblade].ucg[i]=0.0;
			ibm[ibi].bladestructure[iblade].ucg_o[i]=0.0;
			ibm[ibi].bladestructure[iblade].ucg_om1[i]=0.0;

			ibm[ibi].bladestructure[iblade].vcg[i]=0.0;
			ibm[ibi].bladestructure[iblade].vcg_o[i]=0.0;
			ibm[ibi].bladestructure[iblade].vcg_om1[i]=0.0;
		}
         	
	    	PetscPrintf(PETSC_COMM_WORLD, "Interpolate blade property! \n");
	        /* setup blade property on the structure grid nodes */
         	
	        ibm[ibi].bladestructure[iblade].rstart = ibm[ibi].bladestructure[iblade].r_inp[0];
	 	ibm[ibi].bladestructure[iblade].rend = ibm[ibi].bladestructure[iblade].r_inp[ ibm[ibi].bladestructure[iblade].NumberofInp-1];
                
	        ibm[ibi].bladestructure[iblade].dh=(ibm[ibi].bladestructure[iblade].rend-ibm[ibi].bladestructure[iblade].rstart)/(((double)ibm[ibi].bladestructure[iblade].NumberofElements)-1.0);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) {
	        	ibm[ibi].bladestructure[iblade].r[i]=((double)i)*ibm[ibi].bladestructure[iblade].dh;
	 	}
               
	 	vector<double> r0, g0, r1, g1;
               
		r0.resize(ibm[ibi].bladestructure[iblade].NumberofInp); 
		g0.resize(ibm[ibi].bladestructure[iblade].NumberofInp); 
		r1.resize(ibm[ibi].bladestructure[iblade].NumberofElements); 
		g1.resize(ibm[ibi].bladestructure[iblade].NumberofElements); 
    ibm[ibi].bladestructure[iblade].hubrad = ibm[ibi].r_nacelle;

	 	for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) r0.at(i)=ibm[ibi].bladestructure[iblade].r_inp[i];
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) r1.at(i)=ibm[ibi].bladestructure[iblade].r[i];
         	
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].mass_inp[i];
	 	_interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].mass[i]=g1.at(i);

	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].E_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].E[i]=g1.at(i);
         	
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].G_inp[i];
	 	_interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].G[i]=g1.at(i);
         	
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].Iksi_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].Iksi[i]=g1.at(i);
                
	 	for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].Ieta_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	 	for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].Ieta[i]=g1.at(i);
                
	 	for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].Icg_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].Icg[i]=g1.at(i);
                
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].twistpre_inp[i];
	 	_interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].twistpre[i]=g1.at(i);
         	
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].lcg_inp[i];
	 	_interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].lcg[i]=g1.at(i);

	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].lpi_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].lpi[i]=g1.at(i);

	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].Ietaetaksi_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].Ietaetaksi[i]=g1.at(i);
 
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].Ietaksiksi_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].Ietaksiksi[i]=g1.at(i);
 
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofInp;i++) g0.at(i)=ibm[ibi].bladestructure[iblade].J_inp[i];
	        _interpolation_leastsquares(r0, g0, r1, g1, ibm[ibi].bladestructure[iblade].NumberofInp, ibm[ibi].bladestructure[iblade].NumberofElements);
		smoothprofile(g1,Nsmooth);
	        for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) ibm[ibi].bladestructure[iblade].J[i]=g1.at(i);
        	
	}
	}

        if (ti == tistart){
             Init_TurbineStructure(ibm);
        }

	if (!rank) {

		FILE *fd;
    	    	char filen[80]; 
        
		for (ibi=0;ibi<NumberOfTurbines;ibi++) {
		for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {
       
    	    	        sprintf(filen,"%s/bladeproperty_interpolation%3.3d_%3.3d.out" , path, ibi, iblade);
    	    	        fd = fopen(filen, "w"); 
    	    		PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = r mass E G Iksi Ieta Icg twistpre lcg lpi Ietaetaksi Ietaksiksi J \n");
    	    		
			for(i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le %le %le %le %le %le %le %le %le %le %le %le \n", ibm[ibi].bladestructure[iblade].r[i], ibm[ibi].bladestructure[iblade].mass[i], ibm[ibi].bladestructure[iblade].E[i], ibm[ibi].bladestructure[iblade].G[i], ibm[ibi].bladestructure[iblade].Iksi[i], ibm[ibi].bladestructure[iblade].Ieta[i], ibm[ibi].bladestructure[iblade].Icg[i], ibm[ibi].bladestructure[iblade].twistpre[i], ibm[ibi].bladestructure[iblade].lcg[i], ibm[ibi].bladestructure[iblade].lpi[i] , ibm[ibi].bladestructure[iblade].Ietaetaksi[i], ibm[ibi].bladestructure[iblade].Ietaksiksi[i], ibm[ibi].bladestructure[iblade].J[i] );
			}
			fclose(fd);
		}
		}
	}

	MPI_Barrier(PETSC_COMM_WORLD);
 	PetscPrintf(PETSC_COMM_WORLD, "Read saved blade structure files2! \n");
	PetscErrorCode readfiles_bladestructure(IBMNodes *ibm);
	if (restart_turbinestructure) readfiles_bladestructure(ibm);

//	MPI_Barrier(PETSC_COMM_WORLD);
 //	PetscPrintf(PETSC_COMM_WORLD, "Read saved blade structure files2 \n");

 //   	PetscPrintf(PETSC_COMM_SELF, "Compute the derivatives of blade structure property! \n");
	computespatialderivatives0_bladestructure(ibm);

//    	PetscPrintf(PETSC_COMM_SELF, "Finish Initialization of turbine structure model! \n");
	return(0);
}


/** 
 * Compute spatial derivatives of displacement u 

 *	bctype = 0: 	u(0)=dudr(0)=0; 	d2udr2(L)=d3udr3(L)=0
 *			v(0)=dvdr(0)=0;		d2vdr2(L)=d3vdr3(L)=0
 *			theta(0)=0; 		dthetadr(L)=0
 */
PetscErrorCode compute_spatialderivatives_u_bladestructure(Beam *blade, vector<double> &dudr, vector<double> &d2udr2, vector<double> &d3udr3, vector<double> &d4udr4)
{
	int i, j;
	double fm4, fm3, fm2, fm1, f, fp1, fp2, fp3, fp4;

	double h=blade->dh;

	int NElmts = blade->NumberofElements;

	/* first derivatives */
	for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->u[j-1];
		fp1=blade->u[j+1];
		dudr.at(j)=dfdr(fm1, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	f=blade->u[j];
        fm1=blade->u[j-1];
	dudr.at(j) = (f-fm1)/h;

	/* second derivatives */
        for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->u[j-1];
		f=blade->u[j];
		fp1=blade->u[j+1];
		d2udr2.at(j)=d2fdr2(fm1, f, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	d2udr2.at(j)=0.0;

        /* third derivatives */
	for (j=2;j<NElmts-2;j++) 
	{
		fm2=blade->u[j-2];
		fm1=blade->u[j-1];
		f=blade->u[j];
		fp1=blade->u[j+1];
		fp2=blade->u[j+2];
		d3udr3.at(j)=d3fdr3(fm2, fm1, f, fp1, fp2, h);
	}
        
	// BC at 1
	j=1;
        fm1=blade->u[j-1];
        f=blade->u[j];
        fp1=blade->u[j+1];
        fp2=blade->u[j+2];	       
	d3udr3.at(j)=(fm1-0.5*f-fp1+0.5*fp2)/pow(h,3);

	// BC at NElmts-1
	j=NElmts-1;
	d3udr3.at(j)=0.0;

	// BC at NElmts-2
	j=NElmts-2;
	fm2=blade->u[j-2];
	fm1=blade->u[j-1];
	f=blade->u[j];
	d3udr3.at(j)=(-0.5*fm2+fm1-0.3*f)/pow(j,3);
 
	/* fourth derivatives */
	for (j=2;j<NElmts-2;j++) 
	{
                fm2=blade->u[j-2];
                fm1=blade->u[j-1];
                f=blade->u[j];
                fp1=blade->u[j+1];
                fp2=blade->u[j+2];
		d4udr4.at(j)=d4fdr4(fm2, fm1, f, fp1, fp2, h);
	}

	// BC at 1
	j=1;
	fm1=blade->u[j-1];
        f=blade->u[j];
        fp1=blade->u[j+1];
        fp2=blade->u[j+2];        
	d4udr4.at(j)=(-4*fm1+7*f-4*fp1+fp2)/pow(h,4);

	// BC at NElmts-1
	j=NElmts-1;
        fm2=blade->u[j-2];
	fm1=blade->u[j-1];
        f=blade->u[j];
	d4udr4.at(j)=(2*fm2-4*fm1+2*f)/pow(h,4);

	// BC at NElmts-2
	j=NElmts-2;
        fm2=blade->u[j-2];
	fm1=blade->u[j-1];
        f=blade->u[j];
        fp1=blade->u[j+1];
	d4udr4.at(j)=(fm2-4*fm1+5*f-2*fp1)/pow(h,4);

	return(0);
}

/** 
 * Compute spatial derivatives of displacement v 

 *	bctype = 0: 	u(0)=dudr(0)=0; 	d2udr2(L)=d3udr3(L)=0
 *			v(0)=dvdr(0)=0;		d2vdr2(L)=d3vdr3(L)=0
 *			theta(0)=0; 		dthetadr(L)=0
 */
PetscErrorCode compute_spatialderivatives_v_bladestructure(Beam *blade, vector<double> &dvdr, vector<double> &d2vdr2, vector<double> &d3vdr3, vector<double> &d4vdr4)
{
	int i, j;
	double fm4, fm3, fm2, fm1, f, fp1, fp2, fp3, fp4;

	double h=blade->dh;

	int NElmts = blade->NumberofElements;

	/* first derivatives */
	for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->v[j-1];
		fp1=blade->v[j+1];
		dvdr.at(j)=dfdr(fm1, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	f=blade->v[j];
        fm1=blade->v[j-1];
	dvdr.at(j) = (f-fm1)/h;

	/* second derivatives */
        for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->v[j-1];
		f=blade->v[j];
		fp1=blade->v[j+1];
		d2vdr2.at(j)=d2fdr2(fm1, f, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	d2vdr2.at(j)=0.0;

        /* third derivatives */
	for (j=2;j<NElmts-2;j++) 
	{
		fm2=blade->v[j-2];
		fm1=blade->v[j-1];
		f=blade->v[j];
		fp1=blade->v[j+1];
		fp2=blade->v[j+2];
		d3vdr3.at(j)=d3fdr3(fm2, fm1, f, fp1, fp2, h);
	}
        
	// BC at 1
	j=1;
        fm1=blade->v[j-1];
        f=blade->v[j];
        fp1=blade->v[j+1];
        fp2=blade->v[j+2];	       
	d3vdr3.at(j)=(fm1-0.5*f-fp1+0.5*fp2)/pow(h,3);

	// BC at NElmts-1
	j=NElmts-1;
	d3vdr3.at(j)=0.0;

	// BC at NElmts-2
	j=NElmts-2;
	fm2=blade->v[j-2];
	fm1=blade->v[j-1];
	f=blade->v[j];
	d3vdr3.at(j)=(-0.5*fm2+fm1-0.3*f)/pow(j,3);
 
	/* fourth derivatives */
	for (j=2;j<NElmts-2;j++) 
	{
                fm2=blade->v[j-2];
                fm1=blade->v[j-1];
                f=blade->v[j];
                fp1=blade->v[j+1];
                fp2=blade->v[j+2];
		d4vdr4.at(j)=d4fdr4(fm2, fm1, f, fp1, fp2, h);
	}

	// BC at 1
	j=1;
	fm1=blade->v[j-1];
        f=blade->v[j];
        fp1=blade->v[j+1];
        fp2=blade->v[j+2];        
	d4vdr4.at(j)=(-4*fm1+7*f-4*fp1+fp2)/pow(h,4);

	// BC at NElmts-1
	j=NElmts-1;
        fm2=blade->v[j-2];
	fm1=blade->v[j-1];
        f=blade->v[j];
	d4vdr4.at(j)=(2*fm2-4*fm1+2*f)/pow(h,4);

	// BC at NElmts-2
	j=NElmts-2;
        fm2=blade->v[j-2];
	fm1=blade->v[j-1];
        f=blade->v[j];
        fp1=blade->v[j+1];
	d4vdr4.at(j)=(fm2-4*fm1+5*f-2*fp1)/pow(h,4);

	return(0);
}

PetscErrorCode compute_spatialderivatives_dudtdvdt_bladestructure(Beam *blade, vector<double> &ddt_dudr, vector<double> &ddt_dvdr)
{
	int i, j;
	double fm1, f, fp1;

	double h=blade->dh;

	int NElmts = blade->NumberofElements;
	/* first derivatives */
	for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->dudt[j-1];
		fp1=blade->dudt[j+1];
		ddt_dudr.at(j)=dfdr(fm1, fp1, h);

		fm1=blade->dvdt[j-1];
		fp1=blade->dvdt[j+1];
		ddt_dvdr.at(j)=dfdr(fm1, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	fm1=blade->dudt[j-1];
        f=blade->dudt[j];
        ddt_dudr.at(j)=(f-fm1)/h;

        fm1=blade->dvdt[j-1];
        f=blade->dvdt[j];
        ddt_dvdr.at(j)=(f-fm1)/h;

	return(0);
}

/** 
 * Compute spatial derivatives of displacements 

 *	bctype = 0: 	u(0)=dudr(0)=0; 	d2udr2(L)=d3udr3(L)=0
 *			v(0)=dvdr(0)=0;		d2vdr2(L)=d3vdr3(L)=0
 *			theta(0)=0; 		dthetadr(L)=0
 */

PetscErrorCode compute_spatialderivatives_theta_bladestructure(Beam *blade, vector<double> &dthetadr, vector<double> &d2thetadr2)
{
	int i, j;
	double fm1, f, fp1;

	double h=blade->dh;
	int NElmts = blade->NumberofElements;

	/* first derivatives */
	for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->theta[j-1];
		fp1=blade->theta[j+1];
		dthetadr.at(j)=dfdr(fm1, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	dthetadr.at(j)=0.0;

	/* second derivatives */
        for (j=1;j<NElmts-1;j++) 
	{
		fm1=blade->theta[j-1];
		f=blade->theta[j];
		fp1=blade->theta[j+1];
		d2thetadr2.at(j)=d2fdr2(fm1, f, fp1, h);
	}

	// BC at NElmts-1
	j=NElmts-1;
	fm1=blade->theta[j-1];
	f=blade->theta[j];

	d2thetadr2.at(j)=(2*fm1-2*f)/pow(h,2); 

	return(0);
}

/** 
 * Compute temporal derivatives of phi 
 */
PetscErrorCode compute_temporalderivatives_phi_bladestructure(Beam *blade, double &dphidt, double &d2phidt2)
{
	int i, j;
	double fm3, fm2, fm1, f, fp1, fp2, fp3;

	double dt=dt_turbinestructure;

	/* first derivatives */
	fm2=blade->phi_om1;
	fm1=blade->phi_o;
	f=blade->phi;
        PetscPrintf(PETSC_COMM_WORLD, "####### compute_temporalderivatives_phi %le %le \n",f,fm1);
	dphidt=dfdr_backward(fm2, fm1, f, dt);

	/* second derivatives */
	fm3=blade->phi_om2;
	fm2=blade->phi_om1;
	fm1=blade->phi_o;
	f=blade->phi;
	d2phidt2=d2fdr2_backward(fm3, fm2, fm1, f, dt);

        PetscPrintf(PETSC_COMM_WORLD, "####### Angular Velocity and Acceleration %le %le \n",dphidt,d2phidt2);
	return(0);
}

/** 
 * Compute temporal derivatives of beta 
 */
PetscErrorCode compute_temporalderivatives_beta_bladestructure(Beam *blade, double &dbetadt, double &d2betadt2)
{
	int i, j;
	double fm3, fm2, fm1, f, fp1, fp2, fp3;

	double dt=dt_turbinestructure;

	/* first derivatives */
	fm2=blade->beta_om1;
	fm1=blade->beta_o;
	f=blade->beta;
	dbetadt=dfdr_backward(fm2, fm1, f, dt);

	/* second derivatives */
	fm3=blade->beta_om2;
	fm2=blade->beta_om1;
	fm1=blade->beta_o;
	f=blade->beta;
	d2betadt2=d2fdr2_backward(fm3, fm2, fm1, f, dt);

	return(0);
}

/** 
 * Compute temporal derivatives of theta
 */
PetscErrorCode compute_temporalderivatives_theta_bladestructure(Beam *blade, vector<double> &d2thetadt2)
{
	int i, j;
	double fm3, fm2, fm1, f, fp1, fp2, fp3;

	double dt=dt_turbinestructure;
	for (i=0; i<blade->NumberofElements; i++) 
	{
		fm2=blade->dthetadt_om1[i];
		fm1=blade->dthetadt_o[i];
		f=blade->dthetadt[i];

		if (ti==0 && !restart_turbinestructure) 
		{
			d2thetadt2.at(i)=0.0;
		} 
		else if (ti==1 && !restart_turbinestructure) 
		{
			d2thetadt2.at(i)=(f-fm1)/dt;
		} 
		else 
		{
			d2thetadt2.at(i)=dfdr_backward(fm2, fm1, f, dt);
		}
	}
	
	return(0);
}

/** 
 * Compute temporal derivatives of u
 */
PetscErrorCode compute_temporalderivatives_u_bladestructure(Beam *blade, vector<double> &d2udt2)
{
	int i, j;
	double fm3, fm2, fm1, f, fp1, fp2, fp3;

	double dt=dt_turbinestructure;
	for (i=0; i<blade->NumberofElements; i++) 
	{
		fm2=blade->dudt_om1[i];
		fm1=blade->dudt_o[i];
		f=blade->dudt[i];

		if (ti==0 && !restart_turbinestructure) 
		{
			d2udt2.at(i)=0.0;
		} 
		else if (ti==1 && !restart_turbinestructure) 
		{
			d2udt2.at(i)=(f-fm1)/dt;
		} 
		else 
		{
			d2udt2.at(i)=dfdr_backward(fm2, fm1, f, dt);
		}
	}
	
	return(0);
}


/** 
 * Compute temporal derivatives of v
 */
PetscErrorCode compute_temporalderivatives_v_bladestructure(Beam *blade, vector<double> &d2vdt2)
{
	int i, j;
	double fm3, fm2, fm1, f, fp1, fp2, fp3;

	double dt=dt_turbinestructure;
	for (i=0; i<blade->NumberofElements; i++) 
	{
		fm2=blade->dvdt_om1[i];
		fm1=blade->dvdt_o[i];
		f=blade->dvdt[i];

		if (ti==0 && !restart_turbinestructure) 
		{
			d2vdt2.at(i)=0.0;
		} 
		else if (ti==1 && !restart_turbinestructure) 
		{
			d2vdt2.at(i)=(f-fm1)/dt;
		} 
		else 
		{
			d2vdt2.at(i)=dfdr_backward(fm2, fm1, f, dt);
		}
	}
	
	return(0);
}

/**
 * Compute w0 and w1 on page 211 of Kallesoe's paper 
 */
PetscErrorCode computew0w1_bladestructure(Beam *blade, vector<double> &w0, vector<double> &w1, vector<double> dudr, vector<double> dvdr)
{
 	int i, j, k;

	w0.at(blade->NumberofElements-1)=0;
	w1.at(blade->NumberofElements-1)=0;

	for (i=1; i<blade->NumberofElements; i++) {
		double _w0=0;
		double _w1=0;

//		for (j=i; j<blade->NumberofElements-1; j++) {
		for (j=0; j<i; j++) {
			double dlpidr = blade->dlpidr[j];
			double dr = blade->r[j+1]-blade->r[j];

			_w0 += sqrt(1-dlpidr*dlpidr)*dr;
			_w1 += -0.5*sqrt(pow(dudr.at(j),2)+pow(dvdr.at(j),2)+2.0*dlpidr*dudr.at(j))*dr;
		}
		w0.at(i)=_w0+blade->hubrad;
		w1.at(i)=_w1;
	}

	return(0);
}

/** 
 * Compute first derivative of w0
 */
PetscErrorCode compute_dw0dr_bladestructure(Beam *blade, vector<double> &dw0dr, vector<double> w0)
{
	int i, j;
	double fm3, fm2, fm1, f, fp1, fp2, fp3, dfdrm1, d2fdr2p1, d2fdr2p2, d3fdr3p1, d3fdr3p2, d3fdr3p3, d4fdr4p1;

	double h=blade->dh;
	/* first derivatives */
	for (j=1;j<blade->NumberofElements-1;j++) {

		fm1=w0.at(j-1);
		fp1=w0.at(j+1);
		dw0dr.at(j)=dfdr(fm1, fp1, h);
	}

	return(0);
}

/**
 * Compute Fu0 in Eq. (12b)
 */ 

PetscErrorCode computeFu0_bladestructure(Beam *blade, vector<double> &Fu0, vector<double> d2thetadt2)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) {
		double Theta = blade->twistpre[i]+blade->theta[i];
//		double Theta = blade->theta[i];
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];

//		Fu0.at(i)=-mass*d2thetadt2.at(i)*lcg*sin(Theta);
		Fu0.at(i)=-d2thetadt2.at(i)*lcg*sin(Theta);
//    PetscPrintf(PETSC_COMM_WORLD,"####### Compute Fu0: %le %le %le \n", mass, lcg, d2thetadt2.at(i));
	}

	return(0);
}

/**
 * Compute Fv0 in Eq. (12b)
 */
PetscErrorCode computeFv0_bladestructure(Beam *blade, vector<double> &Fv0, vector<double> d2thetadt2)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double Theta = blade->twistpre[i]+blade->theta[i];
//		double Theta = blade->theta[i];
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];

//		Fv0.at(i)=mass*d2thetadt2.at(i)*lcg*cos(Theta);
		Fv0.at(i)=d2thetadt2.at(i)*lcg*cos(Theta);
	}

	return(0);
}

/**
 * Compute ucg and vcg in Eq. (13)
 */
PetscErrorCode computeucgvcg4Fu1Fv1_bladestructure(Beam *blade)
{

 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double thetabar = blade->twistpre[i];
		double theta = blade->theta[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double lcg = blade->lcg[i];

		blade->ucg_om1[i]=blade->ucg_o[i];
		blade->ucg_o[i]=blade->ucg[i];
		blade->ucg[i]=u+lpi+lcg*cos(thetabar)-lcg*theta*sin(thetabar);

		blade->vcg_om1[i]=blade->vcg_o[i];
		blade->vcg_o[i]=blade->vcg[i];
		blade->vcg[i]=v+lcg*sin(thetabar)+lcg*theta*cos(thetabar);
	}

	return(0);
}

/**
 * Compute time derivatives of ucg and vcg in eq. (13)
 */
PetscErrorCode computeucgvcgtimederivative4Fu1Fv1_bladestructure(Beam *blade, vector<double> &ducgdt, vector<double> &dvcgdt)
{

 	int i, j, k;
	double f, fm1, fm2;

	double dt=dt_turbinestructure;
	for (i=0; i<blade->NumberofElements; i++) 
	{
		fm2=blade->ucg_om1[i];
		fm1=blade->ucg_o[i];
		f=blade->ucg[i];
		ducgdt.at(i)=dfdr_backward(fm2, fm1, f, dt);

		fm2=blade->vcg_om1[i];
		fm1=blade->vcg_o[i];
		f=blade->vcg[i];
		dvcgdt.at(i)=dfdr_backward(fm2, fm1, f, dt);
	}

	return(0);
}

/**
 * Compute integral of T1 in eq. (13)
 */
PetscErrorCode computeT1integral_bladestructure(Beam *blade, vector<double> &T1, vector<double> &T1integral, double& dbetadt, double& dphidt)
{
 	int i, j, k;

	vector<double> r1;

	r1.resize(blade->NumberofElements);

	for (i=0;i<blade->NumberofElements;i++) 
	{
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double lcg = blade->lcg[i];

		T1.at(i)=2*mass*dbetadt*dphidt*((u+lpi)*sin(beta)+v*cos(beta)+lcg*sin(thetabar+beta));
	}

	T1integral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{

		double _T1integral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) 
		{
			double dr = r1.at(j+1)-r1.at(j);
			_T1integral+=T1.at(j)*dr;
		}
		T1integral.at(i)=_T1integral;
	}

	return(0);
}

/**
 * Compute integral of dT1dr in eq. (13)
 */
PetscErrorCode computedT1drintegral_bladestructure(Beam *blade, vector<double> &dT1dr, vector<double> &dT1drintegral, vector<double> T1, double& dbetadt, double& dphidt, vector<double> dudr, vector<double> dvdr)
{

 	int i, j, k;
	vector<double> r1;

	r1.resize(blade->NumberofElements);

	for (i=0;i<blade->NumberofElements;i++) {
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		dT1dr.at(i)=2*dmassdr*dbetadt*dphidt*((u+lpi)*sin(beta)+v*cos(beta)+lcg*sin(thetabar+beta))+2*mass*dbetadt*dphidt*((dudr.at(i)+dlpidr)*sin(beta)+dvdr.at(i)*cos(beta)+dlcgdr*sin(thetabar+beta));

	}


	dT1drintegral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) {

		double _dT1drintegral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) {
			double dr = r1[j+1]-r1[j];
			_dT1drintegral+=dT1dr.at(j)*dr;
		}

//		int itip=blade->NumberofElements-1;
//
//		_dT1drintegral += T1.at(itip)*blade->r[itip] - T1.at(i)*blade->r[i];

		dT1drintegral.at(i)=_dT1drintegral;
	}

	return(0);
}

/**
 * Compute Fu1 eq. (13a)
 */
PetscErrorCode computeFu1_bladestructure(Beam *blade, vector<double> &Fu1, double dbetadt, double d2betadt2, double dphidt, vector<double> dudr, vector<double> d2udr2, vector<double> ducgdt, vector<double> dvcgdt, vector<double> T1, vector<double> T1integral, vector<double> dT1dr, vector<double> dT1drintegral, vector<double> dthetadr )
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{

		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];
		double dthetabardr = dthetadr.at(i);

		//blade->coefu_dudr[i] = dT1drintegral.at(i);
		//blade->coefu_d2udr2[i] = T1integral.at(i);

		double A_Fu1 = -d2betadt2*mass*blade->vcg[i];
		double B_Fu1 = -pow(dbetadt,2)*mass*blade->ucg[i];
		double C_Fu1 = -2*dbetadt*mass*dvcgdt.at(i);
		double D_Fu1 = dT1dr.at(i)*lcg*cos(thetabar)+T1.at(i)*dlcgdr*cos(thetabar)-T1.at(i)*lcg*sin(thetabar)*dthetabardr;
		double E_Fu1 = (d2udr2.at(i)+d2lpidr2)*T1integral.at(i)+(dudr.at(i)+dlpidr)*dT1drintegral.at(i);

		Fu1.at(i)=A_Fu1+B_Fu1+C_Fu1+D_Fu1+E_Fu1/* + blade->coefu_dudr[i]*dudr.at(i) + blade->coefu_d2udr2[i]*d2udr2.at(i)*/;

	}

//   WriteForces_Term("Fu1.dat",blade->NumberofElements, blade, Fu1);

	return(0);
}

/**
 * Compute Fv1 eq. (13b)
 */
PetscErrorCode computeFv1_bladestructure(Beam *blade, vector<double> &Fv1, double dbetadt, double d2betadt2, double dphidt, vector<double> dvdr, vector<double> d2vdr2, vector<double> ducgdt, vector<double> dvcgdt, vector<double> T1, vector<double> T1integral, vector<double> dT1dr, vector<double> dT1drintegral, vector<double> dthetadr)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) {

		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];
		double dthetabardr = dthetadr.at(i);
		double ucg = blade->ucg[i];
		double vcg = blade->vcg[i];

//		blade->coefv_dvdr[i] = dT1drintegral.at(i);
//		blade->coefv_d2vdr2[i] = T1integral.at(i);

		double A_Fv1 = d2betadt2*mass*ucg;
		double B_Fv1 = -pow(dbetadt,2)*mass*vcg;
		double C_Fv1 = 2*dbetadt*mass*ducgdt.at(i);
		double D_Fv1 = dT1dr.at(i)*lcg*sin(thetabar)+T1.at(i)*dlcgdr*sin(thetabar)+T1.at(i)*lcg*cos(thetabar)*dthetabardr;
		double E_Fv1 = d2vdr2.at(i)*T1integral.at(i)+dvdr.at(i)*dT1drintegral.at(i);
		//double E_Fv1 = 0.0;

		Fv1.at(i)=A_Fv1+B_Fv1+C_Fv1+D_Fv1+E_Fv1/*+blade->coefv_dvdr[i]*dvdr.at(i) + blade->coefv_d2vdr2[i]*d2vdr2.at(i)*/;
 
    
	}

//    WriteForces_Term("Fv1.dat",blade->NumberofElements, blade, Fv1);
	return(0);
}

/**
 * Compute integral of T2 in eq. (14)
 */
PetscErrorCode computeT2integral_bladestructure(Beam *blade, vector<double> &T2, vector<double> &T2integral, double dbetadt, double dphidt)
{

 	int i, j, k;
	vector<double> r1;

	r1.resize(blade->NumberofElements);

	for (i=0;i<blade->NumberofElements;i++) {
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double u = blade->u[i];
		double dudt = blade->dudt[i];
		double v = blade->v[i];
		double dvdt = blade->dvdt[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double lcg = blade->lcg[i];

		T2.at(i)=2*mass*dphidt*(dudt*cos(beta)-dvdt*sin(beta));
                
//  	        PetscPrintf(PETSC_COMM_WORLD,"####### Compute T2: %le %le %le \n", dudt, dvdt, beta);
	}

	T2integral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{
		double _T2integral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) {
			double dr = r1[j+1]-r1[j];
			_T2integral+=T2.at(j)*dr;
		}
		T2integral.at(i)=_T2integral;
	}

	return(0);
}

/**
 * Compute integral of dT2dr in eq. (14)
 */
PetscErrorCode computedT2drintegral_bladestructure(Beam *blade, vector<double> &dT2dr, vector<double> &dT2drintegral, double dbetadt, double dphidt, vector<double> ddt_dudr, vector<double> ddt_dvdr, vector<double> T2)
{

 	int i, j, k;
	vector<double> r1;

	r1.resize(blade->NumberofElements);

	for (i=0;i<blade->NumberofElements;i++) 
	{
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double u = blade->u[i];
		double dudt = blade->dudt[i];
		double v = blade->v[i];
		double dvdt = blade->dvdt[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		dT2dr.at(i)=2*dmassdr*dphidt*(dudt*cos(beta)-dvdt*sin(beta))+2*mass*dphidt*(ddt_dudr.at(i)*cos(beta)-ddt_dvdr.at(i)*sin(beta));
	}

	dT2drintegral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{
		double _dT2drintegral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) {
			double dr = r1.at(j+1)-r1.at(j);
			_dT2drintegral+=dT2dr.at(j)*dr;
		}

//		int itip=blade->NumberofElements-1;

//		_dT2drintegral += T2.at(itip)*blade->r[itip] - T2.at(i)*blade->r[i];

		dT2drintegral.at(i)=_dT2drintegral;
	}

	return(0);
}

/**
 * Compute integral of T3 in eq. (14)
 */
PetscErrorCode computeT3integral_bladestructure(Beam *blade, vector<double> &T3, vector<double> &T3integral, double dbetadt, double dphidt, vector<double> w0)
{

 	int i, j, k;
	vector<double> r1;

	r1.resize(blade->NumberofElements);

	for (i=0;i<blade->NumberofElements;i++) 
	{
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double u = blade->u[i];
		double dudt = blade->dudt[i];
		double v = blade->v[i];
		double dvdt = blade->dvdt[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double lcg = blade->lcg[i];

		T3.at(i)=pow(dphidt,2)*mass*w0.at(i);

   // PetscPrintf(PETSC_COMM_WORLD,"Compute  T3int %d %e %e %e\n",i,dphidt,w0.at(i),mass);
	}


	T3integral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{
		double _T3integral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) {
			double dr = r1.at(j+1)-r1.at(j);
			_T3integral+=T3.at(j)*dr;
		}
		T3integral.at(i)=_T3integral;
	}

	return(0);
}

/**
 * Compute integral of dT3dr in eq. (14)
 */
PetscErrorCode computedT3drintegral_bladestructure(Beam *blade, vector<double> &dT3dr, vector<double> &dT3drintegral, double dbetadt, double dphidt, vector<double> ddt_dudr, vector<double> ddt_dvdr, vector<double> w0, vector<double> dw0dr, vector<double> T3)
{

 	int i, j, k;
	vector<double> r1;

	r1.resize(blade->NumberofElements);

	for (i=0;i<blade->NumberofElements;i++) 
	{
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double u = blade->u[i];
		double dudt = blade->dudt[i];
		double v = blade->v[i];
		double dvdt = blade->dvdt[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		dT3dr.at(i)=pow(dphidt,2)*dmassdr*w0.at(i)+pow(dphidt,2)*mass*dw0dr.at(i);
	}


	dT3drintegral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{
		double _dT3drintegral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) {
			double dr = r1.at(j+1)-r1.at(j);
			_dT3drintegral+=dT3dr.at(j)*dr;
		}

//		int itip=blade->NumberofElements-1;
//
//		_dT3drintegral += T3.at(itip)*blade->r[itip] - T3.at(i)*blade->r[i];

		dT3drintegral.at(i)=_dT3drintegral;
	}


	return(0);
}

/**
 * Compute ucghat in eq. (14)
 */
PetscErrorCode computeucghat4Fu2Fv2_bladestructure(Beam *blade, vector<double> &ucghat)
{

 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{

		double thetabar = blade->twistpre[i]+blade->theta[i];
		double theta = blade->twistpre[i]+blade->theta[i];
		double mass = blade->mass[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double lcg = blade->lcg[i];

		ucghat.at(i)=(u+lpi)*cos(beta)-v*sin(beta)+lcg*cos(thetabar+beta)-lcg*theta*sin(thetabar+beta);

	}

	return(0);
}

/**
 * Compute Fu2 of eq. (14a)
 */
PetscErrorCode computeFu2_bladestructure(Beam *blade, vector<double> &Fu2, double dphidt, vector<double> w0, vector<double> dw0dr, vector<double> dudr, vector<double> ddt_dudr, vector<double> d2udr2, vector<double> ddt_dvdr, vector<double> ucghat, vector<double> T2, vector<double> dT2dr, vector<double> T2integral, vector<double> dT2drintegral, vector<double> T3, vector<double> dT3dr, vector<double> T3integral, vector<double> dT3drintegral, vector<double> dthetadr)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{

		double thetabar = blade->twistpre[i];
		double theta = blade->theta[i];
		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		double A_Fu2 = -pow(dphidt,2)*mass*ucghat.at(i)*cos(beta);
		double B_Fu2 = -pow(dphidt,2)*(dmassdr*lcg*w0.at(i)*(cos(thetabar)-theta*sin(thetabar))+mass*dlcgdr*w0.at(i)*(cos(thetabar)-theta*sin(thetabar))+mass*lcg*dw0dr.at(i)*(cos(thetabar)-theta*sin(thetabar))+mass*lcg*w0.at(i)*(-sin(thetabar)*dthetadr.at(i)-dthetadr.at(i)*sin(thetabar)-theta*cos(thetabar)*dthetadr.at(i)));
		double C_Fu2 = -(dlcgdr*T2.at(i)+lcg*dT2dr.at(i))*cos(thetabar);
		double D_Fu2 = -2*dphidt*mass*lcg*(ddt_dudr.at(i)*cos(thetabar)+ddt_dvdr.at(i)*sin(thetabar))*cos(beta);

		//double E_Fu2 = (d2udr2.at(i)+d2lpidr2)*(T3integral.at(i)+T2integral.at(i))+(dudr.at(i)+dlpidr)*(dT3drintegral.at(i)+dT2drintegral.at(i));

		double E_Fu2 = -(d2lpidr2)*(T3integral.at(i)+T2integral.at(i))-(dlpidr)*(dT3drintegral.at(i)+dT2drintegral.at(i));

		blade->coefu_d2udr2[i] = 1.*(T3integral.at(i)+T2integral.at(i)); 
		blade->coefu_dudr[i] =  1.*(dT3drintegral.at(i)+dT2drintegral.at(i));

//		double E_Fu2 = (d2lpidr2)*(T3integral.at(i)+T2integral.at(i))+(dlpidr)*(dT3drintegral.at(i)+dT2drintegral.at(i));
//
//
//		blade->coefu_d2udr2[i] = (T3integral.at(i)+T2integral.at(i)); 
//		blade->coefu_dudr[i] =  (dT3drintegral.at(i)+dT2drintegral.at(i));

//		blade->coefu_d2udr2[i] = 0.; 
//		blade->coefu_dudr[i]   = 0.;

		Fu2.at(i)=-(A_Fu2+B_Fu2+C_Fu2+D_Fu2+E_Fu2);


//    PetscPrintf(PETSC_COMM_WORLD,"Compute Fu2 %d %e %e %e\n",i,E_Fu2,T3integral.at(i),T3integral.at(i));
	}

//    PetscPrintf(PETSC_COMM_WORLD,"####### turbin structure dphidt %le %le \n", dphidt,blade->beta);
//      WriteForces_Term("Fu2.dat",blade->NumberofElements, blade, Fu2);
//      WriteForces_Term("T2.dat",blade->NumberofElements, blade, T2integral);
//      WriteForces_Term("T3.dat",blade->NumberofElements, blade, T3integral);
//      WriteForces_Term("dT2dr.dat",blade->NumberofElements, blade, dT2dr);
//      WriteForces_Term("dT3dr.dat",blade->NumberofElements, blade, dT3dr);
//      WriteForces_Term("dT2drintegral.dat",blade->NumberofElements, blade, dT2drintegral);
//      WriteForces_Term("dT3drintegral.dat",blade->NumberofElements, blade, dT3drintegral);
//      WriteForces_Term("Wo.dat",blade->NumberofElements, blade, w0);
//    WriteForces_Term("dlpidr.dat",blade->NumberofElements, blade, dT2drintegral);
//    WriteForces_Term("dlpidr.dat",blade->NumberofElements, blade, ucghat);
	return(0);
}

/**
 * Compute Fv2 of eq. (14b)
 */
PetscErrorCode computeFv2_bladestructure(Beam *blade, vector<double> &Fv2, double dphidt, vector<double> w0, vector<double> dw0dr, vector<double> dvdr, vector<double> ddt_dudr, vector<double> d2vdr2, vector<double> ddt_dvdr, vector<double> ucghat, vector<double> T2, vector<double> dT2dr, vector<double> T2integral, vector<double> dT2drintegral, vector<double> T3, vector<double> dT3dr, vector<double> T3integral, vector<double> dT3drintegral, vector<double> dthetadr)

{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double theta = blade->theta[i];
		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		double A_Fv2 = pow(dphidt,2)*mass*ucghat.at(i)*sin(beta);
		double B_Fv2 = -pow(dphidt,2)*(dmassdr*lcg*w0.at(i)*(sin(thetabar)+theta*cos(thetabar))+mass*dlcgdr*w0.at(i)*(sin(thetabar)+theta*cos(thetabar))+mass*lcg*dw0dr.at(i)*(sin(thetabar)+theta*cos(thetabar))+mass*lcg*w0.at(i)*(cos(thetabar)*dthetadr.at(i)+dthetadr.at(i)*cos(thetabar)+theta*sin(thetabar)*dthetadr.at(i)));
		double C_Fv2 = dlcgdr*T2.at(i)*sin(thetabar)+lcg*dT2dr.at(i)*sin(thetabar);
		double D_Fv2 = 2*dphidt*mass*lcg*(ddt_dudr.at(i)*cos(thetabar)+ddt_dvdr.at(i)*sin(thetabar))*sin(beta);
		//double E_Fv2 = -d2vdr2.at(i)*(T3integral.at(i)+T2integral.at(i))-dvdr.at(i)*(dT3drintegral.at(i)+dT2drintegral.at(i));

		blade->coefv_d2vdr2[i] = -1.*(T3integral.at(i)+T2integral.at(i));
		blade->coefv_dvdr[i] = -1.*(T3integral.at(i)+T2integral.at(i));

//		Fv2.at(i)=A_Fv2+B_Fv2-C_Fv2+D_Fv2+E_Fv2;
		Fv2.at(i)=(A_Fv2+B_Fv2-C_Fv2+D_Fv2);
	}

//    WriteForces_Term("Fv2.dat",blade->NumberofElements, blade, Fv2);
	return(0);
}

/**
 * Compute integral of T4 of eq. (15)
 */
PetscErrorCode computeT4integral_bladestructure(Beam *blade, vector<double> &T4, vector<double> &T4integral)
{

 	int i, j, k;
	vector<double> r1;

	r1.resize(blade->NumberofElements);

	double gravity = -sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0;i<blade->NumberofElements;i++) 
	{
		r1.at(i)=blade->r[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		double mass = blade->mass[i];
		double phi = blade->phi;
		double dphidt = blade->dphidt;
		double u = blade->u[i];
		double dudt = blade->dudt[i];
		double v = blade->v[i];
		double dvdt = blade->dvdt[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double lcg = blade->lcg[i];

		T4.at(i)=mass*gravity*cos(phi);
	}

	T4integral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{

		double _T4integral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) {
			double dr = r1.at(j+1)-r1.at(j);
			_T4integral+=T4.at(j)*dr;
		}
		T4integral.at(i)=_T4integral;
	}


	return(0);
}

/**
 * Compute integral of dT4dr of eq. (15)
 */
PetscErrorCode computedT4drintegral_bladestructure(Beam *blade, vector<double> &dT4dr, vector<double> &dT4drintegral, vector<double> T4)
{

 	int i, j, k;
	vector<double> r1;
	r1.resize(blade->NumberofElements);

	double gravity = -sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0;i<blade->NumberofElements;i++) 
	{
		r1[i]=blade->r[i];

		double dmassdr = blade->dmassdr[i];
		double phi = blade->phi;

		dT4dr.at(i)=dmassdr*gravity*cos(phi);
	}

	dT4drintegral.at(blade->NumberofElements-1)=0;
	for (i=0;i<blade->NumberofElements-1;i++) 
	{
		double _dT4drintegral = 0.0;
		for (j=i;j<blade->NumberofElements-1;j++) 
		{
			double dr = r1[j+1]-r1[j];
			_dT4drintegral+=dT4dr.at(j)*dr;
		}

//		int itip=blade->NumberofElements-1;
//
//		_dT4drintegral += T4.at(itip)*r1.at(itip) - T4.at(i)*r1.at(i);

		dT4drintegral.at(i)=_dT4drintegral;
	}


	return(0);
}

/**
 * Compute Fu3 of eq. (15a) Gravity
 */
PetscErrorCode computeFu3_bladestructure(Beam *blade, vector<double> &Fu3, vector<double> dudr, vector<double> d2udr2, vector<double> dvdr, vector<double> d2vdr2, vector<double> T4integral, vector<double> dT4drintegral, vector<double> dthetadr)
{
 	int i, j, k;

	double gravity = -sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double thetabar = blade->twistpre[i];
		double theta = blade->theta[i];
		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double phi = blade->phi;
		double dphidt = blade->dphidt;
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		double A_Fu3 = mass*gravity*sin(phi)*cos(beta);

		double B_Fu3 = ((dlcgdr*(dlpidr)*cos(thetabar)*cos(beta)+lcg*(d2lpidr2)*cos(thetabar)*cos(beta)))*mass*gravity*sin(phi);
		double C_Fu3 = -(dmassdr*lcg*(cos(thetabar)-theta*sin(thetabar))+mass*dlcgdr*(cos(thetabar)-theta*sin(thetabar))+mass*lcg*(-sin(thetabar)*dthetadr.at(i)-dthetadr.at(i)*sin(thetabar)-theta*cos(thetabar)*dthetadr.at(i)))*gravity*cos(phi);
		double D_Fu3 = (d2lpidr2)*T4integral.at(i)+(dlpidr)*dT4drintegral.at(i);



                blade->coefu_dudr[i] += ( cos(thetabar)*cos(beta)*dlcgdr*mass*gravity*sin(phi) + dT4drintegral.at(i) );
                blade->coefu_d2udr2[i] += (lcg*cos(thetabar)*cos(beta)*mass*gravity*sin(phi) + T4integral.at(i));

                blade->coefu_dvdr[i] = (dlcgdr*sin(thetabar)*cos(beta)*mass*gravity*sin(phi));
                blade->coefu_d2vdr2[i] = (lcg*sin(thetabar)*cos(beta)*mass*gravity*sin(phi));

		Fu3.at(i)=(A_Fu3+B_Fu3+C_Fu3+D_Fu3);

//    PetscPrintf(PETSC_COMM_WORLD, "i:%i A: %e B: %e  C:%e  D: %e E: %e \n", i, A_Fu3 ,B_Fu3,C_Fu3,D_Fu3,dlpidr);
//		Fu3.at(i)=0.0;
//
	}

//    WriteForces_Term("Fu3.dat",blade->NumberofElements, blade, Fu3);
	return(0);
}


/**
 * Compute Fv3 of eq. (15b) Gravity
 */
PetscErrorCode computeFv3_bladestructure(Beam *blade, vector<double> &Fv3, vector<double> dudr, vector<double> d2udr2, vector<double> dvdr, vector<double> d2vdr2, vector<double> T4integral, vector<double> dT4drintegral, vector<double> dthetadr)
{
 	int i, j, k;

	double gravity = -sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double thetabar = blade->twistpre[i];
		double theta = blade->theta[i];
		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double phi = blade->phi;
		double u = blade->u[i];
		double v = blade->v[i];
		double beta = blade->beta;
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];

		double A_Fv3 = -mass*gravity*sin(phi)*sin(beta);
		double B_Fv3 = ((dlcgdr*(dlpidr)*sin(thetabar)*cos(beta)+lcg*(d2lpidr2)*sin(thetabar)*cos(beta)))*mass*gravity*sin(phi);
		double C_Fv3 = -(dmassdr*lcg*(sin(thetabar)-theta*cos(thetabar))+mass*dlcgdr*(sin(thetabar)-theta*cos(thetabar))+mass*lcg*(cos(thetabar)*dthetadr.at(i)-dthetadr.at(i)*cos(thetabar)+theta*sin(thetabar)*dthetadr.at(i)))*gravity*cos(phi);
		double D_Fv3 = 0.;

                
                blade->coefv_dvdr[i] += (-dlcgdr*sin(thetabar)*sin(beta)*mass*gravity*sin(phi)-dT4drintegral.at(i));
                blade->coefv_d2vdr2[i] += (lcg*sin(thetabar)*sin(beta)*mass*gravity*sin(phi)-T4integral.at(i));

                blade->coefv_dudr[i] = (dlcgdr*sin(thetabar)*cos(beta)*mass*gravity*sin(phi));
                blade->coefv_d2udr2[i] = (sin(thetabar)*cos(beta)*mass*gravity*sin(phi));



		Fv3.at(i)=(A_Fv3+B_Fv3+C_Fv3+D_Fv3);

 //   PetscPrintf(PETSC_COMM_WORLD, "i:%i A: %e B: %e  C:%e  D: %e E: %e \n", i, A_Fv3 ,B_Fv3,C_Fv3,D_Fv3,dT4drintegral.at(i));
//		Fv3.at(i)=0.0;
	}

//    WriteForces_Term("Fv3.dat",blade->NumberofElements, blade, Fv3);
	return(0);
}

/**
 * Compute Fu4 of eq. (16a)
 */


/**
 * Compute Fu4 of eq. (16a)
 */
PetscErrorCode computeFu4_bladestructure(Beam *blade, vector<double> &Fu4, vector<double> d2udr2, vector<double> d3udr3, vector<double> d4udr4, vector<double> d2vdr2, vector<double> d3vdr3, vector<double> d4vdr4, vector<double> dthetadr, vector<double> d2thetadr2)
{
 	int i, j, k;

	//double gravity = sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double theta = blade->theta[i];
		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double phi = blade->phi;
		double dphidt = blade->dphidt;
		double u = blade->u[i];
		double v = blade->v[i];
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double d3lpidr3 = blade->d3lpidr3[i];
		double d4lpidr4 = blade->d4lpidr4[i];
		double lcg = blade->lcg[i];
		double dlcgdr = blade->dlcgdr[i];
		double Iksi = blade->Iksi[i];
		double dIksidr = blade->dIksidr[i];
		double d2Iksidr2 = blade->d2Iksidr2[i];
		double Ieta = blade->Ieta[i];
		double dIetadr = blade->dIetadr[i];
		double d2Ietadr2 = blade->d2Ietadr2[i];

		double E = blade->E[i];
		double dEdr = blade->dEdr[i];
		double d2Edr2 = blade->d2Edr2[i];

		double twistpre = blade->twistpre[i] ;
		double dtwistpredr = blade->dtwistpredr[i] ;
		double d2twistpredr2 = blade->d2twistpredr2[i] ;

//		double   twistpre    =0.0; 
//		double dtwistpredr   =0.0; 
//		double d2twistpredr2 =0.0;  


		// diff(E*(Iksi*cos(thetatilde)*cos(thetatilde)+Ieta*sin(thetatilde)*sin(thetatilde))*d2udr2,2)

		// A_Fu4
		blade->coefu_d3udr3[i] = 2*dEdr*(Iksi*pow(cos(twistpre),2) + Ieta*pow(sin(twistpre),2)); 
		blade->coefu_d2udr2[i] += d2Edr2*(Iksi*pow(cos(twistpre),2) + Ieta*pow(sin(twistpre),2));
      		blade->coefu_d4udr4[i] = E*(Iksi*pow(cos(twistpre),2) + Ieta*pow(sin(twistpre),2)); 
		blade->coefu_d3udr3[i] += 2*E*(pow(cos(twistpre),2)*dIksidr + pow(sin(twistpre),2)*dIetadr + 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Ieta - 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Iksi);
		blade->coefu_d2udr2[i] += 2*dEdr*(pow(cos(twistpre),2)*dIksidr + pow(sin(twistpre),2)*dIetadr + 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Ieta - 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Iksi);
		blade->coefu_d2udr2[i] += E*(d2Iksidr2*pow(cos(twistpre),2) + d2Ietadr2*pow(sin(twistpre),2) + 2*pow(dtwistpredr,2)*pow(cos(twistpre),2)*
                                        Ieta - 2*pow(dtwistpredr,2)*pow(cos(twistpre),2)*Iksi - 2*pow(dtwistpredr,2)*pow(sin(twistpre),2)*Ieta +
                                        2*pow(dtwistpredr,2)*pow(sin(twistpre),2)*Iksi + 4*dtwistpredr*cos(twistpre)*sin(twistpre)*dIetadr +
                                        2*cos(twistpre)*sin(twistpre)*d2twistpredr2*Ieta - 4*dtwistpredr*cos(twistpre)*sin(twistpre)*dIksidr -
                                        2*cos(twistpre)*sin(twistpre)*d2twistpredr2*Iksi);

		// B_Fu4
		blade->coefu_d3vdr3[i] = 2*dtwistpredr*pow(sin(twistpre),2)*E*(Ieta - Iksi);
		blade->coefu_d2vdr2[i] = -2*dtwistpredr*pow(cos(twistpre),2)*dEdr*(Ieta - Iksi);
		blade->coefu_d2vdr2[i] += -pow(cos(twistpre),2)*d2twistpredr2*E*(Ieta - Iksi);
		blade->coefu_d3vdr3[i] += -2*dtwistpredr*pow(cos(twistpre),2)*E*(Ieta - Iksi);
		blade->coefu_d2vdr2[i] += 2*dtwistpredr*pow(sin(twistpre),2)*dEdr*(Ieta - Iksi);
		blade->coefu_d2vdr2[i] += pow(sin(twistpre),2)*d2twistpredr2*E*(Ieta - Iksi); 
		blade->coefu_d3vdr3[i] += -2*cos(twistpre)*sin(twistpre)*dEdr*(Ieta - Iksi);
		blade->coefu_d2vdr2[i] += -d2Edr2*cos(twistpre)*sin(twistpre)*(Ieta - Iksi);
		blade->coefu_d4vdr4[i] = -cos(twistpre)*sin(twistpre)*E*(Ieta - Iksi); 
		blade->coefu_d2vdr2[i] += -2*dtwistpredr*pow(cos(twistpre),2)*E*(dIetadr - dIksidr);
		blade->coefu_d2vdr2[i] += 2*dtwistpredr*pow(sin(twistpre),2)*E*(dIetadr - dIksidr);
		blade->coefu_d3vdr3[i] += -2*cos(twistpre)*sin(twistpre)*E*(dIetadr - dIksidr);	
		blade->coefu_d2vdr2[i] += -2*cos(twistpre)*sin(twistpre)*dEdr*(dIetadr - dIksidr);
		blade->coefu_d2vdr2[i] += -cos(twistpre)*sin(twistpre)*E*(d2Ietadr2 - d2Iksidr2);
		blade->coefu_d2vdr2[i] += 4*pow(dtwistpredr,2)*cos(twistpre)*sin(twistpre)*E*(Ieta - Iksi);

		// C_Fu4
		blade->coefu_d2udr2[i] += 2*dthetadr.at(i)*E*(dIetadr - dIksidr)*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += -2*dthetadr.at(i)*E*(dIetadr - dIksidr)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += 2*theta*dEdr*(dIetadr - dIksidr)*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += -2*theta*dEdr*(dIetadr - dIksidr)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += E*theta*(d2Ietadr2 - d2Iksidr2)*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += -E*theta*(d2Ietadr2 - d2Iksidr2)*cos(2*twistpre);
		blade->coefu_d3udr3[i] += 2*dthetadr.at(i)*E*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefu_d3vdr3[i] += -2*dthetadr.at(i)*E*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += 2*dthetadr.at(i)*E*(Ieta - Iksi)*2*dtwistpredr*cos(2*twistpre);
		blade->coefu_d2vdr2[i] += 2*dthetadr.at(i)*E*(Ieta - Iksi)*2*dtwistpredr*sin(2*twistpre);
		blade->coefu_d3udr3[i] += 2*theta*dEdr*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefu_d3vdr3[i] += -2*theta*dEdr*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += 2*theta*dEdr*(Ieta - Iksi)*2*dtwistpredr*cos(2*twistpre);
		blade->coefu_d2vdr2[i] += 2*theta*dEdr*(Ieta - Iksi)*2*dtwistpredr*sin(2*twistpre);
		blade->coefu_d4udr4[i] += E*theta*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefu_d4vdr4[i] +=  -E*theta*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefu_d3udr3[i] += E*theta*(Ieta - Iksi)*4*dtwistpredr*cos(2*twistpre);
		blade->coefu_d2udr2[i] += E*theta*(Ieta - Iksi)*2*d2twistpredr2*cos(2*twistpre);
		blade->coefu_d3vdr3[i] += E*theta*(Ieta - Iksi)*4*dtwistpredr*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += E*theta*(Ieta - Iksi)*2*d2twistpredr2*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += E*theta*(Ieta - Iksi)*4*pow(dtwistpredr,2)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += E*theta*(Ieta - Iksi)*4*pow(dtwistpredr,2)*sin(2*twistpre);
		blade->coefu_d2udr2[i] += 2*dthetadr.at(i)*dEdr*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += -2*dthetadr.at(i)*dEdr*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += d2Edr2*theta*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefu_d2vdr2[i] += -d2Edr2*theta*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefu_d3udr3[i] += 2*E*theta*(dIetadr - dIksidr)*sin(2*twistpre);
		blade->coefu_d3vdr3[i] += -2*E*theta*(dIetadr - dIksidr)*cos(2*twistpre);
		blade->coefu_d2udr2[i] += 2*E*theta*(dIetadr - dIksidr)*2*dtwistpredr*cos(2*twistpre);
		blade->coefu_d2vdr2[i] += 2*E*theta*(dIetadr - dIksidr)*2*dtwistpredr*sin(2*twistpre);

		double A_Fu4 = 0.0;	
		double B_Fu4 = 0.0;
		double C_Fu4 = 2*dthetadr.at(i)*E*(dIetadr - dIksidr)*(cos(twistpre)*sin(twistpre)*d2lpidr2) + 
			2*theta*dEdr*(dIetadr - dIksidr)*(cos(twistpre)*sin(twistpre)*d2lpidr2) + 
			E*theta*(d2Ietadr2 - d2Iksidr2)*(cos(twistpre)*sin(twistpre)*d2lpidr2) +
			2*dthetadr.at(i)*E*(Ieta - Iksi)*( 
			dtwistpredr*pow(cos(twistpre),2)*d2lpidr2 - dtwistpredr*pow(sin(twistpre),2)*d2lpidr2 + 
			d3lpidr3*cos(twistpre)*sin(twistpre)) + 
			2*theta*dEdr*(Ieta - Iksi)*( dtwistpredr*pow(cos(twistpre),2)*d2lpidr2 - dtwistpredr*pow(sin(twistpre),2)*d2lpidr2 +
			d3lpidr3*cos(twistpre)*sin(twistpre)) + 
			E*theta*(Ieta - Iksi)*(
			2*d3lpidr3*dtwistpredr*pow(cos(twistpre),2) + pow(cos(twistpre),2)*d2twistpredr2*d2lpidr2 -
			2*d3lpidr3*dtwistpredr*pow(sin(twistpre),2) - pow(sin(twistpre),2)*d2twistpredr2*d2lpidr2 +
			cos(twistpre)*sin(twistpre)*d4lpidr4 - 4*pow(dtwistpredr,2)*cos(twistpre)*sin(twistpre)*d2lpidr2) +
			2*dthetadr.at(i)*dEdr*(Ieta - Iksi)*(cos(twistpre)*sin(twistpre)*d2lpidr2) + 
			d2Edr2*theta*(Ieta - Iksi)*(cos(twistpre)*sin(twistpre)*d2lpidr2) + 
			d2thetadr2.at(i)*E*(Ieta - Iksi)*(cos(twistpre)*sin(twistpre)*d2lpidr2) +
			2*E*theta*(dIetadr - dIksidr)*(
			dtwistpredr*pow(cos(twistpre),2)*d2lpidr2 - dtwistpredr*pow(sin(twistpre),2)*d2lpidr2 + d3lpidr3*cos(twistpre)*sin(twistpre));

		Fu4.at(i)=A_Fu4+B_Fu4+C_Fu4;
/*+d2udr2.at(i)*blade->coefu_d2udr2[i] + d3udr3.at(i)*blade->coefu_d3udr3[i] + d4udr4.at(i)*blade->coefu_d4udr4[i]+
					    d2vdr2.at(i)*blade->coefu_d2vdr2[i] + d3vdr3.at(i)*blade->coefu_d3vdr3[i] + d4vdr4.at(i)*blade->coefu_d4vdr4[i]*/;
		//Fu4.at(i)=B_Fu4+C_Fu4;
		//Fu4.at(i)=A_Fu4;
//		Fu4.at(i)=0.0;
	}


//        WriteForces_Term("Fu4.dat",blade->NumberofElements, blade, Fu4);
//        WriteForces_Term("Fu4_d2udr2.dat",blade->NumberofElements, blade, d2udr2);
//        WriteForces_Term("Fu4_d3udr3.dat",blade->NumberofElements, blade, d3udr3);
//        WriteForces_Term("Fu4_d4udr4.dat",blade->NumberofElements, blade, d4udr4);
//        WriteForces_Term("Fu4_d2vdr2.dat",blade->NumberofElements, blade, d2vdr2);
//        WriteForces_Term("Fu4_d3vdr3.dat",blade->NumberofElements, blade, d3vdr3);
//        WriteForces_Term("Fu4_dthetadr.dat",blade->NumberofElements, blade, dthetadr);
//        WriteForces_Term("Fu4_d2thetadr2.dat",blade->NumberofElements, blade, d2thetadr2);

	return(0);
}



/**
 * Compute Fv4 of eq. (16b)
 */
PetscErrorCode computeFv4_bladestructure(Beam *blade, vector<double> &Fv4, vector<double> d2udr2, vector<double> d3udr3, vector<double> d4udr4, vector<double> d2vdr2, vector<double> d3vdr3, vector<double> d4vdr4, vector<double> dthetadr, vector<double> d2thetadr2)

{
 	int i, j, k;

//	cout << "testing \n";

	//double gravity = sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double theta = blade->theta[i];
		double mass = blade->mass[i];
		double dmassdr = blade->dmassdr[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double lpi = blade->lpi[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double d3lpidr3 = blade->d3lpidr3[i];
		double d4lpidr4 = blade->d4lpidr4[i];
		double Iksi = blade->Iksi[i];
		double dIksidr = blade->dIksidr[i];
		double d2Iksidr2 = blade->d2Iksidr2[i];
		double Ieta = blade->Ieta[i];
		double dIetadr = blade->dIetadr[i];
		double d2Ietadr2 = blade->d2Ietadr2[i];

		double E = blade->E[i];
		double dEdr = blade->dEdr[i];
		double d2Edr2 = blade->d2Edr2[i];

		double twistpre = blade->twistpre[i];
		double dtwistpredr = blade->dtwistpredr[i];
		double d2twistpredr2 = blade->d2twistpredr2[i];

//		double   twistpre    =0.0; 
//		double dtwistpredr   =0.0; 
//		double d2twistpredr2 =0.0;  

		// A_Fv4
		blade->coefv_d3vdr3[i] = 2*dEdr*(Ieta*pow(cos(twistpre),2) + Iksi*pow(sin(twistpre),2));
		blade->coefv_d2vdr2[i] += d2Edr2*(Ieta*pow(cos(twistpre),2) + Iksi*pow(sin(twistpre),2));
      		blade->coefv_d4vdr4[i] = E*(Ieta*pow(cos(twistpre),2) + Iksi*pow(sin(twistpre),2));
		blade->coefv_d3vdr3[i] += 2*E*(pow(cos(twistpre),2)*dIetadr + pow(sin(twistpre),2)*dIksidr + 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Iksi - 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Ieta);
		blade->coefv_d2vdr2[i] += 2*dEdr*(pow(cos(twistpre),2)*dIetadr + pow(sin(twistpre),2)*dIksidr + 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Iksi - 2*dtwistpredr*cos(twistpre)*sin(twistpre)*Ieta);
                blade->coefv_d2vdr2[i] +=   
                        E*(d2Ietadr2*pow(cos(twistpre),2) + d2Iksidr2*pow(sin(twistpre),2) + 2*pow(dtwistpredr,2)*pow(cos(twistpre),2)*Iksi - 
                        2*pow(dtwistpredr,2)*pow(cos(twistpre),2)*Ieta - 2*pow(dtwistpredr,2)*pow(sin(twistpre),2)*Iksi + 
                        2*pow(dtwistpredr,2)*pow(sin(twistpre),2)*Ieta + 4*dtwistpredr*cos(twistpre)*sin(twistpre)*dIksidr + 
                        2*cos(twistpre)*sin(twistpre)*d2twistpredr2*Iksi - 4*dtwistpredr*cos(twistpre)*sin(twistpre)*dIetadr - 
                        2*cos(twistpre)*sin(twistpre)*d2twistpredr2*Ieta);

		
                // B_Fv4
//                blade->coefv_d3udr3[i] = 2*dtwistpredr*pow(sin(twistpre),2)*E*(Ieta - Iksi);
//                blade->coefv_d2udr2[i] = -2*dtwistpredr*pow(cos(twistpre),2)*dEdr*(Ieta - Iksi);
//                blade->coefv_d2udr2[i] += -pow(cos(twistpre),2)*d2twistpredr2*E*(Ieta - Iksi);
//                blade->coefv_d3udr3[i] += -2*dtwistpredr*pow(cos(twistpre),2)*E*(Ieta - Iksi);
//                blade->coefv_d2udr2[i] += 2*dtwistpredr*pow(sin(twistpre),2)*dEdr*(Ieta - Iksi);
//                blade->coefv_d2udr2[i] += pow(sin(twistpre),2)*d2twistpredr2*E*(Ieta - Iksi);
//                blade->coefv_d3udr3[i] += -2*cos(twistpre)*sin(twistpre)*dEdr*(Ieta - Iksi);
//                blade->coefv_d2udr2[i] += -d2Edr2*cos(twistpre)*sin(twistpre)*(Ieta - Iksi);
//                blade->coefv_d4udr4[i] = -cos(twistpre)*sin(twistpre)*E*(Ieta - Iksi);
//                blade->coefv_d2udr2[i] += -2*dtwistpredr*pow(cos(twistpre),2)*E*(dIetadr - dIksidr);
//                blade->coefv_d2udr2[i] += 2*dtwistpredr*pow(sin(twistpre),2)*E*(dIetadr - dIksidr);
//                blade->coefv_d3udr3[i] += -2*cos(twistpre)*sin(twistpre)*E*(dIetadr - dIksidr);
//                blade->coefv_d2udr2[i] += -2*cos(twistpre)*sin(twistpre)*dEdr*(dIetadr - dIksidr);
//                blade->coefv_d2udr2[i] += -cos(twistpre)*sin(twistpre)*E*(d2Ietadr2 - d2Iksidr2);
//                blade->coefv_d2udr2[i] += 4*pow(dtwistpredr,2)*cos(twistpre)*sin(twistpre)*E*(Ieta - Iksi);

// B_Fv4
                blade->coefv_d3udr3[i] = -2*dtwistpredr*pow(sin(twistpre),2)*E*(Ieta - Iksi);
                blade->coefv_d2udr2[i] = 2*dtwistpredr*pow(cos(twistpre),2)*dEdr*(Ieta - Iksi);
                blade->coefv_d2udr2[i] += pow(cos(twistpre),2)*d2twistpredr2*E*(Ieta - Iksi);
                blade->coefv_d3udr3[i] += 2*dtwistpredr*pow(cos(twistpre),2)*E*(Ieta - Iksi);
                blade->coefv_d2udr2[i] += -2*dtwistpredr*pow(sin(twistpre),2)*dEdr*(Ieta - Iksi);
                blade->coefv_d2udr2[i] += -pow(sin(twistpre),2)*d2twistpredr2*E*(Ieta - Iksi);
                blade->coefv_d3udr3[i] += 2*cos(twistpre)*sin(twistpre)*dEdr*(Ieta - Iksi);
                blade->coefv_d2udr2[i] += d2Edr2*cos(twistpre)*sin(twistpre)*(Ieta - Iksi);
                blade->coefv_d4udr4[i] = cos(twistpre)*sin(twistpre)*E*(Ieta - Iksi);
                blade->coefv_d2udr2[i] += 2*dtwistpredr*pow(cos(twistpre),2)*E*(dIetadr - dIksidr);
                blade->coefv_d2udr2[i] += -2*dtwistpredr*pow(sin(twistpre),2)*E*(dIetadr - dIksidr);
                blade->coefv_d3udr3[i] += 2*cos(twistpre)*sin(twistpre)*E*(dIetadr - dIksidr);
                blade->coefv_d2udr2[i] += 2*cos(twistpre)*sin(twistpre)*dEdr*(dIetadr - dIksidr);
                blade->coefv_d2udr2[i] += cos(twistpre)*sin(twistpre)*E*(d2Ietadr2 - d2Iksidr2);
                blade->coefv_d2udr2[i] += -4*pow(dtwistpredr,2)*cos(twistpre)*sin(twistpre)*E*(Ieta - Iksi);



		// C_Fv4
		blade->coefv_d2udr2[i] += - 2*dthetadr.at(i)*E*cos(2*twistpre)*(dIetadr - dIksidr);
		blade->coefv_d2vdr2[i] += - 2*dthetadr.at(i)*E*sin(2*twistpre)*(dIetadr - dIksidr);
		blade->coefv_d2udr2[i] += - 2*theta*dEdr*cos(2*twistpre)*(dIetadr - dIksidr); 
		blade->coefv_d2vdr2[i] += - 2*theta*dEdr*sin(2*twistpre)*(dIetadr - dIksidr);
		blade->coefv_d3udr3[i] += - 2*dthetadr.at(i)*E*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefv_d3vdr3[i] += - 2*dthetadr.at(i)*E*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefv_d2vdr2[i] += - 2*dthetadr.at(i)*E*(Ieta - Iksi)*2*dtwistpredr*cos(2*twistpre);
		blade->coefv_d2udr2[i] +=   2*dthetadr.at(i)*E*(Ieta - Iksi)*2*dtwistpredr*sin(2*twistpre);
		blade->coefv_d3udr3[i] += - 2*theta*dEdr*(Ieta - Iksi)*cos(2*twistpre);
		blade->coefv_d3vdr3[i] += - 2*theta*dEdr*(Ieta - Iksi)*sin(2*twistpre);
		blade->coefv_d2vdr2[i] += - 2*theta*dEdr*(Ieta - Iksi)*2*dtwistpredr*cos(2*twistpre);
		blade->coefv_d2udr2[i] +=   2*theta*dEdr*(Ieta - Iksi)*2*dtwistpredr*sin(2*twistpre);
		blade->coefv_d2udr2[i] += - E*theta*cos(2*twistpre)*(d2Ietadr2 - d2Iksidr2);
		blade->coefv_d2vdr2[i] += - E*theta*sin(2*twistpre)*(d2Ietadr2 - d2Iksidr2);
		blade->coefv_d3udr3[i] += - 2*E*theta*(dIetadr - dIksidr)*cos(2*twistpre); 
		blade->coefv_d3vdr3[i] += - 2*E*theta*(dIetadr - dIksidr)*sin(2*twistpre); 
		blade->coefv_d2vdr2[i] += - 2*E*theta*(dIetadr - dIksidr)*2*dtwistpredr*cos(2*twistpre);
		blade->coefv_d2udr2[i] +=   2*E*theta*(dIetadr - dIksidr)*2*dtwistpredr*sin(2*twistpre);
		blade->coefv_d4udr4[i] += - E*theta*(Ieta - Iksi)*cos(2*twistpre); 
		blade->coefv_d4vdr4[i] += - E*theta*(Ieta - Iksi)*sin(2*twistpre); 
		blade->coefv_d3vdr3[i] += - E*theta*(Ieta - Iksi)*4*dtwistpredr*cos(2*twistpre);
		blade->coefv_d2vdr2[i] += - E*theta*(Ieta - Iksi)*2*d2twistpredr2*cos(2*twistpre);
		blade->coefv_d3udr3[i] +=   E*theta*(Ieta - Iksi)*4*dtwistpredr*sin(2*twistpre);
		blade->coefv_d2udr2[i] +=   E*theta*(Ieta - Iksi)*2*d2twistpredr2*sin(2*twistpre);
		blade->coefv_d2udr2[i] +=   E*theta*(Ieta - Iksi)*4*pow(dtwistpredr,2)*cos(2*twistpre);
		blade->coefv_d2vdr2[i] +=   E*theta*(Ieta - Iksi)*4*pow(dtwistpredr,2)*sin(2*twistpre);
		blade->coefv_d2udr2[i] += - 2*dthetadr.at(i)*dEdr*cos(2*twistpre)*(Ieta - Iksi);
		blade->coefv_d2vdr2[i] += - 2*dthetadr.at(i)*dEdr*sin(2*twistpre)*(Ieta - Iksi);
		blade->coefv_d2udr2[i] += - d2Edr2*theta*cos(2*twistpre)*(Ieta - Iksi);
		blade->coefv_d2vdr2[i] += - d2Edr2*theta*sin(2*twistpre)*(Ieta - Iksi);
		blade->coefv_d2udr2[i] += - d2thetadr2.at(i)*E*cos(2*twistpre)*(Ieta - Iksi);
		blade->coefv_d2vdr2[i] += - d2thetadr2.at(i)*E*sin(2*twistpre)*(Ieta - Iksi);

//                double A_Fv4 = 0.0; 
//		double B_Fv4 = 0.0;
//		double C_Fv4 = 0.0;	
		double D_Fv4 = - 2*d3lpidr3*dthetadr.at(i)*E*(Ieta*pow(cos(theta),2) + Iksi*pow(sin(theta),2)) - 
			2*d3lpidr3*theta*dEdr*(Ieta*pow(cos(theta),2) +
			Iksi*pow(sin(theta),2)) - 2*dthetadr.at(i)*d2lpidr2*dEdr*(Ieta*pow(cos(theta),2) + Iksi*pow(sin(theta),2)) -
			d2Edr2*d2lpidr2*theta*(Ieta*pow(cos(theta),2) + Iksi*pow(sin(theta),2)) - d4lpidr4*E*theta*(Ieta*pow(cos(theta),2) +
			Iksi*pow(sin(theta),2)) - d2thetadr2.at(i)*E*d2lpidr2*(Ieta*pow(cos(theta),2) + Iksi*pow(sin(theta),2)) -
			2*d3lpidr3*E*theta*(pow(cos(theta),2)*dIetadr + pow(sin(theta),2)*dIksidr - 2*dthetadr.at(i)*cos(theta)*sin(theta)*Ieta +
			2*dthetadr.at(i)*cos(theta)*sin(theta)*Iksi) - 2*dthetadr.at(i)*E*d2lpidr2*(pow(cos(theta),2)*dIetadr + pow(sin(theta),2)*dIksidr -
			2*dthetadr.at(i)*cos(theta)*sin(theta)*Ieta + 2*dthetadr.at(i)*cos(theta)*sin(theta)*Iksi) -
			2*d2lpidr2*theta*dEdr*(pow(cos(theta),2)*dIetadr + pow(sin(theta),2)*dIksidr - 2*dthetadr.at(i)*cos(theta)*sin(theta)*Ieta +
			2*dthetadr.at(i)*cos(theta)*sin(theta)*Iksi) - E*d2lpidr2*theta*(d2Ietadr2*pow(cos(theta),2) + d2Iksidr2*pow(sin(theta),2) -
			2*pow(dthetadr.at(i),2)*pow(cos(theta),2)*Ieta + 2*pow(dthetadr.at(i),2)*pow(cos(theta),2)*Iksi + 
			2*pow(dthetadr.at(i),2)*pow(sin(theta),2)*Ieta -
			2*pow(dthetadr.at(i),2)*pow(sin(theta),2)*Iksi - 4*dthetadr.at(i)*cos(theta)*sin(theta)*dIetadr - 
			2*cos(theta)*sin(theta)*d2thetadr2.at(i)*Ieta +
			4*dthetadr.at(i)*cos(theta)*sin(theta)*dIksidr + 2*cos(theta)*sin(theta)*d2thetadr2.at(i)*Iksi);

		Fv4.at(i)=D_Fv4;/*+d2udr2.at(i)*blade->coefv_d2udr2[i] + d3udr3.at(i)*blade->coefv_d3udr3[i] + d4udr4.at(i)*blade->coefv_d4udr4[i]+ d2vdr2.at(i)*blade->coefv_d2vdr2[i] + d3vdr3.at(i)*blade->coefv_d3vdr3[i] + d4vdr4.at(i)*blade->coefv_d4vdr4[i]*/
		//Fv4.at(i)=B_Fv4+C_Fv4+D_Fv4;
		//Fv4.at(i)=A_Fv4;
		//cout << "FV4: " << i << " " << A_Fv4 << " " << B_Fv4 << " " << C_Fv4 << " " << D_Fv4 << endl;
//		Fv4.at(i)=0.0;
	}


//    WriteForces_Term("Fv4.dat",blade->NumberofElements, blade, Fv4);
	return(0);
}


/**
 * Compute Fu5 of eq. (17)
 */
PetscErrorCode computeFu5_bladestructure(Beam *blade, vector<double> &Fu5, double d2phidt2, vector<double> w0)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double beta = blade->beta;

		Fu5.at(i)=mass*d2phidt2*w0.at(i)*cos(beta);

//    PetscPrintf(PETSC_COMM_WORLD,"Compute Fu5 %d %e %e %e\n",i,w0.at(i),d2phidt2,beta);

	}
//    WriteForces_Term("Fu5.dat",blade->NumberofElements, blade, Fu5);

	return(0);
}


/**
 * Compute Fv5 of eq. (17)
 */
PetscErrorCode computeFv5_bladestructure(Beam *blade, vector<double> &Fv5, double d2phidt2, vector<double> w0)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double beta = blade->beta;

		Fv5.at(i)=-mass*d2phidt2*w0.at(i)*sin(beta);
	}

//    WriteForces_Term("Fv5.dat",blade->NumberofElements, blade, Fv5);
	return(0);
}

/**
 * Compute Ftheta0 of eq. (18)
 */ 
PetscErrorCode computeFtheta0_bladestructure(Beam *blade, vector<double> &Ftheta0, vector<double> d2udt2, vector<double> d2vdt2)
{
 	int i, j, k;


	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];

		Ftheta0.at(i)=-mass*lcg*(d2udt2.at(i)*sin(thetabar)-d2vdt2.at(i)*cos(thetabar));
	}

	return(0);
}

/**
 * Compute Ftheta1 of eq. (19)
 */ 
PetscErrorCode computeFtheta1_bladestructure(Beam *blade, vector<double> &Ftheta1, vector<double> dudr, vector<double> dvdr, vector<double> w0, vector<double> ucghat, double dphidt)
{
 	int i, j, k;


	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];
		double dlpidr = blade->dlpidr[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double beta = blade->beta;
		double phi = blade->phi;
		Ftheta1.at(i)=mass*lcg*pow(dphidt,2)*ucghat.at(i)*sin(thetabar+beta)+mass*lcg*w0.at(i)*pow(phi,2)*(dvdr.at(i)*cos(thetabar)-(dudr.at(i)+dlpidr)*sin(thetabar));
	}

	return(0);
}

/**
 * Compute Ftheta2 of eq. (20)
 */ 
PetscErrorCode computeFtheta2_bladestructure(Beam *blade, vector<double> &Ftheta2, vector<double> w0, double dbetadt, double d2betadt2)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];
		double lpi = blade->lpi[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double Icg = blade->Icg[i];
		double dlpidr = blade->dlpidr[i];
		double dudt = blade->dudt[i];
		double dvdt = blade->dvdt[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double w0 = blade->w0[i];
		double beta = blade->beta;
		double phi = blade->phi;
		Ftheta2.at(i)=(Icg+mass*pow(lcg,2))*d2betadt2-d2betadt2*mass*lcg*(u*cos(thetabar)+v*sin(thetabar))+mass*lcg*pow(dbetadt,2)*((u+lpi)*sin(thetabar)-v*cos(thetabar))+2*mass*lcg*dbetadt*(dudt*cos(thetabar)+dvdt*sin(thetabar));
	}

	return(0);
}

/**
 * Compute Ftheta3 of eq. (21)
 */ 
PetscErrorCode computeFtheta3_bladestructure(Beam *blade, vector<double> &Ftheta3, vector<double>w0, double d2phidt2)
{
 	int i, j, k;

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];
		double lpi = blade->lpi[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double Icg = blade->Icg[i];
		double dlpidr = blade->dlpidr[i];
		double dudt = blade->dudt[i];
		double dvdt = blade->dvdt[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double beta = blade->beta;
		Ftheta3.at(i)=-mass*d2phidt2*w0.at(i)*lcg*sin(thetabar+beta);
	}

	return(0);
}


/**
 * Compute Ftheta4 of eq. (22)
 */
PetscErrorCode computeFtheta4_bladestructure(Beam *blade, vector<double> &Ftheta4, vector<double> dudr, vector<double> dvdr)
{
 	int i, j, k;

	double gravity = sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

//	  PetscPrintf(PETSC_COMM_WORLD,"Send F theta4 to old\n");
//	for (i=0; i<blade->NumberofElements; i++) 
//	{
//		blade->Ftheta4_o[i]=blade->Ftheta4[i];
//	}

//	  PetscPrintf(PETSC_COMM_WORLD,"Compute new Ftheta4 \n");
	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double lcg = blade->lcg[i];
		double lpi = blade->lpi[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double Icg = blade->Icg[i];
		double dlpidr = blade->dlpidr[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double theta = blade->theta[i];
		double beta = blade->beta;
		double phi = blade->phi;
		Ftheta4.at(i)=-lcg*(sin(beta+thetabar)+theta*cos(beta+thetabar))*mass*gravity*sin(phi)+lcg*(dvdr.at(i)*cos(thetabar)-(dudr.at(i)+dlpidr)*sin(thetabar))*mass*gravity*cos(phi);
	}

	return(0);
}

/**
 * Compute Ftheta5 of eq. (23)
 */
PetscErrorCode computeFtheta5_bladestructure(Beam *blade, vector<double> &Ftheta5, vector<double> d2udr2, vector<double> d3udr3, vector<double> d2vdr2, vector<double> d3vdr3, vector<double> dthetadr, vector<double> d2thetadr2)
{
 	int i, j, k;

	double gravity = sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));

	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double E = blade->E[i];
		double dEdr = blade->dEdr[i];
		double Ietaetaksi = blade->Ietaetaksi[i];
		double dIetaetaksidr = blade->dIetaetaksidr[i];
		double Ietaksiksi = blade->Ietaksiksi[i];
		double dIetaksiksidr = blade->dIetaksiksidr[i];
		double Ieta = blade->Ieta[i];
		double Iksi = blade->Iksi[i];

		double lcg = blade->lcg[i];
		double lpi = blade->lpi[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double Icg = blade->Icg[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double theta = blade->theta[i];
		double dthetabardr = dthetadr.at(i);
		double twistpre = blade->twistpre[i];
		double dtwistpredr = blade->dtwistpredr[i];
		double d2twistpredr2 = blade->d2twistpredr2[i];
		double beta = blade->beta;
		double phi = blade->phi;

		//diff(E*Ietaetaksi*(dtwistpredr+dthetadr)*(d2udr2*sin(thetabar)-d2vdr2*cos(theta)))

		double A_Ftheta5 = - (
			E*Ietaetaksi*(dthetadr.at(i) + dtwistpredr)*(d3udr3.at(i)*sin(thetabar) - d3vdr3.at(i)*cos(theta) + dthetabardr*cos(thetabar)*d2udr2.at(i)
			+ dthetadr.at(i)*sin(theta)*d2vdr2.at(i)) - E*dIetaetaksidr*(dthetadr.at(i) + dtwistpredr)*(cos(theta)*d2vdr2.at(i) - sin(thetabar)*d2udr2.at(i)) -
			Ietaetaksi*dEdr*(dthetadr.at(i) + dtwistpredr)*(cos(theta)*d2vdr2.at(i) - sin(thetabar)*d2udr2.at(i)) - E*Ietaetaksi*(cos(theta)*d2vdr2.at(i)
			- sin(thetabar)*d2udr2.at(i))*(d2thetadr2.at(i) + d2twistpredr2)
			);

		// diff(E*Ietaksiksi*(dtwistpredr+dthetadr)*(d2udr2*cos(twistpre)+d2vdr2*sin(twistpre)))

		double B_Ftheta5 = E*Ietaksiksi*(dthetadr.at(i) + dtwistpredr)*(d3udr3.at(i)*cos(twistpre) + d3vdr3.at(i)*sin(twistpre) +
			dtwistpredr*cos(twistpre)*d2vdr2.at(i) - dtwistpredr*sin(twistpre)*d2udr2.at(i)) + E*dIetaksiksidr*(dthetadr.at(i) +
			dtwistpredr)*(cos(twistpre)*d2udr2.at(i) + sin(twistpre)*d2vdr2.at(i)) + Ietaksiksi*dEdr*(dthetadr.at(i) +
			dtwistpredr)*(cos(twistpre)*d2udr2.at(i) + sin(twistpre)*d2vdr2.at(i)) + E*Ietaksiksi*(cos(twistpre)*d2udr2.at(i) +
			sin(twistpre)*d2vdr2.at(i))*(d2thetadr2.at(i) + d2twistpredr2);

		double C_Ftheta5 = -(E*Iksi-E*Ieta)*((pow(d2udr2.at(i),2)-pow(d2vdr2.at(i),2))*cos(twistpre)*sin(twistpre)-d2udr2.at(i)*d2vdr2.at(i)*cos(2*twistpre)); 

		double D_Ftheta5 = -E*Iksi*d2lpidr2*(d2udr2.at(i)*cos(twistpre)+d2vdr2.at(i)*sin(twistpre))*sin(twistpre);
		
		double E_Ftheta5 = E*Ieta*d2lpidr2*(d2udr2.at(i)*sin(twistpre)-d2vdr2.at(i)*cos(twistpre))*cos(twistpre);

		Ftheta5.at(i)=A_Ftheta5+B_Ftheta5+C_Ftheta5+D_Ftheta5+E_Ftheta5;


//                WriteForces_Term("Ftheta5.dat",blade->NumberofElements, blade, Ftheta5);

	}

	return(0);
}

/**
 * Compute Ftheta6 of eq. (24)
 */
PetscErrorCode computeFtheta6_bladestructure(Beam *blade, vector<double> &Ftheta6, vector<double> d2udr2, vector<double> d3udr3, vector<double> dvdr, vector<double> d2vdr2, vector<double> d3vdr3, vector<double> dthetadr, vector<double> d2thetadr2)

{
 	int i, j, k;

	double gravity = sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2));


	for (i=0; i<blade->NumberofElements; i++) 
	{
		double mass = blade->mass[i];
		double E = blade->E[i];
		double dEdr = blade->dEdr[i];
		double Ietaetaksi = blade->Ietaetaksi[i];
		double dIetaetaksidr = blade->dIetaetaksidr[i];
		double Ietaksiksi = blade->Ietaksiksi[i];
		double dIetaksiksidr = blade->dIetaksiksidr[i];
		double Ieta = blade->Ieta[i];
		double Iksi = blade->Iksi[i];

		double lcg = blade->lcg[i];
		double lpi = blade->lpi[i];
		double u = blade->u[i];
		double v = blade->v[i];
		double Icg = blade->Icg[i];
		double dlpidr = blade->dlpidr[i];
		double d2lpidr2 = blade->d2lpidr2[i];
		double d3lpidr3 = blade->d3lpidr3[i];
		double thetabar = blade->twistpre[i]+blade->theta[i];
		double theta = blade->theta[i];
		double twistpre = blade->twistpre[i];
		double dtwistpredr = blade->dtwistpredr[i];
		double d2twistpredr2 = blade->d2twistpredr2[i];
		double beta = blade->beta;
		double phi = blade->phi;
		double J = blade->J[i];
		double dJdr = blade->dJdr[i];
		double G = blade->G[i];
		double dGdr = blade->dGdr[i];
			
		// diff(G*J*(dthetadr+dvdr*(d2udr2+d2lpidr2)))

		Ftheta6.at(i) = -(
			dJdr*G*(dthetadr.at(i) + dvdr.at(i)*(d2lpidr2 + d2udr2.at(i))) + J*dGdr*(dthetadr.at(i) + dvdr.at(i)*(d2lpidr2 + d2udr2.at(i))) +
			G*J*(d2thetadr2.at(i) + d2vdr2.at(i)*(d2lpidr2 + d2udr2.at(i)) + dvdr.at(i)*(d3lpidr3 + d3udr3.at(i))));
	}

//                WriteForces_Term("Ftheta6.dat",blade->NumberofElements, blade, Ftheta6);
	return(0);
}




/**
 * compute the right-hand-side of the direct solver
 */
PetscErrorCode computerhsDS_bladestructure(Beam *blade)
{
	using namespace std;
 	int i, j, k;

	int NElmts=blade->NumberofElements;

	// compute temporal derivatives, phi, beta
	double dbetadt, d2betadt2, dphidt, d2phidt2;
	compute_temporalderivatives_beta_bladestructure(blade, dbetadt, d2betadt2);

   dphidt = blade->dphidt;
   d2phidt2 = blade->d2phidt2;
// This function is deprecated. Angular velocity and acceleration obtained from turbine control
//	compute_temporalderivatives_phi_bladestructure(blade, dphidt, d2phidt2);

	// compute temporal derivatives 
	vector<double> d2udt2, d2vdt2, d2thetadt2; 
	d2udt2.resize(NElmts);
	d2vdt2.resize(NElmts);
	d2thetadt2.resize(NElmts);

	compute_temporalderivatives_u_bladestructure(blade, d2udt2);
	compute_temporalderivatives_v_bladestructure(blade, d2vdt2);

//	PetscPrintf(PETSC_COMM_WORLD,"Compute temporal derivatives theta \n");
	if (torsion_turbinestructure) compute_temporalderivatives_theta_bladestructure(blade, d2thetadt2);

	// compute spatial derivatives of dudt, dvdt
	vector<double> ddt_dudr, ddt_dvdr;
	ddt_dudr.resize(NElmts);
	ddt_dvdr.resize(NElmts);
	compute_spatialderivatives_dudtdvdt_bladestructure(blade, ddt_dudr, ddt_dvdr);

	// compute spatial derivatives 
	vector<double> dudr, d2udr2, d3udr3, d4udr4;
	dudr.resize(NElmts);
	d2udr2.resize(NElmts);
	d3udr3.resize(NElmts);
	d4udr4.resize(NElmts);

	vector<double> dvdr, d2vdr2, d3vdr3, d4vdr4;
	dvdr.resize(NElmts);
	d2vdr2.resize(NElmts);
	d3vdr3.resize(NElmts);
	d4vdr4.resize(NElmts);

	vector<double> dthetadr, d2thetadr2;
	dthetadr.resize(NElmts);
	d2thetadr2.resize(NElmts);

	compute_spatialderivatives_u_bladestructure(blade, dudr, d2udr2, d3udr3, d4udr4);
	compute_spatialderivatives_v_bladestructure(blade, dvdr, d2vdr2, d3vdr3, d4vdr4);

//	PetscPrintf(PETSC_COMM_WORLD,"Compute spatial derivatives theta \n");
	if (torsion_turbinestructure) compute_spatialderivatives_theta_bladestructure(blade, dthetadr, d2thetadr2); 

	// compute w0, w1
	vector<double> w0, w1;
	w0.resize(NElmts);
	w1.resize(NElmts);

	computew0w1_bladestructure(blade, w0, w1, dudr, dvdr);

	// compute dw0dr
	vector<double> dw0dr;
	dw0dr.resize(NElmts);
	compute_dw0dr_bladestructure(blade, dw0dr, w0);

	// compute Fu0
	vector<double> Fu0;
	Fu0.resize(NElmts);
	computeFu0_bladestructure(blade, Fu0, d2thetadt2);


	// compute Fv0
	vector<double> Fv0;
	Fv0.resize(NElmts);
	computeFv0_bladestructure(blade, Fv0, d2thetadt2);

	computeucgvcg4Fu1Fv1_bladestructure(blade);
	
	// compute ducgdt, dvcgdt for Fu1, Fv1
	vector<double> ducgdt, dvcgdt;
	ducgdt.resize(NElmts);
	dvcgdt.resize(NElmts);
	computeucgvcgtimederivative4Fu1Fv1_bladestructure(blade, ducgdt, dvcgdt);

	// compute T1, T1integral
	vector<double> T1, T1integral;
	T1.resize(NElmts);
	T1integral.resize(NElmts);
	computeT1integral_bladestructure(blade, T1, T1integral, dbetadt, dphidt);

	// compute dT1dr, dT1drintegral
	vector<double> dT1dr, dT1drintegral;
	dT1dr.resize(NElmts);
	dT1drintegral.resize(NElmts);
	computedT1drintegral_bladestructure(blade, dT1dr, dT1drintegral, T1, dbetadt, dphidt, dudr, dvdr);

	// compute Fu1
	vector<double> Fu1;
	Fu1.resize(NElmts);
//	computeFu1_bladestructure(blade, Fu1, dbetadt, d2betadt2, dphidt, dudr, d2udr2, ducgdt, dvcgdt, T1, T1integral, dT1dr, dT1drintegral, dthetadr);

	// compute Fv1
	vector<double> Fv1;
	Fv1.resize(NElmts);
//	computeFv1_bladestructure(blade, Fv1, dbetadt, d2betadt2, dphidt, dvdr, d2vdr2, ducgdt, dvcgdt, T1, T1integral, dT1dr, dT1drintegral, dthetadr);

	// compute T2, T2integral
	vector<double> T2, T2integral;
	T2.resize(NElmts);
	T2integral.resize(NElmts);
	computeT2integral_bladestructure(blade, T2, T2integral, dbetadt, dphidt);

	// compute dT2dr, dT2drintegral
	vector<double> dT2dr, dT2drintegral;
	dT2dr.resize(NElmts);
	dT2drintegral.resize(NElmts);
	computedT2drintegral_bladestructure(blade, dT2dr, dT2drintegral, dbetadt, dphidt, ddt_dudr, ddt_dvdr, T2);

	// compute T3, T3integral
	vector<double> T3, T3integral;
	T3.resize(NElmts);
	T3integral.resize(NElmts);
	computeT3integral_bladestructure(blade, T3, T3integral, dbetadt, dphidt, w0);

	// compute T3, T3integral
	// compute dT3dr, dT3drintegral
	vector<double> dT3dr, dT3drintegral;
	dT3dr.resize(NElmts);
       	dT3drintegral.resize(NElmts);	
	computedT3drintegral_bladestructure(blade, dT3dr, dT3drintegral, dbetadt, dphidt, ddt_dudr, ddt_dvdr, w0, dw0dr, T3);

	// compute T3, T3integral
	// compute ucghat for Fu2, Fv2
	vector<double> ucghat;
	ucghat.resize(NElmts);
	computeucghat4Fu2Fv2_bladestructure(blade, ucghat);

  resetMatrixCoeff(blade);

	// compute T3, T3integral
	// compute Fu2
	vector<double> Fu2;
	Fu2.resize(NElmts);
	computeFu2_bladestructure(blade, Fu2, dphidt, w0, dw0dr, dudr, ddt_dudr, d2udr2, ddt_dvdr, ucghat, T2, dT2dr, T2integral, dT2drintegral, T3, dT3dr, T3integral, dT3drintegral, dthetadr);

	// compute T3, T3integral
	// compute Fv2
	vector<double> Fv2;
	Fv2.resize(NElmts);
	computeFv2_bladestructure(blade, Fv2, dphidt, w0, dw0dr, dvdr, ddt_dudr, d2vdr2, ddt_dvdr, ucghat, T2, dT2dr, T2integral, dT2drintegral, T3, dT3dr, T3integral, dT3drintegral, dthetadr);

	// compute T4, T4integral
	vector<double> T4, T4integral;
	T4.resize(NElmts);
	T4integral.resize(NElmts);
	computeT4integral_bladestructure(blade, T4, T4integral);

	// compute dT4dr, dT4drintegral
	vector<double> dT4dr, dT4drintegral;
	dT4dr.resize(NElmts);
	dT4drintegral.resize(NElmts);
	computedT4drintegral_bladestructure(blade, dT4dr, dT4drintegral, T4);

	// compute Fu3
	vector<double> Fu3;
	Fu3.resize(NElmts);
	computeFu3_bladestructure(blade, Fu3, dudr, d2udr2, dvdr, d2vdr2, T4integral, dT4drintegral, dthetadr);

	// compute Fv3
	vector<double> Fv3;
	Fv3.resize(NElmts);
	computeFv3_bladestructure(blade, Fv3, dudr, d2udr2, dvdr, d2vdr2, T4integral, dT4drintegral, dthetadr);

	// compute Fu4
	vector<double> Fu4;
	Fu4.resize(NElmts);
	computeFu4_bladestructure(blade, Fu4, d2udr2, d3udr3, d4udr4, d2vdr2, d3vdr3, d4vdr4, dthetadr, d2thetadr2);

	// compute Fv4
	vector<double> Fv4;
	Fv4.resize(NElmts);
	computeFv4_bladestructure(blade, Fv4, d2udr2, d3udr3, d4udr4, d2vdr2, d3vdr3, d4vdr4, dthetadr, d2thetadr2);

	// compute Fu5
	vector<double> Fu5;
	Fu5.resize(NElmts);
	computeFu5_bladestructure(blade, Fu5, d2phidt2, w0);

	// compute Fv5
	vector<double> Fv5;
	Fv5.resize(NElmts);
	computeFv5_bladestructure(blade, Fv5, d2phidt2, w0);

	// torsion
	vector<double> Ftheta0, Ftheta1, Ftheta2, Ftheta3, Ftheta4, Ftheta5, Ftheta6;

//	PetscPrintf(PETSC_COMM_WORLD,"Compute F theta \n");
	if (torsion_turbinestructure) 
	{
//	  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta0 \n");
		// compute Ftheta0
		Ftheta0.resize(NElmts);
//		computeFtheta0_bladestructure(blade, Ftheta0, d2udt2, d2vdt2);

//	  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta1 \n");
		// compute Ftheta0
		// compute Ftheta1
		Ftheta1.resize(NElmts);
//		computeFtheta1_bladestructure(blade, Ftheta1, dudr, dvdr, w0, ucghat, dphidt);

//	  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta2 \n");
		// compute Ftheta2
		Ftheta2.resize(NElmts);
//		computeFtheta2_bladestructure(blade, Ftheta2, w0, dbetadt, d2betadt2);

//	  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta3 \n");
		// compute Ftheta3
		Ftheta3.resize(NElmts);
//		computeFtheta3_bladestructure(blade, Ftheta3, w0, d2phidt2);

//	  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta4 \n");
		// compute Ftheta4
		Ftheta4.resize(NElmts);
//		computeFtheta4_bladestructure(blade, Ftheta4, dudr, dvdr);


//	  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta5 \n");
		// compute Ftheta5
		Ftheta5.resize(NElmts);
		computeFtheta5_bladestructure(blade, Ftheta5, d2udr2, d3udr3, d2vdr2, d3vdr3, dthetadr, d2thetadr2);

	//  PetscPrintf(PETSC_COMM_WORLD,"Compute F theta6 \n");
		// compute Ftheta6
		Ftheta6.resize(NElmts);
		computeFtheta6_bladestructure(blade, Ftheta6, d2udr2, d3udr3, dvdr, d2vdr2, d3vdr3, dthetadr, d2thetadr2);
	}

//	PetscPrintf(PETSC_COMM_WORLD,"Compute spatial derivatives F theta end\n");

	//cout << "rhsDS\n";
	for (i=0; i<NElmts; i++) 
	{
		double dt=dt_turbinestructure;
		double fac=0.5;

		double _mass = 1.0/(blade->mass[i]+1.e-19);
		double dudt=blade->dudt[i];
		double dudt_o=blade->dudt_o[i];
		blade->rhsdudt[i] = dudt_o 
			-dt*Fu0.at(i)-dt*_mass*(Fu1.at(i)+Fu2.at(i)+Fu3.at(i)+Fu4.at(i)+Fu5.at(i)-blade->fu[i]); 
//			-dt*Fu0.at(i)-dt*_mass*(Fu1.at(i)+Fu2.at(i)+Fu3.at(i)+Fu4.at(i)+Fu5.at(i)); 
//		blade->rhsdudt[i] = dudt_o+ 
//			fac*(-dt*Fu0.at(i)-dt*_mass*(Fu1.at(i)+Fu2.at(i)+Fu3.at(i)+Fu4.at(i)+Fu5.at(i)-blade->fu[i])) + (1-fac)*blade->rhsdudt_o[i]; 
//
                
//	        PetscPrintf(PETSC_COMM_WORLD,"Fu and Fv: %le %le %le\n",Fu3.at(i),Fv3.at(i),blade->fu[i]);
//	        PetscPrintf(PETSC_COMM_WORLD,"Fu and Fv: %le %le %le %le \n",blade->fu[i],blade->fv[i],Fu4.at(i),Fv4.at(i));

		double u_o=blade->u_o[i];
		blade->rhsu[i] = u_o+blade->rhsdudt[i]*dt; 


//		cout << "~~~~ " << blade->rhsdudt[i] << " " << Fu0.at(i) << " " << Fu1.at(i) << " " << Fu2.at(i) << " " << Fu3.at(i) << " " << Fu4.at(i) << " " << Fu5.at(i) << " " << blade->fu[i] << endl;

		double dvdt=blade->dvdt[i];
		double dvdt_o=blade->dvdt_o[i];
		blade->rhsdvdt[i] = dvdt_o 
	    -dt*Fv0.at(i)-dt*_mass*(Fv1.at(i)+Fv2.at(i)+Fv3.at(i)+Fv4.at(i)+Fv5.at(i)-blade->fv[i]); 
//			-dt*Fv0.at(i)-dt*_mass*(Fv1.at(i)+Fv2.at(i)+Fv3.at(i)+Fv4.at(i)+Fv5.at(i)); 
//			-dt*Fv0.at(i)-dt*_mass*(Fv1.at(i)+Fv2.at(i)+Fv3.at(i)+Fv4.at(i)+Fv5.at(i)-blade->fv[i]); 
//		blade->rhsdvdt[i] = dvdt_o+ 
//			fac*(-dt*Fv0.at(i)-dt*_mass*(Fv1.at(i)+Fv2.at(i)+Fv3.at(i)+Fv4.at(i)+Fv5.at(i)-blade->fv[i])) + (1-fac)*blade->rhsdvdt_o[i]; 

		double v_o=blade->v_o[i];
		blade->rhsv[i] = v_o+blade->rhsdvdt[i]*dt; 

//	        PetscPrintf(PETSC_COMM_WORLD,"rhsu and rhsv: %le %le %le %le \n",blade->rhsdudt[i],blade->rhsdvdt[i], blade->rhsu[i], blade->rhsv[i]);

//		cout << "++++ " << blade->rhsdvdt[i] << " " << Fv0.at(i) << " " << Fv1.at(i) << " " << Fv2.at(i) << " " << Fv3.at(i) << " " << Fv4.at(i) << " " << Fv5.at(i) << " " << blade->fv[i] << endl;
		// not changed for DS yet
//	PetscPrintf(PETSC_COMM_WORLD,"Compute RHS theta end\n");
		if (torsion_turbinestructure) 
		{
	                double dthetadt=blade->dthetadt[i];
        	        double dthetadt_o=blade->dthetadt_o[i];
			double _Icg = 1.0/(blade->Icg[i]+1.e-19);
                	blade->rhsdthetadt[i] = dthetadt_o
				-dt*_Icg*(Ftheta0.at(i)+Ftheta1.at(i)+Ftheta2.at(i)+Ftheta3.at(i)+Ftheta4.at(i)+Ftheta5.at(i)+Ftheta6.at(i)); 

			double theta=blade->theta[i];
			double theta_o=blade->theta_o[i];
			blade->rhstheta[i] = theta_o+dt*blade->rhsdthetadt[i]; 
		}
	}

//	int itmp;
//	cin >> itmp;

	i=0;
	blade->rhsdudt[i] = 0.0;
	blade->rhsu[i] = 0.0;
	blade->rhsdvdt[i] = 0.0;
	blade->rhsv[i] = 0.0;
	
	return 0;

}

/**
 * contruct the matrix of the LU solver
 */
PetscErrorCode initmatrixA_bladestructure(Beam *blade)
{
	int i, j, k;

	int NElmts=blade->NumberofElements;

	for (i=0;i<4*NElmts;i++)
	for (j=0;j<4*NElmts;j++)
	{
		blade->A[i][j]=0.0;
	}

}

/**
 * contruct the matrix of the LU solver
 */
PetscErrorCode constructmatrixA_bladestructure(Beam *blade)
{
	using namespace std;
 	int i, j, k;

	int NElmts=blade->NumberofElements;

	for (i=0;i<4*NElmts;i++)
	for (j=0;j<4*NElmts;j++)
	{
		blade->A[i][j]=0.0;
	}

	double dt=dt_turbinestructure;
	double _h=1.0/blade->dh;
	double fac=1.0; // CN scheme

	for (i=1;i<NElmts;i++)
	{
		double Iksi = blade->Iksi[i];
		double Ieta = blade->Ieta[i];
		double E = blade->E[i];
		double twistpre = blade->twistpre[i];
		double _mass = 1.0/(blade->mass[i]+1.e-19);

		double coefu = fac*dt*_mass*E*(Iksi*pow(cos(twistpre),2)+Ieta*pow(sin(twistpre),2))*pow(_h,4);
		double coefv = fac*dt*_mass*E*(Ieta*pow(cos(twistpre),2)+Iksi*pow(sin(twistpre),2))*pow(_h,4);

		double coef1 = fac*dt*_mass;

		double coefu_dudr   = coef1*blade->coefu_dudr[i]*pow(_h,1);	
		double coefu_d2udr2 = coef1*blade->coefu_d2udr2[i]*pow(_h,2);	
		double coefu_d3udr3 = coef1*blade->coefu_d3udr3[i]*pow(_h,3);	
		double coefu_d4udr4 = coef1*blade->coefu_d4udr4[i]*pow(_h,4);	

		double coefu_dvdr   = coef1*blade->coefu_dvdr[i]*pow(_h,1);	
		double coefu_d2vdr2 = coef1*blade->coefu_d2vdr2[i]*pow(_h,2);	
		double coefu_d3vdr3 = coef1*blade->coefu_d3vdr3[i]*pow(_h,3);	
		double coefu_d4vdr4 = coef1*blade->coefu_d4vdr4[i]*pow(_h,4);	

		double coefv_dudr   = coef1*blade->coefv_dudr[i]*pow(_h,1);	
		double coefv_d2udr2 = coef1*blade->coefv_d2udr2[i]*pow(_h,2);	
		double coefv_d3udr3 = coef1*blade->coefv_d3udr3[i]*pow(_h,3);	
		double coefv_d4udr4 = coef1*blade->coefv_d4udr4[i]*pow(_h,4);	

		double coefv_dvdr   = coef1*blade->coefv_dvdr[i]*pow(_h,1);	
		double coefv_d2vdr2 = coef1*blade->coefv_d2vdr2[i]*pow(_h,2);	
		double coefv_d3vdr3 = coef1*blade->coefv_d3vdr3[i]*pow(_h,3);	
		double coefv_d4vdr4 = coef1*blade->coefv_d4vdr4[i]*pow(_h,4);	

//    if (i==NElmts)
//    {
//       coefv_d2vdr2 = 0.;
//       coefv_d3vdr3 = 0.;
//       coefv_d2udr2 = 0.;
//       coefv_d3udr3 = 0.;
//
//       coefu_d2udr2 = 0.;
//       coefu_d3udr3 = 0.;
//       coefu_d2vdr2 = 0.;
//       coefu_d3vdr3 = 0.;
//    }

//	 PetscPrintf(PETSC_COMM_WORLD,"Debug Blade 1: %le %le %le %le \n",blade->coefu_dudr[i],blade->coefu_d2udr2[i],blade->coefu_d3udr3[i],blade->coefu_d4udr4[i]);
	 //PetscPrintf(PETSC_COMM_WORLD,"Debug Blade 2: %le %le %le %le \n",blade->coefu_dvdr[i],blade->coefu_d2vdr2[i],blade->coefu_d3vdr3[i],blade->coefu_d4vdr4[i]);

		/*
		coefu_dudr   = 0.0;	
                coefu_d2udr2 = 0.0;
                coefu_d3udr3 = 0.0;
                coefu_d4udr4 = coefu; 
                               
                coefu_dvdr   = 0.0; 
                coefu_d2vdr2 = 0.0;
                coefu_d3vdr3 = 0.0;
                coefu_d4vdr4 = 0.0;
                               
                coefv_dudr   = 0.0;
                coefv_d2udr2 = 0.0;
                coefv_d3udr3 = 0.0;
                coefv_d4udr4 = 0.0;
                               
                coefv_dvdr   = 0.0;
	        coefv_d2vdr2 = 0.0;
                coefv_d3vdr3 = 0.0;
                coefv_d4vdr4 = coefv;
		*/

		//if (i==8) cout << coefu << "\t" << coefv << endl;
		//if (i==8) cout << fac << " " << dt << " " << _mass<< " "  << E << " " << Iksi << " " << Ieta << " " << twistpre << " " << _h << endl;
		if (i==1) 
		{
			/*** dudt ***/
			// d4udr4
			blade->A[i][i]			+= 1.0;	
			blade->A[i][i+NElmts-1]		+= -4*coefu_d4udr4;
			blade->A[i][i+NElmts]		+= 7*coefu_d4udr4;
			blade->A[i][i+NElmts+1]		+= -4*coefu_d4udr4;
			blade->A[i][i+NElmts+2]		+= 1*coefu_d4udr4;

			// d4vdr4
			blade->A[i][i+3*NElmts-1]	+= -4*coefu_d4vdr4;
			blade->A[i][i+3*NElmts]		+= 7*coefu_d4vdr4;
			blade->A[i][i+3*NElmts+1]	+= -4*coefu_d4vdr4;
			blade->A[i][i+3*NElmts+2]	+= 1*coefu_d4vdr4;

			// d3udr3
			blade->A[i][i+NElmts-1]		+= 1*coefu_d3udr3;
			blade->A[i][i+NElmts]		+= -0.5*coefu_d3udr3;
			blade->A[i][i+NElmts+1]		+= -1*coefu_d3udr3;
			blade->A[i][i+NElmts+2]		+= 0.5*coefu_d3udr3;

			// d3vdr3
			blade->A[i][i+3*NElmts-1]	+= 1*coefu_d3vdr3;
			blade->A[i][i+3*NElmts]		+= -0.5*coefu_d3vdr3;
			blade->A[i][i+3*NElmts+1]	+= -1*coefu_d3vdr3;
			blade->A[i][i+3*NElmts+2]	+= 0.5*coefu_d3vdr3;

			/*** u ***/
			// d4udr4
//			blade->A[i+NElmts][i+NElmts]=1;
//			blade->A[i+NElmts][i]=-dt;
			blade->A[i+NElmts][i+NElmts-1]		+= -4*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts]		+= 1+7*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts+1]		+= -4*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts+2]		+= 1*coefu_d4udr4*dt;

			// d3vdr4
			blade->A[i+NElmts][i+3*NElmts-1]	+= -4*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= 7*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts+1]	+= -4*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts+2]	+= 1*coefu_d4vdr4*dt;

			// d3udr3
			blade->A[i+NElmts][i+NElmts-1]		+= 1*coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts]		+= -0.5*coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts+1]		+= -1*coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts+2]		+= 0.5*coefu_d3udr3*dt;

			// d3vdr3
			blade->A[i+NElmts][i+3*NElmts-1]	+= 1*coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= -0.5*coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts+1]	+= -1*coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts+2]	+= 0.5*coefu_d3vdr3*dt;

			/*** dvdt ***/
			// d4vdr4
			blade->A[i+2*NElmts][i+2*NElmts]		+= 1.0;	
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= -4*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= 7*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+1]	+= -4*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+2]	+= 1*coefv_d4vdr4;

			// d4udr4
			blade->A[i+2*NElmts][i+NElmts-1]		+= -4*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts]			+= 7*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts+1]		+= -4*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts+2]		+= 1*coefv_d4udr4;

			// d3vdr3 
                        blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= 1*coefv_d3vdr3;
                        blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= -0.5*coefv_d3vdr3;
                        blade->A[i+2*NElmts][i+2*NElmts+NElmts+1]	+= -1*coefv_d3vdr3;
                        blade->A[i+2*NElmts][i+2*NElmts+NElmts+2]	+= 0.5*coefv_d3vdr3;	

			// d3udr3
                        blade->A[i+2*NElmts][i+NElmts-1]		+= 1*coefv_d3udr3;
                        blade->A[i+2*NElmts][i+NElmts]			+= -0.5*coefv_d3udr3;
                        blade->A[i+2*NElmts][i+NElmts+1]		+= -1*coefv_d3udr3;
                        blade->A[i+2*NElmts][i+NElmts+2]		+= 0.5*coefv_d3udr3;	

			/*** v ***/
			// d4vdr4
//			blade->A[i+3*NElmts][i+3*NElmts]=1;
//			blade->A[i+3*NElmts][i+2*NElmts]=-dt;
                        blade->A[i+3*NElmts][i+3*NElmts-1] 	+= -4*coefv_d4vdr4*dt;
                        blade->A[i+3*NElmts][i+3*NElmts]	+= 1+7*coefv_d4vdr4*dt;
                        blade->A[i+3*NElmts][i+3*NElmts+1]	+= -4*coefv_d4vdr4*dt;
                        blade->A[i+3*NElmts][i+3*NElmts+2]	+= 1*coefv_d4vdr4*dt;

			// d4udr4
                        blade->A[i+3*NElmts][i+NElmts-1] 	+= -4*coefv_d4udr4*dt;
                        blade->A[i+3*NElmts][i+NElmts]		+= 7*coefv_d4udr4*dt;
                        blade->A[i+3*NElmts][i+NElmts+1]	+= -4*coefv_d4udr4*dt;
                        blade->A[i+3*NElmts][i+NElmts+2]	+= 1*coefv_d4udr4*dt;

			// d3vdr3
                        blade->A[i+3*NElmts][i+3*NElmts-1]	+= 1*coefv_d3vdr3*dt;
                        blade->A[i+3*NElmts][i+3*NElmts]	+= -0.5*coefv_d3vdr3*dt;
                        blade->A[i+3*NElmts][i+3*NElmts+1]	+= -1*coefv_d3vdr3*dt;
                        blade->A[i+3*NElmts][i+3*NElmts+2]	+= 0.5*coefv_d3vdr3*dt;	

			// d3udr3
                        blade->A[i+3*NElmts][i+NElmts-1]	+= 1*coefv_d3udr3*dt;
                        blade->A[i+3*NElmts][i+NElmts]		+= -0.5*coefv_d3udr3*dt;
                        blade->A[i+3*NElmts][i+NElmts+1]	+= -1*coefv_d3udr3*dt;
                        blade->A[i+3*NElmts][i+NElmts+2]	+= 0.5*coefv_d3udr3*dt;	
		}
		else if (i==NElmts-2)
		{
			/** dudt **/
			// d4udr4
			blade->A[i][i]			+= 1.0;	
			blade->A[i][i+NElmts-2]		+= 1*coefu_d4udr4;
			blade->A[i][i+NElmts-1]		+= -4*coefu_d4udr4;
			blade->A[i][i+NElmts]		+= 5*coefu_d4udr4;
			blade->A[i][i+NElmts+1]		+= -2*coefu_d4udr4;

			// d4vdr4
			blade->A[i][i+3*NElmts-2]	+= 1*coefu_d4vdr4;
			blade->A[i][i+3*NElmts-1]	+= -4*coefu_d4vdr4;
			blade->A[i][i+3*NElmts]		+= 5*coefu_d4vdr4;
			blade->A[i][i+3*NElmts+1]	+= -2*coefu_d4vdr4;

			// d3udr3
			blade->A[i][i+NElmts-2] 	+= -0.5*coefu_d3udr3;
			blade->A[i][i+NElmts-1] 	+= 1*coefu_d3udr3;
			blade->A[i][i+NElmts]		+= -0.5*coefu_d3udr3;

			// d3vdr3
			blade->A[i][i+3*NElmts-2] 	+= -0.5*coefu_d3vdr3;
			blade->A[i][i+3*NElmts-1] 	+= 1*coefu_d3vdr3;
			blade->A[i][i+3*NElmts]		+= -0.5*coefu_d3vdr3;

			/*** u ***/
			// d4udr4
//			blade->A[i+NElmts][i+NElmts]=1;
//			blade->A[i+NElmts][i]=-dt;
			blade->A[i+NElmts][i+NElmts-2] 		+= 1*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts-1] 		+= -4*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts]  		+= 1+5*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts+1]		+= -2*coefu_d4udr4*dt;

			// d4vdr4
			blade->A[i+NElmts][i+3*NElmts-2] 	+= 1*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts-1] 	+= -4*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts]  	+= 5*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts+1]	+= -2*coefu_d4vdr4*dt;

			// d3udr3
			blade->A[i+NElmts][i+NElmts-2] 		+= -0.5*coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts-1] 		+= 1*coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts]		+= -0.5*coefu_d3udr3*dt;

			// d3vdr3
			blade->A[i+NElmts][i+3*NElmts-2] 	+= -0.5*coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts-1] 	+= 1*coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= -0.5*coefu_d3vdr3*dt;

			/*** dvdt ***/
			// d4vdr4
			blade->A[i+2*NElmts][i+2*NElmts]		+= 1.0;	
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-2]	+= 1*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= -4*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= 5*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+1]	+= -2*coefv_d4vdr4;
		
			// d4udr4
			blade->A[i+2*NElmts][i+NElmts-2]		+= 1*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts-1]		+= -4*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts]			+= 5*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts+1]		+= -2*coefv_d4udr4;

			// d3vdr3
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-2]	+= -0.5*coefv_d3vdr3;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= 1*coefv_d3vdr3;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= -0.5*coefv_d3vdr3;

			// d3udr3
			blade->A[i+2*NElmts][i+NElmts-2]		+= -0.5*coefv_d3udr3;
			blade->A[i+2*NElmts][i+NElmts-1]		+= 1*coefv_d3udr3;
			blade->A[i+2*NElmts][i+NElmts]			+= -0.5*coefv_d3udr3;

			//
			/*** v ***/
			// d4vdr4
//			blade->A[i+3*NElmts][i+3*NElmts]=1;
//			blade->A[i+3*NElmts][i+2*NElmts]=-dt;
			blade->A[i+3*NElmts][i+3*NElmts-2]	+= 1*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= -4*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts]	+= 1+5*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts+1]	+= -2*coefv_d4vdr4*dt;

			// d4udr4
			blade->A[i+3*NElmts][i+NElmts-2]	+= 1*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts-1]	+= -4*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts]		+= 5*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts+1]	+= -2*coefv_d4udr4*dt;

			// d3vdr3
			blade->A[i+3*NElmts][i+3*NElmts-2]	+= -0.5*coefv_d3vdr3*dt;
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= 1*coefv_d3vdr3*dt;
			blade->A[i+3*NElmts][i+3*NElmts]	+= -0.5*coefv_d3vdr3*dt;

			// d3udr3
			blade->A[i+3*NElmts][i+NElmts-2]	+= -0.5*coefv_d3udr3*dt;
			blade->A[i+3*NElmts][i+NElmts-1]	+= 1*coefv_d3udr3*dt;
			blade->A[i+3*NElmts][i+NElmts]		+= -0.5*coefv_d3udr3*dt;

		}
		else if (i==NElmts-1)
		{
			/*** dudt ***/
			// d4udr4
			blade->A[i][i]			+= 1.0;	
			blade->A[i][i+NElmts-2]		+= 2*coefu_d4udr4;
			blade->A[i][i+NElmts-1]		+= -4*coefu_d4udr4;
			blade->A[i][i+NElmts]		+= 2*coefu_d4udr4;

			// d4vdr4
			blade->A[i][i+3*NElmts-2]	+= 2*coefu_d4vdr4;
			blade->A[i][i+3*NElmts-1]	+= -4*coefu_d4vdr4;
			blade->A[i][i+3*NElmts]		+= 2*coefu_d4vdr4;

			// dudr
			blade->A[i][i+NElmts-1]		+= -1*coefu_dudr;
			blade->A[i][i+NElmts]		+= 1*coefu_dudr;

			// dvdr
			blade->A[i][i+3*NElmts-1]	+= -1*coefu_dvdr;
			blade->A[i][i+3*NElmts]		+= 1*coefu_dvdr;

			/*** u ***/
			// d4udr4
//			blade->A[i+NElmts][i+NElmts]=1;
//			blade->A[i+NElmts][i]=-dt;
			blade->A[i+NElmts][i+NElmts-2] 		+= 2*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts-1]		+= -4*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts]		+= 1+2*coefu_d4udr4*dt;

			// d4vdr4
			blade->A[i+NElmts][i+3*NElmts-2] 	+= 2*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts-1]	+= -4*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= 2*coefu_d4vdr4*dt;

			// dudr
			blade->A[i+NElmts][i+NElmts-1]		+= -1*coefu_dudr*dt;
			blade->A[i+NElmts][i+NElmts]		+= 1*coefu_dudr*dt;

			// dvdr
			blade->A[i+NElmts][i+3*NElmts-1]	+= -1*coefu_dvdr*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= 1*coefu_dvdr*dt;

			/*** dvdt ***/
			// d4vdr4
			blade->A[i+2*NElmts][i+2*NElmts]		+= 1.0;	
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-2]	+= 2*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= -4*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= 2*coefv_d4vdr4;

			// d4udr4
			blade->A[i+2*NElmts][i+NElmts-2]		+= 2*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts-1]		+= -4*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts]			+= 2*coefv_d4udr4;

			// dvdr
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= -1*coefv_dvdr*dt;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= 1*coefv_dvdr*dt;

			// dudr
			blade->A[i+2*NElmts][i+NElmts-1]		+= -1*coefv_dudr*dt;
			blade->A[i+2*NElmts][i+NElmts]			+= 1*coefv_dudr*dt;

			/*** v ***/
			// d4vdr4
//			blade->A[i+3*NElmts][i+3*NElmts]=1;
//			blade->A[i+3*NElmts][i+2*NElmts]=-dt;
			blade->A[i+3*NElmts][i+3*NElmts-2]	+= 2*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= -4*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts]	+= 1+2*coefv_d4vdr4*dt;

			// d4udr4
			blade->A[i+3*NElmts][i+NElmts-2]	+= 2*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts-1]	+= -4*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts]		+= 2*coefv_d4udr4*dt;

			// dvdr
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= -1*coefv_dvdr*dt;
			blade->A[i+3*NElmts][i+3*NElmts]	+= 1*coefv_dvdr*dt;

			// dudr
			blade->A[i+3*NElmts][i+NElmts-1]	+= -1*coefv_dudr*dt;
			blade->A[i+3*NElmts][i+NElmts]		+= 1*coefv_dudr*dt;
		}
		else 
		{
			/*** dudt ***/
			// d4udr4
			blade->A[i][i]			+= 1.0;	
			blade->A[i][i+NElmts-2]		+= coefu_d4udr4;
			blade->A[i][i+NElmts-1]		+= -4*coefu_d4udr4;
			blade->A[i][i+NElmts]		+= 6*coefu_d4udr4;
			blade->A[i][i+NElmts+1]		+= -4*coefu_d4udr4;
			blade->A[i][i+NElmts+2]		+= coefu_d4udr4;

			// d4vdr4
			blade->A[i][i+3*NElmts-2]	+= coefu_d4vdr4;
			blade->A[i][i+3*NElmts-1]	+= -4*coefu_d4vdr4;
			blade->A[i][i+3*NElmts]		+= 6*coefu_d4vdr4;
			blade->A[i][i+3*NElmts+1]	+= -4*coefu_d4vdr4;
			blade->A[i][i+3*NElmts+2]	+= coefu_d4vdr4;

			// d3udr3
			blade->A[i][i+NElmts-2]		+= -0.5*coefu_d3udr3;
			blade->A[i][i+NElmts-1]		+= coefu_d3udr3;
//			blade->A[i][i+NElmts]=coefu_d3udr3;
			blade->A[i][i+NElmts+1]		+= -coefu_d3udr3;
			blade->A[i][i+NElmts+2]		+= 0.5*coefu_d3udr3;

			// d3vdr3
			blade->A[i][i+3*NElmts-2]	+= -0.5*coefu_d3vdr3;
			blade->A[i][i+3*NElmts-1]	+= coefu_d3vdr3;
//			blade->A[i][i+3*NElmts]=coefu_d3udr3;
			blade->A[i][i+3*NElmts+1]	+= -coefu_d3vdr3;
			blade->A[i][i+3*NElmts+2]	+= 0.5*coefu_d3vdr3;

			// d2udr2
			blade->A[i][i+NElmts-1]		+=coefu_d2udr2;
			blade->A[i][i+NElmts]		+=-2*coefu_d2udr2;
			blade->A[i][i+NElmts+1]		+=coefu_d2udr2;

			// d2vdr2
			blade->A[i][i+3*NElmts-1]	+=coefu_d2vdr2;
			blade->A[i][i+3*NElmts]		+=-2*coefu_d2vdr2;
			blade->A[i][i+3*NElmts+1]	+=coefu_d2vdr2;

			/*** u ***/
			// d4udr4
//			blade->A[i+NElmts][i+NElmts]=1;
//			blade->A[i+NElmts][i]=-dt;
			blade->A[i+NElmts][i+NElmts-2]		+= coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts-1]		+= -4*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts]		+= 1+6*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts+1]		+= -4*coefu_d4udr4*dt;
			blade->A[i+NElmts][i+NElmts+2]		+= coefu_d4udr4*dt;

			// d4vdr4
			blade->A[i+NElmts][i+3*NElmts-2]	+= coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts-1]	+= -4*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= 6*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts+1]	+= -4*coefu_d4vdr4*dt;
			blade->A[i+NElmts][i+3*NElmts+2]	+= coefu_d4vdr4*dt;

			// d3udr3
			blade->A[i+NElmts][i+NElmts-2]		+= -0.5*coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts-1]		+= coefu_d3udr3*dt;
//			blade->A[i+NElmts][i+NElmts]		+= coefu_d3udr3;
			blade->A[i+NElmts][i+NElmts+1]		+= -coefu_d3udr3*dt;
			blade->A[i+NElmts][i+NElmts+2]		+= 0.5*coefu_d3udr3*dt;

			// d3vdr3
			blade->A[i+NElmts][i+3*NElmts-2]	+= -0.5*coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts-1]	+= coefu_d3vdr3*dt;
//			blade->A[i+NElmts][i+3*NElmts]		+= coefu_d3udr3;
			blade->A[i+NElmts][i+3*NElmts+1]	+= -coefu_d3vdr3*dt;
			blade->A[i+NElmts][i+3*NElmts+2]	+= 0.5*coefu_d3vdr3*dt;

			// d2udr2
			blade->A[i+NElmts][i+NElmts-1]		+= coefu_d2udr2*dt;
			blade->A[i+NElmts][i+NElmts]		+= -2*coefu_d2udr2*dt;
			blade->A[i+NElmts][i+NElmts+1]		+= coefu_d2udr2*dt;

			// d2vdr2
			blade->A[i+NElmts][i+3*NElmts-1]	+= coefu_d2vdr2*dt;
			blade->A[i+NElmts][i+3*NElmts]		+= -2*coefu_d2vdr2*dt;
			blade->A[i+NElmts][i+3*NElmts+1]	+= coefu_d2vdr2*dt;

			/*** dvdt ***/
			// d4vdr4
			blade->A[i+2*NElmts][i+2*NElmts]		+= 1.0;	
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-2]	+= coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= -4*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= 6*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+1]	+= -4*coefv_d4vdr4;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+2]	+= coefv_d4vdr4;

			// d4udr4
			blade->A[i+2*NElmts][i+NElmts-2]		+= coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts-1]		+= -4*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts]			+= 6*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts+1]		+= -4*coefv_d4udr4;
			blade->A[i+2*NElmts][i+NElmts+2]		+= coefv_d4udr4;

			// d3vdr3
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-2]	+= -0.5*coefv_d3vdr3;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= coefv_d3vdr3;
//			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= coefv_d3vdr3;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+1]	+= -coefv_d3vdr3;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+2]	+= 0.5*coefv_d3vdr3;

			// d3udr3
			blade->A[i+2*NElmts][i+NElmts-2]		+= -0.5*coefv_d3udr3;
			blade->A[i+2*NElmts][i+NElmts-1]		+= coefv_d3udr3;
//			blade->A[i+2*NElmts][i+NElmts]			+= coefv_d3udr3;
			blade->A[i+2*NElmts][i+NElmts+1]		+= -coefv_d3udr3;
			blade->A[i+2*NElmts][i+NElmts+2]		+= 0.5*coefv_d3udr3;

			// d2vdr2
			blade->A[i+2*NElmts][i+2*NElmts+NElmts-1]	+= coefv_d2vdr2;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts]		+= -2*coefv_d2vdr2;
			blade->A[i+2*NElmts][i+2*NElmts+NElmts+1]	+= coefv_d2vdr2;

			// d2udr2
			blade->A[i+2*NElmts][i+NElmts-1]	+= coefv_d2udr2;
			blade->A[i+2*NElmts][i+NElmts]		+= -2*coefv_d2udr2;
			blade->A[i+2*NElmts][i+NElmts+1]	+= coefv_d2udr2;

			/*** v ***/
			// d4vdr4
//			blade->A[i+3*NElmts][i+3*NElmts]=1;
//			blade->A[i+3*NElmts][i+2*NElmts]=-dt;
			blade->A[i+3*NElmts][i+3*NElmts-2]	+= coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= -4*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts]	+= 1+6*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts+1]	+= -4*coefv_d4vdr4*dt;
			blade->A[i+3*NElmts][i+3*NElmts+2]	+= coefv_d4vdr4*dt;

			// d4udr4
			blade->A[i+3*NElmts][i+NElmts-2]	+= coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts-1]	+= -4*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts]		+= 6*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts+1]	+= -4*coefv_d4udr4*dt;
			blade->A[i+3*NElmts][i+NElmts+2]	+= coefv_d4udr4*dt;

			// d3vdr3
			blade->A[i+3*NElmts][i+3*NElmts-2]	+= -0.5*coefv_d3vdr3*dt;
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= coefv_d3vdr3*dt;
//			blade->A[i+3*NElmts][i+3*NElmts]	+= coefv_d3vdr3*dt;
			blade->A[i+3*NElmts][i+3*NElmts+1]	+= -coefv_d3vdr3*dt;
			blade->A[i+3*NElmts][i+3*NElmts+2]	+= 0.5*coefv_d3vdr3*dt;

			// d3udr3
			blade->A[i+3*NElmts][i+NElmts-2]	+= -0.5*coefv_d3udr3*dt;
			blade->A[i+3*NElmts][i+NElmts-1]	+= coefv_d3udr3*dt;
//			blade->A[i+3*NElmts][i+NElmts]		+= coefv_d3udr3*dt;
			blade->A[i+3*NElmts][i+NElmts+1]	+= -coefv_d3udr3*dt;
			blade->A[i+3*NElmts][i+NElmts+2]	+= 0.5*coefv_d3udr3*dt;

			// d2vdr2
			blade->A[i+3*NElmts][i+3*NElmts-1]	+= coefv_d2vdr2*dt;
			blade->A[i+3*NElmts][i+3*NElmts]	+= -2*coefv_d2vdr2*dt;
			blade->A[i+3*NElmts][i+3*NElmts+1]	+= coefv_d2vdr2*dt;

			// d2udr2
			blade->A[i+3*NElmts][i+NElmts-1]	+= coefv_d2udr2*dt;
			blade->A[i+3*NElmts][i+NElmts]		+= -2*coefv_d2udr2*dt;
			blade->A[i+3*NElmts][i+NElmts+1]	+= coefv_d2udr2*dt;
		}
	}

	// BCs
	i = 0;
	blade->A[i][i]=1.0; // dudt
	blade->A[i+NElmts][i+NElmts]=1.0; // u
	blade->A[i+2*NElmts][i+2*NElmts]=1.0; // dvdt
	blade->A[i+3*NElmts][i+3*NElmts]=1.0; // v


	for (i=0;i<4*NElmts;i++)
	{	
//		cout << endl << i << endl;
//		cout << "\t";
		for (j=0;j<4*NElmts;j++)
		{
//			cout << blade->A[i][j] << " "; 
			blade->A_LU[i][j] = blade->A[i][j];

		}
	}
}

/**
 * contruct the rhs for the LU solver
 */
PetscErrorCode constructrhsb_bladestructure(Beam *blade)
{
	using namespace std;
 	int i, j, k;

	int NElmts=blade->NumberofElements;

	for (i=0;i<NElmts;i++)
	{

		if (i==0)
		{
			// dudt
			blade->b[i]=0.0;

			// u
			blade->b[i+NElmts]=0.0;

			// dvdt
			blade->b[i+2*NElmts]=0.0;

			// v
			blade->b[i+3*NElmts]=0.0;

		} else
		{
			// dudt
			blade->b[i]=blade->rhsdudt[i];

			// u
			blade->b[i+NElmts]=blade->rhsu[i];

			// dvdt
			blade->b[i+2*NElmts]=blade->rhsdvdt[i];

			// v
			blade->b[i+3*NElmts]=blade->rhsv[i];
		}
	}


//	for (i=0;i<4*NElmts;i++)
//		blade->bb[i]=blade->b[i];
}

/**
* Solve the structure equations only on rank 0. It should be ok for the actuator line simulations. Need to test. 
*/  
PetscErrorCode solvestructurefunctions_bladestructure(IBMNodes *ibm)
{
 	int i, j, k;
	int ibi, iblade;
	double norm;

//  	PetscPrintf(PETSC_COMM_WORLD,"start solve structure equations\n");
	int my_rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);


	// save previous displacements and velocity 
	for (ibi=0; ibi<NumberOfTurbines; ibi++) 
	for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) 
	{

//		PetscPrintf(PETSC_COMM_WORLD,"%d th turbine %d th blade\n", ibi, iblade);
//		cout << ibi << "th turbine " << iblade << "th blade" << endl;

		int NElmts = ibm[ibi].bladestructure[iblade].NumberofElements;

		for (i=0;i<NElmts;i++) 
		{
	    	        ibm[ibi].bladestructure[iblade].u_om2[i] = ibm[ibi].bladestructure[iblade].u_om1[i];
    		        ibm[ibi].bladestructure[iblade].u_om1[i] = ibm[ibi].bladestructure[iblade].u_o[i];
    	        	ibm[ibi].bladestructure[iblade].u_o[i] = ibm[ibi].bladestructure[iblade].u[i];

	    	        ibm[ibi].bladestructure[iblade].v_om2[i] = ibm[ibi].bladestructure[iblade].v_om1[i];
    		        ibm[ibi].bladestructure[iblade].v_om1[i] = ibm[ibi].bladestructure[iblade].v_o[i];
    	        	ibm[ibi].bladestructure[iblade].v_o[i] = ibm[ibi].bladestructure[iblade].v[i];

	    	        ibm[ibi].bladestructure[iblade].rhsu_om2[i] = ibm[ibi].bladestructure[iblade].rhsu_om1[i];
    		        ibm[ibi].bladestructure[iblade].rhsu_om1[i] = ibm[ibi].bladestructure[iblade].rhsu_o[i];
    	        	ibm[ibi].bladestructure[iblade].rhsu_o[i] = ibm[ibi].bladestructure[iblade].rhsu[i];

	    	        ibm[ibi].bladestructure[iblade].rhsv_om2[i] = ibm[ibi].bladestructure[iblade].rhsv_om1[i];
    		        ibm[ibi].bladestructure[iblade].rhsv_om1[i] = ibm[ibi].bladestructure[iblade].rhsv_o[i];
    	        	ibm[ibi].bladestructure[iblade].rhsv_o[i] = ibm[ibi].bladestructure[iblade].rhsv[i];

			if (torsion_turbinestructure) 
			{ 

			        ibm[ibi].bladestructure[iblade].theta_om2[i] = ibm[ibi].bladestructure[iblade].theta_om1[i];
    		        	ibm[ibi].bladestructure[iblade].theta_om1[i] = ibm[ibi].bladestructure[iblade].theta_o[i];
	    		        ibm[ibi].bladestructure[iblade].theta_o[i] = ibm[ibi].bladestructure[iblade].theta[i];

			        ibm[ibi].bladestructure[iblade].rhstheta_om2[i] = ibm[ibi].bladestructure[iblade].rhstheta_om1[i];
    		        	ibm[ibi].bladestructure[iblade].rhstheta_om1[i] = ibm[ibi].bladestructure[iblade].rhstheta_o[i];
	    		        ibm[ibi].bladestructure[iblade].rhstheta_o[i] = ibm[ibi].bladestructure[iblade].rhstheta[i];
			}

    	        	ibm[ibi].bladestructure[iblade].dudt_om2[i] = ibm[ibi].bladestructure[iblade].dudt_om1[i];
	    	        ibm[ibi].bladestructure[iblade].dudt_om1[i] = ibm[ibi].bladestructure[iblade].dudt_o[i];
    		        ibm[ibi].bladestructure[iblade].dudt_o[i] = ibm[ibi].bladestructure[iblade].dudt[i];

	    	        ibm[ibi].bladestructure[iblade].dvdt_om2[i] = ibm[ibi].bladestructure[iblade].dvdt_om1[i];
    		        ibm[ibi].bladestructure[iblade].dvdt_om1[i] = ibm[ibi].bladestructure[iblade].dvdt_o[i];
    	        	ibm[ibi].bladestructure[iblade].dvdt_o[i] = ibm[ibi].bladestructure[iblade].dvdt[i];

    	        	ibm[ibi].bladestructure[iblade].rhsdudt_om2[i] = ibm[ibi].bladestructure[iblade].rhsdudt_om1[i];
	    	        ibm[ibi].bladestructure[iblade].rhsdudt_om1[i] = ibm[ibi].bladestructure[iblade].rhsdudt_o[i];
    		        ibm[ibi].bladestructure[iblade].rhsdudt_o[i] = ibm[ibi].bladestructure[iblade].rhsdudt[i];

	    	        ibm[ibi].bladestructure[iblade].rhsdvdt_om2[i] = ibm[ibi].bladestructure[iblade].rhsdvdt_om1[i];
    		        ibm[ibi].bladestructure[iblade].rhsdvdt_om1[i] = ibm[ibi].bladestructure[iblade].rhsdvdt_o[i];
    	        	ibm[ibi].bladestructure[iblade].rhsdvdt_o[i] = ibm[ibi].bladestructure[iblade].rhsdvdt[i];


			if (torsion_turbinestructure) 
			{
			        ibm[ibi].bladestructure[iblade].dthetadt_om2[i] = ibm[ibi].bladestructure[iblade].dthetadt_om1[i];
		        	ibm[ibi].bladestructure[iblade].dthetadt_om1[i] = ibm[ibi].bladestructure[iblade].dthetadt_o[i];
	    		        ibm[ibi].bladestructure[iblade].dthetadt_o[i] = ibm[ibi].bladestructure[iblade].dthetadt[i];

			        ibm[ibi].bladestructure[iblade].rhsdthetadt_om2[i] = ibm[ibi].bladestructure[iblade].rhsdthetadt_om1[i];
		        	ibm[ibi].bladestructure[iblade].rhsdthetadt_om1[i] = ibm[ibi].bladestructure[iblade].rhsdthetadt_o[i];
	    		        ibm[ibi].bladestructure[iblade].rhsdthetadt_o[i] = ibm[ibi].bladestructure[iblade].rhsdthetadt[i];
			}
		}

		if (!my_rank) 
		{
			int N = 4*NElmts;
			double d;

			vector<double> u_pre, v_pre, dudt_pre, dvdt_pre;

			u_pre.resize(NElmts);
			v_pre.resize(NElmts);
			dudt_pre.resize(NElmts);
			dvdt_pre.resize(NElmts);

			double res=0, res_max=1.e-6;
			int iter=0, iter_max=2;
//			do 
//			{
				for (i=0;i<NElmts;i++)
				{
					dudt_pre.at(i)	=ibm[ibi].bladestructure[iblade].dudt[i];
					dvdt_pre.at(i)	=ibm[ibi].bladestructure[iblade].dvdt[i];
					u_pre.at(i)	=ibm[ibi].bladestructure[iblade].u[i];
					v_pre.at(i)	=ibm[ibi].bladestructure[iblade].v[i];
				}

//		    PetscPrintf(PETSC_COMM_WORLD,"InitmatrixA \n");
				initmatrixA_bladestructure(&(ibm[ibi].bladestructure[iblade]));

//		    PetscPrintf(PETSC_COMM_WORLD,"Compute RHS DS \n");
				computerhsDS_bladestructure(&(ibm[ibi].bladestructure[iblade]));

//		    PetscPrintf(PETSC_COMM_WORLD,"Construct matrix A \n");
				constructmatrixA_bladestructure(&(ibm[ibi].bladestructure[iblade]));

//                                {
//                                char filen[80];
//                                sprintf(filen, "matrixA1_%3.3d.dat", iblade);
//                                Write_Amatrix(filen,NElmts, iblade, &(ibm[ibi].bladestructure[iblade]), ibm);
//                                }


//		    PetscPrintf(PETSC_COMM_WORLD,"LU Decomposition \n");
				ludcmp((ibm[ibi].bladestructure[iblade].A_LU), N, (ibm[ibi].bladestructure[iblade].indx), &d);

//                                {
//                                char filen[80];
//                                sprintf(filen, "matrixA2_%3.3d.dat", iblade);
//                                Write_Amatrix(filen,NElmts, iblade, &(ibm[ibi].bladestructure[iblade]), ibm);
//                                }

//		    PetscPrintf(PETSC_COMM_WORLD,"Compute RHS \n");
				constructrhsb_bladestructure(&(ibm[ibi].bladestructure[iblade]));

//                                {
//                                char filen[80];
//                                sprintf(filen, "matrixB_%3.3d.dat", iblade);
//                                Write_Bmatrix(filen,NElmts, iblade, &(ibm[ibi].bladestructure[iblade]), ibm);
//                                }
		
				// NOTE: haven't added torsion_turbinestructure

//		    PetscPrintf(PETSC_COMM_WORLD,"LUBack subs \n");
				lubksb((ibm[ibi].bladestructure[iblade].A_LU), N, (ibm[ibi].bladestructure[iblade].indx), (ibm[ibi].bladestructure[iblade].b));

				for (i=0;i<NElmts;i++)
				{
					ibm[ibi].bladestructure[iblade].dudt[i] = ibm[ibi].bladestructure[iblade].b[i];
					ibm[ibi].bladestructure[iblade].u[i] = ibm[ibi].bladestructure[iblade].b[i+NElmts];
					ibm[ibi].bladestructure[iblade].dvdt[i] = ibm[ibi].bladestructure[iblade].b[i+2*NElmts];
					ibm[ibi].bladestructure[iblade].v[i] = ibm[ibi].bladestructure[iblade].b[i+3*NElmts];

//         		                PetscPrintf(PETSC_COMM_WORLD,"var sol   %le %le %le %le \n",ibm[ibi].bladestructure[iblade].b[i],ibm[ibi].bladestructure[iblade].b[i+NElmts],ibm[ibi].bladestructure[iblade].b[i+2*NElmts], ibm[ibi].bladestructure[iblade].b[i+3*NElmts]);
				}


//    		                PetscPrintf(PETSC_COMM_WORLD,"var indx  %d\n",ibm[ibi].bladestructure[iblade].indx);
                                char filen[80];
                                sprintf(filen, "Displacement_%3.3d.dat", iblade);
                                Write_Displacement(filen,NElmts, iblade, &(ibm[ibi].bladestructure[iblade]), ibm);

				for (i=0;i<NElmts;i++)
				{
//					res = fmax(fabs(ibm[ibi].bladestructure[iblade].dudt[i]-dudt_pre.at(i)), res);
//					res = fmax(fabs(ibm[ibi].bladestructure[iblade].dvdt[i]-dvdt_pre.at(i)), res);
//					res = fmax(fabs(ibm[ibi].bladestructure[iblade].u[i]-u_pre.at(i)), res);
//					res = fmax(fabs(ibm[ibi].bladestructure[iblade].v[i]-v_pre.at(i)), res);

					res += pow(fabs(ibm[ibi].bladestructure[iblade].dudt[i]-dudt_pre.at(i)), 2);
					res += pow(fabs(ibm[ibi].bladestructure[iblade].dvdt[i]-dvdt_pre.at(i)), 2);
					res += pow(fabs(ibm[ibi].bladestructure[iblade].u[i]-u_pre.at(i)), 2);
					res += pow(fabs(ibm[ibi].bladestructure[iblade].v[i]-v_pre.at(i)), 2);
				}

				res /= (4*NElmts);

				res = sqrt(res);
				iter++;

//				cout << "\titer " << iter << " residual: " << res << endl;
//			}  while (res>res_max && iter<iter_max);

//		  PetscPrintf(PETSC_COMM_WORLD,"Residuals \n");
//			computeRES_bladestructure(&(ibm[ibi].bladestructure[iblade]));

			/*
			double res[N];
			for (i=0;i<N;i++)
				res[i]=0.0;

			double res_max=0;
			for (i=0;i<N;i++)
			{
				for (j=0;j<N;j++)
				{
					res[i]+=ibm[ibi].bladestructure[iblade].A[i][j]*ibm[ibi].bladestructure[iblade].b[j];	
				}
				res[i]-=ibm[ibi].bladestructure[iblade].bb[i];
				if (fabs(res[i])>res_max) res_max = fabs(res[i]);
//				cout << i << " residual: " << res[i] << endl;
			}

			cout << "The maximum residual (LU): " << res_max << endl;
			*/


    			MPI_Bcast(ibm[ibi].bladestructure[iblade].u, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		        MPI_Bcast(ibm[ibi].bladestructure[iblade].v, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    			if (torsion_turbinestructure) MPI_Bcast(ibm[ibi].bladestructure[iblade].theta, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);

    		        MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		        MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    		        if (torsion_turbinestructure) MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);

		} else 
		{

		        MPI_Bcast(ibm[ibi].bladestructure[iblade].u, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].v, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
		        if (torsion_turbinestructure) MPI_Bcast(ibm[ibi].bladestructure[iblade].theta, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);

		        MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
		        MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
		        if (torsion_turbinestructure) MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
		}

//csantoni
//        if(iblade==1) PetscFinalize();
	}

//  	PetscPrintf(PETSC_COMM_WORLD,"finish solve structure equations: %i \n", ti);
//        if(ti >tistart){
//        PetscFinalize();
//        }


	return(0);
}	

/** 
 * Transfer and reconstruct the computed aerodynamic forces for structure solver 
 */
PetscErrorCode getforcesfromflowsolver_bladestructure(IBMNodes *ibm, FSInfo *fsi)
{
 	int i, j, k;
	int ibi, iblade;

	//double *r0, *g0, *r1, *g1;
	//double *fx, *fy, *fz;

	// compute the angle of the blade with vertical z-axis, pitch
	for (ibi=0; ibi<NumberOfTurbines; ibi++) 
	for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) 
	{
		ibm[ibi].bladestructure[iblade].beta_om2=ibm[ibi].bladestructure[iblade].beta_om1;
		ibm[ibi].bladestructure[iblade].beta_om1=ibm[ibi].bladestructure[iblade].beta_o;
		ibm[ibi].bladestructure[iblade].beta_o=ibm[ibi].bladestructure[iblade].beta;
		ibm[ibi].bladestructure[iblade].beta = -acos(-1.0)*ibm[ibi].pitch[iblade]/180.0; 


		ibm[ibi].bladestructure[iblade].phi_om2 = ibm[ibi].bladestructure[iblade].phi_om1;
		ibm[ibi].bladestructure[iblade].phi_om1 = ibm[ibi].bladestructure[iblade].phi_o;
		ibm[ibi].bladestructure[iblade].phi_o   = ibm[ibi].bladestructure[iblade].phi;

		ibm[ibi].bladestructure[iblade].phi     = abs(fsi[ibi].ang_axis)  + ((double) iblade)*2.*M_PI/3. + M_PI/2 ; 

//    PetscPrintf(PETSC_COMM_WORLD, "#######  Phi: %e %e\n", ibm[ibi].bladestructure[iblade].phi,fsi[ibi].ang_axis);

// This parameter is obtained dimensional at Calc_turbineangvel function
//    ibm[ibi].bladestructure[iblade].d2phidt2  = fsi[ibi].angaccel_axis;
    ibm[ibi].bladestructure[iblade].d2phidt2  = abs(fsi[ibi].angaccel_axis);
//     // Obtained from Calc_turbineangvel function
//    ibm[ibi].bladestructure[iblade].dphidt    = fsi[ibi].angvel_axis*refvel_wt/reflength_wt;
    ibm[ibi].bladestructure[iblade].dphidt    = abs(fsi[ibi].angvel_axis)*refvel_wt/reflength_wt;
//    ibm[ibi].bladestructure[iblade].phi       = fsi[ibi].ang_axis;
	}
	

	vector<double> r0, g0, r1, g1, fx, fy, fz;	
//	FILE *fd;
//    	char filen[80]; 

//	PetscPrintf(PETSC_COMM_WORLD, "start getforcesfromflowsolver_bladestructure  \n");
	for (ibi=0; ibi<NumberOfTurbines; ibi++) {

		

		// the structure intertial system  
		double nY_x, nY_y, nY_z; // axis
		double nZ_x, nZ_y, nZ_z; // vertical
		double nX_x, nX_y, nX_z; // X
		double rr;
				
		rr=sqrt(pow(fsi[ibi].nx_tb,2)+pow(fsi[ibi].ny_tb,2)+pow(fsi[ibi].nz_tb,2))+1.e-19;

		nY_x = fsi[ibi].nx_tb/rr; 
		nY_y = fsi[ibi].ny_tb/rr; 
		nY_z = fsi[ibi].nz_tb/rr; 

		rr=sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2))+1.e-19;

		nZ_x = -gravity_x/rr; 
		nZ_y = -gravity_y/rr; 
		nZ_z = -gravity_z/rr; 

		nX_x=nY_z*nZ_y-nY_y*nZ_z;
		nX_y=nY_x*nZ_z-nY_y*nZ_z;
		nX_z=nY_x*nZ_y-nY_y*nZ_x;

//		PetscPrintf(PETSC_COMM_WORLD, "nX_x, nX_y nX_z: %le %le %le  \n", nX_x,nX_y,nX_z);
//		PetscPrintf(PETSC_COMM_WORLD, "nY_x, nY_y nY_z: %le %le %le  \n", nY_x,nY_y,nY_z);
//		PetscPrintf(PETSC_COMM_WORLD, "nZ_x, nZ_y nZ_z: %le %le %le  \n", nZ_x,nZ_y,nZ_z);

		int NumberofElements_AL = ibm[ibi].n_elmt/ibm[ibi].num_blade;

//		PetscPrintf(PETSC_COMM_WORLD, "Allocate r0, g0  \n");

		//r0=(double *) malloc(NumberofElements_AL*sizeof(double));
		//g0=(double *) malloc(NumberofElements_AL*sizeof(double));

		r0.resize(NumberofElements_AL);
		g0.resize(NumberofElements_AL);
//		PetscPrintf(PETSC_COMM_WORLD, "Done Allocate r0, g0  \n");

		for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {

//			PetscPrintf(PETSC_COMM_WORLD, "iblade=%d Allocate r1, g1, fx, fy, fz  \n", iblade);
			r1.resize(ibm[ibi].bladestructure[iblade].NumberofElements);
			g1.resize(ibm[ibi].bladestructure[iblade].NumberofElements);
			fx.resize(ibm[ibi].bladestructure[iblade].NumberofElements);
			fy.resize(ibm[ibi].bladestructure[iblade].NumberofElements);
			fz.resize(ibm[ibi].bladestructure[iblade].NumberofElements);

//			PetscPrintf(PETSC_COMM_WORLD, "iblade=%d Done Allocate r1, g1, fx, fy, fz  \n", iblade);

			//double Phi = ibm[ibi].bladestructure[iblade].phi;

			// Phi is not from the torque cotroller yet
			int NumberofNodes_AL = ibm[ibi].n_v/ibm[ibi].num_blade;
		
			int ii = iblade * NumberofNodes_AL + 5; 

			double rx = ibm[ibi].x_bp[ii]-fsi[ibi].x_c;
			double ry = ibm[ibi].y_bp[ii]-fsi[ibi].y_c;
			double rz = ibm[ibi].z_bp[ii]-fsi[ibi].z_c;

			double rr = sqrt(rx*rx+ry*ry+rz*rz);	
			rx /= rr; ry /= rr; rz /= rr;

			double gx=nZ_x, gy=nZ_y, gz=nZ_z; // assuming y is the vertical direction in CFD

			//double Phi = acos(rx*gx+ry*gy+rz*gz);
      double Phi=ibm[ibi].bladestructure[iblade].phi;


			double Beta = -acos(-1.0)*ibm[ibi].pitch[iblade]/180.0;

			for (i=0;i<NumberofElements_AL;i++){
				int ii = iblade * NumberofElements_AL + i; 

				double rx = ibm[ibi].cent_x[ii]-fsi[ibi].x_c;
				double ry = ibm[ibi].cent_y[ii]-fsi[ibi].y_c;
				double rz = ibm[ibi].cent_z[ii]-fsi[ibi].z_c;
	        	
				r0.at(i) = sqrt(rx*rx+ry*ry+rz*rz);
			}
        	


			for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) r1.at(i)=ibm[ibi].bladestructure[iblade].r[i]/reflength_wt;

			// scale the force from CFD to real
			double rho_cfd=1.0; // no levelset 	
			double ratio_rho=rho_air/rho_cfd;
			double ratio_V=refvel_wt/refvel_cfd;

			//PetscPrintf(PETSC_COMM_WORLD, "iblade=%d interpolation F_lagr_x 1 \n", iblade);
			// interpolate force to the structure model grid
			for (i=0;i<NumberofElements_AL;i++) 
			{
				int ii = iblade * NumberofElements_AL + i; 
				g0.at(i)=-ibm[ibi].F_lagr_x[ii]*(ratio_rho*pow(ratio_V,2)*pow(reflength_wt,1));
			}

			//PetscPrintf(PETSC_COMM_WORLD, "iblade=%d interpolation F_lagr_x 2 \n", iblade);
			_interpolation_leastsquares(r0, g0, r1, fx, NumberofElements_AL, ibm[ibi].bladestructure[iblade].NumberofElements);

			//PetscPrintf(PETSC_COMM_WORLD, "iblade=%d interpolation F_lagr_y 1 \n", iblade);
			for (i=0;i<NumberofElements_AL;i++) 
			{
				int ii = iblade * NumberofElements_AL + i; 
				g0.at(i)=-ibm[ibi].F_lagr_y[ii]*(ratio_rho*pow(ratio_V,2)*pow(reflength_wt,1));
			}

			//PetscPrintf(PETSC_COMM_WORLD, "iblade=%d interpolation F_lagr_y 2 \n", iblade);
			_interpolation_leastsquares(r0, g0, r1, fy, NumberofElements_AL, ibm[ibi].bladestructure[iblade].NumberofElements);

	//		PetscPrintf(PETSC_COMM_WORLD, "iblade=%d interpolation F_lagr_z 1 \n", iblade);
			for (i=0;i<NumberofElements_AL;i++) 
			{
				int ii = iblade * NumberofElements_AL + i; 
				g0.at(i)=-ibm[ibi].F_lagr_z[ii]*(ratio_rho*pow(ratio_V,2)*pow(reflength_wt,1));
			}

	//		PetscPrintf(PETSC_COMM_WORLD, "iblade=%d interpolation F_lagr_z 2 \n", iblade);
			_interpolation_leastsquares(r0, g0, r1, fz, NumberofElements_AL, ibm[ibi].bladestructure[iblade].NumberofElements);


//			PetscPrintf(PETSC_COMM_WORLD, "iblade=%d projection inertial-->(XYZ) \n", iblade);
			// projection: from inertial frame of the flow to the intertial frame of the structure (X,Y,Z) 
//		int rank=0;
//          FILE *f;
//	  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//		if (!rank) {
//          char filen[80];
//              sprintf(filen, "./Tstruct_%06d_%02d.dat", ti,iblade);
//              f = fopen(filen, "w");
//    }

			for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				double fx_tmp=fx.at(i);	
				double fy_tmp=fy.at(i);	
				double fz_tmp=fz.at(i);	

				fx.at(i)=fx_tmp*nX_x+fy_tmp*nX_y+fz_tmp*nX_z;
				fy.at(i)=fx_tmp*nY_x+fy_tmp*nY_y+fz_tmp*nY_z;
				fz.at(i)=fx_tmp*nZ_x+fy_tmp*nZ_y+fz_tmp*nZ_z;

     // if (!rank) PetscFPrintf(PETSC_COMM_WORLD, f, "%i %le %le %le %le %le %le %le %le %le \n", i, fx_tmp,fy_tmp,fz_tmp, fx.at(i),fy.at(i),fz.at(i),nX_x,nX_y,nX_z);
//			//PetscPrintf(PETSC_COMM_WORLD, "Fx Fy Fz: %i %le %le %le %le %le %le \n", i, fx_tmp,fy_tmp,fz_tmp, fx.at(i),fy.at(i),fz.at(i));

			}

		
		
//			PetscPrintf(PETSC_COMM_WORLD, "iblade=%d projection (XYZ)--->(xhatyhatzhat) \n", iblade);
			// projection: from the structure inertial frame (X,Y,Z) to frame rotating with the hub (xhat,yhat,zhat)
			for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				double fX=fx.at(i);	
				double fY=fy.at(i);	
				double fZ=fz.at(i);	
				double fxhat, fyhat, fzhat;
        double Phi = ibm[ibi].bladestructure[iblade].phi;
//        double pphi=ibm[ibi].bladestructure[0].phi;

				xyzhat_fromXYZ(Phi, fX, fY, fZ, &fxhat, &fyhat, &fzhat);

				fx.at(i)=fxhat;
				fy.at(i)=fyhat;
				fz.at(i)=fzhat;


			}

	
//			PetscPrintf(PETSC_COMM_WORLD, "iblade=%d projection (xhatyhatzhat)---> rotating with pi \n", iblade);
			// projection: from the structure inertial frame (X,Y,Z) to frame rotating with the hub (xhat,yhat,zhat)
			// projection: from frame rotaing with the hub (xhat,yhat,zhat) to frame rotating with pi because of pitch action
			for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				double fxhat=fx.at(i);	
				double fyhat=fy.at(i);	
				double fzhat=fz.at(i);	
				double fx_tmp, fy_tmp, fz_tmp;

				xyz_fromxyzhat(Beta, fxhat, fyhat, fzhat, &fx_tmp, &fy_tmp, &fz_tmp);
  //    if (!rank) PetscFPrintf(PETSC_COMM_WORLD, f, "%i %le %le %le %le %le %le %le\n", i, fxhat,fyhat,fzhat, fx.at(i),fy.at(i),fz.at(i),Phi);

				fx.at(i)=fx_tmp;
				fy.at(i)=fy_tmp;
				fz.at(i)=fz_tmp;


			}
   //   if (!rank) fclose(f);


	//		PetscPrintf(PETSC_COMM_WORLD, "iblade=%d save old \n", iblade);
			// save 
			for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				ibm[ibi].bladestructure[iblade].fu_o[i]=ibm[ibi].bladestructure[iblade].fu[i];
				ibm[ibi].bladestructure[iblade].fv_o[i]=ibm[ibi].bladestructure[iblade].fv[i];
				ibm[ibi].bladestructure[iblade].M_o[i]=ibm[ibi].bladestructure[iblade].M[i];
			}



	//		PetscPrintf(PETSC_COMM_WORLD, "iblade=%d set new \n", iblade);
			for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) {
				ibm[ibi].bladestructure[iblade].fu[i]=fx.at(i);
				ibm[ibi].bladestructure[iblade].fv[i]=fy.at(i);
				ibm[ibi].bladestructure[iblade].M[i]=0;
			}

		}	
	}
//		PetscPrintf(PETSC_COMM_WORLD, "end getforcesfromflowsolver_bladestructure  \n");

	return(0);
}



// transfer the blade displacement from structure solver to flow solver  
PetscErrorCode transferbladedisplacement2flowblade(IBMNodes *ibm, FSInfo *fsi)
{
 	int i, j, k;
	int ibi, iblade;


	using namespace std;

	vector<double> r0, r0_cent, r1, g1, g2;
	vector<double> fx, fy, fz, ftheta;
	vector<double> ux, uy, uz, utheta;

//	MPI_Barrier(PETSC_COMM_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD, "Turbine Loop  \n");

	for (ibi=0; ibi<NumberOfTurbines; ibi++) {

//		cout << "transferbladedisplacement2flowblade: dir \n";
		// the structure intertial system  
		double nY_x, nY_y, nY_z; // axis
		double nZ_x, nZ_y, nZ_z; // vertical
		double nX_x, nX_y, nX_z; // X
		double rr;
				
		rr=sqrt(pow(fsi[ibi].nx_tb,2)+pow(fsi[ibi].ny_tb,2)+pow(fsi[ibi].nz_tb,2))+1.e-19;

//		cout << "transferbladedisplacement2flowblade: dir 1 \n";
		nY_x = fsi[ibi].nx_tb/rr; 
		nY_y = fsi[ibi].ny_tb/rr; 
		nY_z = fsi[ibi].nz_tb/rr; 

		rr=sqrt(pow(gravity_x,2)+pow(gravity_y,2)+pow(gravity_z,2))+1.e-19;

//		cout << "transferbladedisplacement2flowblade: dir 2\n";
		nZ_x = -gravity_x/rr; 
		nZ_y = -gravity_y/rr; 
		nZ_z = -gravity_z/rr; 

//		cout << "transferbladedisplacement2flowblade: dir 3 \n";
		nX_x=nY_y*nZ_z-nY_z*nZ_y;
		nX_y=nY_z*nZ_x-nY_x*nZ_z;
		nX_z=nY_x*nZ_y-nY_y*nZ_x;

//		cout << "transferbladedisplacement2flowblade: dir 4 \n";

//		cout << "transferbladedisplacement2flowblade: allocate mem \n";
		int NumberofNodes_AL = ibm[ibi].n_v/ibm[ibi].num_blade;
		int NumberofElements_AL = ibm[ibi].n_elmt/ibm[ibi].num_blade;

		r0.resize(NumberofNodes_AL);
		r0_cent.resize(NumberofElements_AL);

		fx.resize(NumberofNodes_AL);
		fy.resize(NumberofNodes_AL);
		fz.resize(NumberofNodes_AL);
		if (torsion_turbinestructure) ftheta.resize(NumberofElements_AL);

		ux.resize(NumberofNodes_AL);
		uy.resize(NumberofNodes_AL);
		uz.resize(NumberofNodes_AL);
		if (torsion_turbinestructure) utheta.resize(NumberofElements_AL);


		for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) 
		{

//			cout << "transferbladedisplacement2flowblade: mem, r0, r0_cent \n";

			r1.resize(ibm[ibi].bladestructure[iblade].NumberofElements);
			g1.resize(ibm[ibi].bladestructure[iblade].NumberofElements);
			g2.resize(ibm[ibi].bladestructure[iblade].NumberofElements);

			double Beta = -acos(-1.0)*ibm[ibi].pitch[iblade]/180.0;

			for (i=0;i<NumberofNodes_AL;i++){
				int ii = iblade * NumberofNodes_AL + i; 

				double rx = ibm[ibi].x_bp[ii]-fsi[ibi].x_c;
				double ry = ibm[ibi].y_bp[ii]-fsi[ibi].y_c;
				double rz = ibm[ibi].z_bp[ii]-fsi[ibi].z_c;
	        	
				r0.at(i) = sqrt(rx*rx+ry*ry+rz*rz);
			}
        	
			for (i=0;i<NumberofElements_AL;i++){
				int ii = iblade * NumberofElements_AL + i; 

				double rx = ibm[ibi].cent_x[ii]-fsi[ibi].x_c;
				double ry = ibm[ibi].cent_y[ii]-fsi[ibi].y_c;
				double rz = ibm[ibi].cent_z[ii]-fsi[ibi].z_c;
	        	
				r0_cent.at(i) = sqrt(rx*rx+ry*ry+rz*rz);
			}


			for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) r1.at(i)=ibm[ibi].bladestructure[iblade].r[i]/reflength_wt;

//			cout << "transferbladedisplacement2flowblade: u \n";
			// interpolate displacement to the actuator line grid of flow solver
			// u, dudt
			for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) 
			{
				g1.at(i)=ibm[ibi].bladestructure[iblade].u[i];
				g2.at(i)=ibm[ibi].bladestructure[iblade].dudt[i];
			}

			_interpolation_leastsquares(r1, g1, r0, fx, ibm[ibi].bladestructure[iblade].NumberofElements, NumberofNodes_AL);
			_interpolation_leastsquares(r1, g2, r0, ux, ibm[ibi].bladestructure[iblade].NumberofElements, NumberofNodes_AL);

//			cout << "transferbladedisplacement2flowblade: v \n";
			// v
			for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) 
			{
				g1.at(i)=ibm[ibi].bladestructure[iblade].v[i];
				g2.at(i)=ibm[ibi].bladestructure[iblade].dvdt[i];
			}

			_interpolation_leastsquares(r1, g1, r0, fy, ibm[ibi].bladestructure[iblade].NumberofElements, NumberofNodes_AL);
			_interpolation_leastsquares(r1, g2, r0, uy, ibm[ibi].bladestructure[iblade].NumberofElements, NumberofNodes_AL);

//			cout << "transferbladedisplacement2flowblade: w \n";
			// w
			for (i=0;i<NumberofNodes_AL;i++)
			{
				fz[i]=0;
				uz[i]=0;
			}

//			cout << "transferbladedisplacement2flowblade: theta \n";
			// theta
			if (torsion_turbinestructure) 
			{
				for (i=0;i<ibm[ibi].bladestructure[iblade].NumberofElements;i++) 
				{
					g1.at(i)=ibm[ibi].bladestructure[iblade].theta[i];
					g2.at(i)=ibm[ibi].bladestructure[iblade].dthetadt[i];
				}

				_interpolation_leastsquares(r1, g1, r0_cent, ftheta, ibm[ibi].bladestructure[iblade].NumberofElements, NumberofElements_AL);
				_interpolation_leastsquares(r1, g2, r0_cent, utheta, ibm[ibi].bladestructure[iblade].NumberofElements, NumberofElements_AL);
			}

	
//			cout << "transferbladedisplacement2flowblade: x,y,z->xhat,yhat,zhat \n";
			// projection: from frame rotating with pi because of pitch action to frame rotaing with the hub (xhat,yhat,zhat) 
			for (i=0; i<NumberofNodes_AL; i++) 
			{
				double fxhat, fyhat, fzhat, fx_tmp, fy_tmp, fz_tmp;

				// u
				fx_tmp=fx[i];
				fy_tmp=fy[i];
				fz_tmp=fz[i];

				xyzhat_fromxyz(Beta, fx_tmp, fy_tmp, fz_tmp, &fxhat, &fyhat, &fzhat);

//				cout << fx_tmp << "\t" << fy_tmp << "\t" << fz_tmp << endl;
//				cout << fxhat << "\t" << fyhat << "\t" << fzhat << endl;
//				cout << endl;

				fx[i]=fxhat;	
				fy[i]=fyhat;	
				fz[i]=fzhat;	

				// dudt
				fx_tmp=ux[i];
				fy_tmp=uy[i];
				fz_tmp=uz[i];

				xyzhat_fromxyz(Beta, fx_tmp, fy_tmp, fz_tmp, &fxhat, &fyhat, &fzhat);

				ux[i]=fxhat;	
				uy[i]=fyhat;	
				uz[i]=fzhat;	
			}

//			cout << "transferbladedisplacement2flowblade: xhat,yhat,zhat->X,Y,Z \n";
			// projection: from frame rotating with the hub (xhat,yhat,zhat) to the structure inertial frame (X,Y,Z)
//			double Phi = ibm[ibi].bladestructure[iblade].phi;
			// Phi is not from the torque cotroller yet
			//
			
//	    MPI_Barrier(PETSC_COMM_WORLD);
//	    PetscPrintf(PETSC_COMM_WORLD, "Node Loop  \n");

			int NumberofNodes_AL = ibm[ibi].n_v/ibm[ibi].num_blade;
		
			int ii = iblade * NumberofNodes_AL + 5; 

			double rx = ibm[ibi].x_bp[ii]-fsi[ibi].x_c;
			double ry = ibm[ibi].y_bp[ii]-fsi[ibi].y_c;
			double rz = ibm[ibi].z_bp[ii]-fsi[ibi].z_c;

			double rr = sqrt(rx*rx+ry*ry+rz*rz);	
			rx /= rr; ry /= rr; rz /= rr;

			double gx=nZ_x, gy=nZ_y, gz=nZ_z; // assuming y is the vertical direction in CFD

			double Phi = acos(rx*gx+ry*gy+rz*gz);

			for (i=0; i<NumberofNodes_AL; i++) 
			{
				double fxhat, fyhat, fzhat, fX, fY, fZ;

				// u
				fxhat=fx[i];
				fyhat=fy[i];
				fzhat=fz[i];

				XYZ_fromxyzhat(Phi, fxhat, fyhat, fzhat, &fX, &fY, &fZ);

//				cout << fxhat << "\t" << fyhat << "\t" << fzhat << endl;
//				cout << fX << "\t" << fY << "\t" << fZ << endl;
//				cout << endl;

				fx[i]=fX;	
				fy[i]=fY;	
				fz[i]=fZ;	

				// dudt
				fxhat=ux[i];
				fyhat=uy[i];
				fzhat=uz[i];

				XYZ_fromxyzhat(Phi, fxhat, fyhat, fzhat, &fX, &fY, &fZ);

				ux[i]=fX;	
				uy[i]=fY;	
				uz[i]=fZ;	

			}

//	    MPI_Barrier(PETSC_COMM_WORLD);
//	    PetscPrintf(PETSC_COMM_WORLD, "Node Loop End \n");
				
//			cout << "transferbladedisplacement2flowblade: X,Y,Z->inertial frame \n";
//			cout << nX_x << "\t" << nX_y << "\t" << nX_z << endl;
//			cout << nY_x << "\t" << nY_y << "\t" << nY_z << endl;
//			cout << nZ_x << "\t" << nZ_y << "\t" << nZ_z << endl;
			// projection: from the intertial frame of the structure (X,Y,Z) to inertial frame of the flow
			for (i=0; i<NumberofNodes_AL; i++) 
			{
				// u
				double fX=fx[i];	
				double fY=fy[i];	
				double fZ=fz[i];	

				fx[i]=fX*nX_x+fY*nY_x+fZ*nZ_x;
				fy[i]=fX*nX_y+fY*nY_y+fZ*nZ_y;
				fz[i]=fX*nX_z+fY*nY_z+fZ*nZ_z;

//				cout << fX << "\t" << fY << "\t" << fZ << endl;
//				cout << fx[i] << "\t" << fy[i] << "\t" << fz[i] << endl;
//				cout << endl;

				// dudt
				fX=ux[i];	
				fY=uy[i];	
				fZ=uz[i];	

				ux[i]=fX*nX_x+fY*nY_x+fZ*nZ_x;
				uy[i]=fX*nX_y+fY*nY_y+fZ*nZ_y;
				uz[i]=fX*nX_z+fY*nY_z+fZ*nZ_z;
			}
	
//			cout << "transferbladedisplacement2flowblade: transfer disp \n";
//


//	    MPI_Barrier(PETSC_COMM_WORLD);
//	    PetscPrintf(PETSC_COMM_WORLD, "Number of Nodes reflength %f \n", reflength_wt);
//	    MPI_Barrier(PETSC_COMM_WORLD);

			double ratio_V=refvel_cfd/refvel_wt;

//	    MPI_Barrier(PETSC_COMM_WORLD);
//	    PetscPrintf(PETSC_COMM_WORLD, "iblade nuberof nodes i %i   %i \n", iblade, NumberofNodes_AL );
//	    PetscPrintf(PETSC_COMM_WORLD, "Number of Nodes ratioV %f \n", ratio_V);
//	    MPI_Barrier(PETSC_COMM_WORLD);

			for (i=0; i<NumberofNodes_AL; i++) 
			{
				int ii = iblade * NumberofNodes_AL + i; 

				ibm[ibi].disp_x[ii]=fx[i]/reflength_wt;
				ibm[ibi].disp_y[ii]=fy[i]/reflength_wt;
				ibm[ibi].disp_z[ii]=fz[i]/reflength_wt;

				ibm[ibi].ux_struct[ii]=ux[i]*ratio_V;
				ibm[ibi].uy_struct[ii]=uy[i]*ratio_V;
				ibm[ibi].uz_struct[ii]=uz[i]*ratio_V;
			}


//

			if (torsion_turbinestructure) 
			{ 
				for (i=0; i<NumberofElements_AL; i++) 
				{
					int ii = iblade * NumberofElements_AL + i; 
					ibm[ibi].disp_theta[ii]=ftheta[i];
//	        MPI_Barrier(PETSC_COMM_WORLD);
					ibm[ibi].utheta_struct[ii]=utheta[i];
				}
			}

			// scale the displacement from real to CFD simulation

//			for (i=0; i<NumberofNodes_AL; i++) {
//				int ii = iblade * NumberofNodes_AL + i; 
//				ibm[ibi].disp_x[ii]=fx[i]/reflength_wt;
//				ibm[ibi].disp_y[ii]=fy[i]/reflength_wt;
//				ibm[ibi].disp_z[ii]=fz[i]/reflength_wt;
//				if (torsion_turbinestructure) 
//			}

		} // blades 	

	} // turbines 


	for (ibi=0; ibi<NumberOfTurbines; ibi++) 
	{
		int NElmts = (int)((double)(ibm[ibi].n_elmt/ibm[ibi].num_blade)+0.1);

	        FILE *f;
		char filen[80];  

		for (int nb=0; nb<ibm[ibi].num_blade; nb++)
		{
			int i = (nb+1)*NElmts-1;
        		sprintf(filen, "%s/tipDeformation_%03d_%d_CFD.dat", path,ibi,nb);
	        	f = fopen(filen, "a");

			double x = ibm[ibi].disp_x[i];
			double y = ibm[ibi].disp_y[i];
			double z = ibm[ibi].disp_z[i];
			double u = ibm[ibi].ux_struct[i];
			double v = ibm[ibi].uy_struct[i];
			double w = ibm[ibi].uz_struct[i];
      double pphi=ibm[ibi].bladestructure[nb].phi;
	        	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le \n", ti, pphi, x, y, z, u, v, w);

			fclose(f);
		}
	}

	return(0);
}



PetscErrorCode writefiles_bladestructure(IBMNodes *ibm)
{
 	int i, j, k;
	int ibi, iblade;

   	int rank=0;
    	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


	if (!rank) 
	{

		// write displacements on tips
		for (ibi=0; ibi<NumberOfTurbines; ibi++) 
		{
			for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) 
			{
				int NElmts = ibm[ibi].bladestructure[iblade].NumberofElements;
	        	        FILE *f;
				char filen[80];  
	        	        sprintf(filen, "%s/tipDeformation_%03d_%d.dat", path,ibi,iblade);
			        f = fopen(filen, "a");

				i = NElmts-1;
				double r = ibm[ibi].bladestructure[iblade].r[i];
    		        	double u = ibm[ibi].bladestructure[iblade].u[i];
    		        	double dudt = ibm[ibi].bladestructure[iblade].dudt[i];
    		        	double v = ibm[ibi].bladestructure[iblade].v[i];
    		        	double dvdt = ibm[ibi].bladestructure[iblade].dvdt[i];
				double theta = ibm[ibi].bladestructure[iblade].theta[i];
    			       	double dthetadt = ibm[ibi].bladestructure[iblade].dthetadt[i];
                  double pphi = ibm[ibi].bladestructure[iblade].phi;

                        	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le \n", ti, pphi, r, u, dudt, v, dvdt, theta, dthetadt);

	        		fclose(f);
			} // iblade 
		} // ibi

		if (ti==tistart || ti==tistart+1 || ti % tiout == 0) 
		{
			for (ibi=0; ibi<NumberOfTurbines; ibi++) 
			{
				for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) 
				{
		        	        FILE *f;
	    				char filen[80];  
		        	        sprintf(filen, "%s/bladeDeformation_u%06d_%03d_%d.dat", path,ti,ibi,iblade);
	        		        f = fopen(filen, "w");

					PetscFPrintf(PETSC_COMM_WORLD, f, "VARIABLES = r u u_o u_om1 u_om2 dudt dudt_o dudt_om1 dudt_om2\n");

					for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
					{
						double r = ibm[ibi].bladestructure[iblade].r[i];

    	    		        		double u = ibm[ibi].bladestructure[iblade].u[i];
						double u_o = ibm[ibi].bladestructure[iblade].u_o[i];
						double u_om1 = ibm[ibi].bladestructure[iblade].u_om1[i];
						double u_om2 = ibm[ibi].bladestructure[iblade].u_om2[i];
        
    	    		        		double dudt = ibm[ibi].bladestructure[iblade].dudt[i];
						double dudt_o = ibm[ibi].bladestructure[iblade].dudt_o[i];
						double dudt_om1 = ibm[ibi].bladestructure[iblade].dudt_om1[i];
						double dudt_om2 = ibm[ibi].bladestructure[iblade].dudt_om2[i];
        
        	                		PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le \n", r, u, u_o, u_om1, u_om2, dudt, dudt_o, dudt_om1, dudt_om2);
					}

	                		fclose(f);

		        	        sprintf(filen, "%s/bladeDeformation_v%06d_%03d_%d.dat", path,ti,ibi,iblade);
	        		        f = fopen(filen, "w");

					PetscFPrintf(PETSC_COMM_WORLD, f, "VARIABLES = r v v_o v_om1 v_om2 dvdt dvdt_o dvdt_om1 dvdt_om2\n");

					for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
					{
						double r = ibm[ibi].bladestructure[iblade].r[i];

    	    		        		double v = ibm[ibi].bladestructure[iblade].v[i];
						double v_o = ibm[ibi].bladestructure[iblade].v_o[i];
						double v_om1 = ibm[ibi].bladestructure[iblade].v_om1[i];
						double v_om2 = ibm[ibi].bladestructure[iblade].v_om2[i];
        
    	    		        		double dvdt = ibm[ibi].bladestructure[iblade].dvdt[i];
						double dvdt_o = ibm[ibi].bladestructure[iblade].dvdt_o[i];
						double dvdt_om1 = ibm[ibi].bladestructure[iblade].dvdt_om1[i];
						double dvdt_om2 = ibm[ibi].bladestructure[iblade].dvdt_om2[i];
        
        	                		PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le \n", r, v, v_o, v_om1, v_om2, dvdt, dvdt_o, dvdt_om1, dvdt_om2);
					}

	                		fclose(f);

					if (torsion_turbinestructure) 
					{
			        	        sprintf(filen, "%s/bladeDeformation_theta%06d_%03d_%d.dat", path,ti,ibi,iblade);
		        		        f = fopen(filen, "w");

						PetscFPrintf(PETSC_COMM_WORLD, f, "VARIABLES = r theta theta_o theta_om1 theta_om2 dthetadt dthetadt_o dthetadt_om1 dthetadt_om2 \n");

						for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
						{
							double r = ibm[ibi].bladestructure[iblade].r[i];

    	    		        			double theta = ibm[ibi].bladestructure[iblade].theta[i];
							double theta_o = ibm[ibi].bladestructure[iblade].theta_o[i];
							double theta_om1 = ibm[ibi].bladestructure[iblade].theta_om1[i];
							double theta_om2 = ibm[ibi].bladestructure[iblade].theta_om2[i];
        
    	    			        		double dthetadt = ibm[ibi].bladestructure[iblade].dthetadt[i];
							double dthetadt_o = ibm[ibi].bladestructure[iblade].dthetadt_o[i];
							double dthetadt_om1 = ibm[ibi].bladestructure[iblade].dthetadt_om1[i];
							double dthetadt_om2 = ibm[ibi].bladestructure[iblade].dthetadt_om2[i];
        
        	                			PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le\n", r, theta, theta_o, theta_om1, theta_om2, dthetadt, dthetadt_o, dthetadt_om1, dthetadt_om2);
						}

		                		fclose(f);
					} // torsion

				} // iblade 
			} // turbine 
		} // ti
	} // rank


	return(0);
}

// Read the blade deformation profiles for restart 
PetscErrorCode readfiles_bladestructure(IBMNodes *ibm)
{
 	PetscPrintf(PETSC_COMM_WORLD, "PetscErrorCode readfiles_bladestructure\n");
 	int i, j, k;
	int ibi, iblade;

   	int rank=0;
    	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  	char string[60];

 	PetscPrintf(PETSC_COMM_WORLD, "NUMBER OF TURBINES LOOP\n");
	for (ibi=0; ibi<NumberOfTurbines; ibi++) 
	{
		for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) 
		{
			if (!rank) 
			{

                     	        FILE *fd;
                     		char filen[80];  
                     	        sprintf(filen, "%s/bladeDeformation_u%06d_%03d_%d.dat", path,tistart,ibi,iblade);

 	                   PetscPrintf(PETSC_COMM_WORLD, "Opening file %s \n", filen);

                     	        fd = fopen(filen, "r");

 	                   PetscPrintf(PETSC_COMM_WORLD, "File Open %s \n", filen);
            
      		        	fgets(string, 128, fd);
 	                   PetscPrintf(PETSC_COMM_WORLD, "Header %s \n", string);
            
 	         //          PetscPrintf(PETSC_COMM_WORLD, " Loop through blade number of elements %i \n", ibm[ibi].bladestructure[iblade].NumberofElements);

                     		for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
				{

            				double r, u, u_o, u_om1, u_om2, dudt, dudt_o, dudt_om1, dudt_om2; 


                     fscanf(fd, "%le %le %le %le %le %le %le %le %le \n", &r, &u, &u_o, &u_om1, &u_om2, &dudt, &dudt_o, &dudt_om1, &dudt_om2 );

                             		ibm[ibi].bladestructure[iblade].r[i]=r;
 	                   PetscPrintf(PETSC_COMM_WORLD, " Blade r %i %f\n",i,r);
                          		ibm[ibi].bladestructure[iblade].u[i]=u;
                     			ibm[ibi].bladestructure[iblade].u_o[i]=u_o;
                     			ibm[ibi].bladestructure[iblade].u_om1[i]=u_om1;
                     			ibm[ibi].bladestructure[iblade].u_om2[i]=u_om2;
                                            
                             		ibm[ibi].bladestructure[iblade].dudt[i]=dudt;
                     			ibm[ibi].bladestructure[iblade].dudt_o[i]=dudt_o;
                     			ibm[ibi].bladestructure[iblade].dudt_om1[i]=dudt_om1;
                     			ibm[ibi].bladestructure[iblade].dudt_om2[i]=dudt_om2;
                     		}
            
                     		fclose(fd);


                     	        sprintf(filen, "%s/bladeDeformation_v%06d_%03d_%d.dat", path,tistart,ibi,iblade);
                     	        fd = fopen(filen, "r");
            
      		 //       	fgets(string, 128, f);
            
      		        	fgets(string, 128, fd);
                     		for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
				{
            				double r, v, v_o, v_om1, v_om2, dvdt, dvdt_o, dvdt_om1, dvdt_om2; 
                             		fscanf(fd, "%le %le %le %le %le %le %le %le %le \n", &r, &v, &v_o, &v_om1, &v_om2, &dvdt, &dvdt_o, &dvdt_om1, &dvdt_om2);

                          ibm[ibi].bladestructure[iblade].r[i]=r;
                          ibm[ibi].bladestructure[iblade].v[i]=v;
                     			ibm[ibi].bladestructure[iblade].v_o[i]=v_o;
                     			ibm[ibi].bladestructure[iblade].v_om1[i]=v_om1;
                     			ibm[ibi].bladestructure[iblade].v_om2[i]=v_om2;
                                            
                             		ibm[ibi].bladestructure[iblade].dvdt[i]=dvdt;
                     			ibm[ibi].bladestructure[iblade].dvdt_o[i]=dvdt_o;
                     			ibm[ibi].bladestructure[iblade].dvdt_om1[i]=dvdt_om1;
                     			ibm[ibi].bladestructure[iblade].dvdt_om2[i]=dvdt_om2;
                     		}
            
                     		fclose(fd);

 	                     PetscPrintf(PETSC_COMM_WORLD, "Reading blade Deformation data V\n");

    				MPI_Bcast(ibm[ibi].bladestructure[iblade].r, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);


				if (torsion_turbinestructure) 
				{
                     	        	sprintf(filen, "%s/bladeDeformation_theta%06d_%03d_%d.dat", path,tistart,ibi,iblade);
	                     	        fd = fopen(filen, "r");
            
      			        	fgets(string, 128, fd);
            
                	     		for (i=0; i<ibm[ibi].bladestructure[iblade].NumberofElements; i++) 
					{
            					double r, theta, theta_o, theta_om1, theta_om2, dthetadt, dthetadt_o, dthetadt_om1, dthetadt_om2; 
                             			fscanf(fd, "%le %le %le %le %le %le %le %le %le \n", &r, &theta, &theta_o, &theta_om1, &theta_om2, &dthetadt, &dthetadt_o, &dthetadt_om1, &dthetadt_om2 );

                             			ibm[ibi].bladestructure[iblade].r[i]=r;
	                             		ibm[ibi].bladestructure[iblade].theta[i]=theta;
        	             			ibm[ibi].bladestructure[iblade].theta_o[i]=theta_o;
                	     			ibm[ibi].bladestructure[iblade].theta_om1[i]=theta_om1;
                     				ibm[ibi].bladestructure[iblade].theta_om2[i]=theta_om2;
                                            
                             			ibm[ibi].bladestructure[iblade].dthetadt[i]=dthetadt;
                     				ibm[ibi].bladestructure[iblade].dthetadt_o[i]=dthetadt_o;
	                     			ibm[ibi].bladestructure[iblade].dthetadt_om1[i]=dthetadt_om1;
	                     			ibm[ibi].bladestructure[iblade].dthetadt_om2[i]=dthetadt_om2;
        	             		}
            	
                	     		fclose(fd);

    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
				} // torsion structure 

			} // rank=0
			else 
			{

    				MPI_Bcast(ibm[ibi].bladestructure[iblade].r, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].u_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dudt_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].v_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    				MPI_Bcast(ibm[ibi].bladestructure[iblade].dvdt_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);


				if (torsion_turbinestructure) 
				{
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].theta_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt_o, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt_om1, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
    					MPI_Bcast(ibm[ibi].bladestructure[iblade].dthetadt_om2, ibm[ibi].bladestructure[iblade].NumberofElements, MPIU_REAL, 0, PETSC_COMM_WORLD);
				} // torsion structure 


			}
		} // blades 
	} // turbines 


	return(0);
}

// transfer the force from actuator line to actuator surface 
PetscErrorCode DispProjection_l2s(IBMNodes *ibm_surface, IBMNodes *ibm_line, int NumberOfObjects) 
{
	int ibi, elmt_s, elmt_l, l, node_s;
	vector<double> count;	

  	for (ibi=0; ibi<NumberOfObjects; ibi++) 
	{
		count.resize(ibm_surface[ibi].n_v);
	
		for (node_s=0; node_s<ibm_surface[ibi].n_v; node_s++)
		{
			count.at(node_s) = 1.e-19;
			ibm_surface[ibi].disp_x[node_s] = 0.0;
			ibm_surface[ibi].disp_y[node_s] = 0.0;
			ibm_surface[ibi].disp_z[node_s] = 0.0;

			ibm_surface[ibi].ux_struct[node_s] = 0.0;
			ibm_surface[ibi].uy_struct[node_s] = 0.0;
			ibm_surface[ibi].uz_struct[node_s] = 0.0;
		}

		for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) 
		{
			int nv1_s = ibm_surface[ibi].nv1[elmt_s];
			int nv2_s = ibm_surface[ibi].nv2[elmt_s];
			int nv3_s = ibm_surface[ibi].nv3[elmt_s];

			elmt_l=ibm_surface[ibi].s2l[elmt_s];
			int nv1_l = ibm_line[ibi].nv1[elmt_l];
			int nv2_l = ibm_line[ibi].nv2[elmt_l];

                        ibm_surface[ibi].disp_x[nv1_s] += ibm_line[ibi].disp_x[nv1_l];
                        ibm_surface[ibi].disp_x[nv1_s] += ibm_line[ibi].disp_x[nv2_l];
                        ibm_surface[ibi].disp_y[nv1_s] += ibm_line[ibi].disp_y[nv1_l];
                        ibm_surface[ibi].disp_y[nv1_s] += ibm_line[ibi].disp_y[nv2_l];
                        ibm_surface[ibi].disp_z[nv1_s] += ibm_line[ibi].disp_z[nv1_l];
                        ibm_surface[ibi].disp_z[nv1_s] += ibm_line[ibi].disp_z[nv2_l];

                        ibm_surface[ibi].ux_struct[nv1_s] += ibm_line[ibi].ux_struct[nv1_l];
                        ibm_surface[ibi].ux_struct[nv1_s] += ibm_line[ibi].ux_struct[nv2_l];
                        ibm_surface[ibi].uy_struct[nv1_s] += ibm_line[ibi].uy_struct[nv1_l];
                        ibm_surface[ibi].uy_struct[nv1_s] += ibm_line[ibi].uy_struct[nv2_l];
                        ibm_surface[ibi].uz_struct[nv1_s] += ibm_line[ibi].uz_struct[nv1_l];
                        ibm_surface[ibi].uz_struct[nv1_s] += ibm_line[ibi].uz_struct[nv2_l];

			count.at(nv1_s) += 2;

                        ibm_surface[ibi].disp_x[nv2_s] += ibm_line[ibi].disp_x[nv1_l];
                        ibm_surface[ibi].disp_x[nv2_s] += ibm_line[ibi].disp_x[nv2_l];
                        ibm_surface[ibi].disp_y[nv2_s] += ibm_line[ibi].disp_y[nv1_l];
                        ibm_surface[ibi].disp_y[nv2_s] += ibm_line[ibi].disp_y[nv2_l];
                        ibm_surface[ibi].disp_z[nv2_s] += ibm_line[ibi].disp_z[nv1_l];
                        ibm_surface[ibi].disp_z[nv2_s] += ibm_line[ibi].disp_z[nv2_l];

                        ibm_surface[ibi].ux_struct[nv2_s] += ibm_line[ibi].ux_struct[nv1_l];
                        ibm_surface[ibi].ux_struct[nv2_s] += ibm_line[ibi].ux_struct[nv2_l];
                        ibm_surface[ibi].uy_struct[nv2_s] += ibm_line[ibi].uy_struct[nv1_l];
                        ibm_surface[ibi].uy_struct[nv2_s] += ibm_line[ibi].uy_struct[nv2_l];
                        ibm_surface[ibi].uz_struct[nv2_s] += ibm_line[ibi].uz_struct[nv1_l];
                        ibm_surface[ibi].uz_struct[nv2_s] += ibm_line[ibi].uz_struct[nv2_l];

			count.at(nv2_s) += 2;

                        ibm_surface[ibi].disp_x[nv3_s] += ibm_line[ibi].disp_x[nv1_l];
                        ibm_surface[ibi].disp_x[nv3_s] += ibm_line[ibi].disp_x[nv2_l];
                        ibm_surface[ibi].disp_y[nv3_s] += ibm_line[ibi].disp_y[nv1_l];
                        ibm_surface[ibi].disp_y[nv3_s] += ibm_line[ibi].disp_y[nv2_l];
                        ibm_surface[ibi].disp_z[nv3_s] += ibm_line[ibi].disp_z[nv1_l];
                        ibm_surface[ibi].disp_z[nv3_s] += ibm_line[ibi].disp_z[nv2_l];

                        ibm_surface[ibi].ux_struct[nv3_s] += ibm_line[ibi].ux_struct[nv1_l];
                        ibm_surface[ibi].ux_struct[nv3_s] += ibm_line[ibi].ux_struct[nv2_l];
                        ibm_surface[ibi].uy_struct[nv3_s] += ibm_line[ibi].uy_struct[nv1_l];
                        ibm_surface[ibi].uy_struct[nv3_s] += ibm_line[ibi].uy_struct[nv2_l];
                        ibm_surface[ibi].uz_struct[nv3_s] += ibm_line[ibi].uz_struct[nv1_l];
                        ibm_surface[ibi].uz_struct[nv3_s] += ibm_line[ibi].uz_struct[nv2_l];

			count.at(nv3_s) += 2;
		}

		for (node_s=0; node_s<ibm_surface[ibi].n_v; node_s++)
		{
			ibm_surface[ibi].disp_x[node_s] /= count.at(node_s);
			ibm_surface[ibi].disp_y[node_s] /= count.at(node_s);
			ibm_surface[ibi].disp_z[node_s] /= count.at(node_s);

			ibm_surface[ibi].ux_struct[node_s] /= count.at(node_s);
			ibm_surface[ibi].uy_struct[node_s] /= count.at(node_s);
			ibm_surface[ibi].uz_struct[node_s] /= count.at(node_s);
		}

  	}

	return (0);
}

PetscErrorCode resetMatrixCoeff(Beam *blade)
{
 	int i;

	for (i=0; i<blade->NumberofElements; i++) 
	{
    blade->coefu_dudr[i]=0.;
    blade->coefu_d2udr2[i]=0.;
    blade->coefu_d3udr3[i]=0.;
    blade->coefu_d4udr4[i]=0.;

    blade->coefu_dvdr[i]=0.;
    blade->coefu_d2vdr2[i]=0.;
    blade->coefu_d3vdr3[i]=0.;
    blade->coefu_d4vdr4[i]=0.;

    blade->coefv_dvdr[i]=0.;
    blade->coefv_d2vdr2[i]=0.;
    blade->coefv_d3vdr3[i]=0.;
    blade->coefv_d4vdr4[i]=0.;

    blade->coefv_dudr[i]=0.;
    blade->coefv_d2udr2[i]=0.;
    blade->coefv_d3udr3[i]=0.;
    blade->coefv_d4udr4[i]=0.;

	}

	return(0);
}
