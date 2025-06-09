/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

extern Vec LevelSet, LevelSet0, LevelSet_o;
extern double M(double a, double b);
extern int immersed, NumberOfBodies;
extern int i_periodic, j_periodic, k_periodic;

extern PetscInt ti, tiout;
extern double wave_source_cent_z;

void Initialize_wave(UserCtx *user)
{
	int i;
	int NZMOD, NXMOD;
	int size_wave=1;
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	PetscMalloc(sizeof(WAVEInfo), &user[0].wave_inf);
	if(wave_momentum_source==1 || wave_momentum_source==4){
		if(!rank) { // root processor reads the wave data
			FILE *fd;
			char filen2[256];
			char path_wave[256];
			sprintf(path_wave, "./");
			PetscOptionsGetString(PETSC_NULL,"-path_wave", path_wave, 256, PETSC_NULL);		
			sprintf(filen2, "%s/WAVE_info%06d.dat", path_wave, 0);
			fd = fopen(filen2, "r");
			if(!fd) {
				printf("\n******************* Cannot open %s ! *******************\n", path_wave);
			}
			else {
				printf("\n**Succesfully read wave file in %s/WAVE_info%06d.dat **\n", path_wave, 0);
			}
			double temp;
			fscanf(fd,"%le %d %d %le %le\n", &temp, &user[0].wave_inf[0].NZMOD, &user[0].wave_inf[0].NXMOD, &user[0].wave_inf[0].PEZ, &user[0].wave_inf[0].PEX);
			NZMOD=user[0].wave_inf[0].NZMOD;
			NXMOD=user[0].wave_inf[0].NXMOD;
			printf("NZMOD: %d,  NXMOD: %d\n", NZMOD, NXMOD);
			MPI_Bcast(&NZMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			fclose(fd);
		}
		else if (rank) {	
			MPI_Bcast(&NZMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			user[0].wave_inf[0].NZMOD=NZMOD;
			user[0].wave_inf[0].NXMOD=NXMOD;			
		}	
		size_wave=(NZMOD)*(NXMOD+1)/2;;
	}
	if(wave_momentum_source==2 || wave_momentum_source==3)size_wave=1;//one wave simulation
	if(!rank)printf("Memory allocated for importing WAVE info\n");
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_a);
	PetscMalloc(size_wave*sizeof(int)   , &user[0].wave_inf[0].WAVE_ind);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_theta);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_Kz);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_Kx);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_KK);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_P_0);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_P_1);	
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_src_dist);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_omega);	
	user[0].wave_inf[0].WAVE_ti=wave_start_read;
};

double r_wave(double omega_wave, double omega_peak){
	double r_wave;
	double sigma_wave;
	if (omega_wave<=omega_peak){
		sigma_wave=0.07;
	}else{
		sigma_wave=0.09;
	}
	r_wave=exp(-pow(omega_wave-omega_peak,2)/(2.*sigma_wave*sigma_wave*omega_peak*omega_peak));
	return r_wave;
}



double D_wave(double omega_wave, double wave_angle, double N_wave, double s_wave){
	double D_wave;
	D_wave=N_wave*pow(cos(-wave_angle/2.),2.*s_wave);
	return D_wave;
}

void Initialize_wave_IB(UserCtx *user)
{
	int i,j;
	int NZMOD, NXMOD;
	double temp_kk, temp_kk2;
	PetscMalloc(sizeof(WAVEInfo), &user[0].wave_inf);


	//Set the value for this parameters
	int 		NTHE=10;
	int 		NK=20;	
	int 		n_direction=NTHE;//set to odd number so that it has 0 angle direction 
	int 		n_wave=NK;

	double grav=9.81;//gravitational acceeleration
	double Fetch=7905.;
	double U_10=12.;//velocity at elevvation equal to 10m
	double alpha_wave=0.076*pow(U_10*U_10/Fetch/grav,.22);
	double omega_peak=0.855*grav/U_10;
	omega_peak=22.*pow(grav*grav/U_10/Fetch,.333333);
	double s_wave, mu_omega, DD_wave, omega_wave, wave_angle;
	double omega_max=4*omega_peak;
	double omega_min=0.1*omega_peak;
	double rr_wave;
	double r_wave(double omega_wave, double omega_peak);
	double D_wave(double omega_wave, double wave_angle, double N_wave, double s_wave);
	double S_jonswap, N_wave;
	double var_dir;
	double var_omega;	
	PetscInt	rank;
	int size_wave=(n_direction)*(n_wave);//total number of wave frequencies	
	if(!rank)printf("omega_peak, omega_max, omega_min: %le %le %le \n",omega_peak, omega_max, omega_min);
	

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if(!rank)printf("Memory allocated wave_IB\n");
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_a);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_theta);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_Kz);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_Kx);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_KK);
	PetscMalloc(size_wave*sizeof(double), &user[0].wave_inf[0].WAVE_omega);	
	


	
	//fortran code
	double	U10=12.;
	double	FETCH=7905.;
	double	G=9.8;
	double	PI = acos(-1.0);
	double	TWOPI = 2. * acos(-1.0);
	double	ALPHA_P = 0.076 * pow((pow(U10,2.) / G / FETCH),0.22);
	double	KP = pow((G * pow((pow(22.,3.) / U10 / FETCH),2.)),0.33333);
	double	OMEGAP = 22.*pow((G*G/U10/FETCH),(1./3.));
	double	GAMMA = 3.3	;
	double 	KMIN=.01;
	double	KMAX=15.;//314
	double 	SIGMA;
	double Lz_domain=25.24;
	double Lx_domain=12.62;	
	srand( time(NULL)) ;
	PetscMalloc((NK)*sizeof(double), &user[0].wave_inf[0].WVN_z);
	PetscMalloc((NTHE)*sizeof(double), &user[0].wave_inf[0].WVN_x);
	PetscMalloc(NTHE*(NK)*sizeof(double), &user[0].wave_inf[0].SK2D);	

	int wave_IB_periodic=1;
	for(i=0;i<NK;i++){
			j=-NK+2*i;
			user[0].wave_inf[0].WVN_z[i]=2.*PI/Lz_domain*(double)(i+1);
		//printf("WVN: %f \n",user[0].wave_inf[0].WVN[i]);
	}
	for(i=0;i<NTHE;i++){
			j=-NTHE+2*i;
			user[0].wave_inf[0].WVN_x[i]=2.*PI/Lx_domain*(double)(j+1);
		//printf("WVN: %f \n",user[0].wave_inf[0].WVN[i]);
	}	
	//user[0].wave_inf[0].WVN[NK]=exp(log(KMIN) + (double)(NK+1 - 1) * (log(KMAX) - log(KMIN)) / (double)(NK - 1));	

	for(j=0;j<NTHE;j++){
		//double THETA = -PI / 2.0 + (double)(j +1.0 - 1.0) / (double)NTHE * PI;
		//double THETA_p1 = -PI / 2.0 + (double)(j +2.0 - 1.0) / (double)NTHE * PI;
		for(i=0;i<NK;i++){
			user[0].wave_inf[0].WAVE_Kz			[i+j*(NK)]=user[0].wave_inf[0].WVN_z[i];
			user[0].wave_inf[0].WAVE_Kx			[i+j*(NK)]=user[0].wave_inf[0].WVN_x[j];			
			user[0].wave_inf[0].WAVE_KK[i+j*(NK)]=sqrt(pow(user[0].wave_inf[0].WVN_z[i],2.)+pow(user[0].wave_inf[0].WVN_x[j],2.));
			double THETA=atan(user[0].wave_inf[0].WVN_x[j]/user[0].wave_inf[0].WVN_z[i]);
			double OMEGA    = sqrt(G * user[0].wave_inf[0].WAVE_KK[i+j*(NK)]);
			if (user[0].wave_inf[0].WAVE_KK[i+j*(NK)] <= KP){
				 SIGMA = 0.07;
			}
			else{
				 SIGMA = 0.09		;		
			}				
			user[0].wave_inf[0].SK2D[i+j*(NK)] = ALPHA_P * G*G / pow(OMEGA,5.)
			* exp(-5./4.*pow((OMEGA/OMEGAP),(-4.)))
			* pow(GAMMA,(exp(-pow((OMEGA-OMEGAP),2.)/2./pow(SIGMA,2.)/pow(OMEGAP,2.))))
			* (4.0 / TWOPI) * pow(user[0].wave_inf[0].WVN_x[j]/user[0].wave_inf[0].WAVE_KK[i+j*(NK)],2.);
			user[0].wave_inf[0].SK2D[i+j*(NK)] = user[0].wave_inf[0].SK2D[i+j*(NK)] * pow(G,2.) / 2. / pow(OMEGA,3.);
			if(!rank)printf("SK2D: %f \n",user[0].wave_inf[0].SK2D[i+j*(NK)]);		

			user[0].wave_inf[0].WAVE_omega	[i+j*(NK)]=OMEGA;
			user[0].wave_inf[0].WAVE_a			[i+j*(NK)]=sqrt(2.*user[0].wave_inf[0].SK2D[i+j*(NK)]*TWOPI/Lz_domain*TWOPI/Lx_domain);

			int n = rand() % 20000;
			n -= 10000;
			user[0].wave_inf[0].WAVE_theta	[i+j*(NK)]=((double)n)/10000.*PI;
			if(!rank)printf("A %le, Kz %le, Kx %le, Omega %le, K %le, theta: %f \n",user[0].wave_inf[0].WAVE_a[i+j*(NK)],user[0].wave_inf[0].WAVE_Kz[i+j*(NK)],user[0].wave_inf[0].WAVE_Kx[i+j*(NK)],user[0].wave_inf[0].WAVE_omega[i+j*(NK)],user[0].wave_inf[0].WAVE_KK[i+j*(NK)],user[0].wave_inf[0].WAVE_theta	[i+j*(NK)] );
		}
	}
/*	double xxx=25.;
	double zzz=30.;
	double uxx=0.;
	double uyy=0.;
	double uzz=0.;
	for(i=0;i<size_wave;i++){
		double sinhKH=1./sinh(user[0].wave_inf[0].WAVE_KK[i]*wave_depth);
		double sinalpha=user[0].wave_inf[0].WAVE_Kx[i]/user[0].wave_inf[0].WAVE_KK[i];
		double cosalpha=user[0].wave_inf[0].WAVE_Kz[i]/user[0].wave_inf[0].WAVE_KK[i];
		double wave_aa=user[0].wave_inf[0].WAVE_a[i];
		double coef_0=user[0].wave_inf[0].WAVE_Kz[i]*zzz+user[0].wave_inf[0].WAVE_Kx[i]*xxx-user[0].wave_inf[0].WAVE_omega[i]*.25;
		double elev=wave_aa*cos(coef_0);
		double coef_1=elev*user[0].wave_inf[0].WAVE_omega[i]*cosh(user[0].wave_inf[0].WAVE_KK[i]*(wave_depth+elev))*sinhKH;
		//printf("elev: %f",elev);
		if(coef_1 != coef_1){
			printf("here \n");
			coef_1=0.;
		}
			//int n = rand() % 20000;
			//n -= 10000;
			//user[0].wave_inf[0].WAVE_theta	[i]=((double)n)/10000.*PI;		
		uxx += sinalpha*coef_1;
		uyy += user[0].wave_inf[0].WAVE_omega[i]*wave_aa*sinh(user[0].wave_inf[0].WAVE_KK[i]*(wave_depth+elev))*(sinhKH)*sin(coef_0);
		uzz += cosalpha*coef_1;
		if(!rank)printf("%f %f %f \n",user[0].wave_inf[0].WAVE_KK[i],coef_1,sinhKH);
	}
	*/
	if(!rank)printf("JONSWAP,    kp: %f,   Lp: %f  \n",KP,TWOPI/KP);
/*	
	//Computation of the spreading function D_wave(omega_wave, wave_angle)
	double coef_2=11.5*pow(grav/omega_peak/U_10,2.5);
	var_dir=3.1415926535897932384626433832795/((double)n_direction-1.);
	var_omega=(omega_max-omega_min)/((double)n_wave-1.);

	if(!rank)printf("coef_2, var_dir, var_omega: %le %le %le \n",coef_2, var_dir, var_omega);	

	for(i=0;i<n_direction-1;i++)
	for(j=0;j<n_wave-1;j++){
		wave_angle=3.1415926535897932384626433832795/((double)n_direction-1.)*(double)i-1.5707963267948966192313216916398;
		omega_wave=omega_min+(omega_max-omega_min)/((double)n_wave-1.)*(double)j;
		user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]=omega_wave*omega_wave/grav;
		user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]=user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]/tanh(user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]*wave_depth);
		user[0].wave_inf[0].WAVE_omega[i*(n_wave-1)+j]=omega_wave;
		if(omega_wave<=omega_peak){
			mu_omega=5.;
		}else{
			mu_omega=-2.5;
		}
		//Computation of the JONSWAP spectrum S_jonswap
		rr_wave=r_wave(omega_wave, omega_peak);	
		s_wave=coef_2*pow(omega_wave/omega_peak,mu_omega);
		N_wave=0.5/sqrt(3.1415926535897932384626433832795)*abs(gamma(s_wave+1.00000001)/gamma(s_wave+.50000001));		
		DD_wave=D_wave(omega_wave, wave_angle, N_wave, s_wave);
		S_jonswap=alpha_wave*grav*grav/pow(omega_wave,5.)*exp(-1.25*pow(omega_peak/omega_wave,4.))*pow(3.3,rr_wave);
		//if(!rank)printf("wave_angle, omega_wave, rr_wave: %le %le %le \n",wave_angle, omega_wave, rr_wave);			
		//if(!rank)printf("S_jonswap: %le \n",S_jonswap);	
		//if(!rank)printf("s_wave, N_wave, DD_wave: %le %le %le \n",s_wave, N_wave, DD_wave);			
		user[0].wave_inf[0].WAVE_a[i*(n_wave-1)+j]=sqrt(2.*S_jonswap*DD_wave*var_dir*var_omega);
		user[0].wave_inf[0].WAVE_Kz[i*(n_wave-1)+j]=user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]*cos(wave_angle);
		user[0].wave_inf[0].WAVE_Kx[i*(n_wave-1)+j]=user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]*sin(wave_angle);
		user[0].wave_inf[0].WAVE_theta[i*(n_wave-1)+j]=0.;
		if(i*(n_wave-1)+j==0){
			user[0].wave_inf[0].WAVE_a[i*(n_wave-1)+j]=.05;
			user[0].wave_inf[0].WAVE_Kz[i*(n_wave-1)+j]=2*3.1415/200.;
			user[0].wave_inf[0].WAVE_Kx[i*(n_wave-1)+j]=0.;
			user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]=user[0].wave_inf[0].WAVE_Kz[i*(n_wave-1)+j];
			user[0].wave_inf[0].WAVE_theta[i*(n_wave-1)+j]=0.;
			user[0].wave_inf[0].WAVE_omega[i*(n_wave-1)+j]=sqrt(user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]*grav*tanh(user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j]*wave_depth));
		}
		if(!rank)printf("A %le, Kz %le, Kx %le, Omega %le, K %le \n",user[0].wave_inf[0].WAVE_a[i*(n_wave-1)+j],user[0].wave_inf[0].WAVE_Kz[i*(n_wave-1)+j],user[0].wave_inf[0].WAVE_Kx[i*(n_wave-1)+j],user[0].wave_inf[0].WAVE_omega[i*(n_wave-1)+j],user[0].wave_inf[0].WAVE_KK[i*(n_wave-1)+j] );
	}
*/
	user[0].wave_inf[0].WAVE_num_max=size_wave;	
};


void Initialize_wind(UserCtx *user)
{
	int i,j;
	int NXMOD, NZ;
	double temp;
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if(!wave_momentum_source)PetscMalloc(sizeof(WAVEInfo), &user[0].wave_inf);
	if(!rank) { // root processor reads the wave data
		FILE *fd;
		char filen[256];
		char path_wind[256];
		sprintf(path_wind, "./");
		PetscOptionsGetString(PETSC_NULL,"-path_wind", path_wind, 256, PETSC_NULL);		
		sprintf(filen, "%s/WAVE_wind%06d.dat", path_wind, 0);
		fd = fopen(filen, "r");
		if(!fd) {
			printf("\n******************* Cannot open %s ! *******************\n", path_wind);
		}
		else {
			printf("\n****Succesfully read wind file in %s **\n", path_wind);
		}
		fscanf(fd,"%le %d %d\n", &temp, &NXMOD, &NZ);
		printf("NXMOD: %d,  NZ: %d\n", NXMOD, NZ);
		user[0].wave_inf[0].NZ=NZ;
		MPI_Bcast(&NZ, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		fclose(fd);
	}
	else if (rank) {		
		MPI_Bcast(&NZ, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	}	
	int size_wind=NXMOD*NZ;
	user[0].wave_inf[0].NZ=NZ;
	user[0].wave_inf[0].NXMOD=NXMOD;
	if(!rank)printf("Memory allocated for Importing WIND_info\n");
	//In the far field simulation using Shen Codes x,i is streamwise, y,j is spanwise, and z,k is vertical
	PetscMalloc(size_wind*sizeof(double), &user[0].wave_inf[0].WIND_U );
	PetscMalloc(size_wind*sizeof(double), &user[0].wave_inf[0].WIND_V );
	PetscMalloc(size_wind*sizeof(double), &user[0].wave_inf[0].WIND_W );
	PetscMalloc(size_wind*sizeof(double), &user[0].wave_inf[0].WIND_Y );
	PetscMalloc(size_wind*sizeof(double), &user[0].wave_inf[0].WIND_Z );
	user[0].wave_inf[0].WIND_ti=wind_start_read;
	
	DA da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;	

	user[0].wave_inf[0].WIND_id_x =(int**) malloc ( info.xm * sizeof(int*));
	user[0].wave_inf[0].WIND_id_y =(int**) malloc ( info.xm * sizeof(int*));	
	for(i=0;i<info.xm;i++){
		user[0].wave_inf[0].WIND_id_x[i] = (int*)malloc ( info.ym * sizeof(int));
		user[0].wave_inf[0].WIND_id_y[i] = (int*)malloc ( info.ym * sizeof(int));			
		for(j=0;j<info.ym;j++){
			user[0].wave_inf[0].WIND_id_x[i][j]=0;
			user[0].wave_inf[0].WIND_id_y[i][j]=0;

		}
	}
	//Farfield 
	//z,k=vertical direction y,j=spanwise direction  x,i=streamwise direction  
	//NearField
	//y,j=vertical dir       x,i=span                z,k=stream
	if(!rank)printf("wind initialitzation completed \n");
};

void WIND_vel_interpolate(UserCtx *user, double *u1, double *v1, double *w1, double x1, double y1, int i, int j)
{
	int NXMOD=user[0].wave_inf[0].NXMOD;
	int jj =(i  )+(j  )*NXMOD;
	int jjj=(i+1)+(j  )*NXMOD;
	int jjk=(i  )+(j+1)*NXMOD;	
	double vary_varz=(user[0].wave_inf[0].WIND_Y[jjj]-user[0].wave_inf[0].WIND_Y[jj])*(user[0].wave_inf[0].WIND_Z[jjk]-user[0].wave_inf[0].WIND_Z[jj]);
	double A1=(user[0].wave_inf[0].WIND_Y[jjj]-x1)*(user[0].wave_inf[0].WIND_Z[jjk]-y1);
	double A2=(x1-user[0].wave_inf[0].WIND_Y[jj])*(user[0].wave_inf[0].WIND_Z[jjk]-y1);
	double A3=(user[0].wave_inf[0].WIND_Y[jjj]-x1)*(y1-user[0].wave_inf[0].WIND_Z[jj]);
	double A4=(x1-user[0].wave_inf[0].WIND_Y[jj])*(y1-user[0].wave_inf[0].WIND_Z[jj]);
	*u1=user[0].wave_inf[0].WIND_V[(i  )+(j  )*NXMOD]*A1+
			user[0].wave_inf[0].WIND_V[(i+1)+(j  )*NXMOD]*A2+
			user[0].wave_inf[0].WIND_V[(i  )+(j+1)*NXMOD]*A3+
			user[0].wave_inf[0].WIND_V[(i+1)+(j+1)*NXMOD]*A4;
	*u1=*u1/vary_varz;
	*v1=user[0].wave_inf[0].WIND_W[(i  )+(j  )*NXMOD]*A1+
			user[0].wave_inf[0].WIND_W[(i+1)+(j  )*NXMOD]*A2+
			user[0].wave_inf[0].WIND_W[(i  )+(j+1)*NXMOD]*A3+
			user[0].wave_inf[0].WIND_W[(i+1)+(j+1)*NXMOD]*A4;
	*v1=*v1/vary_varz;
	*w1=user[0].wave_inf[0].WIND_U[(i  )+(j  )*NXMOD]*A1+
			user[0].wave_inf[0].WIND_U[(i+1)+(j  )*NXMOD]*A2+
			user[0].wave_inf[0].WIND_U[(i  )+(j+1)*NXMOD]*A3+
			user[0].wave_inf[0].WIND_U[(i+1)+(j+1)*NXMOD]*A4;
	*w1=*w1/vary_varz;	
	//printf("u1:%f, v1:%f, w1:%f\n",*u1,*v1,*w1);

}

void WAVE_DATA_input(UserCtx *user)
{
/*This function is for setting the water wave field information (amplitude, frequencies, direction angle, ..) to be simulated. 
The main control option is defined by "wave_momentum_source" which can adopt one of the following values:
	1.	Simulation of a multi-frequency wave field importing the wave information from an external 
			file (from far-field solver)
	2.	Simulation of a single frequency wave field with defined at the control file
	3.	Simulation of a single frequency wave field using periodic boundary conditions. Wave defined 
			at control file
	4.	Simulation of a multi-frequency wave field using periodic boundary conditions and importing wave 
			information from external file (far-field).
*/
	int i,l,m;
	int ti_wave;
	double time;
	int WAVE_num_max;
	double temp_kk,temp_kk2;
	double twopi=2.*acos(-1.),onepi=acos(-1.), onepi2=onepi*onepi;
	double eps, eps_number=2.00001,ff;
	//double ro_w=1000, ro_air=1.2;
	double ro_w=rho_water, ro_air=rho_air;
	double grav=fabs(gravity_z)+fabs(gravity_y);//gravity must be defined in one direction only;
	double grav2=grav*grav;
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int NZMOD,NXMOD;
	int WAVE_ti=user[0].wave_inf[0].WAVE_ti;
	if(!rank)printf("Reading WAVE info for WAVE time step: \n",WAVE_ti);	
	int ind_max=0;
	//double L_ref=1000.0;
	if(wave_momentum_source==1){
	
		if(!rank) { // root processor reads the wave data
			FILE *fd;
			char filen[256];
			char path_wave[256];
			sprintf(path_wave, "./");
			PetscOptionsGetString(PETSC_NULL,"-path_wave", path_wave, 256, PETSC_NULL);		
			sprintf(filen, "%s/WAVE_info%06d.dat", path_wave, WAVE_ti);
			fd = fopen(filen, "r");
			if(!fd) {
				printf("\n******************* Cannot open %s ! *******************\n", path_wave);
			}
			else {
				printf("\n**Succesfully read wave file in %s/WAVE_info%06d.dat **\n", path_wave, WAVE_ti);
			}		
			fscanf(fd,"%le %d %d %le %le\n", &time, &user[0].wave_inf[0].NZMOD, &user[0].wave_inf[0].NXMOD, &user[0].wave_inf[0].PEZ, &user[0].wave_inf[0].PEX);
			NZMOD=user[0].wave_inf[0].NZMOD;
			NXMOD=user[0].wave_inf[0].NXMOD;
			
			MPI_Bcast(&time, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NZMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].PEZ, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].PEX, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

			
			WAVE_num_max = (NZMOD)*(NXMOD+1)/2;
			user[0].wave_inf[0].WAVE_num_max = WAVE_num_max;
			
			for(i=0;i<WAVE_num_max;i++){
				fscanf(fd,"%le %le\n",&user[0].wave_inf[0].WAVE_a[i],&user[0].wave_inf[0].WAVE_theta[i]);
				user[0].wave_inf[0].WAVE_a[i]*=wave_wind_reflength;
			}
			MPI_Bcast(&user[0].wave_inf[0].WAVE_a[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WAVE_theta[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			
			printf("\nReading WAVE_info%06d.dat file:\n",WAVE_ti);
			printf("time=%le, NZMOD=%d, NXMOD=%d, PEZ=%le, PEX=%le \n", time, NZMOD, NXMOD, user[0].wave_inf[0].PEZ, user[0].wave_inf[0].PEX);
			for(i=0;i<4;i++){
				printf("%le %le\n",user[0].wave_inf[0].WAVE_a[i],user[0].wave_inf[0].WAVE_theta[i]);
			}
			printf("...     ...            WAVE_num_max=%d  \n",WAVE_num_max);
			printf("%le %le\n",user[0].wave_inf[0].WAVE_a[WAVE_num_max-1],user[0].wave_inf[0].WAVE_theta[WAVE_num_max-1]);
			printf("Reading WAVE_info%06d.dat file completed\n\n",WAVE_ti);
			fclose(fd);
		}
		else if (rank) {
		
			MPI_Bcast(&time, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].NZMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].PEZ, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].PEX, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			NZMOD=user[0].wave_inf[0].NZMOD;
			NXMOD=user[0].wave_inf[0].NXMOD;
			WAVE_num_max = (NZMOD)*(NXMOD+1)/2;
			user[0].wave_inf[0].WAVE_num_max = WAVE_num_max;
			MPI_Bcast(&user[0].wave_inf[0].WAVE_a[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WAVE_theta[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
		
		}
		


		//defining wave numbers in x and y direction
		for (i=0; i<NZMOD/2; i++){
			user[0].wave_inf[0].WAVE_Kz[i]=(1.0+(double)i)*user[0].wave_inf[0].PEZ/wave_wind_reflength;
			user[0].wave_inf[0].WAVE_Kx[i]=0.0;
		}
		for (m=1; m<NXMOD/2+1; m++)
		for (l=0; l<NZMOD/2; l++){
			user[0].wave_inf[0].WAVE_Kz[(2*m-1)*(NZMOD/2)+l]=(1.0+(double)l)*user[0].wave_inf[0].PEZ/wave_wind_reflength;
			user[0].wave_inf[0].WAVE_Kx[(2*m-1)*(NZMOD/2)+l]=((double)m)*user[0].wave_inf[0].PEX/wave_wind_reflength;
			user[0].wave_inf[0].WAVE_Kz[(2*m)*(NZMOD/2)+l]=(1.0+(double)l)*user[0].wave_inf[0].PEZ/wave_wind_reflength;
			user[0].wave_inf[0].WAVE_Kx[(2*m)*(NZMOD/2)+l]=(-(double)m)*user[0].wave_inf[0].PEX/wave_wind_reflength;
		}

		for(i=0;i<WAVE_num_max;i++){
			temp_kk=user[0].wave_inf[0].WAVE_Kz[i]*user[0].wave_inf[0].WAVE_Kz[i]+user[0].wave_inf[0].WAVE_Kx[i]*user[0].wave_inf[0].WAVE_Kx[i];
			temp_kk2=sqrt(temp_kk);
			//user[0].wave_inf[0].WAVE_omega[i]=sqrt(temp_kk2/1.84);//Froud number instead of gravity
			//user[0].wave_inf[0].WAVE_omega[i]=sqrt(temp_kk2*9.81);//Froud number instead of gravity
			user[0].wave_inf[0].WAVE_omega[i]=sqrt(temp_kk2*grav*tanh(temp_kk2*wave_depth));
			eps=twopi/user[0].wave_inf[0].WAVE_Kz[i]/eps_number;
			ff=  onepi2*sin(user[0].wave_inf[0].WAVE_Kz[i]*eps)/(user[0].wave_inf[0].WAVE_Kz[i]*(onepi2-eps*eps*user[0].wave_inf[0].WAVE_Kz[i]*user[0].wave_inf[0].WAVE_Kz[i]));
			user[0].wave_inf[0].WAVE_P_0[i] = user[0].wave_inf[0].WAVE_a[i]*grav2*eps*2.0*ro_w/(user[0].wave_inf[0].WAVE_omega[i]*user[0].wave_inf[0].WAVE_omega[i]*ff*(ro_air+ro_w))*(user[0].wave_inf[0].WAVE_Kz[i])/(temp_kk2);
			user[0].wave_inf[0].WAVE_src_dist[i]	= eps;
				if(user[0].wave_inf[0].WAVE_a[i]>1.e-05*wave_wind_reflength){
					user[0].wave_inf[0].WAVE_ind[ind_max]=i;
					ind_max+=1;
				}
		}


		if(!rank)printf("Number of waves considered: Wave ind max = %d\n",ind_max);
		user[0].wave_inf[0].WAVE_ind_max=ind_max;
		user[0].wave_inf[0].WAVE_ti+=1;
		user[0].wave_inf[0].WAVE_time_since_read=0.;
		if(user[0].wave_inf[0].WAVE_ti==wave_recicle)user[0].wave_inf[0].WAVE_ti=0;//recicle
	}
	if(wave_momentum_source==2){
		i=0;
		double angle_=wave_angle_single;
		double K_wave=wave_K_single;

		double ramp;//ramp up function at wave_ti_start;
		if(ti<wave_ti_start)ramp=0.;
		if(ti>=wave_ti_start || ti<wave_ti_start+200)ramp=1./200.*(double)(ti-wave_ti_start);
		if(ti>=wave_ti_start+200)ramp=1.;
		user[0].wave_inf[0].WAVE_Kz[i]=K_wave*cos(angle_);
		user[0].wave_inf[0].WAVE_Kx[i]=K_wave*sin(angle_);
		user[0].wave_inf[0].WAVE_a[i]=wave_a_single;
		user[0].wave_inf[0].WAVE_theta[i]=0.;
		user[0].wave_inf[0].WAVE_num_max=1;
		temp_kk=user[0].wave_inf[0].WAVE_Kz[i]*user[0].wave_inf[0].WAVE_Kz[i]+user[0].wave_inf[0].WAVE_Kx[i]*user[0].wave_inf[0].WAVE_Kx[i];
		temp_kk2=sqrt(temp_kk);
		user[0].wave_inf[0].WAVE_omega[i]=sqrt(temp_kk2*grav*tanh(temp_kk2*wave_depth));
		eps=twopi/user[0].wave_inf[0].WAVE_Kz[i]/eps_number;
		//eps=twopi/K_wave/eps_number*cos(angle_);
		ff=  onepi2*sin(user[0].wave_inf[0].WAVE_Kz[i]*eps)/(user[0].wave_inf[0].WAVE_Kz[i]*(onepi2-eps*eps*user[0].wave_inf[0].WAVE_Kz[i]*user[0].wave_inf[0].WAVE_Kz[i]));
		user[0].wave_inf[0].WAVE_P_0[i] = ramp*user[0].wave_inf[0].WAVE_a[i]*grav2*eps*2.*ro_w/(user[0].wave_inf[0].WAVE_omega[i]*user[0].wave_inf[0].WAVE_omega[i]*ff*(ro_air+ro_w))*(user[0].wave_inf[0].WAVE_Kz[i])/(temp_kk2);
		user[0].wave_inf[0].WAVE_src_dist[i] = eps;
		user[0].wave_inf[0].WAVE_ind_max=1;
		user[0].wave_inf[0].WAVE_ind[0]=0;
		if (rank==0)printf("Information for wave generation: ti %d, wave_ti_start %d, ramp function %f\n", ti, wave_ti_start, ramp);
	}
	if(wave_momentum_source==3){
		i=0;
		double angle_=wave_angle_single;
		double K_wave=wave_K_single;
		user[0].wave_inf[0].WAVE_Kz[i]=K_wave*cos(angle_);
		user[0].wave_inf[0].WAVE_Kx[i]=K_wave*sin(angle_);			
		user[0].wave_inf[0].WAVE_a[i]=wave_a_single;
		user[0].wave_inf[0].WAVE_theta[i]=0.;			
		user[0].wave_inf[0].WAVE_num_max=1;
		temp_kk=user[0].wave_inf[0].WAVE_Kz[i]*user[0].wave_inf[0].WAVE_Kz[i]+user[0].wave_inf[0].WAVE_Kx[i]*user[0].wave_inf[0].WAVE_Kx[i];
		temp_kk2=sqrt(temp_kk);
		user[0].wave_inf[0].WAVE_omega[i]=sqrt(temp_kk2*grav*tanh(temp_kk2*wave_depth));
		eps=twopi/user[0].wave_inf[0].WAVE_Kz[i]/eps_number;
		//eps=twopi/K_wave/eps_number*cos(angle_);
		
		double ramp, ramp2;//ramp up function at wave_ti_start;		
		double time=user->dt*ti;
		double delta_t=.08*twopi/user[0].wave_inf[0].WAVE_omega[i];
		ramp = 0.5/delta_t*(1.+cos(.5*twopi/delta_t*(time-delta_t)));
		if (time-delta_t>delta_t)ramp=0.;
		ramp2= 0.5/delta_t*(1.+cos(.5*twopi/delta_t*(time-delta_t-.25*twopi/user[0].wave_inf[0].WAVE_omega[i])))	;
		if (time-.25*twopi/user[0].wave_inf[0].WAVE_omega[i]<0.)ramp2=0.;
		if (time-delta_t-.25*twopi/user[0].wave_inf[0].WAVE_omega[i]>delta_t)ramp2=0.;

		user[0].wave_inf[0].WAVE_P_0[i] = -ramp *user[0].wave_inf[0].WAVE_a[i]/user[0].wave_inf[0].WAVE_omega[i]*21.5053;
		user[0].wave_inf[0].WAVE_P_1[i] = -ramp2*user[0].wave_inf[0].WAVE_a[i]/user[0].wave_inf[0].WAVE_omega[i]*21.5053;
		user[0].wave_inf[0].WAVE_src_dist[i]	= 5000.;
		user[0].wave_inf[0].WAVE_ind_max=1;		
		user[0].wave_inf[0].WAVE_ind[0]=0;		
		if (rank==0)printf("Information for wave generation: ti %d, wave_ti_start %d, ramp function %f\n", ti, wave_ti_start, ramp);
	}
	if(wave_momentum_source==4){
		if(ti==tistart) {
			if(!rank) { // root processor reads the wave data
				FILE *fd;
				char filen[80];
				//sprintf(filen, "../wave/WAVE_info%06d.dat",WAVE_ti);
				sprintf(filen, "WAVE_info%06d.dat",0);
				fd = fopen(filen, "r");
				if(!fd)printf("Error reading wave input file WAVE_info\n");
				
				fscanf(fd,"%le %d %d %le %le\n", &time, &user[0].wave_inf[0].NZMOD, &user[0].wave_inf[0].NXMOD, &user[0].wave_inf[0].PEZ, &user[0].wave_inf[0].PEX);
				NZMOD=user[0].wave_inf[0].NZMOD;
				NXMOD=user[0].wave_inf[0].NXMOD;
				
				MPI_Bcast(&time, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&NZMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].PEZ, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].PEX, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

				
				WAVE_num_max = (NZMOD)*(NXMOD+1)/2;
				user[0].wave_inf[0].WAVE_num_max = WAVE_num_max;
				
				for(i=0;i<WAVE_num_max;i++){
					fscanf(fd,"%le %le\n",&user[0].wave_inf[0].WAVE_a[i],&user[0].wave_inf[0].WAVE_theta[i]);
					user[0].wave_inf[0].WAVE_a[i]*=wave_wind_reflength;
				}
				MPI_Bcast(&user[0].wave_inf[0].WAVE_a[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].WAVE_theta[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
				
				printf("\nReading WAVE_info%06d.dat file:\n",WAVE_ti);
				printf("time=%le, NZMOD=%d, NXMOD=%d, PEZ=%le, PEX=%le \n", time, NZMOD, NXMOD, user[0].wave_inf[0].PEZ, user[0].wave_inf[0].PEX);
				for(i=0;i<4;i++){
					printf("%le %le\n",user[0].wave_inf[0].WAVE_a[i],user[0].wave_inf[0].WAVE_theta[i]);
				}
				printf("...     ...            WAVE_num_max=%d  \n",WAVE_num_max);
				printf("%le %le\n",user[0].wave_inf[0].WAVE_a[WAVE_num_max-1],user[0].wave_inf[0].WAVE_theta[WAVE_num_max-1]);
				printf("Reading WAVE_info%06d.dat file completed\n\n",WAVE_ti);
				fclose(fd);
			}
			else if (rank) {
				MPI_Bcast(&time, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].NZMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].PEZ, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].PEX, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
				NZMOD=user[0].wave_inf[0].NZMOD;
				NXMOD=user[0].wave_inf[0].NXMOD;
				WAVE_num_max = (NZMOD)*(NXMOD+1)/2;
				user[0].wave_inf[0].WAVE_num_max = WAVE_num_max;
				MPI_Bcast(&user[0].wave_inf[0].WAVE_a[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
				MPI_Bcast(&user[0].wave_inf[0].WAVE_theta[0], WAVE_num_max, MPIU_REAL, 0, PETSC_COMM_WORLD);	
			}

			//defining wave numbers in x and y direction
			for (i=0; i<NZMOD/2; i++){
				user[0].wave_inf[0].WAVE_Kz[i]=-(1.0+(double)i)*user[0].wave_inf[0].PEZ/wave_wind_reflength;
				user[0].wave_inf[0].WAVE_Kx[i]=0.0;
			}
			for (m=1; m<NXMOD/2+1; m++)
			for (l=0; l<NZMOD/2; l++){
				user[0].wave_inf[0].WAVE_Kz[(2*m-1)*(NZMOD/2)+l]=-(1.0+(double)l)*user[0].wave_inf[0].PEZ/wave_wind_reflength;
				user[0].wave_inf[0].WAVE_Kx[(2*m-1)*(NZMOD/2)+l]=-((double)m)*user[0].wave_inf[0].PEX/wave_wind_reflength;
				user[0].wave_inf[0].WAVE_Kz[(2*m)*(NZMOD/2)+l]=-(1.0+(double)l)*user[0].wave_inf[0].PEZ/wave_wind_reflength;
				user[0].wave_inf[0].WAVE_Kx[(2*m)*(NZMOD/2)+l]=-(-(double)m)*user[0].wave_inf[0].PEX/wave_wind_reflength;			
			}
			for(i=0;i<WAVE_num_max;i++){
				if(user[0].wave_inf[0].WAVE_a[i]>1.e-05*wave_wind_reflength){
					user[0].wave_inf[0].WAVE_ind[ind_max]=i;
					ind_max+=1;
				}
			}
			user[0].wave_inf[0].WAVE_ind_max=ind_max;
		}
		
		int ii,jj;
		for(jj=0;jj<user[0].wave_inf[0].WAVE_ind_max;jj++){
			ii=user[0].wave_inf[0].WAVE_ind[jj];
			temp_kk=user[0].wave_inf[0].WAVE_Kz[ii]*user[0].wave_inf[0].WAVE_Kz[ii]+user[0].wave_inf[0].WAVE_Kx[ii]*user[0].wave_inf[0].WAVE_Kx[ii];
			temp_kk2=sqrt(temp_kk);
			user[0].wave_inf[0].WAVE_omega[ii]=sqrt(temp_kk2*grav*tanh(temp_kk2*wave_depth));
			double ramp1, ramp2;//ramp up function at wave_ti_start;		
			double time=user->dt*ti;
			double delta_t=50*user->dt;//.08*twopi/user[0].wave_inf[0].WAVE_omega[ii];
			ramp1 = 0.5/delta_t*(1.+cos(.5*twopi/delta_t*(time-delta_t)));
			if (time-delta_t>delta_t)ramp1=0.;
			ramp2= 0.5/delta_t*(1.+cos(.5*twopi/delta_t*(time-delta_t-.25*twopi/user[0].wave_inf[0].WAVE_omega[ii])))	;
			if (time-.25*twopi/user[0].wave_inf[0].WAVE_omega[ii]<0.)ramp2=0.;
			if (time-delta_t-.25*twopi/user[0].wave_inf[0].WAVE_omega[ii]>delta_t)ramp2=0.;
			double temp_2=user[0].wave_inf[0].WAVE_a[ii]/user[0].wave_inf[0].WAVE_omega[ii]*21.5053;
			user[0].wave_inf[0].WAVE_P_0[ii] = -ramp1*temp_2;
			user[0].wave_inf[0].WAVE_P_1[ii] = -ramp2*temp_2;		
		}
		

	}	
	
};

void WIND_DATA_input(UserCtx *user)
{
  int i,j,k,l,m;
  int WIND_ti=user[0].wave_inf[0].WIND_ti;
	int NXMOD,NZ,i_max;
	double temp,temp222;
	double time;
	double x1,y1,z1;
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);	
	if(!rank)printf("Reading WIND info for WAVE time step: %d\n",WIND_ti);
	
	//locate near field coordinates in the far field grid
	Cmpnts ***cent;	
	DA da = user->da, fda = user->fda;
	DALocalInfo	info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;	
	
	DAVecGetArray(fda, user->lCent, &cent);		

	if(air_flow_levelset==2){
	
		if(!rank) { // root processor reads the wave data
			FILE *fd;
			char filen[256];
			char path_wind[256];
			sprintf(path_wind, "./");
			PetscOptionsGetString(PETSC_NULL,"-path_wind", path_wind, 256, PETSC_NULL);		
			sprintf(filen, "%s/WAVE_wind%06d.dat", path_wind, WIND_ti);
			fd = fopen(filen, "r");
			if(!fd) {
				printf("\n******************* Cannot open %s ! *******************\n", path_wind);
			}
			else {
				printf("\n**Succesfully read wind file in %s/WAVE_wind%06d.dat **\n", path_wind, WIND_ti);
			}			
			fscanf(fd,"%le %d %d\n", &time, &NXMOD, &NZ);
			i_max= NXMOD*NZ;
			user[0].wave_inf[0].NZ=NZ;
			k=0;
			for(i=0;i<i_max;i++){
				fscanf(fd,"%le %le %le %le %le %le\n",&temp,&user[0].wave_inf[0].WIND_Y[i],&user[0].wave_inf[0].WIND_Z[i],&user[0].wave_inf[0].WIND_U[i],&user[0].wave_inf[0].WIND_V[i],&user[0].wave_inf[0].WIND_W[i]);
				user[0].wave_inf[0].WIND_Y[i]+=wave_wind_yshift;
			}		

			MPI_Bcast(&time, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NZ, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_Y[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_Z[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_U[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_V[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_W[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);			

			printf("\nReading WAVE_wind%06d.dat file:\n",WIND_ti);
			printf("time=%le, NXMOD=%d, NZ=%d\n", time, NXMOD, NZ);

			
			for(i=0;i<2;i++){
				printf("%le %le %le %le %le\n",user[0].wave_inf[0].WIND_Y[i],user[0].wave_inf[0].WIND_Z[i],user[0].wave_inf[0].WIND_U[i],user[0].wave_inf[0].WIND_V[i],user[0].wave_inf[0].WIND_W[i]);
			}
			printf("...     ...            i_max=%d  \n",i_max);
			printf("%le %le %le %le %le\n",user[0].wave_inf[0].WIND_Y[NXMOD-1],user[0].wave_inf[0].WIND_Z[NZ-1],user[0].wave_inf[0].WIND_U[i_max-1],user[0].wave_inf[0].WIND_V[i_max-1],user[0].wave_inf[0].WIND_W[i_max-1]);
			printf("Reading WAVE_wind%06d.dat file completed\n\n",WIND_ti);
			fclose(fd);
		}
		else if (rank) {
				
			MPI_Bcast(&time, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NZ, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&NXMOD, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			i_max= NXMOD*NZ;
			MPI_Bcast(&user[0].wave_inf[0].WIND_Y[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_Z[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_U[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_V[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);
			MPI_Bcast(&user[0].wave_inf[0].WIND_W[0], i_max, MPIU_REAL, 0, PETSC_COMM_WORLD);	
		}

		user[0].wave_inf[0].WIND_ti+=1;
		if(user[0].wave_inf[0].WIND_ti==wind_recicle)user[0].wave_inf[0].WIND_ti=0;//recicle
		user[0].wave_inf[0].WIND_time_since_read=0.;
		//double 	scale_vel=9.1439;
		//double  scale_L=1000.;
		for(i=0;i<i_max;i++){
			user[0].wave_inf[0].WIND_U[i]*=wave_wind_refvel;
			user[0].wave_inf[0].WIND_V[i]*=wave_wind_refvel;
			user[0].wave_inf[0].WIND_W[i]*=wave_wind_refvel;
			user[0].wave_inf[0].WIND_Y[i]*=wave_wind_reflength;
			user[0].wave_inf[0].WIND_Z[i]*=wave_wind_reflength;			
		}
		if(zs==zs){
			k=lzs;
			//map near field nodes into farfield nodes 
			for (i=lxs; i<lxe; i++)
			for (j=lys; j<lye; j++){
				z1=cent[k][j][i].z;
				y1=cent[k][j][i].y;
				x1=cent[k][j][i].x;
				
				//Use simetry out of X boundary
				if(x1>user[0].wave_inf[0].WIND_Y[NXMOD-1]){
					x1=2*user[0].wave_inf[0].WIND_Y[NXMOD-1]-x1;
				}				
				if(x1<user[0].wave_inf[0].WIND_Y[0]){
					x1=2*user[0].wave_inf[0].WIND_Y[0]-x1;
				}
				
				for (l=0;l<NXMOD-1;l++){
					if(x1>=user[0].wave_inf[0].WIND_Y[l] && x1<user[0].wave_inf[0].WIND_Y[l+1]){
						user[0].wave_inf[0].WIND_id_x[i-lxs][j-lys]=l;
						break;
					}
					if(l==(NXMOD-2))printf("ERROR: Point outside of the Farfield x domain!!!!!!!!\n");		
				}
				//printf("i:%d,j:%d,id:%d,x:%f,y:%f\n",i,j,l,x1,y1);

				if(y1<user[0].wave_inf[0].WIND_Z[l]){user[0].wave_inf[0].WIND_id_y[i-lxs][j-lys]=0;}
				else{		
					for (m=0;m<NZ-1;m++){					
						if(y1>=user[0].wave_inf[0].WIND_Z[l+m*NXMOD] && y1<user[0].wave_inf[0].WIND_Z[l+(m+1)*NXMOD]){
							user[0].wave_inf[0].WIND_id_y[i-lxs][j-lys]=m;
							break;
						}
						if(m==(NZ-2))printf("i,j,k:%d,%d   y1=%f   ERROR: Point outside of the Farfield vertical domain!!!!!!!!\n",i,j,y1);
					}
				}	
				//printf("i:%d,j:%d,idx:%d,idy:%d,x:%f,y:%f\n",i,j,user[0].wave_inf[0].WIND_id_x[i][j],user[0].wave_inf[0].WIND_id_y[i][j],x1,y1);
			}	
		}
		user[0].wave_inf[0].NXMOD=NXMOD;
		//printf("myid:%d, xs:%d, xe:%d, mx:%d, ys:%d, ye:%d, my:%d, zs:%d, ze:%d, mz:%d",rank,xs,xe,info.xm,ys,ye,info.ym,zs,ze,info.zm);
	}

	
	
	DAVecRestoreArray(fda, user->lCent, &cent);
}

void WAVE_Formfuction2(UserCtx *user)
{
	PetscReal ts,te,cput;
	PetscGetTime(&ts);
	int ti_wave;
	double time=user->dt*ti;	
	int WAVE_num_max;
	double temp_kk,temp_kk2;
	double twopi=2.*acos(-1.),onepi=acos(-1.);
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int ii,jj;
	int NZMOD,NXMOD;
	int WAVE_ti=user[0].wave_inf[0].WAVE_ti;
	if(!rank)printf("Forming RHS for WAVE sourve\n");	

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscReal ***level,***aj,***nvert;
	Cmpnts ***cent;
	UserCtx *user2=&user[1];
	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts ***WAVE_fp;

	DA		da = user->da, fda = user->fda;	
	
	DAGetLocalInfo(user->da, &info);
	mx = info.mx; my = info.my; mz = info.mz;
	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	//VecSet(user->WAVE_fp, 0);
	DAVecGetArray(user->fda, user->lCent, &cent);
	DAVecGetArray(user->da, user->lLevelset, &level);
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda, user->WAVE_fp, &WAVE_fp);
	int numprint=user[0].wave_inf[0].WAVE_ind_max;
	if (numprint>=40)numprint=40;
	
	for(jj=0;jj<numprint;jj++){
		ii=user[0].wave_inf[0].WAVE_ind[jj];
		if(!rank) printf("A:%f, Kz:%f, Kx:%f, Omega:%f, Theta:%f \n",user[0].wave_inf[0].WAVE_a[ii],user[0].wave_inf[0].WAVE_Kz[ii],user[0].wave_inf[0].WAVE_Kx[ii],user[0].wave_inf[0].WAVE_omega[ii],user[0].wave_inf[0].WAVE_theta[ii]);
	}
	double level_thick=2.0*dthick;
	double wave_source_distance_z=1.0;
	if(wave_momentum_source==1) {
		wave_source_distance_z=12.;//maximum wave length/2
		if(floating_turbine_case){
			wave_source_distance_z=200.;//maximum wave length/2
		}
		if(!rank) printf("WAVE Time since read:%f \n",user[0].wave_inf[0].WAVE_time_since_read);

	}
	if(wave_momentum_source==2) wave_source_distance_z=user[0].wave_inf[0].WAVE_src_dist[0]*1.2;	
	if(wave_momentum_source==3) wave_source_distance_z=user[0].wave_inf[0].WAVE_src_dist[0]*1.2;	
	if(wave_momentum_source==4) wave_source_distance_z=1.e5;	
	if(!rank) printf("WAVE Source distance:%f \n",wave_source_distance_z);
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		//variables
		double x = cent[k][j][i].x;	
		//double y = cent[k][j][i].y;				
		double z = cent[k][j][i].z;
		//source at x=0
		if (fabs(level[k][j][i])<(level_thick*1.2) && z>(wave_source_cent_z - wave_source_distance_z) && z<(wave_source_cent_z + wave_source_distance_z)){
			double dldc, dlde, dldz;
			double dl_dx, dl_dy, dl_dz;
			double ajc = aj[k][j][i];
			double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
			double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
			double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
			Compute_dscalar_center (i, j, k, mx, my, mz, level, nvert, &dldc, &dlde, &dldz );
			Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0,  zet1,  zet2,  ajc, dldc, dlde, dldz, &dl_dx, &dl_dy, &dl_dz);					
			double grad_= sqrt(dl_dx*dl_dx+dl_dy*dl_dy+dl_dz*dl_dz);
			dl_dx=dl_dx/grad_;
			dl_dy=dl_dy/grad_;
			dl_dz=dl_dz/grad_;

			double epscos_1=dl_dx*epscos(level[k][j][i],level_thick);
			double epscos_2=dl_dy*epscos(level[k][j][i],level_thick);
			double epscos_3=dl_dz*epscos(level[k][j][i],level_thick);				

			int ii,jj;
			double temp_11=0.0;
			for(jj=0;jj<user[0].wave_inf[0].WAVE_ind_max;jj++){
				ii=user[0].wave_inf[0].WAVE_ind[jj];
				if(wave_momentum_source==1)	{
					temp_11+=user[0].wave_inf[0].WAVE_P_0[ii]*sin(user[0].wave_inf[0].WAVE_omega[ii]*user[0].wave_inf[0].WAVE_time_since_read-user[0].wave_inf[0].WAVE_Kx[ii]*cent[k][j][i].x-user[0].wave_inf[0].WAVE_theta[ii])*epscos((cent[k][j][i].z - wave_source_cent_z),user[0].wave_inf[0].WAVE_src_dist[ii]);
				}
				else if(wave_momentum_source==2) {
					temp_11+=user[0].wave_inf[0].WAVE_P_0[ii]*sin(user[0].wave_inf[0].WAVE_omega[ii]*time-user[0].wave_inf[0].WAVE_Kx[ii]*cent[k][j][i].x-user[0].wave_inf[0].WAVE_theta[ii])*epscos((cent[k][j][i].z - wave_source_cent_z),user[0].wave_inf[0].WAVE_src_dist[ii]);
				}
				else if(wave_momentum_source==3 || wave_momentum_source==4 ) {
					temp_11+=	user[0].wave_inf[0].WAVE_P_0[ii]*sin(user[0].wave_inf[0].WAVE_Kz[ii]*(cent[k][j][i].z - wave_source_cent_z)+user[0].wave_inf[0].WAVE_Kx[ii]*cent[k][j][i].x-user[0].wave_inf[0].WAVE_theta[ii]);
					temp_11+=	user[0].wave_inf[0].WAVE_P_1[ii]*cos(user[0].wave_inf[0].WAVE_Kz[ii]*(cent[k][j][i].z - wave_source_cent_z)+user[0].wave_inf[0].WAVE_Kx[ii]*cent[k][j][i].x-user[0].wave_inf[0].WAVE_theta[ii]);
				}
			}
			WAVE_fp[k][j][i].x=temp_11*epscos_1;
			WAVE_fp[k][j][i].y=temp_11*epscos_2;
			WAVE_fp[k][j][i].z=temp_11*epscos_3;
		}
		else{
			WAVE_fp[k][j][i].x=.0;
			WAVE_fp[k][j][i].y=.0;
			WAVE_fp[k][j][i].z=.0;				
		}
	}	
	
	DAVecRestoreArray(fda, user->WAVE_fp, &WAVE_fp);
	DAGlobalToLocalBegin(user->fda, user->WAVE_fp, INSERT_VALUES, user->lWAVE_fp);
	DAGlobalToLocalEnd(user->fda, user->WAVE_fp, INSERT_VALUES, user->lWAVE_fp);
	
	DAVecRestoreArray(user->fda, user->lCent, &cent);
	DAVecRestoreArray(user->da, user->lLevelset, &level);
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);	
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lNvert, &nvert);

	PetscGetTime(&te);
	cput=te-ts;
	double cput2;
	PetscGlobalMax(&cput,&cput2,PETSC_COMM_WORLD);
	if(!rank)printf("END Forming RHS for WAVE sourve\n");	
	if(!rank)printf("Time elapsed for forming RHS for wave source: %.2e(s)\n\n", cput2);
};

void WAVE_SL_Formfuction2(UserCtx *user)
{
	/*
	*This Function forms the RHS term in the momentum equation corresponding to the sponge layer 
	*method as descrived in Choi and Yoon 2009 (Coastal Engineering 56)
	*It has the following options that can be controled from the file "control.dat":
		*-wave_sponge_layer 0   In this case the fuction is not in use
		*-wave_sponge_layer 1   Sponge layers used at the z boundaries
		*-wave_sponge_layer 2	  Sponge layers used at both the z and x boundaries
		*-wave_sponge_zs 3.0 	  Lenght of the two sponge layers of the z boundaries
		*-wave_sponge_z01 -10.	  Starting z coordinate of the fisrt sponge layer of the z boundary
		*-wave_sponge_z02 12.0	  Starting z coordinate of the second sponge layer of the z boundary
		*-wave_sponge_xs 3.0		  Lenght of the two sponge layers of the x boundaries
		*-wave_sponge_x01 -10.0	Starting x coordinate of the fisrt sponge layer of the x boundary
		*-wave_sponge_x02 7.0		Starting x coordinateof the second sponge layer of the x boundary
	*Note that wave_momentum_source must be active. Otherwise WAVE_fp won't be reinitialized to zero 
	*and its value will incorrecly increase at every time step. 
	*/
	
	PetscReal ts,te,cput;
	PetscGetTime(&ts);
	int ti_wave;
	double time;
	int WAVE_num_max;
	double temp_kk,temp_kk2;
	double twopi=2.*acos(-1.),onepi=acos(-1.);
	double eps, eps_number=2.00001,ff;
	double ro_w=1000, ro_air=1.2;
	//sponge layer
	double C1=200.;
	double C2=1.*ro_w;
	double ns=10.;		
	double wave_source_heigh=0.000;

	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if(!rank)printf("Forming RHS for WAVE Sponge Layer\n");	

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	PetscReal ***level,***aj,***nvert;
	Cmpnts ***cent;
	UserCtx *user2=&user[1];
	Cmpnts	***csi, ***eta, ***zet;
	Cmpnts	***ucat;
	Cmpnts ***WAVE_fp;
	PetscReal ***rho;

	DA		da = user->da, fda = user->fda;	
	DAGetLocalInfo(user->da, &info);
	mx = info.mx; my = info.my; mz = info.mz;
	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	DAVecGetArray(user->fda, user->lCent, &cent);
	DAVecGetArray(user->da, user->lLevelset, &level);
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(fda, user->WAVE_fp, &WAVE_fp);
	DAVecGetArray(da, user->lDensity, &rho);
	DAVecGetArray(fda, user->Ucat,  &ucat);

	double level_thick=2.0*dthick; // Ali added May 1, 2016

        double min_x=0., min_y=0., min_z=0.;

        if(solitary_wave){  //reinitialize the WAVE_fp to zero for solitary waves
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
			WAVE_fp[k][j][i].x = 0.;
			WAVE_fp[k][j][i].y = 0.;
			WAVE_fp[k][j][i].z = 0.;
                                }
                         }

	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {

		//cell center coordinates
		double z = cent[k][j][i].z;
		double y = cent[k][j][i].y;
		double x = cent[k][j][i].x;
	
		if (wave_sponge_layer==1 || wave_sponge_layer==2) {
		//xmin boundary sponge layer
		if (fabs(level[k][j][i])<(level_thick*1.2) && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) { // Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;		
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((wave_sponge_zs-(z-wave_sponge_z01))/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}		
		//xmax boundary sponge layer
		if (fabs(level[k][j][i])<(level_thick*1.2) && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) { //Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((z-wave_sponge_z02)/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}
                }
		if (wave_sponge_layer==2) {
			//ymin boundary sponge layer
		        if (fabs(level[k][j][i])<(level_thick*1.2) && x<wave_sponge_x01+wave_sponge_xs && x>wave_sponge_x01) { // Ali added May 1, 2016
			//if (y<wave_source_heigh && x<wave_sponge_x01+wave_sponge_xs && x>wave_sponge_x01) {  //Ale deactivated May 1, 2016
				double u=ucat[k][j][i].x;  
				double v=ucat[k][j][i].y;    
				double w=ucat[k][j][i].z; 
				if (ti ==0) {
					w=0.,v=0.,u=0.;
				}
				double coef_1=(exp(powf((wave_sponge_xs-(x-wave_sponge_x01))/(wave_sponge_xs),ns))-1.0)/(exp(1.0)-1.0);
				WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
				WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
				WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
			}		
			//ymax boundary sponge layer
		        if (fabs(level[k][j][i])<(level_thick*1.2) && x<wave_sponge_x02+wave_sponge_xs && x>wave_sponge_x02) { // Ali added May 1, 2016
			//if (y<wave_source_heigh && x<wave_sponge_x02+wave_sponge_xs && x>wave_sponge_x02) { // Ali deactivated May 1, 2016
				double u=ucat[k][j][i].x;  
				double v=ucat[k][j][i].y;    
				double w=ucat[k][j][i].z;
				if (ti ==0) {
					w=0.,v=0.,u=0.;
				}
				double coef_1=(exp(pow((x-wave_sponge_x02)/(wave_sponge_xs),ns))-1.0)/(exp(1.0)-1.0);
				//WAVE_fp stores the value at the cell centers
				WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
				WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
				WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;					
			}										
		}																	
		if (wave_sponge_layer==3) { //Sponge in the water and air
		//xmin boundary sponge layer
		if (z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) { // Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;		
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((wave_sponge_zs-(z-wave_sponge_z01))/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}		
		//xmax boundary sponge layer
		if (z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) { //Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((z-wave_sponge_z02)/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}
                }
		if (wave_sponge_layer==4) { //Sponge in the water
		//xmin boundary sponge layer
		if (level[k][j][i] >= 0. && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) { // Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;		
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((wave_sponge_zs-(z-wave_sponge_z01))/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}		
		//xmax boundary sponge layer
		if (level[k][j][i] >= 0. && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) { //Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((z-wave_sponge_z02)/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}
                }
		if (wave_sponge_layer==5) { //Sponge in the water depth 
		//xmin boundary sponge layer
		if (level[k][j][i] > (-dthick*1.2) && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z01+wave_sponge_zs && z>wave_sponge_z01) { // Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;		
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((wave_sponge_zs-(z-wave_sponge_z01))/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}		
		//xmax boundary sponge layer
		if (level[k][j][i] > (-dthick*1.2) && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) {  // Ali added May 1, 2016
		//if (y<wave_source_heigh && z<wave_sponge_z02+wave_sponge_zs && z>wave_sponge_z02) { //Ali deactivated May 1, 2016
			double u=ucat[k][j][i].x;  
			double v=ucat[k][j][i].y;    
			double w=ucat[k][j][i].z;
			if (ti ==0) {
				w=0.,v=0.,u=0.;
			}
			double coef_1=(exp(powf((z-wave_sponge_z02)/(wave_sponge_zs),ns))-1.0)/(exp(1.0)-1.0);
			WAVE_fp[k][j][i].x-=(C1*u+C2*u*fabs(u))*coef_1;
			WAVE_fp[k][j][i].y-=(C1*v+C2*v*fabs(v))*coef_1;
			WAVE_fp[k][j][i].z-=(C1*w+C2*w*fabs(w))*coef_1;
		}
                }
        min_x=PetscMax(min_x,WAVE_fp[k][j][i].x);
        min_y=PetscMax(min_y,WAVE_fp[k][j][i].y);
        min_z=PetscMax(min_z,WAVE_fp[k][j][i].z);
	}
	DAVecRestoreArray(fda, user->WAVE_fp, &WAVE_fp);
	DAGlobalToLocalBegin(user->fda, user->WAVE_fp, INSERT_VALUES, user->lWAVE_fp);
	DAGlobalToLocalEnd(user->fda, user->WAVE_fp, INSERT_VALUES, user->lWAVE_fp);
	DAVecRestoreArray(user->fda, user->lCent, &cent);
	DAVecRestoreArray(user->da, user->lLevelset, &level);
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);	
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lDensity, &rho);
	DAVecRestoreArray(fda, user->Ucat,  &ucat);
	PetscGetTime(&te);
	cput=te-ts;
	double cput2;
	PetscGlobalMax(&cput,&cput2,PETSC_COMM_WORLD);
	if(!rank)printf("END Forming RHS for WAVE Sponge Layer sourve\n");	
	if(!rank)printf("Time elapsed for forming RHS for wave Sponge Layer source: %.2e(s)\n\n", cput2);
	if(!rank)printf("min_x: %e, min_y: %e, min_z: %e\n", min_x, min_y, min_z);
};

PetscErrorCode Solitary_wave_inlet_velocity_profile_Boussinesq(UserCtx *user, double *uin, double *vin, double *win, double x, double y, double z, double win_inlet_flux, int jj, double time, int ti)
{
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        double h = inlet_y - inlet_bed_elevation; 
        double y0 = y - inlet_bed_elevation; // y0 is measured from any elevation and it is the same as coor[][][].y in the code. The y0 is measured from the bed elevation at the inlet vertically upward
        PetscReal H = solitary_wave_amplitude; 
        double g = 9.81; // gravitational acceleration
        double C = sqrt(g*(h+H));  //celcerity of the solitary wave
        //double inlet_z_for_solitary_wave = z;  //inlet z where the wave enters. It must be equal to the z that comes from bcs.c: Check this out
       
     
	//if(!rank)printf("z comming from bcs.c for inlet z: %le(s)\n\n", z);

        int tisw = ti - ti_start_solitary_wave;
        double time_sw = user->dt * (double)tisw;
           double eps_Hh = H / h;
           double eps_yh = y0 / h;
           double Lambda_sw = 4.12 * h / sqrt(eps_Hh); 
           double z_shift = Lambda_sw * (-1.0) + inlet_z_for_solitary_wave;
           //double z_shift = Lambda_sw * (-1.5);
           double solitary_wave_period = Lambda_sw / C;

           double Eta_sw = H * pow((1. / cosh( sqrt(3. * eps_Hh / 4.) * (z - z_shift - C * time_sw) / h)), 2.);
          double eps_eH = Eta_sw / H;
          double eps_eh = Eta_sw / h;
          double W_sw = C * ((eps_Hh + 3. * eps_Hh * eps_Hh * (1. / 6. - (1. / 2.) * eps_yh * eps_yh ) ) * eps_eH - eps_Hh * eps_Hh * (7. / 9. - (9. / 4.) * eps_yh * eps_yh) * eps_eH * eps_eH);
          double V_sw = C * sqrt(3. * eps_Hh) * eps_yh * eps_eh * tanh (sqrt(3. * eps_Hh / 4.) * (z - z_shift - C * time_sw)) * (1. + 0.5 * eps_Hh * (1. - 7.* eps_eH - eps_yh * eps_yh * (1. - 3. * eps_eH)));
           

           //Grimshaw 3rd order Solitary Wave for H>           
           if(eps_Hh >= 0.25) {
	     if(!rank && ti==ti_start_solitary_wave)printf("Grimshaw 3rd-order Solitary wave\n");
             double e1 = eps_Hh;
             double e2 = eps_Hh * eps_Hh;
             double e3 = eps_Hh * eps_Hh * eps_Hh;
             double e4 = eps_Hh * eps_Hh * eps_Hh * eps_Hh;
             double e6 = eps_Hh * eps_Hh * eps_Hh * eps_Hh * eps_Hh * eps_Hh;

             double a = sqrt(3. * e1 /4.) * (1. - 5. * e1 / 8. + 71. * e2 / 128.);
             C = sqrt(g * h) * sqrt (1. + e2 - 1. * e2 / 20. - 3. * e3 / 70.); 
             double X = z - z_shift - C * time_sw;
             
             double s1 = 1. / cosh (a * X / h);              
             double s2 = s1 * s1;
             double s3 = s2 * s1;
             double s4 = s3 * s1;
             double s5 = s4 * s1;
             double s6 = s5 * s1;

             double q1 = tanh (a * X / h);              
             double q2 = q1 * q1;
              
             W_sw = sqrt(g*h) * (e1 * s2 - e2 * (-s2/4. + s4 + (y0/h)*(y0/h)*(3. * s2/2. - 9. *s4/4.)) - e3 * (19. * s2/40. + s4/5. - 6.*s6/5. + (y0/h)*(y0/h) * (-3. *s2/2.-15.*s4/4.+15.*s6/2.)+
                    (y0/h)*(y0/h)*(y0/h)*(y0/h)*(-3. * s2/8. + 45.*s4/16. - 45.*s6/16.)));
             
             V_sw = sqrt(g*h) *  (sqrt(3. *e1) *(y0/h) * q1 * (-e1*s2 + e2 *(3.*s2/8. + 2.*s4 +(y0/h)*(y0/h)*(s2/2. -3.*s4/2.)) + e3*(49.*s2/640. -17. *s4/29. -18.*s6/5. + (y0/h)*(y0/h)*(-13.*s2/6.                    -25.*s4/16. + 15.*s6/2.) + (y0/h)*(y0/h)*(y0/h)*(y0/h)*(-3.*s2/40. + 9.*s4/8.-27.*s6/16.))));

           } else {if(!rank && ti==ti_start_solitary_wave)printf("Boussinesq Solitary wave\n");}

             
           *uin = 0; 
           *vin = V_sw;
           *win = W_sw + win_inlet_flux;
            
return ( 0 );
}

PetscErrorCode Solitary_wave_inlet_elevation_profile_Boussinesq(UserCtx *user, double *free_surface_elevation, double *free_surface_elevation_to, double x, double y, double z, double z_to, double time, int ti)
{
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        double h = inlet_y - inlet_bed_elevation; 
        double y0 = y - inlet_bed_elevation; // y0 is measured from any elevation and it is the same as coor[][][].y in the code. The y0 is measured from the bed elevation at the inlet vertically upward
        PetscReal H = solitary_wave_amplitude; 
        double g = 9.81; // gravitational acceleration
        double C = sqrt(g*(h+H));  //celcerity of the solitary wave
        //double inlet_z_for_solitary_wave = z;  //inlet z where the wave enters. It must be equal to the z that comes from bcs.c: Check this out
        //double inlet_z_for_solitary_wave_to = -5.0;//z_to;  //inlet z where the wave enters. It must be equal to the z that comes from bcs.c: Check this out
            
	//if(!rank)printf("z comming from level.c for inlet z: %le(s)\n\n", z);
        int tisw = ti - ti_start_solitary_wave;
        double time_sw = user->dt * (double)tisw;

           double eps_Hh = H / h;
           double eps_yh = y0 / h;
           double Lambda_sw = 4.12 * h / sqrt(eps_Hh); 
           double z_shift = Lambda_sw * (-1.0) + inlet_z_for_solitary_wave;
           double solitary_wave_period = Lambda_sw / C;

           double Eta_sw = H * pow((1. / cosh( sqrt(3. * eps_Hh / 4.) * (z - z_shift - C * time_sw) / h)),2.);
           double Eta_sw_to = H * pow((1. / cosh( sqrt(3. * eps_Hh / 4.) * (z_to - z_shift - C * time_sw) / h)),2.);
           
           //Grimshaw 3rd order Solitary Wave for H>           
           if(eps_Hh >= 0.25) {
	     if(!rank && ti==ti_start_solitary_wave)printf("Grimshaw 3rd-order Solitary wave\n");
             double e1 = eps_Hh;
             double e2 = eps_Hh * eps_Hh;
             double e3 = eps_Hh * eps_Hh * eps_Hh;
             double e4 = eps_Hh * eps_Hh * eps_Hh * eps_Hh;
             double e6 = eps_Hh * eps_Hh * eps_Hh * eps_Hh * eps_Hh * eps_Hh;

             double a = sqrt(3. * e1 /4.) * (1. - 5. * e1 / 8. + 71. * e2 / 128.);
             C = sqrt(g * h) * sqrt (1. + e2 - 1. * e2 / 20. - 3. * e3 / 70.); 
             double X = z - z_shift - C * time_sw;
             double X_to = z_to - z_shift - C * time_sw;
             
             double s1 = 1. / cosh (a * X / h);              
             double s2 = s1 * s1;
             double s3 = s2 * s1;
             double s4 = s3 * s1;
             double s5 = s4 * s1;
             double s6 = s5 * s1;

             double s1_to = 1. / cosh (a * X_to / h);              
             double s2_to = s1_to * s1_to;
             double s3_to = s2_to * s1_to;
             double s4_to = s3_to * s1_to;
             double s5_to = s4_to * s1_to;
             double s6_to = s5_to * s1_to;

             double q1 = tanh (a * X / h);              
             double q2 = q1 * q1;
              
             double q1_to = tanh (a * X_to / h);              
             double q2_to = q1_to * q1_to;

             Eta_sw = h * (e1 * s2 - 3. * e2 * s2 * q2 /4. + e3 * (5. * s2 * q2 /8. - 101. * s4 * q2 / 58.));
             Eta_sw_to = h * (e1 * s2_to - 3. * e2 * s2_to * q2_to /4. + e3 * (5. * s2_to * q2_to /8. - 101. *                                 s4_to * q2_to / 58.));
              } else {if(!rank && ti==ti_start_solitary_wave) printf("Boussinesq Solitary wave\n");}

           *free_surface_elevation = inlet_y + Eta_sw;
           *free_surface_elevation_to = inlet_y + Eta_sw_to;
            
return ( 0 );
}


PetscErrorCode Linear_wave_single_inlet_velocity_profile (UserCtx *user, double *uin, double *vin, double *win, double x, double y, double z, double win_inlet_flux, int jj, double time, int ti)
{
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        double PI = acos(-1);
        double h = inlet_y - inlet_bed_elevation; 
        double y0 = y - inlet_bed_elevation; // y0 is measured from any elevation and it is the same as coor[][][].y in the code. The y0 is measured from the bed elevation at the inlet vertically upward
        double Y = y - inlet_y; // Y is the y from still water surface before wave comes in. Y is measured from still-initial-water-level. 
                                // Y for points along vertical that are below below original flat free-surface (still initial-water-level) is negative and above it positive.  
        PetscReal H = linear_wave_single_amplitude; 
        double g = 9.81; // gravitational acceleration
            
	//if(!rank)printf("z comming from level.c for inlet z: %le(s)\n\n", z);
        int ti_lws = ti - ti_start_linear_wave_single;
        double time_lws = user->dt * (double)ti_lws;

           double Lambda_lws = 2.* PI / linear_wave_single_number; 
           double z_shift = acos(0)/linear_wave_single_number - inlet_z_for_linear_wave_single;
           double omega = sqrt( g * linear_wave_single_number * tanh (linear_wave_single_number * h));
           double linear_wave_single_period = 2. * PI / omega;
           double linear_wave_single_speed = Lambda_lws / linear_wave_single_period;

           double Eta_lws    = H * cos(linear_wave_single_number * (z + z_shift) - omega * time_lws);
           //double Eta_lws_to = H * cos(linear_wave_single_number * (z_to - z_shift) - omega * time_lws);


            double W_lws = ((H*omega) * cosh (linear_wave_single_number * (h + Y))/sinh (linear_wave_single_number * h)) * cos(linear_wave_single_number * (z + z_shift) - omega * time_lws);
            double V_lws = ((H*omega) * sinh (linear_wave_single_number * (h + Y))/sinh (linear_wave_single_number * h)) * sin(linear_wave_single_number * (z + z_shift) - omega * time_lws);
            //double W_lws = ((H * g * linear_wave_single_number / omega) * cosh (linear_wave_single_number * (h + Y))/sinh (linear_wave_single_number * h)) * cos(linear_wave_single_number * (z + z_shift) - omega * time_lws);
            //double V_lws = ((H * g * linear_wave_single_number / omega) * sinh (linear_wave_single_number * (h + Y))/sinh (linear_wave_single_number * h)) * sin(linear_wave_single_number * (z + z_shift) - omega * time_lws);
           
           *uin = 0; 
           *vin = V_lws;
           *win = W_lws + win_inlet_flux;
            
return ( 0 );
}


PetscErrorCode Linear_wave_single_inlet_elevation_profile (UserCtx *user, double *free_surface_elevation, double *free_surface_elevation_to, double x, double y, double z, double z_to, double time, int ti)
{
	PetscInt	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        double PI = acos(-1);
        double h = inlet_y - inlet_bed_elevation; 
        double y0 = y - inlet_bed_elevation; // y0 is measured from any elevation and it is the same as coor[][][].y in the code. The y0 is measured from the bed elevation at the inlet vertically upward
        PetscReal H = linear_wave_single_amplitude; 
        double g = 9.81; // gravitational acceleration
            
		//if(!rank)printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n");
		//if(!rank)printf("z comming from level.c for inlet z: %le and z_to from bcs.c: %le \n\n", z, z_to);
        int ti_lws = ti - ti_start_linear_wave_single;
        double time_lws = user->dt * (double)ti_lws;

           double Lambda_lws = 2.* PI / linear_wave_single_number; 
           double z_shift = acos(0)/linear_wave_single_number - inlet_z_for_linear_wave_single;
           double omega = sqrt( g * linear_wave_single_number * tanh (linear_wave_single_number * h));
           double linear_wave_single_period = 2. * PI / omega;
           double linear_wave_single_speed = Lambda_lws / linear_wave_single_period;

           double Eta_lws    = H * cos(linear_wave_single_number * (z + z_shift) - omega * time_lws);
           double Eta_lws_to = H * cos(linear_wave_single_number * (z_to + z_shift) - omega * time_lws);
           
           *free_surface_elevation = inlet_y + Eta_lws;
           *free_surface_elevation_to = inlet_y + Eta_lws_to;

return ( 0 );
}
