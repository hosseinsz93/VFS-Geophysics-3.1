#include "variables.h"
#include <algorithm>
#include <math.h>

using namespace std;

static double R0;
static double rho_cfd;
static double T_wt;
static double T_cfd;

static double ratio_rho;
static double ratio_L;
static double ratio_T;
static double ratio_V;

static double region;

const int nYw=50;
double sNx[nYw];
double sNy[nYw];
double sNz[nYw];

//static double ang_axis_old;

PetscErrorCode IndividualPitchControl(IBMNodes *ibm, FSInfo *fsi,PetscInt ibi, PetscReal dt_wt);
PetscErrorCode ReadYawInput(PetscInt ti, PetscReal& nx, PetscReal& ny, PetscReal& nz);

PetscErrorCode TurbineTorqueControl_Output(FSInfo *fsi, IBMNodes *ibm)
{

        int rank, i, ibi;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//        PetscBarrier(PETSC_NULL);
        if (!rank) {
                for (ibi=0; ibi<NumberOfTurbines; ibi++) {
                        FILE *f;
                        char filen[80];
                        sprintf(filen, "TurbineTorqueControl%6.6d_%3.3d.dat", (int)ti, ibi);
                        f = fopen(filen, "w");
                        PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le \n", fsi[ibi].Torque_aero, fsi[ibi].Torque_generator, fsi[ibi].ang_axis, fsi[ibi].angvel_axis, ibm[ibi].Tipspeedratio, ibm[ibi].pitch[0], ibm[ibi].pitch[1], ibm[ibi].pitch[2], fsi[ibi].angvel_axis_err, fsi[ibi].angvel_axis_err_integral, fsi[ibi].betacos_IPC, fsi[ibi].betasin_IPC);
                        //PetscFPrintf(PETSC_COMM_WORLD, f, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e \n", fsi[ibi].Torque_aero, fsi[ibi].Torque_generator, fsi[ibi].ang_axis_old, fsi[ibi].angvel_axis, ibm[ibi].Tipspeedratio, ibm[ibi].pitch[0], ibm[ibi].pitch[1], ibm[ibi].pitch[2], fsi[ibi].angvel_axis_err, fsi[ibi].angvel_axis_err_integral, fsi[ibi].betacos_IPC, fsi[ibi].betasin_IPC, fsi[ibi].K_Gen, fsi[ibi].GeneratorTransSpeed,fsi[ibi].Alpha_OT);
                        fclose(f);
                }
        }
        return(0);
}

PetscErrorCode TurbineTorqueControl_Input(FSInfo *fsi, IBMNodes *ibm)
{
        int  rank, i, ibi;
//        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        for (ibi=0; ibi<NumberOfTurbines; ibi++) {
                        FILE *f;
                        char filen[80];
                        double angvel, TSR;
                        sprintf(filen, "TurbineTorqueControl%6.6d_%3.3d.dat", (int)ti, ibi);
                        f = fopen(filen, "r");
                        if (f)
                        {
                                PetscPrintf(PETSC_COMM_WORLD, "Read %s for turbine restart file!! \n", filen);
 //                               fscanf( f, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e \n", fsi[ibi].Torque_aero, fsi[ibi].Torque_generator, fsi[ibi].ang_axis, fsi[ibi].angvel_axis, ibm[ibi].Tipspeedratio, ibm[ibi].pitch[0], ibm[ibi].pitch[1], ibm[ibi].pitch[2], fsi[ibi].angvel_axis_err, fsi[ibi].angvel_axis_err_integral, fsi[ibi].betacos_IPC, fsi[ibi].betasin_IPC, fsi[ibi].K_Gen, fsi[ibi].GeneratorTransSpeed,fsi[ibi].Alpha_OT);
 //                               fscanf(f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n", &(fsi[ibi].Torque_aero), &(fsi[ibi].Torque_generator), &(fsi[ibi].ang_axis), &(angvel), &(TSR), &(ibm[ibi].pitch[0]), &(ibm[ibi].pitch[1]), &(ibm[ibi].pitch[2]), &(fsi[ibi].angvel_axis_err), &(fsi[ibi].angvel_axis_err_integral), &(fsi[ibi].betacos_IPC), &(fsi[ibi].betasin_IPC), &(fsi[ibi].K_Gen), &(fsi[ibi].GeneratorTransSpeed),&(fsi[ibi].Alpha_OT));
                                  fscanf(f, "%le %le %le %le %le %le %le %le %le %le %le %le \n", &(fsi[ibi].Torque_aero), &(fsi[ibi].Torque_generator), &(fsi[ibi].ang_axis), &(angvel), &(TSR), &(ibm[ibi].pitch[0]), &(ibm[ibi].pitch[1]), &(ibm[ibi].pitch[2]), &(fsi[ibi].angvel_axis_err), &(fsi[ibi].angvel_axis_err_integral), &(fsi[ibi].betacos_IPC), &(fsi[ibi].betasin_IPC));                
                        } else
                                PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s !! \n", filen);

                        fclose(f);
                        if (!FixTipSpeedRatio) ibm[ibi].Tipspeedratio=TSR;

                        if (fixturbineangvel) fsi[ibi].angvel_axis=fsi[ibi].angvel_fixed;
                        else fsi[ibi].angvel_axis=angvel;

        }


  return(0);
}


PetscErrorCode TurbineTorqueControlInitialization(FSInfo *fsi, IBMNodes *ibm)
{

	PetscInt ibi;
        R0=fsi[0].r_rotor/reflength_wt;
        rho_cfd=1.0;                     // no levelset
        T_wt=fsi[0].r_rotor/refvel_wt;   // the 0 turbine cannot be in the other turbines' wakes.
        T_cfd=R0/refvel_cfd;             // the 0 turbine cannot be in the other turbines' wakes.

        ratio_rho=rho_air/rho_cfd;
        ratio_L=reflength_wt;
        ratio_T=T_wt/T_cfd;
        ratio_V=refvel_wt/refvel_cfd;
 

        //PetscBarrier(PETSC_NULL);
        PetscPrintf(PETSC_COMM_WORLD, "###################################\n");
        PetscPrintf(PETSC_COMM_WORLD, "  Turbine Control Initialization\n");
        PetscPrintf(PETSC_COMM_WORLD, "###################################\n");
        //PetscBarrier(PETSC_NULL);

        if (rstart_turbinerotation) {
        	PetscPrintf(PETSC_COMM_WORLD, "Read turbine rotation restart file \n");
        	TurbineTorqueControl_Input(fsi, ibm);

	        for (ibi=0; ibi<NumberOfTurbines; ibi++) {
			PetscPrintf(PETSC_COMM_WORLD, "angvel %f ang %f \n", fsi[ibi].angvel_axis, fsi[ibi].ang_axis);
		}

	} 
	else
	{
        	PetscPrintf(PETSC_COMM_WORLD, "Initializing Turbine \n");
        	for (ibi=0; ibi<NumberOfTurbines; ibi++)
        	{
        	        fsi[ibi].angvel_axis_err_integral = 0.0;
        	        fsi[ibi].angvel_axis_err = 0.0;

        	        fsi[ibi].betacos_IPC = 0.0;
        	        fsi[ibi].betasin_IPC = 0.0;
                        fsi[ibi].Torque_aero = 0.0;
                        fsi[ibi].Torque_generator = 0.0;
        	        // assuming three blades
        	        ibm[ibi].pitch_IPC[0] = 0.0;
        	        ibm[ibi].pitch_IPC[1] = 0.0;
        	        ibm[ibi].pitch_IPC[2] = 0.0;

        	        ibm[ibi].pitch[0] = 0.0;
        	        ibm[ibi].pitch[1] = 0.0;
        	        ibm[ibi].pitch[2] = 0.0;
        	}

        	for (ibi=0; ibi<NumberOfTurbines; ibi++) {
        	        fsi[ibi].ang_axis=0.0;
                  fsi[ibi].ang_axis_old=0.0;
        	        double R=fsi[ibi].r_rotor/reflength_wt;
        	        double TSR;
        	        if (FixTipSpeedRatio) {
        	                TSR=ibm[ibi].Tipspeedratio;
        	                fsi[ibi].angvel_axis=TSR*ibm[ibi].U_ref/R;
                          PetscPrintf(PETSC_COMM_WORLD, "FixTSR= %f Uref= %f R= %f angl_vel= %f\n", TSR, ibm[ibi].U_ref, R,  fsi[ibi].angvel_axis);
        	        } else if (turbinetorquecontrol) {
        	                double angvel_desired = fsi[ibi].GeneratorSpeed_desired/fsi[ibi].GearBoxRatio;
//Dviding by 3 of optimal to start at a lowe speed
        	                //fsi[ibi].angvel_axis=angvel_desired * ratio_T/3.0; //Temporarily changed to 3.4.18 code
                            fsi[ibi].angvel_axis=angvel_desired * ratio_T;
//        	                PetscPrintf(PETSC_COMM_WORLD, "Nondim AngvelAxis %f \n", fsi[ibi].angvel_axis);
//        	                PetscPrintf(PETSC_COMM_WORLD, "RatioT %f \n", ratio_T);

                          if(fsi[ibi].TPCntrl==2){

                           if(fsi[ibi].K_Gen<=0){
                            fsi[ibi].K_Gen = fabs(0.5*M_PI*ratio_rho*pow(fsi[ibi].r_rotor,5)*fsi[ibi].CP_max/\
                                                        pow(fsi[ibi].TSR_max,3)/pow(fsi[ibi].GearBoxRatio,3));
                           }                      

                            fsi[ibi].Alpha_OT = fabs(fsi[ibi].Torque_generator_max/\
                                                (fabs(fsi[ibi].GeneratorSpeed_desired)-fabs(fsi[ibi].GeneratorSyncSpeed)));

                            fsi[ibi].GeneratorTransSpeed = (fsi[ibi].Alpha_OT -\
                                   sqrt(fsi[ibi].Alpha_OT*(fsi[ibi].Alpha_OT-4.*fsi[ibi].K_Gen*fabs(fsi[ibi].GeneratorSyncSpeed))))/\
                                            (2.*fsi[ibi].K_Gen);
        	               PetscPrintf(PETSC_COMM_WORLD, "K_Generator %f %f %f %f %f %f\n", fsi[ibi].K_Gen,fsi[ibi].TSR_max,fsi[ibi].r_rotor,fsi[ibi].CP_max,fsi[ibi].GearBoxRatio,ratio_rho);
                          }

        	        } else if (fixturbineangvel) {
        	                fsi[ibi].angvel_axis=fsi[ibi].angvel_fixed;
        	                ibm[ibi].Tipspeedratio = fsi[ibi].angvel_axis*R/ibm[ibi].U_ref;
        	        }
        	        PetscPrintf(PETSC_COMM_WORLD, "fixturbineangvel = %i and angvel of turbne at first time step %f \n", fixturbineangvel, fsi[ibi].angvel_axis);
        	}
	}
     
	if ( (ti) % tiout == 0) TurbineTorqueControl_Output(fsi, ibm);

	return(0);
}


PetscErrorCode Calc_GenTorque_Pitch_C96(IBMNodes *ibm, FSInfo *fsi, PetscReal dt, PetscReal& TorqueOut, PetscInt ibi)
{
        // Clipper Liberty C96 Turbine control system
        // The kt and kp constants are particular of this turbine 
        // Must consider the units of Kp [degrees/(rad/s)] and Kt [N m /(rad/s)]
        
        const double rpm2radps = M_PI/30.;
        const double kp_prop = -0.031053366009632  / rpm2radps; //K_pitch,proportional
        const double kp_int  = -0.003879670183856  / rpm2radps;  //K_pitch,integral
        const double kp_ratio = -kp_prop/kp_int;  
        const double kt_prop = -88.471633656702878 / rpm2radps; //K_torque,proportional
        const double kt_int  = -11.053254552452582 / rpm2radps; //K_torque,integral
        const double pitch_max = 180;

        
        double sign = fsi[ibi].GeneratorSpeed_desired/fabs(fsi[ibi].GeneratorSpeed_desired);
        // Computing the absolute to make sure they are all positive regardless of the sign in the input file
        double angvel_wt = fabs( fsi[ibi].angvel_axis / ratio_T );
        double angvel_rated = fabs( fsi[ibi].GeneratorSpeed_desired);
//        double Torque_rated = fabs( fsi[ibi].Torque_generator_max );
        double dt_wt = dt * ratio_T;
        double K=fabs(0.5*M_PI*ratio_rho*pow(fsi[ibi].r_rotor,5)*fsi[ibi].CP_max/pow(fsi[ibi].TSR_max,3)/pow(fsi[ibi].GearBoxRatio,3));
//        double K=fabs(0.5*M_PI*pow(fsi[ibi].r_rotor,5)*fsi[ibi].CP_max/pow(fsi[ibi].TSR_max,3)/pow(fsi[ibi].GearBoxRatio,3));
        
        double angvel_gen = angvel_wt*fsi[ibi].GearBoxRatio;

	      fsi[ibi].angvel_axis_err =  angvel_rated-angvel_gen;
        fsi[ibi].angvel_axis_err_integral += fsi[ibi].angvel_axis_err*dt_wt; //This is an integral

        double PropErr = kp_prop*fsi[ibi].angvel_axis_err;
        double IntErr  = kp_int *fsi[ibi].angvel_axis_err_integral;

        double PitchComm = 0.0;

        
        if( (PropErr+IntErr) < 0. )
        {
//          PetscPrintf(PETSC_COMM_WORLD, " Region 2 \n");
          fsi[ibi].angvel_axis_err_integral = kp_ratio*fsi[ibi].angvel_axis_err; // this line assumes that kp_prop/kp_int = kt_prop/kt_int
          TorqueOut = K*pow(angvel_gen,2);
          //PitchComm = 1.0; //Hossein temporarily activate
          region = 2;      //Hossein temporarily deactivate
          PitchComm = 0.0;
        }
        else
        {
//          PetscPrintf(PETSC_COMM_WORLD, " Region 2.5-3 \n");
//          PitchComm = 1.0 + PropErr + IntErr; //Hossein temporarily activate
          PitchComm =  PropErr + IntErr;    //Hossein temporarily deactivate
          TorqueOut = K*pow(angvel_gen,2)+ kt_prop*fsi[ibi].angvel_axis_err + kt_int*fsi[ibi].angvel_axis_err_integral;
          region = 3;
        }

        //Hossein added from 3.4.18 code
        /*if(abs(TorqueOut)>2.350e+04){
                TorqueOut = 2.350e+04;
        }*/

        if(abs(TorqueOut)>fsi[ibi].Torque_generator_max){
                TorqueOut = fsi[ibi].Torque_generator_max;
        }

//        PitchComm = min(PitchComm,pitch_max);
//        Due to the data obtained from the Eolos turbine seems that the pitch command is negative (System of reference is different)
        ibm[ibi].pitch[0] = PitchComm;
        ibm[ibi].pitch[1] = PitchComm;
        ibm[ibi].pitch[2] = PitchComm;

        TorqueOut *= sign;

	return(0);
}


PetscErrorCode Calc_GenTorque_SWiFT(IBMNodes *ibm, FSInfo *fsi, PetscReal dt, PetscReal& TorqueOut, PetscInt ibi)
{
        // Jonkman Turbine control system
        
        double sign = fsi[ibi].GeneratorSpeed_desired/fabs(fsi[ibi].GeneratorSpeed_desired);
        // Computing the absolute to make sure they are all positive regardless of the sign in the input file
        double angvel_wt = fabs( fsi[ibi].angvel_axis / ratio_T );
        double angvel_rated = fabs( fsi[ibi].GeneratorSpeed_desired);
        double torque_rated = fabs(fsi[ibi].Torque_generator_max);
//        double Torque_rated = fabs( fsi[ibi].Torque_generator_max );
        double dt_wt = dt * ratio_T;
        double angvel_gen = angvel_wt*fsi[ibi].GearBoxRatio;

        double omega_trns = fabs(fsi[ibi].GeneratorTransSpeed);
        double K = fsi[ibi].K_Gen;
        double alpha = fsi[ibi].Alpha_OT;
        double omega_sync = fsi[ibi].GeneratorSyncSpeed;
//        double region = 0;

//        PetscPrintf(PETSC_COMM_WORLD,"####### SWIFT TURBINE CONTROLLER ########\n"); 
//        PetscPrintf(PETSC_COMM_WORLD,"angvel_gen = %le angvel_rated = %le angvel_trns= %le \n",angvel_gen,angvel_rated,omega_trns); 

        if(angvel_gen>=angvel_rated){ // Region 3
          region=3;
          TorqueOut = torque_rated;
        } else if (angvel_gen<omega_trns) { // Region 2
          TorqueOut = K*pow(angvel_gen,2);
          region=2;
        } else { // Region 2.5
          TorqueOut = alpha*(angvel_gen-omega_sync);
          region=2.5;
        }
//        Due to the data obtained from the Eolos turbine seems that the pitch command is negative (System of reference is different)
        ibm[ibi].pitch[0] = 0.0;
        ibm[ibi].pitch[1] = 0.0;
        ibm[ibi].pitch[2] = 0.0;

        TorqueOut *= sign;
//        PetscPrintf(PETSC_COMM_WORLD,"Region = %f \n",region); 

	return(0);
}


PetscErrorCode Calc_turbineangvel(PetscReal dt, IBMNodes *ibm, FSInfo *fsi)
{
        PetscInt      l, ibi;
        PetscReal F_z, P, r, rx, ry, rz;	

        double TorqueAeroOut;
        double angvelout;
        
        
	for (ibi=0; ibi<NumberOfTurbines; ibi++) {

        	double R=fsi[ibi].r_rotor/reflength_wt;
        	if (FixTipSpeedRatio){ 
            fsi[ibi].angvel_axis=ibm[ibi].Tipspeedratio*ibm[ibi].U_ref/R;
			      double TorqueAero=fsi[ibi].Torque_aero*ratio_rho*pow(ratio_V, 2)*pow(ratio_L,3);
            TorqueAeroOut = TorqueAero;
            angvelout = fsi[ibi].angvel_axis;
			//PetscPrintf(PETSC_COMM_WORLD,"Check: For ibi%i, TSR=%f, U_ref=%f and computed AngVel=%f\n", ibi, ibm[ibi].Tipspeedratio, ibm[ibi].U_ref,  fsi[ibi].angvel_axis);
           }
        	else if (fixturbineangvel) {
        	        fsi[ibi].angvel_axis=fsi[ibi].angvel_fixed;
        	        ibm[ibi].Tipspeedratio = fsi[ibi].angvel_axis*R/ibm[ibi].U_ref;
			            double TorqueAero=fsi[ibi].Torque_aero*ratio_rho*pow(ratio_V, 2)*pow(ratio_L,3);
                  TorqueAeroOut = TorqueAero;
                  
        	} 
		else if (turbinetorquecontrol) {

			double J=fsi[ibi].J_rotation;
			double Torque_c, Torque_a;

			double dt_wt=dt * ratio_T; //dimesional time

      double TorqueGen = 0.;

      ibm[ibi].pitch[0] = 0.0;
      ibm[ibi].pitch[1] = 0.0;
      ibm[ibi].pitch[2] = 0.0;

      if (fsi[ibi].TPCntrl==1){
        Calc_GenTorque_Pitch_C96(ibm, fsi, dt, TorqueGen, ibi);
      }else if(fsi[ibi].TPCntrl==2){
        Calc_GenTorque_SWiFT(ibm, fsi, dt, TorqueGen, ibi);
      }else{
         PetscPrintf(PETSC_COMM_WORLD,"#### WARNING NO TURBINE CONTROLLER SPECIFIED\n"); 
      }

			if (turbineindividualpitchcontrol){
        IndividualPitchControl(ibm, fsi, ibi, dt_wt);
      }

			fsi[ibi].Torque_generator = TorqueGen; 
			double TorqueAero=fsi[ibi].Torque_aero*ratio_rho*pow(ratio_V, 2)*pow(ratio_L,3);
      TorqueAeroOut = TorqueAero;

//   //csantoni 
//        {int rank;
//        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//        printf("Torque_Aereo %le %le %le %le %i\n",fsi[ibi].Torque_aero,ratio_rho,ratio_V,ratio_L,rank);
//        }


      double angvel_0=fsi[ibi].angvel_axis / ratio_T;

      fsi[ibi].angaccel_axis = (TorqueAero - TorqueGen*fsi[ibi].GearBoxRatio)/J; // Hossein added from structure code

      double angvel_1=angvel_0 + dt_wt*(TorqueAero - TorqueGen*fsi[ibi].GearBoxRatio)/J;
//                        double angvel_1=angvel_0 + dt_wt*(TorqueAero - TorqueGen)/J;

      angvelout = angvel_0;

      fsi[ibi].angvel_axis=angvel_1 * ratio_T;
      ibm[ibi].Tipspeedratio = fsi[ibi].angvel_axis*R/ibm[ibi].U_ref;

    	double err = (angvel_1-angvel_0)*J/dt_wt-(TorqueAero-TorqueGen);
//                        PetscPrintf(PETSC_COMM_WORLD, "err of solving Torque equation  %le inertial %le Torque %le \n", err, 
//                                        (angvel_1-angvel_0)*J/dt_wt, TorqueAero-TorqueGen);
		 }
   }
      

        for (ibi=0; ibi<NumberOfTurbines; ibi++) {
                fsi[ibi].ang_axis_old=fsi[ibi].ang_axis;
                double angvel_1=fsi[ibi].angvel_axis;
                fsi[ibi].ang_axis+=angvel_1*dt;
        }


       if ( (ti) % tiout == 0) TurbineTorqueControl_Output(fsi, ibm);

// Control system test file
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if(!rank){
                int ibi=0;
                char filen[80];
                FILE *TestFile;
	        sprintf(filen, "TurbineTorqueTest.dat");
//	        filen= "TurbineTorqueTest.dat";
	        TestFile = fopen(filen, "a");
	        if (TestFile){
	        	PetscFPrintf(PETSC_COMM_WORLD, TestFile, "%le %le %le %le %le %le %le %le %le %le\n", TorqueAeroOut, fsi[ibi].Torque_generator, fsi[ibi].angvel_axis, ibm[ibi].Tipspeedratio, ibm[ibi].pitch[0], angvelout, region,ibm[ibi].U_ref,fsi[ibi].ang_axis, fsi[ibi].ang_axis_old);
	        }
		fclose(TestFile);
	}
			
	return(0);
}


PetscErrorCode IndividualPitchControl(IBMNodes *ibm, FSInfo *fsi, PetscInt ibi, PetscReal dt_wt)
{

  double ang1 = fsi[ibi].ang_axis;
  double ang2 = fsi[ibi].ang_axis + 2.*M_PI/3.0; //Hossein changed plus to minus like 3.4.18 code
  double ang3 = fsi[ibi].ang_axis + 4.*M_PI/3.0; //Hossein changed plus to minus like 3.4.18 code
  
  double Mcos = (2.0/3.0)*(cos(ang1)*fsi[ibi].Moment_bladebending[0] +
                    cos(ang2)*fsi[ibi].Moment_bladebending[1] +
                    cos(ang3)*fsi[ibi].Moment_bladebending[2]);
  
  double Msin = (2.0/3.0)*(sin(ang1)*fsi[ibi].Moment_bladebending[0] +
                    sin(ang2)*fsi[ibi].Moment_bladebending[1] +
                    sin(ang3)*fsi[ibi].Moment_bladebending[2]);
  
  Mcos = Mcos*ratio_rho*pow(ratio_V, 2)*pow(ratio_L,3);
  Msin = Msin*ratio_rho*pow(ratio_V, 2)*pow(ratio_L,3);
  
  fsi[ibi].betacos_IPC = fsi[ibi].betacos_IPC + fsi[ibi].Ki_IPC * dt_wt * Mcos;
  fsi[ibi].betasin_IPC = fsi[ibi].betasin_IPC + fsi[ibi].Ki_IPC * dt_wt * Msin;
  
  ibm[ibi].pitch_IPC[0] = cos(ang1)*fsi[ibi].betacos_IPC + sin(ang1)*fsi[ibi].betasin_IPC;        
  ibm[ibi].pitch_IPC[1] = cos(ang2)*fsi[ibi].betacos_IPC + sin(ang2)*fsi[ibi].betasin_IPC;        
  ibm[ibi].pitch_IPC[2] = cos(ang3)*fsi[ibi].betacos_IPC + sin(ang3)*fsi[ibi].betasin_IPC;        
  
  ibm[ibi].pitch[0] += ibm[ibi].pitch_IPC[0];
  ibm[ibi].pitch[1] += ibm[ibi].pitch_IPC[1];
  ibm[ibi].pitch[2] += ibm[ibi].pitch_IPC[2];

	return(0);
}




PetscErrorCode Yaw_Control(PetscInt ti, PetscInt ibi, FSInfo *fsi_rotncll)
{
// This function applies the different Yaw control schemes

//      PetscPrintf(PETSC_COMM_WORLD, "############################\n");
//      PetscPrintf(PETSC_COMM_WORLD, "Turbine Yaw vector Init \n");
//
//      PetscBarrier(PETSC_NULL);

      PetscReal nx,ny,nz;
      ReadYawInput(ti,nx,ny,nz);

//      PetscBarrier(PETSC_NULL);
//      PetscPrintf(PETSC_COMM_WORLD, "############################\n");
//      PetscPrintf(PETSC_COMM_WORLD, "Yaw Out: %f %f %f %i\n",nx,ny,nz,ti);
  
      double rr=sqrt(nx*nx+ny*ny+nz*nz);

      fsi_rotncll[ibi].nx_tbo = fsi_rotncll[ibi].nx_tb;
      fsi_rotncll[ibi].ny_tbo = fsi_rotncll[ibi].ny_tb;
      fsi_rotncll[ibi].nz_tbo = fsi_rotncll[ibi].nz_tb;

      fsi_rotncll[ibi].nx_tb = nx/rr;
      fsi_rotncll[ibi].ny_tb = ny/rr;
      fsi_rotncll[ibi].nz_tb = nz/rr;

      PetscBarrier(PETSC_NULL);
      PetscPrintf(PETSC_COMM_WORLD, "############################\n");
//      PetscPrintf(PETSC_COMM_WORLD, "Turbine Yaw vector: %f %f %f %i\n",nx,ny,nz,ti);
      PetscPrintf(PETSC_COMM_WORLD, "Turbine Yaw vector2: %f %f %f %i\n",fsi_rotncll[ibi].nx_tb,fsi_rotncll[ibi].ny_tb,fsi_rotncll[ibi].nz_tb,ti);
  

  return(0);
}

PetscErrorCode ReadYawInput(PetscInt ti, PetscReal& nx, PetscReal& ny, PetscReal& nz)
{
// This function reads the ReadYaw Input and create a buffer of data to minimize 
// kernel calls to disk
  int indx = ti%nYw;

//  PetscPrintf(PETSC_COMM_WORLD, "############################\n");
//  PetscPrintf(PETSC_COMM_WORLD, " Indx %i\n",indx);
  
  if(indx!=0){
    nx=sNx[indx]; ny=sNy[indx]; nz=sNz[indx];
    return(0);
  } else {

    int nfskp = ti/nYw; 

    FILE *fYaw;
    char filen[80];
    double angvel, TSR;
    double dummy;
    int idummy;

    sprintf(filen, "yaw_prescrb.inp");
    fYaw = fopen(filen, "r");

//    PetscBarrier(PETSC_NULL);
//    PetscPrintf(PETSC_COMM_WORLD, "############################\n");
//    PetscPrintf(PETSC_COMM_WORLD, " Reading file %s\n",filen);
//    PetscPrintf(PETSC_COMM_WORLD, " nfskp %i nYw %i\n",nfskp, nYw);

    if (fYaw){
       
       for (int ii=0; ii<nfskp*nYw; ii++){
        fscanf(fYaw, "%i %le %le %le", &idummy,&dummy,&dummy,&dummy);
       }

       for (int ii=0; ii<nYw; ii++){
        fscanf(fYaw, "%i %le %le %le", &idummy,&sNx[ii],&sNy[ii],&sNz[ii]);
//        fscanf(fYaw, "%i %le %le %le", &idummy,&dummy,&dummy,&dummy);
//        PetscPrintf(PETSC_COMM_WORLD, " Reading true %i %f %f %f\n",ii,sNx[ii],sNy[ii],sNz[ii]);
       }

    }else{
      PetscPrintf(PETSC_COMM_WORLD, "ERROR: yaw_prescrb.inp NOT FOUND\n");
    }
    fclose(fYaw);
    nx=sNx[0]; ny=sNy[0]; nz=sNz[0];
//    PetscPrintf(PETSC_COMM_WORLD, " Yaw read Out %f %f %f\n",nx,ny,nz);
    return(0);
  }

}
