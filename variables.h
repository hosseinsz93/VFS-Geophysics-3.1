#ifndef _VARIABLES_H_
#define _VARIABLES_H_

#include "petscvec.h"
#include "petscda.h"
#include "petscksp.h"
#include "petscsnes.h"

#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_IJ_mv.h"

#include <assert.h>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

// Seokkoo
#if defined (__CPLUSPLUS__) || defined (__cplusplus)

#include <vector>
#include <algorithm>
#include <assert.h>

#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

template < typename T > T sqr(T const & a) { return a*a; }

typedef struct {
PetscScalar   Txx, Txy, Txz, Tyx, Tyy, Tyz, Tzx, Tzy, Tzz;
} Tensor;

struct Cmpnts
{
    Cmpnts() { }
    Cmpnts(Cmpnts const & a) : x(a.x), y(a.y), z(a.z) { }
    Cmpnts(double _x, double _y, double _z) : x(_x), y(_y), z(_z) { }

    double dist(Cmpnts const & b) const
    {
      return sqrt(sqr(x-b.x) + sqr(y-b.y) + sqr(z-b.z));
    }

    PetscScalar x, y, z;
};

typedef struct {
	PetscScalar	x, y;
} Cmpnts2;

typedef struct {
	PetscReal	t, f;
} FlowWave;

typedef struct {
	PetscReal	x, y;
} Cpt2D;

typedef struct {
	Vec	Ubcs; // An ugly hack, waste of memory
} BCS;

//Hossein
// add (Toni)
//for wave_momentum_source
typedef struct {
	double *WAVE_a, *WAVE_theta, *WAVE_Kz, *WAVE_Kx, *WAVE_KK, *WAVE_P_0, *WAVE_P_1, *WAVE_src_dist, *WAVE_omega;
	double time, PEZ, PEX, WAVE_dt, WAVE_time_since_read, WIND_time_since_read;
	int  NZMOD, NXMOD, WAVE_ti,WAVE_num_max, NZ;
	int *WAVE_ind, WAVE_ind_max;
	double *WIND_U, *WIND_V, *WIND_W;
	double *WIND_Y, *WIND_Z;
	int  WIND_ti;
	int **WIND_id_x,**WIND_id_y;
	double *SK2D;
	double *WVN_z;
	double *WVN_x;
} WAVEInfo;
// end (Toni)

/*
typedef struct {
	Cmpnts csi, eta, zet, ucont, ucat;
	double aj, p;
	int process, nvert;
} Neighbor_Element;	// seokkoo

typedef struct {
	Cmpnts csi, eta, zet, ucont, ucat;
	double aj, p;
	int nvert;
	struct Neighbor_Element *neighbor[6];
	struct Element *next;
} Element;	// seokkoo
*/
typedef struct {
        PetscInt c1,c2,c3;
} Cell2Cell;

typedef struct {
        PetscInt c[100];
} Node2Cell;

typedef struct {
        PetscInt cs[200];
} Cell2Cells;

typedef struct {
  PetscInt	i1, j1, k1;
  PetscInt	i2, j2, k2;
  PetscInt	i3, j3, k3;
  PetscReal	cr1, cr2, cr3; // coefficients
  PetscReal	d_i; // distance to interception point on grid cells
  PetscInt	imode; // interception mode

  PetscReal     ucomp, vcomp, wcomp, utaub; // velocity components and bed shear stress for each IB node
  PetscReal     CC; // Concentration above the bed cell
  PetscReal     fl_den; // Fluid Density above the bed cell
  PetscReal     fl_vis; // Fluid kinematic viscostiy above the bed cell

  PetscReal     w_v; // half depth-averaged vertical velocity
  PetscReal     k_v; // half depth-averaged tke
  PetscInt	ni, nj, nk;	// fluid cell
  PetscReal	d_s; // shortest distance to solid surfaces
  Cmpnts	pmin;
  PetscInt	cell; // shortest distance surface cell
  PetscReal	cs1, cs2, cs3;

  PetscInt	i11, j11, k11;
  PetscInt	i22, j22, k22;
  PetscInt	i33, j33, k33;
  PetscReal	cr11, cr22, cr33; // coefficients
  PetscReal	d_ii; // distance to interception point on grid cells
  PetscInt	iimode; // interception mode
  PetscReal	cs11, cs22, cs33;

  PetscInt	ii1, jj1, kk1;
  PetscInt	ii2, jj2, kk2;
  PetscInt	ii3, jj3, kk3;
  PetscReal	ct1, ct2, ct3; // coefficients
  //PetscReal	d_s; // distance to interception point on grid cells
  PetscInt	smode; // interception mode

  PetscInt	ii11, jj11, kk11;
  PetscInt	ii22, jj22, kk22;
  PetscInt	ii33, jj33, kk33;
  PetscReal	ct11, ct22, ct33; // coefficients
  PetscReal	d_ss; // distance to interception point on grid cells
  PetscInt	ssmode; // interception mode
  
  Tensor        Shear;
  PetscReal     Cs_x, Cs_y, Cs_z;
  
/*   PetscInt      bi1, bj1, bk1; */
/*   PetscInt      bi2, bj2, bk2; */
/*   PetscInt      bi3, bj3, bk3; */
/*   PetscInt      bi4, bj4, bk4; */

/*   PetscReal     bcr1,bcr2,bcr3,bcr4; */
} IBMInfo;

//Hossein added from turbine structure
/* Turbine structure */
/* The grid nodes may be different from the one employed in the flow simulation */
typedef struct {
	int		NumberofInp; // Number of elemets

	PetscReal	*r_inp;
	PetscReal	rstart, rend;
	PetscReal	*mass_inp;	// mass per length  	
	PetscReal	*E_inp, *G_inp;	// Young's Modulus
	PetscReal 	*Iksi_inp, *Ieta_inp, *Icg_inp;	// moments of inertia
	PetscReal	*twistpre_inp;
	PetscReal 	*lcg_inp; 	// distance between center of gravity and center of elastic axis
	PetscReal	*lpi_inp;  // distance from center of elastic axis to center of pitch axis 	
	PetscReal	*Ietaetaksi_inp, *Ietaksiksi_inp;
	PetscReal	*J_inp;


	int		NumberofElements; // Number of elemets

	PetscReal	*Ietaetaksi, *Ietaksiksi;
	PetscReal	*dIetaetaksidr, *dIetaksiksidr;
	PetscReal	*r;
	PetscReal	*mass, *dmassdr, *d2massdr2;	// mass per length  	
	PetscReal	*E, *dEdr, *d2Edr2, *G, *dGdr, *d2Gdr2;	// Young's Modulus
	PetscReal 	*Iksi, *dIksidr, *d2Iksidr2, *Ieta, *dIetadr, *d2Ietadr2, *Icg, *dIcgdr, *d2Icgdr2;	// moments of inertia
	PetscReal	*twistpre, *dtwistpredr, *d2twistpredr2;
	PetscReal 	*lcg, *dlcgdr, *d2lcgdr2; 	// distance between center of gravity and center of elastic axis
	PetscReal	*lpi, *dlpidr, *d2lpidr2, *d3lpidr3, *d4lpidr4;  // distance from center of elastic axis to center of pitch axis 	
	PetscReal	*J, *dJdr;
	

	PetscReal	phi, phi_o, phi_om1, phi_om2; // angular position of rotor
	PetscReal	dphidt, dphidt_o, dphidt_om1; // angular velocity of rotor
	PetscReal	d2phidt2, d2phidt2_o, d2phidt2_om1; // angular acceleration of rotor

	PetscReal	beta, beta_o, beta_om1, beta_om2; // pitch angle of blade 
	PetscReal	dbetadt, dbetadt_o, dbetadt_om1; // pitch velocity of blade 
	PetscReal	d2betadt2, d2betadt2_o, d2betadt2_om1; // pitch acceleration of blade 

	PetscReal	*w0, *w1;
	PetscReal	*dw0dr;
	PetscReal	*u, *v, *theta; // bending displacements and tortional displacement at time n+1 
	PetscReal	*u_o, *v_o, *theta_o; // bending displacements and tortional displacement at time n
	PetscReal	*u_om1, *v_om1, *theta_om1; // bending displacements and tortional displacement at time n-1
	PetscReal	*u_om2, *v_om2, *theta_om2; // bending displacements and tortional displacement at time n-2
	PetscReal	*dudt, *dvdt, *dthetadt; // bending and tortional velocities
	PetscReal	*dudt_o, *dvdt_o, *dthetadt_o; // bending and tortional velocities at time n
	PetscReal	*dudt_om1, *dvdt_om1, *dthetadt_om1; // bending and tortional velocities at time n-1
	PetscReal	*dudt_om2, *dvdt_om2, *dthetadt_om2; // bending and tortional velocities at time n-2
	PetscReal	*d2udt2, *d2vdt2, *d2thetadt2; // bending and tortional accelerations

	PetscReal	*err_u, *err_v, *err_theta, *err_dudt, *err_dvdt, *err_dthetadt;

	PetscReal	*ucg, *vcg, *ucghat; 
	PetscReal	*ucg_o, *vcg_o, *ucghat_o; 
	PetscReal	*ucg_om1, *vcg_om1, *ucghat_om1; 
	PetscReal	*ducgdt, *dvcgdt, *ducghatdt; 

	PetscReal	*fu, *fv, *ftheta, *M; 	// external forces and twisting moment. fw is in the radial direction at time n+1
	PetscReal	*fu_o, *fv_o, *ftheta_o, *M_o; 	// external forces and twisting moment. fw is in the radial direction at time n

	/* right-hand-side for bending */
	PetscReal	*Fu0,   *Fu1,   *Fu2,   *Fu3,   *Fu4,   *Fu5,   *Fu6;	
	PetscReal	*Fu0_o, *Fu1_o, *Fu2_o, *Fu3_o, *Fu4_o, *Fu5_o, *Fu6_o;

	PetscReal	*Fv0,   *Fv1,   *Fv2,   *Fv3,   *Fv4,   *Fv5,   *Fv6;
	PetscReal	*Fv0_o, *Fv1_o, *Fv2_o, *Fv3_o, *Fv4_o, *Fv5_o, *Fv6_o;

	/* right-hand-side for twisting */
	PetscReal	*Ftheta0,   *Ftheta1,   *Ftheta2,   *Ftheta3,   *Ftheta4,   *Ftheta5,   *Ftheta6;	
	PetscReal	*Ftheta0_o, *Ftheta1_o, *Ftheta2_o, *Ftheta3_o, *Ftheta4_o, *Ftheta5_o, *Ftheta6_o;

	/* right-hand-side sum for displacements and temporal derivative of displacements */
	PetscReal	*rhsu,     *rhsdudt;
	PetscReal	*rhsv,     *rhsdvdt;
	PetscReal	*rhstheta, *rhsdthetadt;

	PetscReal	*rhsu_o,     *rhsdudt_o;
	PetscReal	*rhsv_o,     *rhsdvdt_o;
	PetscReal	*rhstheta_o, *rhsdthetadt_o;

	PetscReal	*rhsu_om1,     *rhsdudt_om1;
	PetscReal	*rhsv_om1,     *rhsdvdt_om1;
	PetscReal	*rhstheta_om1, *rhsdthetadt_om1;

	PetscReal	*rhsu_om2,     *rhsdudt_om2;
	PetscReal	*rhsv_om2,     *rhsdvdt_om2;
	PetscReal	*rhstheta_om2, *rhsdthetadt_om2;

        PetscReal       *resu,     *resdudt;
        PetscReal       *resv,     *resdvdt;
        PetscReal       *restheta,     *resdthetadt;

	PetscReal	*coefu_dudr, *coefu_d2udr2, *coefu_d3udr3, *coefu_d4udr4;
	PetscReal	*coefu_dvdr, *coefu_d2vdr2, *coefu_d3vdr3, *coefu_d4vdr4;
	PetscReal	*coefu_dthetadr, *coefu_d2thetadr2;

	PetscReal	*coefv_dudr, *coefv_d2udr2, *coefv_d3udr3, *coefv_d4udr4;
	PetscReal	*coefv_dvdr, *coefv_d2vdr2, *coefv_d3vdr3, *coefv_d4vdr4;
	PetscReal	*coefv_dthetadr, *coefv_d2thetadr2;

	PetscReal	*coeftheta_dudr, *coeftheta_d2udr2, *coeftheta_d3udr3, *coeftheta_d4udr4;
	PetscReal	*coeftheta_dvdr, *coeftheta_d2vdr2, *coeftheta_d3vdr3, *coeftheta_d4vdr4;
	PetscReal	*coeftheta_dthetadr, *coeftheta_d2thetadr2;

	PetscReal	*coefdudt_ddt_dudr, *coefdudt_ddt_dvdr;
	PetscReal	*coefdvdt_ddt_dudr, *coefdvdt_ddt_dvdr;

	double 		**A, **A_LU, *b, *bb;

	int 		*indx;
	/* BCs */
	double		u_bc0, dudr_bc0, d2udr2_bc0, d3udr3_bc0, d4udr4_bc0;
	double		u_bc1, dudr_bc1, d2udr2_bc1, d3udr3_bc1, d4udr4_bc1;

	double		v_bc0, dvdr_bc0, d2vdr2_bc0, d3vdr3_bc0, d4vdr4_bc0;
	double		v_bc1, dvdr_bc1, d2vdr2_bc1, d3vdr3_bc1, d4vdr4_bc1;

	double		theta_bc0, dthetadr_bc0, d2thetadr2_bc0, d3thetadr3_bc0, d4thetadr4_bc0;
	double		theta_bc1, dthetadr_bc1, d2thetadr2_bc1, d3thetadr3_bc1, d4thetadr4_bc1;

	int		bctype;

	/* grid spacing */
	double		dh;

	/* size of time step*/
	double 		dt;

  /*  */
  double    hubrad;

} Beam;

typedef struct {
        /* for actuator line model */
        int        	num_AC, num_CD, num_CL, num_Aerodym; // num_AC: number of input points for angle of pitch and chord length of blade
        PetscReal       *r_AC, *angle_twistInp, *chord_bladeInp;  // pitch angel and chord length from input
        PetscReal       *ang_CD, *CDInp, *ang_CL, *CLInp, *CMInp, *ang_Aerodym;
        PetscReal       r_beg, r_end; // This type airfoil begins at r=r_beg and 
} FOIL;

typedef struct {
	PetscInt	nbnumber;
	PetscInt	n_v, n_elmt;	// number of vertices and number of elements
	PetscInt	my_n_v, my_n_elmt;	// seokkoo, my proc
	PetscInt	*nv1, *nv2, *nv3;	// Node index of each triangle
	PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
	PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes
	PetscReal	*x_bp0, *y_bp0, *z_bp0;
        PetscReal       *x_bp_i, *y_bp_i, *z_bp_i;
	PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
	Cmpnts		*u, *uold, *urm1, *u_pre;
  	Cmpnts        	*af, *af_old, *af_pre; //Acceleration

	PetscReal     	*x_bp_pre, *y_bp_pre, *z_bp_pre;	
	PetscReal     	*error;
	PetscReal     	*error_force;

	// Aitken
	PetscReal   	Aitken_factor;
	PetscReal   	Aitken_factor_force;

	// shear stress projection
	PetscReal       *Txx,*Txy, *Txz, *Tyy, *Tzz, *Tyz;
  	PetscReal       *Traction_x, *Traction_y, *Traction_z;
  	PetscReal       *Traction_x_pre, *Traction_y_pre, *Traction_z_pre;
  	PetscReal       *Area_x, *Area_y, *Area_z;
  	PetscReal      	*SumP_Tx,*SumP_Ty,*SumP_Tz,*AbsP_tau;
 	PetscReal     	*SumM_Tx,*SumM_Ty,*SumM_Tz,*AbsM_tau;

	// for thin body
  	PetscReal 	*Traction_x_positiveside;
  	PetscReal 	*Traction_y_positiveside;
  	PetscReal 	*Traction_z_positiveside;

  	PetscReal 	*Traction_x_negativeside;
  	PetscReal 	*Traction_y_negativeside;
  	PetscReal 	*Traction_z_negativeside;

	PetscReal     	*dA ;         // area of an element
	PetscReal     	*nt_x, *nt_y, *nt_z; //tangent dir
	PetscReal     	*ns_x, *ns_y, *ns_z; //azimuthal dir

	PetscReal     	*cent_x,*cent_y,*cent_z;
	PetscReal	*disp_x, *disp_y, *disp_z, *disp_theta;

	//Hossein added from turbine structure
	// element center displacement
	PetscReal	*ux_struct, *uy_struct, *uz_struct, *utheta_struct;


    // Added 1/10/2010 ali
    Cmpnts        *Bvel;
    PetscReal     *Bvel_u, *Bvel_v, *Bvel_w;
    PetscReal     *vx, *vy, *vz;
    PetscReal     *deltz_p_us, *deltz_p_ds, *A_us, *A_ds;
    PetscReal     *BShS, *BCShS, *Shvel, *Shvel_v, *BShS_v;
    PetscReal     *SCont;
    PetscReal     *C;
    PetscReal     *fluid_den, *fluid_vis;
    PetscReal     *w_ave, *k_ave;
    PetscReal     *max_bed_angle, *SedFlux;
    PetscInt      *Rigidity, *Mobility;
    PetscReal     *Delz, *Dely, *dzbp,*dVol,*cent_z_old, *cent_y_old, *cent_zl, *cent_yl, *cent_z_AVE, *cent_y_AVE, *atke_old, *z_bp_l, *y_bp_l;
    PetscReal      Deltab;
    PetscReal      sed_density;
    PetscReal     *elmt_depth, *scc, *sbb;
    PetscReal     *u_grad_x,*v_grad_x,*c_grad_x;//,*z_grad_x;
    PetscReal     *u_grad_y,*v_grad_y,*c_grad_y;//,*z_grad_y;
	PetscReal     *c_grad_z,*w_grad_z,*w_grad_x,*u_grad_z; 
    PetscReal     *u_max,*v_max,*w_max,*c_max;//,*z_max;
    PetscReal     *netflux_old;
	PetscReal     *Dflux; //Hossein
	PetscReal     *Eflux; //Hossein
	PetscReal     *qflux; //Hossein
    PetscReal     *F_flux_12,*F_flux_13,*F_flux_23;
    PetscReal     dtime[1000000],time_bedchange[1000000];
    

    Cell2Cell  *c2c;
    Node2Cell  *n2c;
    Cell2Cells  *c2cs;
	// for radius check
	Cmpnts *qvec;
	PetscReal *radvec; 
	
	//seokkoo
	
	PetscInt *_nv1, *_nv2, *_nv3;	// for rank 0
	PetscInt *count;
	PetscInt *local2global_elmt;
	PetscInt total_n_elmt, total_n_v;
	
	PetscReal *_x_bp, *_y_bp, *_z_bp;	// for rank 0
	PetscReal *shear, *mean_shear;
	PetscReal *reynolds_stress1;
	PetscReal *reynolds_stress2;
	PetscReal *reynolds_stress3;
	PetscReal *pressure;
	PetscInt  fixedNodes;	//Not Rotating
	PetscReal 	*Tmprt_lagr, *Ftmprt_lagr, *tmprt;
	PetscReal	volume, inertiamoment_x, inertiamoment_y, inertiamoment_z;
	PetscReal 	*F_lagr_x, *F_lagr_y, *F_lagr_z; // force at the IB surface points (lagrange points)
	PetscReal 	*dF_lagr_x, *dF_lagr_y, *dF_lagr_z; 
	PetscReal 	*U_lagr_x, *U_lagr_y, *U_lagr_z; // estimated fluid velocity at Lagrange points 
	PetscReal 	*Urelmag; // relative incoming velocity for actuator model 

	Cmpnts		*Urel; // vector of the relative incoming velocity
	Cmpnts		*Uinduced; // vector of the induced velocity
	Cmpnts		*circulation; // circulation vector on the blade
	Cmpnts		*liftdirection; // direction of the lift 
	Cmpnts		*rotationdirection; // direction of the lift 



	Cmpnts		*Urel_mean; // vector of the relative incoming velocity
	Cmpnts		*Uinduced_mean; // vector of the induced velocity
	Cmpnts		*circulation_mean; // circulation vector on the blade
	Cmpnts		*liftdirection_mean; // direction of the lift 

        int        	*i_min, *i_max, *j_min, *j_max, *k_min, *k_max;
        /* ACL */
        PetscReal       *angle_attack, *angle_twist, *chord_blade; // twist angle and chord length at each grid point
        PetscReal       *CD, *CL;
	PetscReal	*pitch;  
	PetscReal	*pitch_old;  
	PetscReal	*pitch_IPC;  
	PetscReal	pitch_min;

	PetscReal 	U_ref;

	

	PetscReal	*dhx, *dhy, *dhz;


	PetscReal 	friction_factor, pressure_factor;

	PetscReal	*frictionfactor;
	
	PetscReal	dh;

	PetscReal	indf_axis, Tipspeedratio, CT, indf_tangent;


        PetscReal 	*Fr_mean, *Fa_mean, *Ft_mean; // force at the IB surface points (lagrange points)
        PetscReal 	*Ur_mean, *Ua_mean, *Ut_mean;
	PetscReal 	*AOA_mean, *Urelmag_mean;


        PetscReal 	*AOAAOA_mean, *FFa_mean, *FFt_mean; 


	PetscReal	*centIP_x, *centIP_y, *centIP_z;

        PetscInt        	*iIP_min, *iIP_max, *jIP_min, *jIP_max, *kIP_min, *kIP_max;
        PetscReal 	*U_IPlagr_x, *U_IPlagr_y, *U_IPlagr_z, *P_IPlagr;

	PetscReal	*dh_IP;

	PetscReal	*Nut_lagr, *Shear_lagr_x, *Shear_lagr_y, *Shear_lagr_z;
	PetscReal	*ShearDesired_lagr_x, *ShearDesired_lagr_y, *ShearDesired_lagr_z;
	PetscReal	*UU_lagr_x, *UU_lagr_y, *UU_lagr_z;

	PetscReal	*random_color;

	PetscInt		num_cf;
	PetscReal	r_in[200], cf_in[200]; // Readed friction coefficient from a file  HARD CODING

	PetscReal	forcecoeffx_actuator, forcecoeffy_actuator, forcecoeffz_actuator, Uref_actuator, Lref_actuator;


	PetscInt		num_circulationInp;
	PetscReal	*circulationInp, *r_circulationInp; // simulation with circulation specified along the blade  
	PetscReal	*circulationSpecified;

	PetscInt		num_blade, num_foiltype;

	Beam		*bladestructure;

	FOIL		*acl;

	PetscReal	r_nearhubinflowcorrection;

	PetscReal	coeff_SG;

	PetscReal	r_nacelle;

	PetscInt 		*node_boundary;
	PetscInt 		*elmt_boundary;
		
	// xiaolei add SEDI
		
    	PetscInt  		*color;
    	PetscInt 		*s2l; // actuator line element index for each actuator surface element 
	
	// end add SEDI

	// sediment end
	//
	// fsi_IBM
	 PetscReal      F_x, F_y, F_z;

  	// Total surface area
     	PetscReal      	Area;
} IBMNodes;

typedef struct node {
	PetscInt Node;
	struct node *next;
} node;

typedef struct list{
	node *head;
} LIST;


typedef struct list_node {
	PetscInt	index;
	struct list_node *next;
} Node_List;

typedef struct IBMListNode {
	IBMInfo ibm_intp;
	struct IBMListNode* next;
} IBMListNode;

typedef struct IBMList {
	IBMListNode *head;
} IBMList;


typedef struct UserCtx {
	
	DA da;	/* Data structure for scalars (include the grid geometry informaion, to obtain the grid information, use DAGetCoordinates) */
	DA fda;	// Data Structure for vectors
	DA fda2;	// Data Structure for vectors with 2 variables
	DALocalInfo info;

	Vec	Cent;	// Coordinates of cell centers
	Vec	Csi, Eta, Zet, Aj;
	Vec	ICsi, IEta, IZet, IAj;
	Vec	JCsi, JEta, JZet, JAj;
	Vec	KCsi, KEta, KZet, KAj;
	//Vec	/*Area,lArea,*/ Volume, lVolume;
	Vec	Ucont;	// Contravariant velocity components
	Vec	Ucat;	// Cartesian velocity components
	Vec	Ucat_o;
	Vec	Ucont_o, Ucont_rm1, Rhs, dUcont, pUcont;

	Cmpnts **ucat_plane;		// for pseudo periodic BC
	Cmpnts **ucat_plane_is;		// for pseudo periodic BC
	Cmpnts **ucat_plane_js;		// for pseudo periodic BC
	Cmpnts **ucat_plane_ks;		// for pseudo periodic BC
	
  
	//Vec	lUcont_c;
	Vec	RHS_o;//, RHS_rm1;		// for AM2 scheme
	Vec	dP;
//	Vec	RHS_IRK[6];	// for implicit RK
	int	current_IRK_stage; 
	//Vec	Ucont_rm2;//Conv_o, Visc_o;//, Conv_rm1;// 	// seokkoo

	Vec Fp;
	Vec	Div1, Div2, Div3;
	Vec	Visc1, Visc2, Visc3;
  
	Vec	K_Omega, lK_Omega, K_Omega_o, lK_Omega_o;//, K_Omega_rm1;
	Vec	Conc, lConc, Conc_o, lConc_o;
	Vec	FCurrent, lFCurrent;
	Vec	Conc_Javg;
	Vec	Distance;
	Vec	lNu_t, lF1;		// eddy viscosity, f1 constant for SST model
	Vec	lUcat_old;
  
	Vec	P, P_o;
	Vec	Phi;
	Vec	GridSpace;
	Vec	Nvert;
	Vec	Nvert_o;
	Vec 	Nvert_o_fixed;
	Vec   Itfc, ItfcP;
	BCS	Bcs;

  Vec cnpyNvert;
  Vec lcnpyF;
  Vec cnpyF;
  Vec cnpyCdy;

	int	rhs_count;
	Vec	rhsD;
  //Vec Density;
// seokkoo, Dirichlet BC RHS
	Vec	Levelset, Levelset_o, lLevelset, lDensity, lMu;	// seokkoo, density of fluids
	Vec	lST;	//surface tension
	Vec	lCs;
	Vec	lUstar;
	Vec	Gid, Gidm;				// seokkoo, Global ID for local array containg ghost points
	Vec	Phi2, B2, Ucont2;
	
	int	local_Phi2_size, p_global_begin;
	int	reduced_p_size;
  FILE *fp_inflow_u, *fp_inflow_ko, *fp_inflow_l, *fp_inflow_t, *fp_inflow_t1, *fp_inflow_t2, *fp_inflow_t3;
  
	double shear_stress[6], shear_force_avg[6];
	double mean_k_flux, mean_k_area;
  double *local_bulk_velocity; 
	//std::vector<double> k_area;
	double *k_area;
  //double ustar_avg, Re_tau_avg;
	double ustar_now[6];
	
	//Hossein
	// add (Toni)
	//for wave_momentum_source
	WAVEInfo *wave_inf;
	Vec	WAVE_fp,lWAVE_fp;	
// end (Toni)
	
	double lA_cyl, A_cyl;
	double lA_cyl_x, A_cyl_x;
	double lA_cyl_z, A_cyl_z;
	double lFvx_cyl, lFvz_cyl, Fvx_cyl, Fvz_cyl;
	double lFpx_cyl, lFpz_cyl, Fpx_cyl, Fpz_cyl;
  
	SNES snes;
	Vec rhs;
	Mat J;
	
	Vec	lUcont, lUcat, lP, lPhi;
	Vec	lUcont_o, lUcont_rm1;//, ldUcont;
	Vec	lCsi, lEta, lZet, lAj;
	Vec	lICsi, lIEta, lIZet, lIAj;
	Vec	lJCsi, lJEta, lJZet, lJAj;
	Vec	lKCsi, lKEta, lKZet, lKAj;
	Vec	lGridSpace;
	Vec	lNvert, lNvert_o, lNFace;
	Vec	lCent;
	Vec   lItfc, lItfcP;
	Vec	inletU;
	Vec	nhostU, nhostP;
	Vec	DUold;
	Vec	Forcing;
	Vec	Ucont_MG;

	Vec	Dt;

	AO	ao;

	PetscReal	ren;	// Reynolds number
	PetscReal	dt; 	// time step
	PetscInt	sed_tio, SimpleCellCheck, osl_inlet_sediment_flux; 	// tio out for sediment
	PetscReal	dts; 	// time step for sediment transport
	PetscReal	st;	// Strouhal number
	PetscReal	Ri;	// Richardson number
  
	PetscInt        Paraview;
	// added for new poisson solver
	PetscReal     FluxIntpSum;
        //PetscReal     u_in[500][500], v_in[500][500],w_in[500][500];
	// Added Iman
	Vec           psuedot;
	PetscReal     cfl, vnn;

	PetscReal	r[101], tin[101], uinr[101][1001];

	PetscReal *itfchostx, *itfchosty, *itfchostz;
	PetscReal FluxInSum, FluxOutSum;

	PetscErrorCode aotopetsc;
	PetscTruth assignedA;

	PetscInt _this;
	PetscInt *idx_from;

	PetscInt bctype[6];
	PetscInt itfcptsnumber;
	PetscInt *itfcI, *itfcJ, *itfcK;
	PetscInt *itfchostI, *itfchostJ, *itfchostK, *itfchostB;
	PetscInt	IM, JM, KM; // dimensions of grid

	PetscInt ibmnumber;
	IBMInfo  *ibm_intp;
	Mat	A, C;
	KSP   ksp;

	IBMNodes *ibm;

	DA	*da_f, *da_c;
	struct UserCtx *user_f, *user_c;
	Vec	*lNvert_c;

	Vec	B;
	Vec	Rhsp, X, R;
  
	Mat	MR, MP;
	MatNullSpace nullsp;

  
  /* Variables for multi-nullspace case */
	PetscInt *KSKE;
	PetscTruth multinullspace;

	IBMList *ibmlist;

	PetscInt thislevel, mglevels;

	PetscTruth *isc, *jsc, *ksc;
  
	FlowWave *inflow;
	PetscInt number_flowwave;
  
	
  #if defined (__CPLUSPLUS__) || defined (__cplusplus)
	/*
	int map_set;
	int solver_set;
	bool has_singletone;
	Epetra_Map *Map;
	Epetra_CrsMatrix *Mat_P;
	Epetra_Vector *Vec_P, *Vec_RHS;
	
	Epetra_LinearProblem *Poisson_Problem;
	ML_Epetra::MultiLevelPreconditioner *MLPrec;
	Epetra_CrsSingletonFilter *SingletonFilter;
	Epetra_LinearProblem *Reduced_Problem;
	*/
  #endif
	Vec Ucat_sum;		// sum of u, v, w over time
	Vec Ucat_cross_sum;	// sum of uv, vw, wu
	Vec Ucat_square_sum;	// sum of uu, vv, ww
        Vec P_sum;		// sum of p
        Vec Conc_sum;		// sum of concentration
        Vec Conc_cross_sum;	// sum of cu, cv, cw
        Vec Level_sum;		// sum of level
        Vec Level_square_sum;		// sum of level*level
        Vec lUstar_sum, lUstar_;// sum and instantanious ustar
	Vec K_sum; // sum of Kinetic energy in RANS
	Vec Nut_sum; // sum of eddy viscosity
	Vec P_square_sum;
	/*Vec P_cross_sum;*/
	
	Vec Udp_sum; //size 1; u*dpdx + v*dpdy + w*dpdz
	Vec dU2_sum; //size 3; (dui_dx)^2 + (dui_dy)^2 + (dui_dz)^2;  ... (3*3=9 components)
	Vec UUU_sum; //size 3; (u^2+v^2+w^2)*ui
	Vec tauS_sum; // size 1; sum of tau_ij Sij = 2 nu_t Sij Sij; tau_ij = +2 nu_t Sij
	
	Vec Vort_sum;
	Vec Vort_square_sum;
	
	Vec lSx, lSy, lSz, lS;
	Vec lLM, lMM, lNM;
	
	//int **free_surface_j;		// free_surf[i][k] = (top j cell)
	double **free_surface_location;	// y or z position[i][k]; only for rank=0
	//double **free_surface_p;	// Dirichlet pressure BC
	//double **vof;			// volume fraction at surface
	//int **bottom_surface_j;		// bottom_surf[i][k] = (bottom j cell)
	//double **bottom_surface_y;	// y position
	Cmpnts **ucont_plane;		// for pseudo periodic BC
	Cmpnts2 **komega_plane;
	Cmpnts2 **komega_plane_is;
	Cmpnts2 **komega_plane_js;
	Cmpnts2 **komega_plane_ks;

	double **tmprt_plane;
        double **tmprt_plane_is;
        double **tmprt_plane_js;
        double **tmprt_plane_ks;
        double **tmprt1_plane;
        double **tmprt1_plane_is;
        double **tmprt1_plane_js;
        double **tmprt1_plane_ks;
        double **tmprt2_plane;
        double **tmprt2_plane_is;
        double **tmprt2_plane_js;
        double **tmprt2_plane_ks;
        double **tmprt3_plane;
        double **tmprt3_plane_is;
        double **tmprt3_plane_js;
        double **tmprt3_plane_ks;

			// levelset inflow
        double **level_plane;
        double **level_plane_is;
        double **level_plane_js;
        double **level_plane_ks;
	
        Vec F_eul;              // Force on the background grids
        Vec FIB_eul;              // Force on the background grids
        Vec lF_eul;
        Vec lFIB_eul;
        Vec Nut_eul;              // distributed nut on the background grid nodes
        Vec lNut_eul;

        Vec Tmprt, lTmprt, Tmprt_o, lTmprt_o, Tmprt_rm1, lTmprt_rm1; 
        Vec FTmprt, lFTmprt;             //buoyancy force

		// more than one scalars 
        Vec Tmprt1, lTmprt1, Tmprt1_o, lTmprt1_o, Tmprt1_rm1, lTmprt1_rm1;  
        Vec Tmprt2, lTmprt2, Tmprt2_o, lTmprt2_o, Tmprt2_rm1, lTmprt2_rm1;  
        Vec Tmprt3, lTmprt3, Tmprt3_o, lTmprt3_o, Tmprt3_rm1, lTmprt3_rm1;  

        Vec Pr_t, lPr_t;

        PetscInt bctype_tmprt[6];
	Vec Curv;
	Vec Grad_abs;
	Vec Heaviside;
	Vec lNvert_4diffusedinterface;
  	PetscReal    rho_fluid;


} UserCtx;

typedef struct {
	UserCtx *user;
	PetscInt thislevel;
} MGCtx;

typedef struct {
	PetscInt mglevels;
	PetscInt thislevel;

	PetscTruth  isc, jsc, ksc;
	MGCtx *mgctx;
} UserMG;

typedef struct {
	PetscReal	P;    //Press on the surface elmt
	PetscInt		n_P; //number of Press Pts on the elmt
	PetscReal	Tow_ws, Tow_wt; //wall shear stress of the elmt
	PetscReal	Tow_wn; // normal stress 

	PetscInt		Clsnbpt_i,Clsnbpt_j,Clsnbpt_k; //Closest Near Bndry Pt to each surf elmt 
	PetscInt		icell,jcell,kcell;
	PetscInt		FoundAroundcell;
	PetscInt		Need3rdPoint;
	//PetscInt      Aroundcellnum;
} SurfElmtInfo;

typedef struct {
	PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
	PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
	PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
	PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
	PetscReal    F_x_old,F_y_old,F_z_old; //Forces & Area
	PetscReal    F_x_real,F_y_real,F_z_real; //Forces & Area
	PetscReal    M_x,M_y,M_z; // Moments
	PetscReal    M_x_old,M_y_old,M_z_old; //Forces & Area
	PetscReal    M_x_real,M_y_real,M_z_real; //Forces & Area
	PetscReal    M_x_rm2,M_y_rm2,M_z_rm2; //Forces & Area
	PetscReal    M_x_rm3,M_y_rm3,M_z_rm3; //Forces & Area
	PetscReal    x_c,y_c,z_c; // center of rotation
	PetscReal    CMx_c, CMy_c, CMz_c; // center of mass
	PetscInt     rot_dir;
	PetscReal    XYangle; 
	PetscReal    fixed_ang_vel; //for forced rotation
	PetscReal    Mdpdn_x, Mdpdn_y,Mdpdn_z;
	PetscReal    Mdpdn_x_old, Mdpdn_y_old,Mdpdn_z_old;

	// Aitkin's iteration
	PetscReal    dS[6],dS_o[6],atk,atk_o;

	// for force calculation
	SurfElmtInfo  *elmtinfo;
	IBMInfo       *fsi_intp;

	//Max & Min of ibm domain where forces are calculated
	PetscReal    Max_xbc,Min_xbc;
	PetscReal    Max_ybc,Min_ybc;
	PetscReal    Max_zbc,Min_zbc;

	// CV bndry
	PetscInt     CV_ys,CV_ye,CV_zs,CV_ze;

  PetscInt TPCntrl;   // Turbine torque and collective pitch controller   1:Clipper_C96   2:Jonkman et al. (no pitch controller)

  PetscInt YawCntrl;   // Turbine Yaw Controller   0:Fixed   1:Prescribed (lookup table)

  PetscInt yawtistrt; // time step at which the turbine start the yawing
	
       PetscReal    nx_tb, ny_tb, nz_tb; // direction vector of rotor axis rotor_model
        PetscReal    nx_tbo, ny_tbo, nz_tbo; // direction vector of rotor axis rotor_model

        PetscReal    angvel_z, angvel_x, angvel_y, angvel_axis, angvel_axis_err, angvel_axis_err_old, angvel_axis_err_integral, angvel_axis_err_derivative;

	PetscReal    K_proportional, K_integral, K_derivative;

	PetscReal    Kratio_torque, angvel_axis_err_relax, angaccel_axis;  
	PetscReal    *Moment_bladebending;  // the bending moment for each blade, for individual pitch control  

	PetscReal    Ki_IPC; // gain for individual pitch control

	PetscReal    betacos_IPC, betasin_IPC;

	PetscReal    WindSpeed_rated, WindSpeed_cutin, WindSpeed_cutout;

        PetscReal    x_c0,y_c0,z_c0;

	PetscReal    Torque_generator, J_rotation, CP_max, TSR_max, r_rotor, Torque_aero, ang_axis, angvel_fixed, Force_axis, Torque_generator_max, GeneratorSpeed_desired, GearBoxRatio;  

  PetscReal    ang_axis_old;

  PetscReal    K_Gen, Alpha_OT; // Generator constant
 
  PetscReal    GeneratorSyncSpeed, GeneratorTransSpeed; //GeneratorSyncSpeed = GeneratorSpeed_desired/(1+slip_percentage)

	int rotate_alongaxis;

	PetscReal    xnacelle_upstreamend, ynacelle_upstreamend, znacelle_upstreamend;

}	FSInfo;

typedef struct {
        PetscReal    x_c, y_c, z_c; // center of mass
        PetscReal    x_r,y_r,z_r; // center of rotation
        PetscInt     rot_dir; //axis of rotation, 0=x, 1=y, 2=z, 3=z in xy plane
        PetscReal    XYangle; //skew or yaw about z axis
        PetscReal    fixed_ang_vel; //for forced rotation

} 	FSISetup;

//#define COEF_TIME_ACCURACY 1.5
extern PetscReal COEF_TIME_ACCURACY;

#if defined (__CPLUSPLUS__) || defined (__cplusplus)
void Contra2Cart(UserCtx *user);
void Contra2Cart_2(UserCtx *user);
void Contra2Contra(UserCtx *user, Vec lUcont_center);
void Contra2Cart_single(Cmpnts &csi, Cmpnts &eta, Cmpnts &zet, Cmpnts &ucont, Cmpnts *ucat);

void DestroyIBMList(IBMList *ilist);
void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp);
void InitIBMList(IBMList *ilist);
void insertnode(LIST *ilist, PetscInt Node);
void destroy(LIST *ilist);
PetscErrorCode Blank_Interface(UserCtx *user);
PetscInt intsect_triangle(PetscReal orig[3], PetscReal dir[3], PetscReal vert0[3], PetscReal vert1[3], PetscReal vert2[3], PetscReal *t, PetscReal *u, PetscReal *v);
PetscTruth ISLineTriangleIntp(Cmpnts p1, Cmpnts p2, IBMNodes *ibm, PetscInt ln_v);
PetscInt ISPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, PetscReal nfx, PetscReal nfy, PetscReal nfz);
PetscErrorCode Dis_P_Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, PetscReal *d);
PetscErrorCode triangle_intp2(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, IBMInfo *ibminfo);
PetscErrorCode triangle_intp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, IBMInfo *ibminfo);
PetscErrorCode triangle_intpp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, IBMInfo *ibminfo);
void initlist(LIST *ilist);
PetscErrorCode InflowWaveFormRead(UserCtx *user);
PetscErrorCode FormMetrics(UserCtx *user);
PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k, PetscInt *ih, PetscInt *jh, PetscInt *kh, UserCtx *user);
PetscInt lidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
PetscErrorCode MG_Initial(UserMG *usermg, IBMNodes *ibm);
//PetscErrorCode ibm_read_ucd(IBMNodes *ibm);
PetscErrorCode ibm_read_ucd(IBMNodes *ibm, PetscInt ibi, PetscReal CMx_c, PetscReal CMy_c, PetscReal CMz_c, PetscReal x_r, PetscReal y_r, PetscReal z_r, PetscReal XYangle);
	// in ibm_io.c
PetscErrorCode FsiInitialize(PetscInt n_elmt, FSInfo *fsi,PetscInt ibi);
PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ti);
PetscErrorCode Read_Rotate_FSI_Parameter_Input(FSISetup *fsimov);
PetscErrorCode Elmt_Move_FSI_ROT(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, PetscInt ibi);
PetscErrorCode Elmt_Move_FSI_TRANS(UserCtx *user, FSInfo *FSinfo, IBMNodes *ibm); //Hossein
PetscErrorCode ibm_surface_out(IBMNodes *ibm, PetscInt ti, PetscInt ibi);
PetscErrorCode ibm_search_advanced(UserCtx *user, IBMNodes *ibm, PetscInt ibi);
PetscErrorCode ibm_interpolation_advanced(UserCtx *user);
PetscErrorCode fluxin(UserCtx *user);
PetscErrorCode Struc_Solver(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi, PetscInt itr_sc, PetscInt tistart, PetscTruth *DoSCLoop);

//Hossein added from turbine structure
/**
 * turbine structure model
 */ 

extern PetscInt turbinestructuremodel;
extern PetscInt torsion_turbinestructure;
extern PetscInt restart_turbinestructure;
PetscErrorCode init_bladestructure(IBMNodes *ibm);
PetscErrorCode solvestructurefunctions_bladestructure(IBMNodes *ibm);

PetscErrorCode solvestructurefunctions1euler_bladestructure(IBMNodes *ibm);

PetscErrorCode transferbladedisplacement2flowblade(IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode getforcesfromflowsolver_bladestructure(IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode DispProjection_l2s(IBMNodes *ibm_surface, IBMNodes *ibm_line, int NumberOfObjects); 
PetscErrorCode writefiles_bladestructure(IBMNodes *ibm);
// add end (xiaolei)

//ali added 2 June 2011
PetscErrorCode Scour(UserCtx * user, IBMNodes  * ibm, PetscInt tistart, PetscInt ti, PetscInt itr_sc);
PetscErrorCode ibm_intp_pj_centroid(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode Projecting_new(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode flow_variables_ref_level_new(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode flow_variables_ref_level(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode Projecting(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode Finalizing_Projecting(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode sediment_variables_projection(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode bed_concentration(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode bed_concentration_levelset(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode sediment_flux(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode sediment_flux_y_direction(UserCtx * user, IBMNodes  * ibm); //Hossein
PetscErrorCode outlet_sediment_flux_bend(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode outlet_sediment_flux_bend_y_direction(UserCtx * user, IBMNodes  * ibm); //Hossein
PetscErrorCode inlet_sediment_flux_contra(UserCtx * user, IBMNodes * ibm);
PetscErrorCode inlet_sediment_flux_contra_TexWash(UserCtx * user, IBMNodes * ibm);
PetscErrorCode outlet_sediment_flux_contra(UserCtx * user, IBMNodes  * ibm);
PetscErrorCode outlet_sediment_flux_contra_y_direction(UserCtx * user, IBMNodes  * ibm); //Hossein
PetscErrorCode check_correct_new_elev(UserCtx * user, IBMNodes  * ibm, PetscInt ti);
PetscErrorCode	avalanche_first_sweep( UserCtx * user, IBMNodes * ibm, PetscInt ti);
PetscErrorCode	avalanche_second_sweep( UserCtx * user, IBMNodes * ibm, PetscInt ti);
PetscErrorCode read_bed_cell_depth(IBMNodes  * ibm, PetscInt tistart);
PetscErrorCode Smoothing_Shear_vel( UserCtx * user, IBMNodes * ibm);
 
PetscErrorCode recomputing_geometry(UserCtx * user, IBMNodes * ibm, PetscInt tistart, PetscInt ti, PetscInt itr_sc, PetscInt aval_loop, PetscInt avalanche_check_number);

PetscErrorCode ibm_surface_VTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti);

// end of added
PetscErrorCode Flow_Solver(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi, PetscInt itr_sc, IBMNodes *wtm, FSInfo *fsi_wt, IBMNodes *ibm_ACD, IBMNodes *ibm_acl2ref, FSInfo *fsi_acl2ref, IBMNodes *ibm_nacelle, FSInfo *fsi_nacelle); //Hossein added PetscInt itr_sc
PetscErrorCode MG_Finalize(UserMG *usermg);
PetscErrorCode GhostNodeVelocity(UserCtx *user);
PetscErrorCode InflowFlux(UserCtx *user) ;
PetscErrorCode OutflowFlux(UserCtx *user);
PetscErrorCode FormBCS(UserCtx *user, FSInfo *fsi);
PetscErrorCode Block_Interface_U(UserCtx *user);
PetscErrorCode ComputeRHS(UserCtx *user, PetscInt istage);
//PetscErrorCode FormFunction1(Vec Ucont, Vec Rhs, UserCtx *user);
//PetscErrorCode FormFunction1_FVM(Vec Ucont, Vec Rhs, UserCtx *user);
PetscErrorCode Spectral(UserCtx *user);
PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv);
PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc);
PetscErrorCode OutflowVelocity(UserCtx *user, Vec Ucont);
PetscErrorCode Find_fsi_2nd_interp_Coeff(PetscInt i, PetscInt j, PetscInt k, PetscInt elmt, Cmpnts p, IBMInfo *ibminfo,UserCtx *user, IBMNodes *ibm);
PetscErrorCode Find_fsi_2nd_interp_Coeff2(PetscInt i, PetscInt j, PetscInt k, PetscInt elmt, Cmpnts p, IBMInfo *ibminfo, UserCtx *user, IBMNodes *ibm);
PetscReal detmnt(PetscReal a[3][3]);
PetscErrorCode MyFieldRestriction(UserCtx *user);
PetscErrorCode Calc_forces_SI(FSInfo *FSinfo,UserCtx *user, IBMNodes *ibm,PetscInt ti, PetscInt ibi, PetscInt bi);
PetscErrorCode Calc_FSI_pos_SC(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime, PetscReal Re);
PetscErrorCode SwingCylinder(FSInfo *fsi, IBMNodes *ibm);
PetscErrorCode Calc_FSI_Ang(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime,PetscInt ibi, UserCtx *user);
PetscErrorCode Calc_FSI_Ang_intg(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt,PetscInt itrSC, PetscInt ibi, UserCtx *user);
PetscErrorCode FSI_DATA_output(FSInfo *FSinf, PetscInt ti);
PetscErrorCode FSI_DATA_Output(FSInfo *FSinfo, PetscInt ibi);
PetscErrorCode CollisionDetectionOfCylinders(FSInfo *fsi, PetscInt NumberOfBodies);
PetscErrorCode MyNvertRestriction(UserCtx *user_h, UserCtx *user_c);
PetscErrorCode ImplicitMomentumSolver(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImplicitMomentumSolver1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImpRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode PoissonSolver_MG(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonSolver_MG_original(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode UpdatePressure(UserCtx *user);
PetscErrorCode Projection(UserCtx *user);
PetscErrorCode Divergence(UserCtx *user);
PetscErrorCode Ucont_P_Binary_Output(UserCtx *user);
PetscErrorCode Ucat_Binary_Output(UserCtx *user);
PetscErrorCode SetInitialGuessToOne(UserCtx *user);
PetscErrorCode VolumeFlux(UserCtx *user, Vec lUcor, PetscReal *ibm_Flux, PetscReal *ibm_Area);
PetscErrorCode triangle_intp3D(double x, double y, double z, 
						double x1, double y1, double z1, 
						double x2, double y2, double z2,
						double x3, double y3, double z3, IBMInfo *ibminfo);
PetscErrorCode RungeKutta_Advection(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, double *dudc, double *dvdc, double *dwdc, double *dude, double *dvde, double *dwde, double *dudz, double *dvdz, double *dwdz);
void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz);
void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz );
void Compute_dscalar_dxyz ( double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
							double dkdc, double dkde, double dkdz, double *dk_dx, double *dk_dy, double *dk_dz);

void Conv_Diff_IC(UserCtx *user);
void Conv_Diff_BC(UserCtx *user);
void Conv_Diff_Function_RHS(UserCtx *user);

void Compute_Distance_Function(UserCtx *user);
void K_Omega_IC(UserCtx *user);
void Compute_Smagorinsky_Constant_1(UserCtx *user, Vec Ucont, Vec Ucat);
void Compute_eddy_viscosity_LES(UserCtx *user);
PetscErrorCode PoissonLHSNew(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo);
void Convert_Phi2_Phi(UserCtx *user);
void PoissonSolver_Hypre(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonRHS2(UserCtx *user, Vec B);

void wall_function (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_freesurface (double nu, double sb, Cmpnts Ua, Cmpnts Ub, PetscReal *ustar, double nx, double ny, double nz); //Hossein
void wall_function_slip (double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_roughness (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_roughness_loglaw (double nu, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void wall_function_roughness_a (UserCtx *user, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
//void wall_function_roughness_a_levelset (UserCtx *user, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz, double fluid_density, double fluid_visosity);

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);
void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, double nx, double ny, double nz);

void Create_Hypre_Solver();
void Create_Hypre_P_Matrix_Vector(UserCtx *user);
void MatHYPRE_IJMatrixCopy(Mat v,HYPRE_IJMatrix &ij);
void Petsc_to_Hypre_Vector(Vec A, HYPRE_IJVector &B, int i_lower);
void Hypre_to_Petsc_Vector(HYPRE_IJVector &B, Vec A, int i_lower);
void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3]);
void Calculate_normal_and_area(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3], double *Ai, double *Aj, double *Ak);

void Compute_du_i (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz);
				
void Compute_du_j (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz);

void Compute_du_k (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
				double *dudc, double *dvdc, double *dwdc, 
				double *dude, double *dvde, double *dwde,
				double *dudz, double *dvdz, double *dwdz);
				
void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz );

void Compute_dscalar_i (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz );
void Compute_dscalar_j (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz );
void Compute_dscalar_k (int i, int j, int k, int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz );

void Calculate_dxdydz(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz);
void AxByC ( double a, Cmpnts &X, double b, Cmpnts &Y, Cmpnts *C);
void AxC ( double a, Cmpnts &X, Cmpnts *C);
double PPM(double WW, double W, double E, double EE, double a);
void Compute_Q(UserCtx *user, Vec Q);
void Levelset_BC(UserCtx *user);
void Distance_Function_RHS (UserCtx *user, Vec Levelset_RHS, int wall_distance);
void Init_Levelset_Vectors(UserCtx *user);
void Destroy_LevelSet_Vectors(UserCtx *user);
void Grad_phi_criterion_Dennis (UserCtx *user); //<--------------- DENNIS ADD
void Reinit_Levelset(UserCtx *user);
void Levelset_Function_IC(UserCtx *user);
void Advect_Levelset(UserCtx *user);

double time_coeff();
void IB_BC_Ucat(UserCtx *user);
void Calc_ShearStress(UserCtx *user);
void write_shear_stress_ibm();
void Calc_k_Flux(UserCtx *user);//, double *Flux, double *Area);

void Compute_Density(UserCtx *user);
double find_utau_Cabot(double nu,  double u, double y, double guess, double dpdn);
double find_utau_Cabot_roughness(double nu, double u, double y, double guess, double dpdn, double ks);
double u_Cabot(double nu, double y, double utau, double dpdn);
double u_Cabot_roughness(double nu, double y, double utau, double dpdn, double ks_plus);
PetscErrorCode Formfunction_2(UserCtx *user, Vec Rhs, double scale);
PetscErrorCode FormFunction_SNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr);;
void initial_guess_for_snes(UserCtx *user, Vec U);
void Pressure_Gradient(UserCtx *user, Vec dP);
double Calc_Minimum_dt (UserCtx *user);
void write_data (UserCtx *user);
void Set ( Cmpnts *A, double a );
double H (double p, double dx);	// Heaviside function
double dH (double p, double dx);	// its derivative
double mean ( double A, double B );
//double mean0 ( double A, double B );
double sign(double a);
double sign1(double a, double dx);
double mod_sign(double d0, double grad, double e);
void get_weight ( int i, int j, int k, int mx, int my, int mz, PetscReal ***aj, PetscReal ***nvert, double nv, double weight[3][3][3]);
double integrate_testfilter(double val[3][3][3], double w[3][3][3]);
void Compute_Surface_Tension(UserCtx *user);
void Compute_dlevel_center_levelset (int i, int j, int k,  int mx, int my, int mz, double sgn, int wall_distance, PetscReal ***level, PetscReal ***nvert, double *dldc, double *dlde, double *dldz);
void Mass_Correction_Levelset (UserCtx *user, Vec D);
double weno3(double f0, double f1, double f2, double f3, double wavespeed);
double eno2(double f0, double f1, double f2, double f3, double a);
void Calc_Inlet_Area(UserCtx *user);
void pre_integrate ();
int file_exist(char *str);
double Upwind(double W, double E, double a);
void Keep_Constant_Flux(UserCtx *user);
//PetscErrorCode VolumeFlux2(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg);
PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg);
void Add_IBMFlux_to_Outlet(UserCtx *user, PetscReal ibm_Flux);
PetscErrorCode ibm_surface_out_with_pressure(IBMNodes *ibm, PetscInt ibi);
void Initialize_free_surface_location_vector(UserCtx *user);
void Calc_free_surface_location(UserCtx *user);
void IB_BC(UserCtx *user);
//extern void SetDirichletPressure(UserCtx *user);

//Hossein
//Solitary waves
extern PetscErrorCode Solitary_wave_inlet_velocity_profile_Boussinesq(UserCtx *user, double *uin, double *vin, double *win, double x, double y, double z, double win_inlet_flux, int jj, double time, int ti);             //Hossein changed coordinates
extern PetscErrorCode Solitary_wave_inlet_elevation_profile_Boussinesq(UserCtx *user, double *free_surface_elevation, double *free_surface_elevation_to, double x, double y, double z, double z_to, double time, int ti); //Hossein changed coordinates

//Linear wave single
extern PetscErrorCode Linear_wave_single_inlet_velocity_profile (UserCtx *user, double *uin, double *vin, double *win, double x, double y, double z, double win_inlet_flux, int jj, double time, int ti);
extern PetscErrorCode Linear_wave_single_inlet_elevation_profile (UserCtx *user, double *free_surface_elevation, double *free_surface_elevation_to, double x, double y, double z, double z_to, double time, int ti);

extern HYPRE_Solver pcg_solver_p, precon_p;
extern HYPRE_IJMatrix Ap;
extern HYPRE_ParCSRMatrix par_Ap;
extern HYPRE_IJVector Vec_p, Vec_p_rhs;
extern HYPRE_ParVector par_Vec_p, par_Vec_p_rhs;

extern double imp_free_tol;
extern IBMNodes	*ibm_ptr;
extern FSInfo   	*fsi_ptr;
extern FSInfo		*fsi_ptr;
extern UserCtx	*user_ptr;
extern PetscInt user_m, user_n, user_p;
extern double mean_pressure_gradient;
extern int testfilter_ik, testfilter_1d;
extern int i_periodic, j_periodic, k_periodic;
extern int ii_periodic, jj_periodic, kk_periodic;
extern int periodic;
extern int poisson_it, amg_agg;
extern double amg_thresh;
extern int i_homo_filter, j_homo_filter, k_homo_filter;
extern int clark, my_rank;
extern PetscInt inletprofile;
extern PetscInt inletCase;
extern PetscInt inletprofile_tmprt;
extern double inlet_flux;
extern int tistart;
extern int les, wallfunction, rans, lowRe, slipbody, central, second_order, weno;
extern int bed_roughness, STRONG_COUPLING;
extern PetscInt sediment, input_ib_depth, conv_diff, SuspendedParticles, mobile_bed, projection_method, periodic_morpho, density_current, sandwave, Barchans, smooth_shear, non_dimensional, English, sediment_influx;
extern PetscReal k_ss, w_s, Nseg, U_Bulk, deltab, d50, FlowDepth, porosity, sbbb, Cs_, Angle_repose,Background_Conc,Inlet_Conc;
extern double inlet_sediment_flux, sediment_thickness, initial_bed_elevation, max_mobilebed_z, min_mobilebed_x,min_mobilebed_y, x_outlet, y_outlet,x_inlet,x_limit_inlet, y_inlet, y_limit_inlet, mm_l, bb_l, mm_l_in, bb_l_in, cell_size, cell_depth;
extern PetscTruth sediment_thickness_flag, XOutlet, YOutlet, XInlet,X_Limit_Inlet, YInlet, Y_Limit_Inlet;
extern PetscInt LiveBed, LOutlet, LInlet;
extern PetscInt no_ibm_search;
extern PetscInt effestive_bed_shear_stress, aval_loop, sand_slide;
extern PetscInt RigidBed;
extern PetscInt zero_grad;
extern PetscInt ti, tiout, mixed, tistart, inflow_recycle_perioid;
extern PetscInt itr_sc;
extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed, invicid;
extern PetscInt TwoD, mixed, averaging;
extern char path[256], gridfile[256];
extern  double dx_min, di_min, dj_min, dk_min;
extern  double di_max, dj_max, dk_max;

extern PetscReal max_cs;//, Fr;
extern IBMNodes	*ibm_ptr;
extern  FSInfo        *fsi_ptr;
extern PetscInt implicit, TwoD, initialzero;
extern PetscInt movefsi, rotatefsi;
extern PetscTruth dpdz_set;
extern int ib_bctype[128];
extern int dynamic_freq;
extern int laplacian, rotational, skew, tiout_ufield, tiend_ufield;
extern int levelset;
extern int dam_break,k_gate;
extern PetscInt flag_level_set_conv;
extern double rho_water, rho_air;
extern double mu_water, mu_air;
extern double dthick;
extern double seudo_dt;
extern PetscTruth dthick_set;
extern int inviscid, surface_tension, poisson;
extern double gravity_x, gravity_y, gravity_z;
extern double inlet_y, outlet_y, inlet_z, outlet_z, inlet_x, outlet_x;
extern int fix_outlet, fix_inlet, fix_level;
//extern int y_direction; //Hossein

extern PetscReal FluxInSum, FluxOutSum;
extern PetscReal FluxInSum_gas, FluxOutSum_gas;
extern PetscTruth inlet_z_flag, inlet_x_flag, inlet_y_flag; //Hossein added all flags
extern int ibm_search;
extern PetscInt  inlet_buffer_k;
extern PetscTruth rough_set;
extern double roughness_size;

extern int save_inflow, save_inflow_period, save_inflow_minus, pseudo_periodic, inletprofile, read_inflow_period;
extern void save_inflow_section(UserCtx *user);
extern void read_inflow_section(UserCtx *user);
extern void store_k_CrossSection(UserCtx *user);
extern PetscErrorCode Calc_Moments(FSInfo *FSinfo, IBMNodes *ibm, SurfElmtInfo *elmtinfo, PetscReal Re, PetscInt ti);
extern PetscErrorCode ibm_Surf_stress(UserCtx *user, IBMNodes *ibm, SurfElmtInfo *elmtinfo, PetscInt ibi);

extern double dthick_const, dt_inflow;

extern int save_ksection[1000];
extern int save_jsection[1000];
extern int save_isection[1000];
extern int nsave_ksection;
extern int nsave_jsection;
extern int nsave_isection;

extern int ucat_plane_allocated;
extern int ucat_plane_allocated_is;
extern int ucat_plane_allocated_js;
extern int ucat_plane_allocated_ks;

extern int ucont_plane_allocated;

extern char path_inflow[256];
extern char path_plane_save[256];

extern PetscReal scale_velocity;

extern int levelset_it;
extern double levelset_tau;	
extern int cross_diffusion;
extern int rotdir;
extern double St_exp;
extern double angvel;
extern int ti_lastsave;
extern double x_r, y_r, z_r;


//Hossein
// add (Toni)
	//for levelset
extern int levelset_it;
extern double levelset_tau;	
extern double level_in_height;
extern int level_in;
extern int levelset_weno;
	//variables for FSI
extern int fsi_6dof;
extern double body_mass;
extern double body_inertia_x, body_inertia_y, body_inertia_z;
extern double body_alpha_rot_x, body_alpha_rot_y, body_alpha_rot_z;
extern double body_alpha_lin_x, body_alpha_lin_y, body_alpha_lin_z;
extern double body_beta_rot_x, body_beta_rot_y, body_beta_rot_z;
extern double body_beta_lin_x, body_beta_lin_y, body_beta_lin_z;
extern int forced_motion, fall_cyll_case;
extern double angle_x0,angle_y0,angle_z0;
	//for solitary waves
extern PetscInt  solitary_wave;
extern PetscInt  ti_start_solitary_wave;
extern PetscInt  ti_restart_solitary_wave;
extern PetscReal inlet_bed_elevation;
extern PetscReal inlet_z_for_solitary_wave;
extern PetscReal solitary_wave_amplitude;
	//for linear wave single
extern PetscInt  ti_start_linear_wave_single;
extern PetscInt  ti_restart_linear_wave_single;
extern PetscReal inlet_z_for_linear_wave_single;
extern PetscReal linear_wave_single_amplitude;
extern PetscReal linear_wave_single_number;
	//for wave_momentum_source
extern int wave_momentum_source, wave_sponge_layer;
extern int wave_IB;
extern double wave_angle_single, wave_K_single, wave_depth, wave_a_single, wave_source_cent_z;
extern double  wave_sponge_zs, wave_sponge_z01, wave_sponge_z02, wave_sponge_xs, wave_sponge_x01, wave_sponge_x02;
double epscos (double x, double eps);
	//for air_flow_levelset
extern int air_flow_levelset, air_flow_levelset_periodic;
extern int wave_average_k;
extern int wind_skip, wind_start_read, wind_recicle, wave_skip, wave_start_read, wave_recicle;
extern int wave_k_ave, wave_i_ave, wave_ti_startave, wave_ti_start;
extern int freesurface_wallmodel, viscosity_wallmodel;
extern double channel_height;
extern double wave_wind_reflength;
extern double wave_wind_refvel;
extern double wave_wind_yshift;
extern int floating_turbine_case;
// end (Toni)


extern PetscErrorCode Export_lagrdata(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension);
extern PetscErrorCode Calc_Ftmprt_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_Ftmprt_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df);
extern PetscErrorCode Calc_Tmprt_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_F_lagr_specified(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_F_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_F_lagr_nacelle(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_F_lagr_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Comput_actualshear_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Comput_modelcoef_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Comput_desiredshear_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_F_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df);
extern PetscErrorCode Calc_FIB_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df);
extern PetscErrorCode Comput_nut_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_Nut_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df);
extern PetscErrorCode Calc_forces_rotor(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi, char fname[80], int NumberOfObjects);
extern PetscErrorCode Calc_forces_actuator(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi);
extern PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm, int NumberOfObjects); 
extern PetscErrorCode Pre_process_all(UserCtx *user, IBMNodes *ibm, int NumberOfObjects);
extern PetscErrorCode ACL_read(IBMNodes *ibm, int ibi, FSInfo *fsi, double reflength);
extern PetscErrorCode Calc_F_lagr_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Export_ForceOnBlade(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Calc_U_IPlagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Coordinates_IP(IBMNodes *ibm, int NumberOfObjects);
extern PetscErrorCode Pre_process_IP(UserCtx *user, IBMNodes *ibm, int NumberOfObjects);

extern PetscErrorCode ColorIB_tmp(IBMNodes *ibm, int NumberOfObjects); // Hossein added from 3.4 code temporarily to solve the issue of rstart rotation

extern PetscErrorCode airfoil_ACL(IBMNodes *ibm,  FSInfo *fsi, int NumberOfObjects);

extern PetscErrorCode Calc_forces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi);
extern PetscErrorCode Calc_eulforces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi);
extern PetscErrorCode Calc_F_lagr_noslip(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfTurbines); // xyang


extern PetscErrorCode Uref_ACL(UserCtx *user, IBMNodes *ibm_ACL, IBMNodes *ibm_ACD, FSInfo *fsi_wt, int NumberOfObjects);
void RHS_Tmprt(UserCtx *user, Vec Tmprt_RHS);

void Tmprt_IC(UserCtx *user);

void Tmprt_BC(UserCtx *user);

extern PetscErrorCode FormFunction_Tmprt(SNES snes, Vec Tmprt, Vec Rhs, void *ptr);
void Solve_Tmprt(UserCtx *user);

void Force_Tmprt(UserCtx *user);

extern PetscErrorCode TECIOOut_rhs1(UserCtx *user, Vec Rhs);
extern PetscErrorCode Add_fluc(UserCtx *user);
extern PetscErrorCode Add_fluc_tmprt(UserCtx *user);
extern PetscErrorCode UpdateXYZ_MoveFrame(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode Trans1DUp_XYZ(IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern PetscErrorCode refAL_Rot(FSInfo *FSinfo, IBMNodes *ibm, IBMNodes *ibm_ref, int ibi);
extern PetscErrorCode rotor_Rot(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension);
extern PetscErrorCode bladepitch_Rot(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects);

extern void rotate_xyz (double ti, double dt, double angvel, double x_c, double y_c, double z_c, double x_bp0, double y_bp0, double z_bp0, double *x_bp, double *y_bp, double *z_bp, double *rot_angle);

extern PetscInt NumberOfTurbines; //xyang
extern PetscInt NumberOfNacelle; //xyang
extern PetscInt NumNacellePerLoc;
extern PetscInt i_periodicIB, j_periodicIB, k_periodicIB;
extern PetscInt maxiteraction_rotormodel;  
extern PetscInt NonUniform_ADModel;  
extern PetscInt rotor_model;  
extern PetscInt TheSameObject;  
extern PetscInt UlagrFromSurface;  
extern PetscInt nacelle_model;  
extern PetscReal indf_ax; // xyang 12-16-2010
extern PetscReal reflength_wt;  
extern PetscReal reflength_nacelle;  
extern PetscReal refvel_wt;  
extern PetscReal refvel_cfd;  
extern PetscInt deltafunc;
extern PetscReal halfwidth_dfunc;
extern PetscReal r_nacelle;
extern PetscReal loc_refvel;

extern PetscInt inflow_levelset;

extern PetscInt temperature; 
extern PetscInt temperature1; 
extern PetscInt temperature2; 
extern PetscInt temperature3; 

extern PetscInt les_prt;
extern PetscReal u_frame, v_frame, w_frame;
extern PetscInt MoveFrame;
extern PetscInt ii_periodicWT, jj_periodicWT, kk_periodicWT; // periodic WT, a row/column of ghost wind turbines needs to be added
extern PetscReal Sx_WT, Sy_WT, Sz_WT;
extern PetscInt Nx_WT, Ny_WT, Nz_WT;
PetscErrorCode disk_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi, int OneDmZ, char fname[80], double reflength);
extern PetscErrorCode deltafunc_test( );


extern PetscErrorCode Connectivity_ib( UserCtx * user, IBMNodes * ibm );
PetscErrorCode surface_read(int fileType, IBMNodes *ibm, int ibi, FSInfo *fsi, char fname[80], double reflength);
PetscErrorCode surface_write(IBMNodes *ibm, int ibi, char fname[80]);
void matchLineColorToSurface(int ibi, IBMNodes * sufibm, IBMNodes * linibm);
PetscErrorCode nacelleYaw_IB(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
PetscErrorCode rotorYaw(IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode calc_s2l(IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects); 
PetscErrorCode ForceProjection_l2s(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects);
PetscErrorCode UlagrProjection_s2l(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects); 


extern PetscInt FixTipSpeedRatio;
extern PetscInt fixturbineangvel;
extern PetscInt rstart_turbinerotation;
extern PetscInt turbinetorquecontrol;
extern PetscInt turbineindividualpitchcontrol;
extern PetscReal particle_dens;
extern PetscInt Shen_AL;
extern PetscInt correction3D_CH;
extern PetscInt correction3D_CL;
extern PetscInt smoothforce_AL;
extern PetscInt Shen1_AL;
extern PetscReal correction_ALShen1;
extern PetscReal c0_CL;
extern PetscReal c1_CH;
extern PetscReal c2_CH;
extern PetscReal c3_CH;
extern PetscReal a_shen;
extern PetscReal b_shen;
extern PetscReal c_shen;
extern PetscInt Prandtl_AL;
extern PetscReal refangle_AL;
extern PetscErrorCode TurbineTorqueControlInitialization(FSInfo *fsi, IBMNodes *ibm);
extern PetscErrorCode Calc_turbineangvel(PetscReal dt, IBMNodes *ibm, FSInfo *fsi);
extern PetscErrorCode Yaw_Control(PetscInt ti, PetscInt ibi, FSInfo *fsi_rotncll);
extern PetscReal count_AL;
extern PetscInt forcewidthfixed; 
extern PetscInt forcedistrwidth_surfacescale; 
extern PetscReal dhi_fixed; 
extern PetscReal dhj_fixed; 
extern PetscReal dhk_fixed; 
extern PetscInt correction3D_DS;
extern PetscReal a_DS;
extern PetscReal b_DS;
extern PetscReal d_DS;
extern PetscReal CD_0;
extern PetscReal AOA_0;
extern PetscReal tipcorrectionratio_Fa;
extern PetscInt cf_nacelle_fromfile;
extern PetscInt specifycirculation;
extern PetscReal coeff_SG;
extern PetscInt turbinestructuremodel;
extern PetscInt torsion_turbinestructure;
PetscErrorCode calc_ibm_volumeFlux(IBMNodes *ibm, PetscReal delti, PetscReal *VolumeFlux);
PetscErrorCode TurbineTorqueControl_Output(FSInfo *fsi, IBMNodes *ibm);
PetscErrorCode Save_IBdelta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, char fname[80], int NumberOfObjects);
PetscErrorCode Update_Nodes_IBdelta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Calc_motion_IBdelta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Calc_intgFTorque_IBdelta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Pre_process_ibdelta(UserCtx *user, IBMNodes *ibm);

/** Hossein added from turbine structure 
 * blade strucutre public variables 
 */
extern PetscReal dt_turbinestructure;

/*
static double Butcher_ESDIRK2_A[4][4] = {
		{ 0., 0., 0., 0. },
		{ 1767732205903./4055673282236., 1767732205903./4055673282236., 0., 0. },
		{ 2746238789719./10658868560708., -640167445237./6845629431997., 1767732205903./4055673282236., 0. },
		{ 1471266399579./7840856788654., -4482444167858./7529755066697., 11266239266428./11593286722821., 1767732205903./4055673282236. },
};

static double Butcher_ESDIRK2_B[4] = { 1471266399579./7840856788654., -4482444167858./7529755066697., 11266239266428./11593286722821., 1767732205903./4055673282236. };
	
static double Butcher_ESDIRK4_A[6][6] = {
		{ 0., 0., 0., 0., 0., 0. },
		{ 1./4., 1./4., 0., 0., 0., 0. } , 
		{ 8611./62500., -1743./31250., 1./4., 0., 0., 0. }, 
		{ 5012029./34652500., -654441./2922500., 174375./388108., 1./4., 0., 0. },
		{ 15267082809./155376265600., -71443401./120774400., 730878875./902184768., 2285395./8070912., 1./4., 0.},
		{ 82889./524892., 0., 15625./83664., 69875./102672., -2260./8211., 1./4. },
};

static double Butcher_ARK4_A[6][6] = {
		{ 0., 0., 0., 0., 0., 0. },
		{ 1./2., 0., 0., 0., 0., 0. },
		{ 13861./62500., 6889./62500., 0., 0., 0., 0. },
		{ -116923316275./2393684061468., -2731218467317./15368042101831., 9408046702089./11113171139209., 0., 0., 0. },
		{ -451086348788./2902428689909., -2682348792572./7519795681897., 12662868775082./11960479115383., 3355817975965./11060851509271., 0., 0. },
		{ 647845179188./3216320057751., 73281519250./8382639484533., 552539513391./3454668386233., 3354512671639./8306763924573., 4040./17871., 0.},
};

static double Butcher_ESDIRK4_B[6] = { 82889./524892., 0., 15625./83664., 69875./102672., -2260./8211., 1./4. };
*/
#endif
#endif
