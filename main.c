#include "variables.h"
#include "canopy.h"
static char help[] = "Testing programming!";
PetscReal COEF_TIME_ACCURACY=1.5;
PetscInt ti,tistart=0;
PetscReal    Flux_in = 4.104388e-04, angle = 0;
PetscInt tiout = 10;
PetscInt sed_tio = 100;
PetscInt SimpleCellCheck = 1;
PetscInt osl_inlet_sediment_flux = 0;
PetscInt sediment_influx = 0;
PetscInt block_number;
PetscReal FluxInSum, FluxOutSum;
PetscReal FluxInSum_gas, FluxOutSum_gas;    // seokkoo
PetscInt immersed = 0;
PetscInt inviscid = 0;
PetscInt movefsi = 0, rotatefsi=0;
PetscInt sediment=0, input_ib_depth=0;
PetscInt English=0;
PetscInt bed_roughness=0;
PetscInt Barchans=0;
PetscInt smooth_shear=0;
PetscReal sediment_thickness=1.0;
PetscTruth sediment_thickness_flag=PETSC_FALSE;
PetscReal initial_bed_elevation=0.;
PetscReal max_mobilebed_z=0.9;
PetscReal min_mobilebed_x=-1.e8;
PetscReal min_mobilebed_y=-1.e8;
PetscReal x_outlet=1000.0;
PetscReal y_outlet=1000.0;
PetscReal x_inlet=1000.0;
PetscReal y_inlet=1000.0;
PetscReal x_limit_inlet=1.e8;
PetscReal y_limit_inlet=1.e8;
PetscReal mm_l=0.0;
PetscReal bb_l=0.0;
PetscReal mm_l_in=0.0;
PetscReal bb_l_in=0.0;
PetscReal cell_size=0.25;
PetscReal cell_depth=0.25;
PetscTruth XOutlet=PETSC_FALSE;
PetscTruth YOutlet=PETSC_FALSE;
PetscTruth XInlet=PETSC_FALSE;
PetscTruth YInlet=PETSC_FALSE;
PetscTruth X_Limit_Inlet=PETSC_FALSE;
PetscTruth Y_Limit_Inlet=PETSC_FALSE;
PetscInt LOutlet=0;
PetscInt LInlet=0;
PetscInt non_dimensional=1;
PetscInt projection_method=1;
PetscReal k_ss=0.001;
PetscInt implicit = 0;
PetscInt imp_MAX_IT = 50; 
PetscInt radi=10;
PetscInt inletprofile=1;
PetscInt inletCase=1;
PetscInt inletprofile_tmprt=1;
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
PetscInt  mg_MAX_IT=30, mg_idx=1, mg_preItr=1, mg_poItr=1;
PetscReal imp_atol=1e-7, imp_rtol=1e-4, imp_stol=1.e-8;
PetscInt TwoD = 0;
PetscInt STRONG_COUPLING=0;
PetscInt rstart_fsi=0;
PetscInt cop=0, regime=1; // 1 escape regime --- 2 cruise regime
PetscInt fish=0;
PetscReal St_exp=0.5,wavelength=0.95;
PetscInt MHV=0;
PetscReal max_angle = -54.*3.1415926/180.;// for MHV; min=0
PetscInt thin=0;
PetscInt dgf_z=0,dgf_y=1,dgf_x=0;
PetscInt dgf_az=0,dgf_ay=0,dgf_ax=1 ;
PetscInt NumberOfBodies=1;
PetscReal L_dim;
PetscInt averaging=0;
PetscInt binary_input=0;
PetscInt xyz_input=0;
PetscInt les=0;
PetscInt inlet_buffer_k=1;
PetscInt wallfunction=0;
PetscInt slipbody=-1000;
PetscInt central=0, second_order=0, weno = 0;
PetscInt initialzero=0;
PetscInt freesurface=0;
PetscInt rans=0, lowRe=0;
PetscInt cross_diffusion=1;
PetscInt surface_tension=0;
PetscInt conv_diff=0;
PetscInt no_ibm_search=0;
PetscInt sandwave=0;
PetscInt density_current=0;
PetscInt SuspendedParticles=0;
PetscInt mobile_bed=0;
PetscReal w_s=0.0;
PetscReal Background_Conc=0.0;
PetscReal Inlet_Conc=0.0;
PetscReal porosity=0.45;
PetscReal sbbb=0.02;
PetscReal Cs_=1.5;
PetscReal Angle_repose=45.;
PetscReal U_Bulk=1.0;
PetscInt LiveBed = 0;
PetscReal deltab = 0.001;
PetscReal sed_density = 2650.0;
PetscReal FlowDepth = 1.0;
PetscReal d50 = 0.0001;
PetscInt effective_bed_shear_stress = 0;
PetscInt aval_loop = 0;
PetscInt sand_slide = 0;
PetscInt RigidBed = 0;
PetscInt zero_grad = 0;
PetscInt periodic_morpho=0;
PetscReal Nseg=23;
PetscInt Paraview=0;
PetscInt y_direction=0;

PetscInt    flag_level_set_conv=0; // <---- DENNIS ADD

//Hossein
    // Solitary waves
PetscInt  solitary_wave = 0;
PetscInt  ti_start_solitary_wave = 0;
PetscInt  ti_restart_solitary_wave = 0;
PetscReal inlet_bed_elevation = 0.5;
PetscReal inlet_z_for_solitary_wave = 0.0;
PetscReal solitary_wave_amplitude = .3;

    // Linear wave single
PetscInt  ti_start_linear_wave_single = 0;
PetscInt  ti_restart_linear_wave_single = 0;
PetscReal inlet_z_for_linear_wave_single = 0.0;
PetscReal linear_wave_single_amplitude = 0.2;
PetscReal linear_wave_single_number = 2.0;
// add (Toni)
	//Variables for levelset
double level_in_height=0.;	
int level_in=0;
int levelset_it=10;
double levelset_tau=0.01;
int levelset_weno=0;
	//Variables for wave_momentum_source
int wave_momentum_source=0;
int wave_sponge_layer=0;
double wave_source_cent_z=0.; //Ali added May 1, 2016
double wave_sponge_zs=0.;//length of the sponge layer
double wave_sponge_z01=-10000.;//start of the sponge layer at the left x boundary
double wave_sponge_z02=10000.;//start of the sponge layer at the right x boundary	
double wave_sponge_xs=0.;//length of the sponge layer
double wave_sponge_x01=-100000.;//start of the sponge layer at the left x boundary
double wave_sponge_x02=100000.;//start of the sponge layer at the right x boundary	
double wave_angle_single=0.;
double wave_K_single=1.0;
double wave_depth=1.0;
double wave_a_single=0.0;	
int wave_ti_start=0;
int wave_IB=0;
	//End variables for wave_momentum_source
	//Variables for air_flow_levelset
int air_flow_levelset=0;
int air_flow_levelset_periodic=0;
int wave_average_k=0;
int wave_skip=0;
int wind_skip=0;
int wind_start_read=0;
int wave_start_read=0;
int wind_recicle=10000;
int wave_recicle=10000;
double wave_wind_reflength=1.0;
double wave_wind_refvel=1.0;
double wave_wind_yshift=0.0;
int floating_turbine_case=0;
int wave_k_ave=0, wave_i_ave=0;
int wave_ti_startave=1000;
int freesurface_wallmodel=0;
int viscosity_wallmodel=0;
double channel_height=1.0;
	//End variables for air_flow_levelset
	//variables for FSI

double dt_inflow;

int levelset=0;
int dam_break=0;
int k_gate=0;
int fix_level=0;
int laplacian=0;
int qcr=0;
int poisson=1;
int amg_agg=1;
double amg_thresh=0.75;
int periodic=0;
int i_periodic=0;
int ii_periodic=0;
int j_periodic=0;
int jj_periodic=0;
int k_periodic=0;
int kk_periodic=0;
int pseudo_periodic=0;
double inlet_flux=-1;
double inlet_sediment_flux=0.;
int delete_previous_file=0;
int mixed=0;
int clark=0;
int vorticity=0;
int initial_perturbation=0;
int initial_perturbation1=0;
int skew=0;
int dynamic_freq=1;
int my_rank;
int ib_bctype[128];
char path[256], gridfile[256];
int i_proc=PETSC_DECIDE, j_proc=PETSC_DECIDE, k_proc=PETSC_DECIDE;
double imp_free_tol=1.e-4;
double poisson_tol=5.e-9;    // relative tolerance
double les_eps=1.e-7;
double mean_pressure_gradient=0;    // relative tolerance
PetscReal max_cs=0.5;
PetscTruth dpdz_set=PETSC_FALSE;
int save_inflow=0;
int save_inflow_period=100;
int save_inflow_minus=0;
PetscInt read_inflow_period=100;

//Add-Hossein for savekplanes
int save_ksection[1000];
int nsave_ksection=0;
int save_jsection[1000];
int nsave_jsection=0;
int save_isection[1000];
int nsave_isection=0;

int ucat_plane_allocated;
int ucat_plane_allocated_is=0;
int ucat_plane_allocated_js=0;
int ucat_plane_allocated_ks=0;

PetscReal scale_velocity=1;

int ucont_plane_allocated=0;
int ti_lastsave=0;
int localstep=1;
int inflow_recycle_perioid=20000;
int save_memory=1;
int ibm_search=0;
PetscTruth rough_set=PETSC_FALSE;
double roughness_size=0.0;
int save_point[3000][10];
double save_coor[6000][10];
int nsave_points=0;
int save_point_level[3000][10];
int nsave_points_level=0;
int testfilter_ik=0;
int testfilter_1d=0;
int i_homo_filter=0;
int j_homo_filter=0;
int k_homo_filter=0;
int poisson_it=10;
int tiout_ufield = -1, tiend_ufield = 10000000;
double dx_min, di_min, dj_min, dk_min;
double di_max, dj_max, dk_max;
double rho_water=1000., rho_air=1.204;    // bubble
double rho_fluid=1000.;
double dthick=1.5;
double seudo_dt=0.01;
PetscTruth dthick_set=PETSC_FALSE;
double mu_water=1.e-3, mu_air=1.78e-5;
double angvel=3.141592;

double gravity_x=0, gravity_y=0, gravity_z=0;
double inlet_x=0, outlet_x=0; //Hossein
double inlet_y=0, outlet_y=0; //Hossein
double inlet_z=0, outlet_z=0; //Hossein
int fix_outlet=0, fix_inlet=0;
int rotdir=2; // 0: rotate around the x-axis, 1:y-axis, 2:z-axis
double x_r=0, y_r=0, z_r=0; // center of rotation of rfsi

char path_inflow[256];
char path_plane_save[256];

PetscInt inflow_levelset=0;

/** Hossein added from turbine structure
 * Blade structure options 
 */

PetscReal 	dt_turbinestructure;

PetscInt 	turbinestructuremodel=0;
PetscInt 	torsion_turbinestructure=0;
PetscInt 	restart_turbinestructure=0;

PetscTruth    rstart_flg=PETSC_FALSE;
PetscTruth    inlet_y_flag=PETSC_FALSE; //Hossein
PetscTruth    inlet_z_flag=PETSC_FALSE; //Hossein


IBMNodes    *ibm_ptr;
PetscInt specifycirculation=0;
PetscInt cf_nacelle_fromfile=0;
PetscReal tipcorrectionratio_Fa=1.0;
PetscReal r_nearhubinflowcorrection=-1.0; // a dimensionless value, using the reference velocity at upstream locations as local inflow velocity for the blade for r<r_nearhubinflowcorrection
PetscInt forcewidthfixed = 0; // fix the width for force distribution
PetscInt forcedistrwidth_surfacescale = 0; // use the surface mesh scale for distributing the force 
PetscReal dhi_fixed, dhj_fixed, dhk_fixed;

PetscInt temperature = 0; 
PetscInt temperature1 = 0; 
PetscInt temperature2 = 0; 
PetscInt temperature3 = 0;

PetscInt maxiteraction_rotormodel=1;
PetscInt NonUniform_ADModel = 0;  
PetscInt NumberOfTurbines=1; 
PetscInt NumberOfActuators=1; 
PetscInt NumberOfNacelle=1; 
PetscInt NumNacellePerLoc = 1;
PetscInt i_periodicIB = 0, j_periodicIB = 0, k_periodicIB = 0;
PetscInt TheSameObject = 1;
PetscInt nacelle_model = 0;  
PetscInt rotor_model = 0;  // 1: actuator disk model 2, 3: actuator line model w Issues, 5: actuator surface, 6: actuator line 
PetscInt UlagrFromSurface = 1;  // the lagrangian velocity for actuator surface model rotor_model=5 is from the lagrangian velocity on the actuator surface
PetscReal indf_ax = 0.25; 
PetscReal reflength_wt = 1.0;  
PetscReal reflength_nacelle = 1.0;  
PetscReal refvel_wt = 1.0;  
PetscReal refvel_cfd = 1.0;  
PetscInt deltafunc = 10;
PetscReal halfwidth_dfunc = 4.0;
PetscReal loc_refvel = 1;
PetscInt les_prt = 0;
PetscReal u_frame, v_frame, w_frame;
PetscInt MoveFrame = 0;
PetscInt ii_periodicWT=0, jj_periodicWT=0, kk_periodicWT=0; // periodic WT, a row/column of ghost wind turbines needs to be added
PetscReal Sx_WT=1.0, Sy_WT=1.0, Sz_WT=1.0;
PetscInt Nx_WT, Ny_WT, Nz_WT;
FSInfo  *fsi_ptr;
UserCtx    *user_ptr;
PetscInt FixTipSpeedRatio=1;
PetscInt fixturbineangvel=0;
PetscInt rstart_turbinerotation=0;
PetscInt turbinetorquecontrol=0;
PetscInt turbineindividualpitchcontrol=0;
PetscInt Shen_AL=0;
PetscReal c0_CL=1;
PetscReal c1_CH=2.2;
PetscReal c2_CH=1;
PetscReal c3_CH=4;
PetscInt correction3D_CH=0;
PetscInt correction3D_CL=0;
PetscInt Shen1_AL=0;
PetscReal correction_ALShen1=1;
PetscReal a_shen=0.125;
PetscReal b_shen=21;
PetscReal c_shen=0.1;
PetscInt Prandtl_AL=0;
PetscInt smoothforce_AL=0;

PetscReal refangle_AL=0.0;
PetscReal count_AL;
PetscInt correction3D_DS=0;
PetscReal a_DS=1.0;
PetscReal b_DS=1.0;
PetscReal d_DS=1.0;

PetscReal CD_0=0.0;
PetscReal AOA_0=0.0;
PetscInt cnpy=0;
PetscInt NumberOfCnpy=0;
PetscReal cnpyHeightDim = 0.;
PetscReal cnpyDispZ = 0;
PetscReal cnpyHeightMlt = 0;
//Meric
PetscReal cnpyCd = 0;

    int file_exist(char *str)
    {
        int r=0;

        if(!my_rank) {
            FILE *fp=fopen(str, "r");
            if(!fp) {
                r=0;
                printf("\nFILE !!! %s does not exist !!!\n", str); 
            }
            else {
                fclose(fp);
                r=1;
            }
        }

        MPI_Bcast(&r, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        return r;
    }

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
    PetscViewer    viewer;
    char filen[90];
    
    PetscInt N;
    VecGetSize(user->Ucont, &N);
    PetscPrintf(PETSC_COMM_WORLD, "PPP - Reading in Restart Data Files %d\n", N);

    sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Ucont));
    PetscViewerDestroy(viewer);

    PetscBarrier(PETSC_NULL);

    PetscOptionsClearValue("-vecload_block_size");

    sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->P));
    PetscViewerDestroy(viewer);

    sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Nvert_o));
    PetscViewerDestroy(viewer);

    sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);    // Seokkoo Kang
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Ucat));
    PetscViewerDestroy(viewer);

    if(!immersed) {
        VecSet(user->Nvert, 0.);
        VecSet(user->Nvert_o, 0.);
    }

    DAGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
    DAGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
    
    DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

    VecCopy(user->Ucont, user->Ucont_o);
    DAGlobalToLocalBegin(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
    DAGlobalToLocalEnd(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);

    DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);

    DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);
    DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat_old);

    DAGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);
    DAGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);

    if(averaging) {    
        sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
        FILE *fp=fopen(filen, "r");

        VecSet(user->Ucat_sum, 0);
        VecSet(user->Ucat_cross_sum, 0);
        VecSet(user->Ucat_square_sum, 0);
        VecSet(user->P_sum, 0);

        if(conv_diff){
            VecSet(user->Conc_sum, 0.);
            VecSet(user->Conc_cross_sum, 0.);
        }

        if(levelset){
            VecSet(user->Level_sum, 0.);
            VecSet(user->Level_square_sum, 0.);
        }

        if(sandwave){
            VecSet(user->lUstar_sum, 0.);
            VecSet(user->lUstar_, 0.);
        }

        if(les) {
            VecSet(user->Nut_sum, 0.);
        }

        if(rans) {
            VecSet(user->K_sum, 0.);
        }

        if(averaging>=2) {
            VecSet(user->P_square_sum, 0);
            //VecSet(user->P_cross_sum, 0);
        }
        if(averaging>=3) {
            if(les) {
                VecSet(user->tauS_sum, 0);
            }
            VecSet(user->Udp_sum, 0);
            VecSet(user->dU2_sum, 0);
            VecSet(user->UUU_sum, 0);
            VecSet(user->Vort_sum, 0);
            VecSet(user->Vort_square_sum, 0);
        }

        if(fp==NULL) {
            PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and contiues the computation ... ***\n\n", filen);
        }
        else {
            fclose(fp);
            PetscBarrier(PETSC_NULL);
            sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, user->Ucat_sum);
            PetscViewerDestroy(viewer);

            sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, user->Ucat_cross_sum);
            PetscViewerDestroy(viewer);

            sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, (user->Ucat_square_sum));
            PetscViewerDestroy(viewer);

            sprintf(filen, "%s/sp_%06d_%1d.dat", path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, user->P_sum);
            PetscViewerDestroy(viewer);

            if(sandwave) {
                sprintf(filen, "%s/sustar_%06d_%1d.dat", path, ti, user->_this);
                if( file_exist(filen) ) {
                    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                    VecLoadIntoVector(viewer, user->lUstar_sum);
                    PetscViewerDestroy(viewer);
                }
            }

            if(conv_diff) {
                sprintf(filen, "%s/sconc_%06d_%1d.dat", path, ti, user->_this);
                if( file_exist(filen) ) {
                    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                    VecLoadIntoVector(viewer, user->Conc_sum);
                    PetscViewerDestroy(viewer);
                }

                sprintf(filen, "%s/sconc_cross_sum_%06d_%1d.dat", path, ti, user->_this);
                if( file_exist(filen) ) {
                    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                    VecLoadIntoVector(viewer, user->Conc_cross_sum);
                    PetscViewerDestroy(viewer);
                }
            }

            if(levelset) {
                sprintf(filen, "%s/slevel_%06d_%1d.dat", path, ti, user->_this);
                if( file_exist(filen) ) {
                    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                    VecLoadIntoVector(viewer, user->Level_sum);
                    PetscViewerDestroy(viewer);
                }
                sprintf(filen, "%s/slevel2_%06d_%1d.dat", path, ti, user->_this);
                if( file_exist(filen) ) {
                    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                    VecLoadIntoVector(viewer, user->Level_square_sum);
                    PetscViewerDestroy(viewer);
                }
            }

            if(les) {
                sprintf(filen, "%s/snut_%06d_%1d.dat", path, ti, user->_this);
                if( file_exist(filen) ) {
                    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                    VecLoadIntoVector(viewer, user->Nut_sum);
                    PetscViewerDestroy(viewer);
                }
            }

            if(rans) {
                sprintf(filen, "%s/sk_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->K_sum);
                PetscViewerDestroy(viewer);
            }

            if(averaging>=2) {
                sprintf(filen, "%s/sp2_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->P_square_sum);
                PetscViewerDestroy(viewer);
            }

            if(averaging>=3) {
                if(les) {
                    sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti, user->_this);
                    if( file_exist(filen) ) {
                        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                        VecLoadIntoVector(viewer, user->tauS_sum);
                        PetscViewerDestroy(viewer);
                    }
                    else PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s !\n", filen);
                }

                sprintf(filen, "%s/su3_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->Udp_sum);
                PetscViewerDestroy(viewer);

                sprintf(filen, "%s/su4_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->dU2_sum);
                PetscViewerDestroy(viewer);

                sprintf(filen, "%s/su5_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->UUU_sum);
                PetscViewerDestroy(viewer);

                sprintf(filen, "%s/svo_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->Vort_sum);
                PetscViewerDestroy(viewer);
                
                sprintf(filen, "%s/svo2_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer, user->Vort_square_sum);
                PetscViewerDestroy(viewer);
            }
            
            PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
        }
    }

    if(levelset) {
        
        sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti, user->_this);
        FILE *fp=fopen(filen, "r");

        if(fp==NULL) {
            PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, terminates ... ***\n\n", filen);
            PetscFinalize();
            exit(0);
        }
        else {
            fclose(fp);
        
            PetscBarrier(PETSC_NULL);
            
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, user->Levelset);
            PetscViewerDestroy(viewer);
            
            VecCopy (user->Levelset, user->Levelset_o);
            DAGlobalToLocalBegin(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
            DAGlobalToLocalEnd(user->da, user->Levelset, INSERT_VALUES, user->lLevelset);
        }
        
    }

    if(les) {
        Vec Cs;
        VecDuplicate(user->P, &Cs);
        
        sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti, user->_this);
        FILE *fp=fopen(filen, "r");

        if(fp==NULL) {
            PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting Cs to 0 and contiues the computation ... ***\n\n", filen);
            VecSet(Cs, 0);
        }
        else {
            fclose(fp);
            PetscBarrier(PETSC_NULL);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, Cs);
            PetscViewerDestroy(viewer);
        }

        DAGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
        DAGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);

        VecDestroy(Cs);
    }

    if(rans) {
        // K-Omega
        sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti, user->_this);
        if( file_exist(filen) ) {
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer,user->K_Omega);
            PetscViewerDestroy(viewer);
        }
        else {
            PetscPrintf(PETSC_COMM_WORLD, "\nInitializing K-omega ... \n\n");
            //PetscBarrier(PETSC_NULL);
            //PetscPrintf(PETSC_COMM_WORLD, "Debug1\n");
            K_Omega_IC(user);
            //PetscBarrier(PETSC_NULL);
            //PetscPrintf(PETSC_COMM_WORLD, "Debug2\n");
            VecSet(user->lNu_t, user->ren);
        }

        VecCopy(user->K_Omega, user->K_Omega_o);

        DAGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
        DAGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);

        DAGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
        DAGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);

        if(rans==3) {
            sprintf(filen, "%s/Distance_%1d.dat", path, user->_this);
            if( file_exist(filen) ) {
                PetscPrintf(PETSC_COMM_WORLD, "Reading %s !\n", filen);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
                VecLoadIntoVector(viewer,user->Distance);
                PetscViewerDestroy(viewer);
            }
            else {
            PetscPrintf(PETSC_COMM_WORLD, "File %s does not exist. Recalculating distance function !\n", filen);
            }
        }
    }

    if(conv_diff) {
        sprintf(filen, "%s/cfield%06d_%1d.dat", path, ti, user->_this);
        if( file_exist(filen) ) {
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
            VecLoadIntoVector(viewer, user->Conc);
            PetscViewerDestroy(viewer);
        }
        
        else {
            PetscPrintf(PETSC_COMM_WORLD,"\n\n Cannot open %s, terminates ... So Will Do Initialization...\n\n", filen);
            PetscPrintf(PETSC_COMM_WORLD, "\nInitializing Conv_Diffusion ... \n\n");
            Conv_Diff_IC(user);
        }    

            VecCopy (user->Conc, user->Conc_o);
            DAGlobalToLocalBegin(user->da, user->Conc, INSERT_VALUES, user->lConc);
            DAGlobalToLocalEnd(user->da, user->Conc, INSERT_VALUES, user->lConc);
            DAGlobalToLocalBegin(user->da, user->Conc_o, INSERT_VALUES, user->lConc_o);
            DAGlobalToLocalEnd(user->da, user->Conc_o, INSERT_VALUES, user->lConc_o);
        }

    return 0;
}

PetscErrorCode Ucat_Binary_Output(UserCtx *user)
{
    PetscViewer    viewer;
    char filen[80];
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat, viewer);
    PetscViewerDestroy(viewer);
    sprintf(filen, "%s/ufield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);

    return 0;
}

int delete_count=0;

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user)
{
    PetscViewer    viewer;
    char filen[80];
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucont, viewer);
    PetscViewerDestroy(viewer);
    sprintf(filen, "%s/vfield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
    PetscBarrier(PETSC_NULL);

    sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat, viewer);
    PetscViewerDestroy(viewer);
    sprintf(filen, "%s/ufield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
    PetscBarrier(PETSC_NULL);

    sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->P, viewer);
    PetscViewerDestroy(viewer);
    sprintf(filen, "%s/pfield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
    PetscBarrier(PETSC_NULL);
    
    if(qcr) {
        Vec Q;
        VecDuplicate(user->P, &Q);
        Compute_Q(user,  Q);
        
        sprintf(filen, "%s/qfield%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(Q, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/qfield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
        
        VecDestroy(Q);
        PetscBarrier(PETSC_NULL);
    }

    sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Nvert, viewer);
    PetscViewerDestroy(viewer);
    sprintf(filen, "%s/nvfield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
    PetscBarrier(PETSC_NULL);

    if(averaging && ti!=0) {
        sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->Ucat_sum, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/su0_%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
        PetscBarrier(PETSC_NULL);

        sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->Ucat_cross_sum, viewer);
        PetscViewerDestroy(viewer);  
        sprintf(filen, "%s/su1_%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
        PetscBarrier(PETSC_NULL);
        
        sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->Ucat_square_sum, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/su2_%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
        PetscBarrier(PETSC_NULL);

        sprintf(filen, "%s/sp_%06d_%1d.dat",path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->P_sum, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/sp_%06d_%1d.dat.info",path, ti, user->_this);    if(!rank) unlink(filen);
        
        if(les) {
            sprintf(filen, "%s/snut_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Nut_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/snut_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }

        if(rans) {
            sprintf(filen, "%s/sk_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->K_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/sk_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }

        if(sandwave) {
            sprintf(filen, "%s/sustar_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->lUstar_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/sustar_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }
        if(conv_diff) {
            sprintf(filen, "%s/sconc_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Conc_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/sconc_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);

            sprintf(filen, "%s/sconc_cross_sum_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Conc_cross_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/sconc_cross_sum_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }

        if(levelset) {
            sprintf(filen, "%s/slevel_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Level_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/slevel_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
            sprintf(filen, "%s/slevel2_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Level_square_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/slevel2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }

        if(averaging>=2) {
            sprintf(filen, "%s/sp2_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->P_square_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/sp2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }

        if(averaging>=3) {
            if(les) {
                sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti, user->_this);
                PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
                VecView(user->tauS_sum, viewer);
                PetscViewerDestroy(viewer);
                sprintf(filen, "%s/stauS_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);
            }

            sprintf(filen, "%s/su3_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Udp_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/su3_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

            sprintf(filen, "%s/su4_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->dU2_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/su4_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

            sprintf(filen, "%s/su5_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->UUU_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/su5_%06d_%1d.dat.info",path, ti, user->_this);if(!rank) unlink(filen);

            sprintf(filen, "%s/svo_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Vort_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/svo_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
            
            sprintf(filen, "%s/svo2_%06d_%1d.dat",path, ti, user->_this);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
            VecView(user->Vort_square_sum, viewer);
            PetscViewerDestroy(viewer);
            sprintf(filen, "%s/svo2_%06d_%1d.dat.info",path, ti, user->_this);        if(!rank) unlink(filen);
        }
        
        PetscBarrier(PETSC_NULL);
    }
    
    if(levelset) {
        sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti, user->_this);
        
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->Levelset, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/lfield%06d_%1d.dat.info",path, ti, user->_this);    if(!rank) unlink(filen);
    }
    
    if(les) {
        Vec Cs;
        
        VecDuplicate(user->P, &Cs);
        DALocalToGlobal(user->da, user->lCs, INSERT_VALUES, Cs);
        
        sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(Cs, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/cs_%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
        
        PetscBarrier(PETSC_NULL);
        VecDestroy(Cs);
    }
    
    if(rans) {
        sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->K_Omega, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/kfield%06d_%1d.dat.info", path, ti, user->_this);    if(!rank) unlink(filen);
        
        PetscBarrier(PETSC_NULL);
    }

    if(sandwave) {
        sprintf(filen, "%s/ustarfield%06d_%1d.dat", path, ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->lUstar_, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/ustarfield%06d_%1d.dat.info",path, ti, user->_this);    if(!rank) unlink(filen);
    }
    

    if(conv_diff) {
        sprintf(filen, "%s/cfield%06d_%1d.dat", path, ti, user->_this);
        
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
        VecView(user->Conc, viewer);
        PetscViewerDestroy(viewer);
        sprintf(filen, "%s/cfield%06d_%1d.dat.info",path, ti, user->_this);    if(!rank) unlink(filen);
    }

    if(!rank && delete_previous_file && delete_count++>=2 && ti-tiout*2!=0) {
        sprintf(filen, "%s/vfield%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
        
        if(!(tiout_ufield>0 && ti == (ti/tiout_ufield) * tiout_ufield && ti<=tiend_ufield)) {
            sprintf(filen, "%s/ufield%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
        }
        sprintf(filen, "%s/pfield%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
        sprintf(filen, "%s/nvfield%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
        if(averaging) {
            sprintf(filen, "%s/sp_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
            sprintf(filen, "%s/su0_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
            sprintf(filen, "%s/su1_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
            sprintf(filen, "%s/su2_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
            if(averaging>=2) {
              //sprintf(filen, "%s/sp1_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
                sprintf(filen, "%s/sp2_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
            }
            if(averaging>=3) {
                sprintf(filen, "%s/su3_%06d_%1d.dat", path, ti-tiout*2, user->_this);   if(!rank) unlink(filen);
                sprintf(filen, "%s/su4_%06d_%1d.dat", path, ti-tiout*2, user->_this);   if(!rank) unlink(filen);
                sprintf(filen, "%s/su5_%06d_%1d.dat", path, ti-tiout*2, user->_this);   if(!rank) unlink(filen);
                sprintf(filen, "%s/svo_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
                sprintf(filen, "%s/svo2_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
            }
            if(les) {
                sprintf(filen, "%s/stauS_%06d_%1d.dat", path, ti-tiout*2, user->_this); if(!rank) unlink(filen);
                sprintf(filen, "%s/snut_%06d_%1d.dat", path, ti-tiout*2, user->_this); if(!rank) unlink(filen);
            }
        }
        if(les>=2) {
            sprintf(filen, "%s/cs_%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
        }

        if(rans) {
            sprintf(filen, "%s/kfield%06d_%1d.dat", path, ti-tiout*2, user->_this);    if(!rank) unlink(filen);
        }

        if(levelset) {
            sprintf(filen, "%s/lfield%06d_%1d.dat", path, ti-tiout*2, user->_this); if(!rank) unlink(filen);
        }
    }
    PetscBarrier(PETSC_NULL);

    return 0;
}

PetscErrorCode Divergence(UserCtx *user)
{
    DA        da = user->da, fda = user->fda;
    DALocalInfo    info = user->info;

    PetscInt    xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt    zs = info.zs, ze = info.zs + info.zm;
    PetscInt    mx = info.mx, my = info.my, mz = info.mz;

    PetscInt    lxs, lys, lzs, lxe, lye, lze;
    PetscInt    i, j, k;

    Vec    Div;
    PetscReal    ***div, ***aj, ***nvert;
    Cmpnts    ***ucont;
    PetscReal    maxdiv;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    DAVecGetArray(fda,user->lUcont, &ucont);
    DAVecGetArray(da, user->lAj, &aj);
    VecDuplicate(user->P, &Div);
    DAVecGetArray(da, Div, &div);
    DAVecGetArray(da, user->lNvert, &nvert);
    for (k=lzs; k<lze; k++) {
        for (j=lys; j<lye; j++) {
            for (i=lxs; i<lxe; i++) {
    
                maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x + ucont[k][j][i].y - ucont[k][j-1][i].y + ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
                if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] + nvert[k][j+1][i] + nvert[k][j-1][i] + nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
                div[k][j][i] = maxdiv;
            }
        }
    }

    if (zs==0) {
        k=0;
        for (j=ys; j<ye; j++) {
            for (i=xs; i<xe; i++) {
                div[k][j][i] = 0.;
            }
        }
    }

    if (ze == mz) {
        k=mz-1;
        for (j=ys; j<ye; j++) {
            for (i=xs; i<xe; i++) {
                div[k][j][i] = 0.;
            }
        }
    }

    if (xs==0) {
        i=0;
        for (k=zs; k<ze; k++) {
            for (j=ys; j<ye; j++) {
                div[k][j][i] = 0.;
            }
        }
    }

    if (xe==mx) {
        i=mx-1;
        for (k=zs; k<ze; k++) {
            for (j=ys; j<ye; j++) {
                div[k][j][i] = 0;
            }
        }
    }

    if (ys==0) {
        j=0;
        for (k=zs; k<ze; k++) {
            for (i=xs; i<xe; i++) {
                div[k][j][i] = 0.;
            }
        }
    }

    if (ye==my) {
        j=my-1;
        for (k=zs; k<ze; k++) {
            for (i=xs; i<xe; i++) {
                div[k][j][i] = 0.;
            }
        }
    }
    DAVecRestoreArray(da, Div, &div);
    VecMax(Div, &i, &maxdiv);
    PetscPrintf(PETSC_COMM_WORLD, "Maxdiv %d %d %e\n", ti, i, maxdiv);
    PetscInt mi;

    for (k=zs; k<ze; k++) {
        for (j=ys; j<ye; j++) {
            for (mi=xs; mi<xe; mi++) {
                if (lidx(mi,j,k,user) ==i) {
                    PetscPrintf(PETSC_COMM_SELF, "MMa %d %d %d\n", mi,j, k);
                }
            }
        }
    }

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (!rank) {
        FILE *f;
        char filen[80];
        sprintf(filen, "%s/Converge_dU", path);
        f = fopen(filen, "a");
        PetscFPrintf(PETSC_COMM_WORLD, f, " Maxdiv=%.2e\n", maxdiv);
        fclose(f);
    }

    DAVecRestoreArray(da, user->lNvert, &nvert);
    DAVecRestoreArray(fda, user->lUcont, &ucont);
    DAVecRestoreArray(da, user->lAj, &aj);
    VecDestroy(Div);
    return(0);
}

void write_data(UserCtx *user)
{
    DA    da = user->da, fda = user->fda;
    DALocalInfo    info = user->info;

    PetscInt    xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt    zs = info.zs, ze = info.zs + info.zm;
    PetscInt    mx = info.mx, my = info.my, mz = info.mz;

    PetscInt    lxs, lys, lzs, lxe, lye, lze;
    PetscInt    i, j, k;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    double lvol=0, vol=0;

    PetscReal ***p, ***aj, ***nvert, ***rho, ***level;
    Cmpnts    ***ucat, ***ucont, ***csi, ***eta, ***zet, ***cent;
  
    if(levelset) {
        DAVecGetArray(da, user->lDensity, &rho);
        DAVecGetArray(da, user->lLevelset, &level);
    }
    DAVecGetArray(da,user->lAj, &aj);
    DAVecGetArray(da, user->lP, &p);
    DAVecGetArray(da, user->lNvert, &nvert);
    DAVecGetArray(fda,user->lCsi, &csi);
    DAVecGetArray(fda,user->lEta, &eta);
    DAVecGetArray(fda,user->lZet, &zet);
    DAVecGetArray(fda,user->lUcat, &ucat);
    DAVecGetArray(fda,user->lUcont, &ucont);
    DAVecGetArray(fda,user->lCent, &cent);

    if (ti==tistart && nsave_points>0) {



		int *iclose, *jclose, *kclose;

		int *sum_iclose, *sum_jclose, *sum_kclose;



		iclose= (int *) malloc(nsave_points*sizeof(int));

		jclose= (int *) malloc(nsave_points*sizeof(int));

		kclose= (int *) malloc(nsave_points*sizeof(int));



		sum_iclose= (int *) malloc(nsave_points*sizeof(int));

		sum_jclose= (int *) malloc(nsave_points*sizeof(int));

		sum_kclose= (int *) malloc(nsave_points*sizeof(int));



		for(int m=0; m<nsave_points; m++) {

			iclose[m]=0;

			jclose[m]=0;

			kclose[m]=0;			

		}



	  	PetscPrintf(PETSC_COMM_WORLD, "save %d flow points\n\n", nsave_points);

		int count[nsave_points], sum_count[nsave_points];

		for(int m=0; m<nsave_points; m++) {

			double XX=save_coor[m][0];

			double YY=save_coor[m][1];

			double ZZ=save_coor[m][2];



			double dis, dis_min;

			dis_min=1.e9;



			for (k=zs; k<ze; k++)

			for (j=ys; j<ye; j++)

			for (i=xs; i<xe; i++) {

				dis=pow(XX-cent[k][j][i].x,2)+pow(YY-cent[k][j][i].y,2)+pow(ZZ-cent[k][j][i].z,2);



				if (dis<dis_min) {

					dis_min=dis;

					iclose[m]=i;

					jclose[m]=j;

					kclose[m]=k;

				}

			}





	

			double dmin_global;

			MPI_Allreduce (&dis_min, &dmin_global, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);

			count[m]=1;

			double diff=fabs(dis_min-dmin_global);

			if (diff>1.e-9) {

				count[m]=0;

				iclose[m]=0; jclose[m]=0; kclose[m]=0;

			}



		}



	  	PetscBarrier(PETSC_NULL);



		MPI_Allreduce (&iclose[0], &sum_iclose[0], nsave_points, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&jclose[0], &sum_jclose[0], nsave_points, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&kclose[0], &sum_kclose[0], nsave_points, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

		MPI_Allreduce (&count[0], &sum_count[0], nsave_points, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);





	        for (int m=0; m<nsave_points; m++) {

			save_point[m][0]=sum_iclose[m]/sum_count[m];

			save_point[m][1]=sum_jclose[m]/sum_count[m];

			save_point[m][2]=sum_kclose[m]/sum_count[m];

	  		PetscPrintf(PETSC_COMM_WORLD, "save flow points at x=%le y=%le z=%le\n", save_coor[m][0],  save_coor[m][1],  save_coor[m][2] );

	  		PetscPrintf(PETSC_COMM_WORLD, "save flow points at i=%d j=%d k=%d\n", save_point[m][0],  save_point[m][1],  save_point[m][2] );

		}





		free(iclose);

		free(jclose);

		free(kclose);



		free(sum_iclose);

		free(sum_jclose);

		free(sum_kclose);



	}
    
    // for sloshing recording
    int ci=mx/2, ck=(mz-3)/2 + 1;
    double lz_sloshing=-100, z_sloshing;
    
    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++) {
        if(nvert[k][j][i]>0.1) continue;
        if(levelset) {
            lvol += rho[k][j][i] / aj[k][j][i];
            if ( i==ci && k==ck) {
                if( level[k][j][i]>=0 && level[k][j+1][i]<0 ) {    // water surface is above my cell center
                    lz_sloshing = cent[k][j][i].z + level[k][j][i];
                }
                else if( level[k][j][i]<0 && level[k][j-1][i]>=0 ) {    // water surface is below my cell center
                    lz_sloshing = cent[k][j][i].z + level[k][j][i];
                }
            }
        }

        for(int m=0; m<nsave_points_level; m++) {
            if(i==save_point_level[m][0] && j==save_point_level[m][1] && k==save_point_level[m][2]) {
                FILE *f;
                char filen[80];

                sprintf(filen, "%s/level_%04d_%04d_%04d_dt_%g.dat", path, i, j, k, user->dt);
                f = fopen(filen, "a");
                if(ti==tistart) fprintf(f, "\n");

                fprintf(f, "%d %.7e %.7e %.7e %.7e\n", ti, cent[k][j][i].x, cent[k][j][i].y, cent[k][j][i].z, level[k][j][i]);
                fclose(f);

                break;
            }
        }

        for(int m=0; m<nsave_points; m++) {
            if(i==save_point[m][0] && j==save_point[m][1] && k==save_point[m][2]) {
                double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
                double dpdc, dpde, dpdz;
                double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
                double dp_dx, dp_dy, dp_dz;
                double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
                double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
                double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
                double ajc = aj[k][j][i];
                double Ai = sqrt ( csi0*csi0 + csi1*csi1 + csi2*csi2 );
                double Aj = sqrt ( eta0*eta0 + eta1*eta1 + eta2*eta2 );
                double Ak = sqrt ( zet0*zet0 + zet1*zet1 + zet2*zet2 );
                double U = 0.5*(ucont[k][j][i].x+ucont[k][j][i-1].x) / Ai;
                double V = 0.5*(ucont[k][j][i].y+ucont[k][j-1][i].y) / Aj;
                double W = 0.5*(ucont[k][j][i].z+ucont[k-1][j][i].z) / Ak;
            
                Compute_du_center ( i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
                Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

                Compute_dscalar_center ( i, j, k, mx, my, mz, p, nvert, &dpdc, &dpde, &dpdz );
                Compute_dscalar_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx, &dp_dy, &dp_dz );

                double vort_x = dw_dy - dv_dz,  vort_y = du_dz - dw_dx, vort_z = dv_dx - du_dy;

                FILE *f;
                char filen[80];

                sprintf(filen, "%s/Flow0_%.2e_%.2e_%.2e_dt_%g.dat", path, save_coor[m][0], save_coor[m][1], save_coor[m][2], user->dt);
                f = fopen(filen, "a");
                if(ti==tistart) fprintf(f, "\n");
                fprintf(f, "%d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n", ti, ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z, p[k][j][i], U, V, W, vort_x, vort_y, vort_z);
                fclose(f);
                
                sprintf(filen, "%s/Flow1_%.2e_%.2e_%.2e_dt_%g.dat", path, save_coor[m][0], save_coor[m][1], save_coor[m][2], user->dt);
                f = fopen(filen, "a");
                if(ti==tistart) fprintf(f, "\n");
                fprintf(f, "%d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n", ti, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dp_dx, dp_dy, dp_dz);
                fclose(f);
                
                break;
            }
        }
    }

    if(levelset) {
        PetscBarrier(PETSC_NULL);
        PetscGlobalSum(&lvol, &vol, PETSC_COMM_WORLD);
        PetscGlobalMax(&lz_sloshing, &z_sloshing, PETSC_COMM_WORLD);
        if(!my_rank) {
            char filen[256];
            sprintf(filen, "%s/mass.dat", path);
            FILE *fp = fopen(filen, "a");
            fprintf(fp, "%d %.10e\n", ti, vol );
            fclose(fp);

            double _time = (ti+1)*user->dt;
            double a = 0.05, b = 2.0, d = 1.0, g = 1., x = 1.0;
            double k2 = 2*M_PI/b, k4 = 4*M_PI/b;
            double w2 = sqrt ( k2 * g * tanh (k2*d) ), w4 = sqrt ( k4 * g * tanh (k4*d) );

            double z_sloshing_exact = a * cos(w2*_time) * cos(k2*x);
            z_sloshing_exact  += 0.125/g * ( 2*pow(w2*a,2.) * cos(2.*w2*_time) + pow(a/w2,2.) * ( pow(k2*g,2.) + pow(w2,4.) ) - pow(a/w2,2.) * ( pow(k2*g,2.) + 3.*pow(w2,4.) ) * cos(w4*_time) );
            
            sprintf(filen, "%s/sloshing.dat", path);
            fp = fopen(filen, "a");
            fprintf(fp, "%e %.7e %.7e %.7e\n", _time, z_sloshing,  -1.0+z_sloshing, z_sloshing_exact);
            fclose(fp);
        }
    }

    PetscBarrier(PETSC_NULL);
    if(user->bctype[0]==11 && ti) {
        PetscGlobalSum(&user->lA_cyl, &user->A_cyl, PETSC_COMM_WORLD);
        PetscGlobalSum(&user->lA_cyl_x, &user->A_cyl_x, PETSC_COMM_WORLD);
        PetscGlobalSum(&user->lA_cyl_z, &user->A_cyl_z, PETSC_COMM_WORLD);
        
        PetscGlobalSum(&user->lFpx_cyl, &user->Fpx_cyl, PETSC_COMM_WORLD);
        PetscGlobalSum(&user->lFpz_cyl, &user->Fpz_cyl, PETSC_COMM_WORLD);
        PetscGlobalSum(&user->lFvx_cyl, &user->Fvx_cyl, PETSC_COMM_WORLD);
        PetscGlobalSum(&user->lFvz_cyl, &user->Fvz_cyl, PETSC_COMM_WORLD);

        double Cpx=user->Fpx_cyl/user->A_cyl_x*4., Cpz=user->Fpz_cyl/user->A_cyl_z*4.;
        double Cvx=user->Fvx_cyl/user->A_cyl_x*4., Cvz=user->Fvz_cyl/user->A_cyl_z*4.;
        double Cdx=Cpx+Cvx, Cdz=Cpz+Cvz;

        if(!my_rank) {
            char filen[256];
            sprintf(filen, "%s/Force_Cylinder.dat", path);
            FILE *fp = fopen(filen, "a");
            if(ti==tistart) fprintf(fp, "\n\n***\n\n");
            fprintf(fp, "%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", ti, Cdx, 0., Cdz, Cpx, 0., Cpz, user->A_cyl_x, 0., user->A_cyl_z, user->A_cyl );
            fclose(fp);
        }
    }
    PetscBarrier(PETSC_NULL);
    
    if(/*inletprofile==13 &&*/ ti) {
        
        if(!my_rank) {
            char filen[256];
            sprintf(filen, "%s/shear_velocity.dat", path);
            FILE *fp = fopen(filen, "a");
            //if(ti==tistart) fprintf(fp, "\n\n***\n\n");
            fprintf(fp, "%d ", ti);
            fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e\n", 
                user->ustar_now[0],user->ustar_now[1],user->ustar_now[2],user->ustar_now[3],user->ustar_now[4],user->ustar_now[5]);
            fclose(fp);
        }
    }

    if(levelset) {
        DAVecRestoreArray(da, user->lDensity, &rho);
        DAVecRestoreArray(da, user->lLevelset, &level);
    }
    DAVecRestoreArray(da,user->lAj, &aj);
    DAVecRestoreArray(da, user->lP, &p);
    DAVecRestoreArray(da, user->lNvert, &nvert);
    DAVecRestoreArray(fda,user->lCsi, &csi);
    DAVecRestoreArray(fda,user->lEta, &eta);
    DAVecRestoreArray(fda,user->lZet, &zet);
    DAVecRestoreArray(fda,user->lUcat, &ucat);
    DAVecRestoreArray(fda,user->lUcont, &ucont);
    DAVecRestoreArray(fda,user->lCent, &cent);
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{

    Vec    ResidualT;
    UserCtx    *user;

    PetscErrorCode    ierr;
    PetscInt    i, bi, ibi;

    PetscReal    norm;

    IBMNodes    *ibm, *ibm0, *ibm1;
    IBMInfo    *ibm_intp;

    // Added for fsi
    FSInfo        *fsi;
    PetscTruth    DoSCLoop;
    PetscInt      itr_sc;
    FSISetup      *fsimov;
    PetscInt level;
    UserMG usermg;

    PetscTruth  flg;

    IBMNodes    *wtm;   
    FSInfo      *fsi_wt;

    IBMNodes    *ibm_ACD;

    IBMNodes    *ibm_acl2ref; 
    FSInfo      *fsi_acl2ref; 

    IBMNodes    *ibm_nacelle; 
    FSInfo      *fsi_nacelle; 

    PetscInitialize(&argc, &argv, (char *)0, help);
    PetscBarrier(PETSC_NULL);
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
  

    srand( time(NULL)) ;// Seokkoo Kang

    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
    PetscPrintf(PETSC_COMM_WORLD, "\nControl.dat file error code is %i\n", ierr); //Debug by Kevin Flora

    if(PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE)==0)

    PetscPrintf(PETSC_COMM_WORLD, "\nControl.dat file was read\n"); //Debug by Kevin Flora
        else
    PetscPrintf(PETSC_COMM_WORLD, "\nControl.dat file was not read\n"); 

    PetscOptionsGetInt(PETSC_NULL, "-tio", &tiout, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-tiou", &tiout_ufield, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-tieu", &tiend_ufield, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-paraview", &Paraview, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-imm", &immersed, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-inv", &inviscid, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &rstart_flg);
    PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-imp_MAX_IT", &imp_MAX_IT, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-fsi", &movefsi, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-rfsi", &rotatefsi, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-radi", &radi, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-inlet", &inletprofile, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-case", &inletCase, PETSC_NULL); //Hossein-1/21/2025 (Lehigh-SBU project) 
   	PetscOptionsGetInt(PETSC_NULL, "-inlet_tmprt", &inletprofile_tmprt, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-str", &STRONG_COUPLING, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-rs_fsi", &rstart_fsi, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-cop", &cop, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);

    	// Hossein added from turbine structure
	
  	PetscOptionsGetInt(PETSC_NULL, "-turbinestructuremodel", &turbinestructuremodel, PETSC_NULL);
  	PetscOptionsGetInt(PETSC_NULL, "-torsion_turbinestructure", &torsion_turbinestructure, PETSC_NULL);
  	PetscOptionsGetInt(PETSC_NULL, "-restart_turbinestructure", &restart_turbinestructure, PETSC_NULL);

//  	PetscOptionsGetReal(PETSC_NULL, "-coeff_SG", &coeff_SG, PETSC_NULL);

	// add end (xiaolei)
    
    PetscOptionsGetInt(PETSC_NULL, "-sediment", &sediment, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-sediment_influx", &sediment_influx, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-English", &English, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-sed_tiout", &sed_tio, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-simple_cell_check", &SimpleCellCheck, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-osl_inlet_sediment_flux", &osl_inlet_sediment_flux, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-Barchans", &Barchans, PETSC_NULL); 
    PetscOptionsGetInt(PETSC_NULL, "-smooth_shear", &smooth_shear, PETSC_NULL); 
    PetscOptionsGetReal(PETSC_NULL, "-sediment_thickness", &sediment_thickness, &sediment_thickness_flag); //ali
    PetscOptionsGetReal(PETSC_NULL, "-x_outlet", &x_outlet, &XOutlet); //ali
    PetscOptionsGetReal(PETSC_NULL, "-y_outlet", &y_outlet, &YOutlet); //ali
    PetscOptionsGetReal(PETSC_NULL, "-x_inlet", &x_inlet, &XInlet); //ali
    PetscOptionsGetReal(PETSC_NULL, "-y_inlet", &y_inlet, &YInlet); //ali
    PetscOptionsGetReal(PETSC_NULL, "-x_limit_inlet", &x_limit_inlet, &X_Limit_Inlet); //ali
    PetscOptionsGetReal(PETSC_NULL, "-y_limit_inlet", &y_limit_inlet, &Y_Limit_Inlet); //ali
    PetscOptionsGetInt(PETSC_NULL, "-line_outlet", &LOutlet, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-line_inlet", &LInlet, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-mm_l", &mm_l, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-bb_l", &bb_l, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-mm_l_in", &mm_l_in, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-bb_l_in", &bb_l_in, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-cell_size", &cell_size, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-cell_depth", &cell_depth, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-non_dimensional", &non_dimensional, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-initial_bed_elevation", &initial_bed_elevation, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-max_mobilebed_z", &max_mobilebed_z, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-min_mobilebed_x", &min_mobilebed_x, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-min_mobilebed_y", &min_mobilebed_y, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-periodic_morpho", &periodic_morpho, PETSC_NULL); //ali
    PetscOptionsGetReal(PETSC_NULL, "-Nseg", &Nseg, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-livebed", &LiveBed, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-effective_bed_shear_stress", &effective_bed_shear_stress, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-aval_loop_num", &aval_loop, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-y_direction", &y_direction, PETSC_NULL); //Hossein
    PetscOptionsGetInt(PETSC_NULL, "-sand_slide", &sand_slide, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-rigidbed", &RigidBed, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-zero_grad", &zero_grad, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-projection_method", &projection_method, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-bed_roughness", &bed_roughness, PETSC_NULL);    // ali: 0 or 1
    PetscOptionsGetReal(PETSC_NULL, "-k_ss", &k_ss, PETSC_NULL);    // ali
    PetscOptionsGetInt(PETSC_NULL, "-ib_depth", &input_ib_depth, PETSC_NULL); //ali
    PetscOptionsGetInt(PETSC_NULL, "-mhv", &MHV, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-reg", &regime, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-twoD", &TwoD, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-thin", &thin, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dgf_z", &dgf_z, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dgf_y", &dgf_y, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dgf_x", &dgf_x, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dgf_az", &dgf_az, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dgf_ay", &dgf_ay, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dgf_ax", &dgf_ax, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-body", &NumberOfBodies, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL);    // Seokkoo Kang: if 1 do averaging; always begin with -rstart 0
    PetscOptionsGetInt(PETSC_NULL, "-binary", &binary_input, PETSC_NULL);    // Seokkoo Kang: if 1 binary PLOT3D file, if 0 ascii.
    PetscOptionsGetInt(PETSC_NULL, "-xyz", &xyz_input, PETSC_NULL);            // Seokkoo Kang: if 1 text xyz format, useful for very big (>1GB) Cartesian grid
    PetscOptionsGetInt(PETSC_NULL, "-les", &les, PETSC_NULL);                // Seokkoo Kang: if 1 Smagorinsky with Cs=0.1, if 2 Dynamic model
    PetscOptionsGetInt(PETSC_NULL, "-inlet_buffer_k", &inlet_buffer_k, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-wallfunction", &wallfunction, PETSC_NULL);    // Seokkoo Kang: 1 or 2
    PetscOptionsGetInt(PETSC_NULL, "-slipbody", &slipbody, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-central", &central, PETSC_NULL);//central differencing
    PetscOptionsGetInt(PETSC_NULL, "-second_order", &second_order, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-weno", &weno, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-initialzero", &initialzero, PETSC_NULL);    // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-freesurface", &freesurface, PETSC_NULL);    // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);            // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-conv_diff", &conv_diff, PETSC_NULL);        // Ali Khosro
    PetscOptionsGetInt(PETSC_NULL, "-no_ibm_search", &no_ibm_search, PETSC_NULL);        // Ali Khosro
    PetscOptionsGetInt(PETSC_NULL, "-sandwave", &sandwave, PETSC_NULL);        // Ali Khosro
    PetscOptionsGetInt(PETSC_NULL, "-density_current", &density_current, PETSC_NULL);        // Ali Khosro
    PetscOptionsGetInt(PETSC_NULL, "-SuspendedParticles", &SuspendedParticles, PETSC_NULL);    // Ali Khosro
    PetscOptionsGetInt(PETSC_NULL, "-mobile_bed", &mobile_bed, PETSC_NULL);        // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-w_s", &w_s, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-background_conc", &Background_Conc, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-inlet_conc", &Inlet_Conc, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-porosity", &porosity, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-sbb", &sbbb, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-Cs", &Cs_, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-Angle_repose", &Angle_repose, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-u_bulk", &U_Bulk, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-d50", &d50, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-deltab", &deltab, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-sed_density", &sed_density, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetReal(PETSC_NULL, "-FlowDepth", &FlowDepth, PETSC_NULL);            // Ali Khosro
    PetscOptionsGetInt(PETSC_NULL, "-cross_diffusion", &cross_diffusion, PETSC_NULL);                     // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-lowRe", &lowRe, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-stension", &surface_tension, PETSC_NULL);            // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-delete", &delete_previous_file, PETSC_NULL);    // Seokkoo Kang: delete previous time step's saved file for saving disk space
    PetscOptionsGetInt(PETSC_NULL, "-mixed", &mixed, PETSC_NULL);            // Seokkoo Kang: mixed model option for LES
    PetscOptionsGetInt(PETSC_NULL, "-clark", &clark, PETSC_NULL);            // Seokkoo Kang: mixed model option for LES
    PetscOptionsGetInt(PETSC_NULL, "-vorticity", &vorticity, PETSC_NULL);            // Seokkoo Kang: vorticity form for viscous terms
    PetscOptionsGetInt(PETSC_NULL, "-pseudo", &pseudo_periodic, PETSC_NULL);    // Seokkoo Kang: pseudo periodic BC in k-direction for genenration of inflow condition
    PetscOptionsGetInt(PETSC_NULL, "-levelset", &levelset, PETSC_NULL);     // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-fix_level", &fix_level, PETSC_NULL);     // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-levelset_it", &levelset_it, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-levelset_tau", &levelset_tau, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-dam_break", &dam_break, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-k_gate", &k_gate, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-rotdir", &rotdir, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, PETSC_NULL);    
    PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, PETSC_NULL);    
    PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, PETSC_NULL);    
    PetscOptionsGetInt(PETSC_NULL, "-laplacian", &laplacian, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-qcrout", &qcr, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-ii_periodic", &ii_periodic, PETSC_NULL);    
    PetscOptionsGetInt(PETSC_NULL, "-jj_periodic", &jj_periodic, PETSC_NULL);    
    PetscOptionsGetInt(PETSC_NULL, "-kk_periodic", &kk_periodic, PETSC_NULL);    
    
    periodic = i_periodic+j_periodic+k_periodic+ii_periodic+jj_periodic+kk_periodic;

    PetscOptionsGetReal(PETSC_NULL, "-scale_velocity", &scale_velocity, PETSC_NULL);
    
    PetscOptionsGetInt(PETSC_NULL, "-perturb", &initial_perturbation, PETSC_NULL);    // Seokkoo Kang: give a random perturbation for initial condition
    PetscOptionsGetInt(PETSC_NULL, "-perturb1", &initial_perturbation1, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-skew", &skew, PETSC_NULL);                // Seokkoo Kang: skew symmetric form of advection term
    PetscOptionsGetInt(PETSC_NULL, "-dynamic_freq", &dynamic_freq, PETSC_NULL);        // Seokkoo Kang: LES dynamic compute frequency 
    if(dynamic_freq<1) dynamic_freq=1;
    
    PetscOptionsGetInt(PETSC_NULL, "-save_inflow", &save_inflow, PETSC_NULL);        // Seokkoo Kang: save infow BC to files; should be used in conjunction wiht -pseudo 1
    PetscOptionsGetInt(PETSC_NULL, "-save_inflow_period", &save_inflow_period, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-save_inflow_minus", &save_inflow_minus, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-ti_lastsave", &ti_lastsave, PETSC_NULL);

    read_inflow_period = save_inflow_period;
	PetscOptionsGetInt(PETSC_NULL, "-read_inflow_period", &read_inflow_period, PETSC_NULL);
    
    PetscOptionsGetInt(PETSC_NULL, "-localstep", &localstep, PETSC_NULL);        // Seokkoo Kang: localstep ( explict + implicit momentum solver )
    PetscOptionsGetInt(PETSC_NULL, "-recycle", &inflow_recycle_perioid, PETSC_NULL);    // Seokkoo Kang, set recycling period of the inflow data
    PetscOptionsGetInt(PETSC_NULL, "-save_memory", &save_memory, PETSC_NULL);    // Seokkoo Kang, save_memory
    PetscOptionsGetInt(PETSC_NULL, "-ibm_search", &ibm_search, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-ip", &i_proc, PETSC_NULL);            // Seokkoo Kang: number of processors in i direction
    PetscOptionsGetInt(PETSC_NULL, "-jp", &j_proc, PETSC_NULL);            // Seokkoo Kang: number of processors in j direction
    PetscOptionsGetInt(PETSC_NULL, "-kp", &k_proc, PETSC_NULL);            // Seokkoo Kang: number of processors in k direction
    PetscOptionsGetInt(PETSC_NULL, "-poisson", &poisson, PETSC_NULL);     // Seokkoo Kang
    PetscOptionsGetInt(PETSC_NULL, "-amg_agg", &amg_agg, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-testfilter_ik", &testfilter_ik, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-testfilter_1d", &testfilter_1d, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-poisson_it", &poisson_it, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-i_homo_filter", &i_homo_filter, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-j_homo_filter", &j_homo_filter, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-k_homo_filter", &k_homo_filter, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-fix_outlet", &fix_outlet, PETSC_NULL);

    PetscOptionsGetInt(PETSC_NULL, "-fix_inlet", &fix_inlet, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-cf_nacelle_fromfile", &cf_nacelle_fromfile, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-tipcorrectionratio_Fa", &tipcorrectionratio_Fa, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-forcewidthfixed", &forcewidthfixed, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetInt(PETSC_NULL, "-forcedistrwidth_surfacescale", &forcedistrwidth_surfacescale, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetReal(PETSC_NULL, "-dhi_fixed", &dhi_fixed, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetReal(PETSC_NULL, "-dhj_fixed", &dhj_fixed, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetReal(PETSC_NULL, "-dhk_fixed", &dhk_fixed, PETSC_NULL); // xyang 6-3-2011

    //Hossein
          //Variables for Solitary waves
	PetscOptionsGetInt(PETSC_NULL, "-solitary_wave", &solitary_wave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ti_start_solitary_wave", &ti_start_solitary_wave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ti_restart_solitary_wave", &ti_restart_solitary_wave, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-solitary_wave_amplitude", &solitary_wave_amplitude, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-inlet_bed_elevation", &inlet_bed_elevation, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-inlet_z_for_solitary_wave", &inlet_z_for_solitary_wave, PETSC_NULL);
	
        //Variables for Linear wave single
	PetscOptionsGetInt(PETSC_NULL, "-ti_start_linear_wave_single", &ti_start_linear_wave_single, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-ti_restart_linear_wave_single", &ti_restart_linear_wave_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-linear_wave_single_amplitude", &linear_wave_single_amplitude, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-linear_wave_single_number", &linear_wave_single_number, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-inlet_z_for_linear_wave_single", &inlet_z_for_linear_wave_single, PETSC_NULL);

        //Variables for wave_momentum_source
	PetscOptionsGetInt(PETSC_NULL, "-wave_momentum_source", &wave_momentum_source, PETSC_NULL); //For momentum source 1 in x direction and 2 both directions
	PetscOptionsGetInt(PETSC_NULL, "-wave_sponge_layer", &wave_sponge_layer, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_source_cent_z", &wave_source_cent_z, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_zs", &wave_sponge_zs, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_z01", &wave_sponge_z01, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_z02", &wave_sponge_z02, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_xs", &wave_sponge_xs, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_x01", &wave_sponge_x01, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_sponge_x02", &wave_sponge_x02, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_angle_single", &wave_angle_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_K_single", &wave_K_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_depth", &wave_depth, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_a_single", &wave_a_single, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_wind_reflength", &wave_wind_reflength, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_wind_refvel", &wave_wind_refvel, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-wave_wind_yshift", &wave_wind_yshift, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_IB", &wave_IB, PETSC_NULL);

    PetscOptionsGetInt(PETSC_NULL, "-inflow_levelset", &inflow_levelset, PETSC_NULL);  // xyang

	//for air_flow_levelset
	PetscOptionsGetInt(PETSC_NULL, "-air_flow_levelset", &air_flow_levelset, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-air_flow_levelset_periodic", &air_flow_levelset_periodic, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_average_k", &wave_average_k, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_skip", &wave_skip, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_ti_start", &wave_ti_start, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wind_skip", &wind_skip, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wind_start_read", &wind_start_read, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_start_read", &wave_start_read, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wind_recicle", &wind_recicle, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_recicle", &wave_recicle, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_k_ave", &wave_k_ave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-wave_i_ave", &wave_i_ave, PETSC_NULL);	
	PetscOptionsGetInt(PETSC_NULL, "-wave_ti_startave", &wave_ti_startave, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-freesurface_wallmodel", &freesurface_wallmodel, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-viscosity_wallmodel", &viscosity_wallmodel, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL, "-channel_height", &channel_height, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-floating_turbine_case", &floating_turbine_case, PETSC_NULL);
	//End variables for wave_momentum_source
	// End (Toni)	

    if(forcewidthfixed) {
        PetscPrintf(PETSC_COMM_WORLD, "\n Actuator model: the width for force distribution is fixed at dhi=%le dhj=%le dhk=%le \n\n", dhi_fixed, dhj_fixed, dhj_fixed);
    }

    PetscOptionsGetInt(PETSC_NULL, "-NonUniform_ADModel", &NonUniform_ADModel, PETSC_NULL); // xyang 12-7-2010
    PetscOptionsGetInt(PETSC_NULL, "-maxiteraction_rotormodel", &maxiteraction_rotormodel, PETSC_NULL); // xyang 12-7-2010

    PetscOptionsGetInt(PETSC_NULL, "-rotor_modeled", &rotor_model, PETSC_NULL); // xyang 12-7-2010
    PetscOptionsGetInt(PETSC_NULL, "-TheSameObject", &TheSameObject, PETSC_NULL); // xyang 12-7-2010
    PetscOptionsGetInt(PETSC_NULL, "-UlagrFromSurface", &UlagrFromSurface, PETSC_NULL); //

    PetscOptionsGetInt(PETSC_NULL, "-nacelle_model", &nacelle_model, PETSC_NULL); // xyang 12-7-2010
    PetscOptionsGetReal(PETSC_NULL, "-indf_ax", &indf_ax, PETSC_NULL); // xyang 12-16-2010

    PetscOptionsGetReal(PETSC_NULL, "-rho_fluid", &rho_fluid, PETSC_NULL); // xyang 3-18-2011

    PetscOptionsGetInt(PETSC_NULL, "-NumNacellePerLoc", &NumNacellePerLoc, PETSC_NULL); // xyang 3-18-2011
    PetscOptionsGetReal(PETSC_NULL, "-reflength_wt", &reflength_wt, PETSC_NULL); // xyang 12-16-2010
    PetscOptionsGetReal(PETSC_NULL, "-reflength_nacelle", &reflength_nacelle, PETSC_NULL); // xyang 12-16-2010
    PetscOptionsGetReal(PETSC_NULL, "-refvel_wt", &refvel_wt, PETSC_NULL); // xyang 12-16-2010
    PetscOptionsGetReal(PETSC_NULL, "-refvel_cfd", &refvel_cfd, PETSC_NULL); // xyang 12-16-2010
    PetscOptionsGetReal(PETSC_NULL, "-loc_refvel", &loc_refvel, PETSC_NULL); // xyang 2-12-2012
    PetscOptionsGetInt(PETSC_NULL, "-deltafunc", &deltafunc, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetReal(PETSC_NULL, "-halfwidth_dfunc", &halfwidth_dfunc, PETSC_NULL); 

    PetscOptionsGetReal(PETSC_NULL, "-halfwidth_dfunc", &halfwidth_dfunc, PETSC_NULL); 
    PetscOptionsGetInt(PETSC_NULL, "-i_periodicIB", &i_periodicIB, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetInt(PETSC_NULL, "-j_periodicIB", &j_periodicIB, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetInt(PETSC_NULL, "-k_periodicIB", &k_periodicIB, PETSC_NULL); // xyang 6-3-2011
    PetscOptionsGetInt(PETSC_NULL, "-les_prt", &les_prt, PETSC_NULL); // xyang 11-03-2011
    PetscOptionsGetInt(PETSC_NULL, "-MoveFrame", &MoveFrame, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-u_frame", &u_frame, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-v_frame", &v_frame, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-w_frame", &w_frame, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-ii_periodicWT", &ii_periodicWT, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-jj_periodicWT", &jj_periodicWT, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-kk_periodicWT", &kk_periodicWT, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-Nx_WT", &Nx_WT, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-Ny_WT", &Ny_WT, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-Nz_WT", &Nz_WT, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-Sx_WT", &Sx_WT, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-Sy_WT", &Sy_WT, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-Sz_WT", &Sz_WT, PETSC_NULL);

    PetscOptionsGetInt(PETSC_NULL, "-turbine", &NumberOfTurbines, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-NumberOfNacelle", &NumberOfNacelle, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-FixTipSpeedRatio", &FixTipSpeedRatio, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-rstart_turbinerotation", &rstart_turbinerotation, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-turbinetorquecontrol", &turbinetorquecontrol, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-turbineindividualpitchcontrol", &turbineindividualpitchcontrol, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-fixturbineangvel", &fixturbineangvel, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-Shen_AL", &Shen_AL, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-correction3D_CH", &correction3D_CH, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-correction3D_CL", &correction3D_CL, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-smoothforce_AL", &smoothforce_AL, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-Shen1_AL", &Shen1_AL, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-correction_ALShen1", &correction_ALShen1, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-c0_CL", &c0_CL, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-c1_CH", &c1_CH, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-c2_CH", &c2_CH, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-c3_CH", &c3_CH, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-a_shen", &a_shen, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-b_shen", &b_shen, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-c_shen", &c_shen, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-Prandtl_AL", &Prandtl_AL, PETSC_NULL);  // xyang
    PetscOptionsGetInt(PETSC_NULL, "-correction3D_DS", &correction3D_DS, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-a_DS", &a_DS, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-b_DS", &b_DS, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-d_DS", &d_DS, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-CD_0", &CD_0, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-AOA_0", &AOA_0, PETSC_NULL);  // xyang
    PetscOptionsGetReal(PETSC_NULL, "-refangle_AL", &refangle_AL, PETSC_NULL);  // xyang
      PetscOptionsGetInt(PETSC_NULL, "-specifycirculation", &specifycirculation, PETSC_NULL);
      //PetscOptionsGetInt(PETSC_NULL, "-turbinestructuremodel", &turbinestructuremodel, PETSC_NULL);
     // PetscOptionsGetInt(PETSC_NULL, "-torsion_turbinestructure", &torsion_turbinestructure, PETSC_NULL);
    
// Canopy
    cnpy =0;
    NumberOfCnpy = 0;
    cnpyHeightDim = 20.;
    //Meric
    cnpyCd = 1.0;
    PetscOptionsGetInt(PETSC_NULL, "-canopy", &cnpy, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-NumberOfCanopies", &NumberOfCnpy, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-CanopyHeightDim", &cnpyHeightDim, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-CanopyHeightMlt", &cnpyHeightMlt, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-CanopyDispZ", &cnpyDispZ, PETSC_NULL);
    // Meric
    PetscOptionsGetReal(PETSC_NULL, "-CanopyCd", &cnpyCd, PETSC_NULL);

    if(movefsi || rotatefsi ) save_memory=0;
    if(!rotatefsi && !movefsi) rstart_fsi=0;

    sprintf(path, ".");
    PetscOptionsGetString(PETSC_NULL,"-path", path, 256, PETSC_NULL);        //  Seokkoo Kang: path for saving output; grid.dat, bcs,dat should be put there, but control.dat should exist in the current directory where job is submitted

    sprintf(gridfile, "grid.dat");
    PetscOptionsGetString(PETSC_NULL,"-grid", gridfile, 256, PETSC_NULL);    //  Seokkoo Kang: the name of the grid file other than grid.dat if you want

    sprintf(path_inflow, "./inflow");
        PetscOptionsGetString(PETSC_NULL,"-path_inflow", path_inflow, 256, PETSC_NULL);         //  Xiaolei Yang: path for inflow field
    
        sprintf(path_plane_save, "./plane");
        PetscOptionsGetString(PETSC_NULL,"-path_plane_save", path_plane_save, 256, PETSC_NULL);         //  Xiaolei Yang: path for inflow field

    int len=strlen(path);
    if(path[len-1]=='/') path[len-1]=0;

    extern void read_grid();
    PetscPrintf(PETSC_COMM_WORLD, "\n**************** Finished reading grid.dat file***********************\n");

    char saveoption_file[400];
    sprintf(saveoption_file, "%s/savepoints", path);
    FILE *fp=fopen(saveoption_file, "r");

    // Savepoints based on indices	

    // if(fp!=NULL) {

       // i=0;

       // do {

       //     fscanf(fp, "%d %d %d\n", &save_point[i][0], &save_point[i][1], &save_point[i][2]);

       //     i++;

       // } while(!feof(fp));

       // nsave_points=i;

       // fclose(fp);

    //}


        // Savepoints based on coordinates
      if(fp!=NULL) {

			i=0;

			do {

				fscanf(fp, "%le %le %le\n", &save_coor[i][0], &save_coor[i][1], &save_coor[i][2]); // xiaolei

				i++;

			} while(!feof(fp));

			nsave_points=i;

			fclose(fp);

		}

    char saveoption_file_level[400];
    sprintf(saveoption_file_level, "%s/savepoints_level", path);
    FILE *fpl=fopen(saveoption_file_level, "r");

    if(fpl!=NULL) {
        i=0;
        do {
            fscanf(fpl, "%d %d %d\n", &save_point_level[i][0], &save_point_level[i][1], &save_point_level[i][2]);
            i++;
        } while(!feof(fpl));
        nsave_points_level=i;
        fclose(fpl);
    }
    
    //Add-Hossein for savekplanes
    	/*if(save_inflow) {
		char save_ksection_file[400];
		sprintf(save_ksection_file, "%s/savekplanes", path);
		FILE *fp=fopen(save_ksection_file, "r");

		if(fp!=NULL) {
			int i=0;
			do {
				fscanf(fp, "%d\n", &save_ksection[i]);
				i++;
			} while(!feof(fp));
			nsave_ksection=i;
			fclose(fp);
		}

		char save_jsection_file[400];
		sprintf(save_jsection_file, "%s/savejplanes", path);
		fp=fopen(save_jsection_file, "r");

		if(fp!=NULL) {
			int i=0;
			do {
				fscanf(fp, "%d\n", &save_jsection[i]);
				i++;
			} while(!feof(fp));
			nsave_jsection=i;
			fclose(fp);
                        PetscPrintf(PETSC_COMM_WORLD, "Save J Planes: %d\n", nsave_jsection);
                }
		char save_isection_file[400];
		sprintf(save_isection_file, "%s/saveiplanes", path);
		fp=fopen(save_isection_file, "r");

		if(fp!=NULL) {
			int i=0;
			do {
				fscanf(fp, "%d\n", &save_isection[i]);
				i++;
			} while(!feof(fp));
			nsave_isection=i;
			fclose(fp);
                        PetscPrintf(PETSC_COMM_WORLD, "Save I Planes: %d\n", nsave_isection);
                }

	}*/

    if(TwoD) PetscPrintf(PETSC_COMM_WORLD, "\n\n!!! 2D computation !!! \n\n");
    if(i_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nI-Periodic\n");
    if(ii_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nII-Periodic\n");
    if(j_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nJ-Periodic \n");
    if(jj_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nJJ-Periodic \n");
    if(k_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nK-Periodic \n");
    if(kk_periodic) PetscPrintf(PETSC_COMM_WORLD, "\nKK-Periodic \n");

    PetscOptionsGetReal(PETSC_NULL, "-max_cs", &max_cs, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-flux", &inlet_flux, PETSC_NULL);            // Seokkoo Kang: the amount of inlet flux, if not set mean bulk velocity is set to 1
    PetscPrintf(PETSC_COMM_WORLD, "\nControl.dat value for flux = %f", *(&inlet_flux), "\n");
    PetscOptionsGetReal(PETSC_NULL, "-imp_tol", &imp_free_tol, PETSC_NULL);        // Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
    PetscOptionsGetReal(PETSC_NULL, "-poisson_tol", &poisson_tol, PETSC_NULL);        // Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
    PetscOptionsGetReal(PETSC_NULL, "-les_eps", &les_eps, PETSC_NULL);        // Seokkoo Kang: small value for preventing very large Cs values in les>1
    PetscOptionsGetReal(PETSC_NULL, "-dpdz", &mean_pressure_gradient, &dpdz_set);        // Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
    PetscOptionsGetReal(PETSC_NULL, "-roughness", &roughness_size, &rough_set);    // Seokkoo Kang: roughness_size
    PetscOptionsGetReal(PETSC_NULL, "-amg_thresh", &amg_thresh, PETSC_NULL);

    PetscOptionsGetReal(PETSC_NULL, "-rho0", &rho_water, PETSC_NULL);    
    PetscOptionsGetReal(PETSC_NULL, "-rho1", &rho_air, PETSC_NULL);        
    PetscOptionsGetReal(PETSC_NULL, "-mu0", &mu_water, PETSC_NULL);    
    PetscOptionsGetReal(PETSC_NULL, "-mu1", &mu_air, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-dthick", &dthick, &dthick_set);
    PetscOptionsGetReal(PETSC_NULL, "-seudo_dt", &seudo_dt, PETSC_NULL);

    if(English){
        rho_water=62.37;
        rho_air=0.0765;
        mu_water=2.1*0.00001;
        mu_air=3.79*0.0000001;
    }

    PetscPrintf(PETSC_COMM_WORLD, "\nrho0=%f, rho1=%f, mu0=%f, mu1=%f\n", rho_water, rho_air, mu_water, mu_air);

    PetscOptionsGetReal(PETSC_NULL, "-gx", &gravity_x, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-gy", &gravity_y, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-gz", &gravity_z, PETSC_NULL);

    PetscOptionsGetReal(PETSC_NULL, "-x_r", &(x_r), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-y_r", &(y_r), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-z_r", &(z_r), PETSC_NULL);

    PetscOptionsGetReal(PETSC_NULL, "-inlet_y", &inlet_y, &inlet_y_flag); //Hossein 
    PetscOptionsGetReal(PETSC_NULL, "-outlet_y", &outlet_y, PETSC_NULL);  //Hossein

    PetscOptionsGetReal(PETSC_NULL, "-inlet_z", &inlet_z, &inlet_z_flag);
    PetscOptionsGetReal(PETSC_NULL, "-outlet_z", &outlet_z, PETSC_NULL);

    PetscOptionsGetReal(PETSC_NULL, "-x_c", &(CMx_c), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-y_c", &(CMy_c), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-z_c", &(CMz_c), PETSC_NULL);

    PetscOptionsGetReal(PETSC_NULL, "-imp_atol", &(imp_atol), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-imp_rtol", &(imp_rtol), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-imp_stol", &(imp_stol), PETSC_NULL);

    PetscOptionsGetReal(PETSC_NULL, "-angvel", &angvel, PETSC_NULL);

    if (fish) {
        PetscOptionsGetReal(PETSC_NULL, "-St_exp", &(St_exp), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, "-wlngth", &(wavelength), PETSC_NULL);
    }

    PetscPrintf(PETSC_COMM_WORLD, "tiout %d %le %le thin %d!\n",tiout, imp_atol,imp_rtol,thin);

    if (MHV) 
        L_dim=1./25.4;//.005;
    else
    L_dim=1.0;//0.0254;

    if (MHV) NumberOfBodies=3;
    if (immersed) {
        PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
        PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
        ibm_ptr = ibm;
        fsi_ptr = fsi;
    }

    if (rotor_model == 1)  {
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_wt[i].angvel_x=0;
            fsi_wt[i].angvel_y=0;
            fsi_wt[i].angvel_z=0;
            fsi_wt[i].angvel_axis=0.0;
        }
    }

    if (rotor_model == 2)  {
            PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
            PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
            for (i=0;i<NumberOfTurbines;i++) {
                fsi_wt[i].angvel_x=0;
                fsi_wt[i].angvel_y=0;
                fsi_wt[i].angvel_z=0;
                fsi_wt[i].angvel_axis=0.0;
            }
            PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
     }

    if (rotor_model == 3)  {
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_wt[i].angvel_x=0;
            fsi_wt[i].angvel_y=0;
            fsi_wt[i].angvel_z=0;
            fsi_wt[i].angvel_axis=0.0;
        }
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
    }

    if (rotor_model == 4)  {
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_wt[i].angvel_x=0;
            fsi_wt[i].angvel_y=0;
            fsi_wt[i].angvel_z=0;
            fsi_wt[i].angvel_axis=0.0;
        }
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
    }

    if (rotor_model == 5)  {
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_wt[i].angvel_x=0;
            fsi_wt[i].angvel_y=0;
            fsi_wt[i].angvel_z=0;
            fsi_wt[i].angvel_axis=0.0;
        }

        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_acl2ref);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_acl2ref);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_acl2ref[i].angvel_x=0;
            fsi_acl2ref[i].angvel_y=0;
            fsi_acl2ref[i].angvel_z=0;
            fsi_acl2ref[i].angvel_axis=0.0;
        }

        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
    }

    if (rotor_model == 6)  {
        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &wtm);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_wt);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_wt[i].angvel_x=0;
            fsi_wt[i].angvel_y=0;
            fsi_wt[i].angvel_z=0;
            fsi_wt[i].angvel_axis=0.0;
        }

        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_acl2ref);
        PetscMalloc(NumberOfTurbines*sizeof(FSInfo), &fsi_acl2ref);
        for (i=0;i<NumberOfTurbines;i++) {
            fsi_acl2ref[i].angvel_x=0;
            fsi_acl2ref[i].angvel_y=0;
            fsi_acl2ref[i].angvel_z=0;
            fsi_acl2ref[i].angvel_axis=0.0;
        }

        PetscMalloc(NumberOfTurbines*sizeof(IBMNodes), &ibm_ACD);
    }

    if (nacelle_model) {

        PetscMalloc(NumberOfNacelle*sizeof(IBMNodes), &ibm_nacelle);
            PetscMalloc(NumberOfNacelle*sizeof(FSInfo), &fsi_nacelle);
            for (i=0;i<NumberOfNacelle;i++) {
                fsi_nacelle[i].angvel_x=0;
                fsi_nacelle[i].angvel_y=0;
                fsi_nacelle[i].angvel_z=0;
                fsi_nacelle[i].angvel_axis=0.0;
                ibm_nacelle[i].coeff_SG=1.0;
            }
    }

    MG_Initial(&usermg, ibm);

    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;

    VecDuplicate(user->lP, &user->lUstar);


    //Hossein
    // add (Toni)
	//Initialitzation of wave_IB
	if(wave_IB==1){
		extern void Initialize_wave_IB(UserCtx *user);
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		for (bi=0; bi<block_number; bi++)Initialize_wave_IB(&user[bi]);
		//Initialize_wave_IB(&user[0]);	//Ali turned on May 2016
	}

	//Initialitzation of wave_momentum_source
	if(wave_momentum_source){
		extern void Initialize_wave(UserCtx *user);
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		for (bi=0; bi<block_number; bi++)Initialize_wave(&user[bi]);
		//Initialize_wave(&user[0]); // Ali turned on May 2016	
	}
	//End Initialitzation of wave_momentum_source	
	if(air_flow_levelset==2){
		extern void Initialize_wind(UserCtx *user);
		extern void WIND_DATA_input(UserCtx *user);				
		level = usermg.mglevels-1;
		user = usermg.mgctx[level].user;
		for (bi=0; bi<block_number; bi++)Initialize_wind(&user[bi]); 
		for (bi=0; bi<block_number; bi++)WIND_DATA_input(&user[bi]); 	
		//Initialize_wind(&user[0]);	//Ali turned on May 2016
		//WIND_DATA_input(&user[0]);	// Ali turned on May 2016
	}		
// End (Toni)


    if(averaging) {    // Seokkoo Kang
        level = usermg.mglevels-1;
        user = usermg.mgctx[level].user;
        VecDuplicate(user->Ucat, &user->Ucat_sum);
        VecDuplicate(user->Ucat, &user->Ucat_cross_sum);
        VecDuplicate(user->Ucat, &user->Ucat_square_sum);
        VecSet(user->Ucat_sum,0);
        VecSet(user->Ucat_cross_sum,0);
        VecSet(user->Ucat_square_sum,0);

        VecDuplicate(user->P, &user->P_sum);
        VecSet(user->P_sum,0);

        if(rans) {
            VecDuplicate(user->P, &user->K_sum);
            VecSet(user->K_sum, 0.);
        }

        if(sandwave) {
            VecDuplicate(user->P, &user->lUstar_sum);
            VecSet(user->lUstar_sum, 0.);
            VecDuplicate(user->P, &user->lUstar_);
            VecSet(user->lUstar_, 0.);
        }

        if(conv_diff) {
            VecDuplicate(user->P, &user->Conc_sum);
            VecSet(user->Conc_sum, 0.);
            VecDuplicate(user->Ucat, &user->Conc_cross_sum);
            VecSet(user->Conc_cross_sum, 0.);
        }

        if(levelset) {
            VecDuplicate(user->P, &user->Level_sum);
            VecDuplicate(user->P, &user->Level_square_sum);
            VecSet(user->Level_sum, 0.);
            VecSet(user->Level_square_sum, 0.);
        }

        if(les) {
            VecDuplicate(user->P, &user->Nut_sum);
            VecSet(user->Nut_sum, 0.);
        }

        if(averaging>=2) {
            VecDuplicate(user->P, &user->P_square_sum);
            VecSet(user->P_square_sum,0);
        }

        if(averaging>=3) {
            if(les) {
                VecDuplicate(user->P, &user->tauS_sum);
                VecSet(user->tauS_sum, 0);
            }

            VecDuplicate(user->P, &user->Udp_sum);
            VecSet(user->Udp_sum, 0);

            VecDuplicate(user->Ucont, &user->dU2_sum);
            VecSet(user->dU2_sum, 0);

            VecDuplicate(user->Ucont, &user->UUU_sum);
            VecSet(user->UUU_sum, 0);

            VecDuplicate(user->Ucont, &user->Vort_sum);
            VecSet (user->Vort_sum, 0);

            VecDuplicate(user->Ucont, &user->Vort_square_sum);
            VecSet (user->Vort_square_sum, 0);
        }
    }


    #ifdef DIRICHLET
    if(freesurface) {
        extern void Initialize_free_surface(UserCtx *user);
        level = usermg.mglevels-1;
        user = usermg.mgctx[level].user;
        for (bi=0; bi<block_number; bi++)  Initialize_free_surface(&user[bi]);
    }
    else {
        PetscPrintf(PETSC_COMM_WORLD, "\n**************** Warning !! Freesurface option not set *********************************!\n");
    }
  #endif
  
    if (immersed) {
        level = usermg.mglevels-1;
        user = usermg.mgctx[level].user;
        for (bi=0; bi<block_number; bi++) {
            PetscMalloc(NumberOfBodies*sizeof(IBMList), &(user[bi].ibmlist));
            for (ibi=0;ibi<NumberOfBodies;ibi++) {
                InitIBMList(&(user[bi].ibmlist[ibi]));
            }
        }

        if (MHV) {
            i=0;
            // read casing
            CMz_c=0.;//1.105;
            // read valves
            CMz_c=4.49+0.31;
            CMy_c=.0;
            L_dim=1./28.;

            for (ibi=1; ibi<NumberOfBodies; ibi++) {
                if (ibi==2) CMy_c=-CMy_c;
                PetscPrintf(PETSC_COMM_WORLD, "Ibm read MHV!\n");
                FsiInitialize(0, &fsi[ibi], ibi);
            }
          
            fsi[1].y_c = -0.0878; fsi[1].z_c = 4.143;//4.21;
            fsi[2].y_c =  0.0878; fsi[2].z_c = 4.143;//4.21;

            fsi[1].S_ang_n[0]= max_angle; fsi[1].S_ang_r[0]= max_angle; fsi[1].S_ang_r[0]= max_angle;  fsi[1].S_ang_rm1[0]= max_angle;
            fsi[2].S_ang_n[0]= -max_angle; fsi[2].S_ang_r[0]= -max_angle; fsi[2].S_ang_r[0]= -max_angle;  fsi[2].S_ang_rm1[0]= -max_angle;

            for (ibi=1; ibi<NumberOfBodies; ibi++) {    
                Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi],0.,ibi);
            }
            PetscBarrier(PETSC_NULL);
            for (ibi=0; ibi<NumberOfBodies; ibi++) {
                ibm_surface_out(&ibm[ibi], 0, ibi);
            }
        } 
        else{
        //KFlora - call the rot_fsi.dat file to read in parameters of rotation
            PetscMalloc(NumberOfBodies*sizeof(FSISetup), &(fsimov));
            if(rotatefsi){
                Read_Rotate_FSI_Parameter_Input(&fsimov[0]);
                PetscBarrier(PETSC_NULL);    
            }

            PetscPrintf(PETSC_COMM_WORLD, "\n**************** Reading IB data for %i IB bodies\n", NumberOfBodies);
            PetscReal setupAngle = 0.0;
            for (i=0;i<NumberOfBodies;i++) {
                ibm[i].fixedNodes = 1;//set default for ibms to fixed

                if (rotatefsi){
                    CMx_c = fsimov[i].x_c;
                    CMy_c = fsimov[i].y_c;
                    CMz_c = fsimov[i].z_c;
                    x_r = fsimov[i].x_r;
                    y_r = fsimov[i].y_r;
                    z_r = fsimov[i].z_r;\
                    setupAngle = fsimov[i].XYangle;
                    if(fsimov[i].fixed_ang_vel > 0.000001 || fsimov[i].fixed_ang_vel < -0.000001){
                        ibm[i].fixedNodes = 0;
                    }
                    PetscPrintf(PETSC_COMM_WORLD, "FixedNodes for IBM%i = %i \n", i,ibm[i].fixedNodes);
                }
                PetscPrintf(PETSC_COMM_WORLD, "Now Read UCD%i called from main.c\n", i);
                ibm_read_ucd(&ibm[i], i, CMx_c, CMy_c, CMz_c, x_r, y_r, z_r, setupAngle);
                PetscBarrier(PETSC_NULL);
                FsiInitialize(0, &fsi[i], i);
                
                PetscPrintf(PETSC_COMM_WORLD, "fsimov= %i, Rot Dir=%i, XY_Angle=%lf, Ang Vel=%lf,  x_r=%lf, y_r=%lf, z_r=%lf, x_c=%lf, y_c=%lf, z_c=%lf\n",\
                    i, fsimov[i].rot_dir, fsimov[i].XYangle, fsimov[i].fixed_ang_vel, fsimov[i].x_r, fsimov[i].y_r, fsimov[i].z_r, fsimov[i].x_c, fsimov[i].y_c, fsimov[i].z_c);

                if (rotatefsi){
                    fsi[i].rot_dir = fsimov[i].rot_dir;
                    if(fsimov[i].rot_dir==0) fsi[i].S_ang_n[1] = fsimov[i].fixed_ang_vel; //x axis
                    else if(fsimov[i].rot_dir==1) fsi[i].S_ang_n[3] = fsimov[i].fixed_ang_vel;//y axis
                    else if(fsimov[i].rot_dir==2) fsi[i].S_ang_n[5] = fsimov[i].fixed_ang_vel;//z axis
                    else if (fsimov[i].rot_dir==3)fsi[i].S_ang_n[1] = fsimov[i].fixed_ang_vel;//skewed in XY plane
                    else PetscPrintf(PETSC_COMM_WORLD, "Rotation direction missing in fsi_input.dat file\n");
                
                    fsi[i].fixed_ang_vel = fsimov[i].fixed_ang_vel;    
                    fsi[i].x_c = fsimov[i].x_c + fsimov[i].x_r;
                    fsi[i].y_c = fsimov[i].y_r + fsimov[i].y_c;
                    fsi[i].z_c = fsimov[i].z_r + fsimov[i].z_c;
                    PetscPrintf(PETSC_COMM_WORLD, "fsi= %i, rot_dir = %i, angl Velocity[1]=%lf\n", i, fsi[i].rot_dir, fsi[i].S_ang_n[1]);
                    fsi[i].XYangle = fsimov[i].XYangle;
                }
            }
        }
    }
  
    if (immersed) {
        ti = 0;
        if (rstart_flg) ti = tistart;
    }

    if (rotor_model) {
        PetscReal cl = 1.;
        PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);
 
        if (!my_rank) {
            FILE *fd;
            char str[256];
            sprintf(str, "%s/Turbine.inp", path);
            fd = fopen(str, "r");
            if(!fd) PetscPrintf(PETSC_COMM_WORLD, "cannot open %s !\n", str),exit(0);
            for (ibi=0;ibi<NumberOfTurbines;ibi++) {
                char string[22];
                fgets(string, 22, fd);
                fscanf(fd, "%i %i %*s %*s", &(wtm[ibi].num_blade), &(wtm[ibi].num_foiltype)); 

                PetscPrintf(PETSC_COMM_WORLD, "Turbine.inp file - Number of Blades: %i Foil Type:%i \n", wtm[ibi].num_blade,wtm[ibi].num_foiltype );

                fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].x_c), &(fsi_wt[ibi].y_c), &(fsi_wt[ibi].z_c)); 
                fscanf(fd, "%le %le %*s %*s",&(fsi_wt[ibi].J_rotation) , &(fsi_wt[ibi].r_rotor) );
                fscanf(fd, "%le %le %*s %*s",&(wtm[ibi].r_nearhubinflowcorrection), &(wtm[ibi].r_nacelle));

                // Controller information
                fscanf(fd, "%i %*s",&(fsi_wt[ibi].TPCntrl)); 

                if(fsi_wt[ibi].TPCntrl==1){ // Clipper C96 Turbine controller
                    PetscPrintf(PETSC_COMM_WORLD, "Clipper Turbine Controller \n");
                    fscanf(fd, "%le %le %*s %*s",&(wtm[ibi].indf_axis), &(wtm[ibi].Tipspeedratio)); 
                    fscanf(fd, "%le %*s",&(wtm[ibi].CT)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].CP_max),&(fsi_wt[ibi].TSR_max), &(fsi_wt[ibi].angvel_fixed)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].Torque_generator_max), &(fsi_wt[ibi].GeneratorSpeed_desired), &(fsi_wt[ibi].GearBoxRatio)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].K_proportional), &(fsi_wt[ibi].K_integral), &(fsi_wt[ibi].K_derivative)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].Ki_IPC), &(fsi_wt[ibi].angvel_axis_err_relax), &(fsi_wt[ibi].Kratio_torque));
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].WindSpeed_rated), &(fsi_wt[ibi].WindSpeed_cutin), &(fsi_wt[ibi].WindSpeed_cutout));
                }
                else if(fsi_wt[ibi].TPCntrl==2){ //Jonkman controller
                    PetscPrintf(PETSC_COMM_WORLD, "Jonkman Controller \n");
                    fscanf(fd, "%le %le %*s %*s",&(fsi_wt[ibi].CP_max),&(fsi_wt[ibi].TSR_max) ); 
                    PetscPrintf(PETSC_COMM_WORLD, "#####CP TSR : %le %le \n", fsi_wt[ibi].CP_max,fsi_wt[ibi].TSR_max );
                    fscanf(fd, "%le %*s", &(fsi_wt[ibi].K_Gen) ); 
                    PetscPrintf(PETSC_COMM_WORLD, "##### K_gen : %le \n", fsi_wt[ibi].K_Gen );
                    fscanf(fd, "%le %le %*s %*s ", &(fsi_wt[ibi].GeneratorSpeed_desired),&(fsi_wt[ibi].GeneratorSyncSpeed)); 
                    PetscPrintf(PETSC_COMM_WORLD, "##### GenSPeedDes , SyncSpeed: %le \n", fsi_wt[ibi].GeneratorSpeed_desired,fsi_wt[ibi].GeneratorSyncSpeed );
                    fscanf(fd, "%le %le %*s %*s",&(fsi_wt[ibi].Torque_generator_max),  &(fsi_wt[ibi].GearBoxRatio)); 
                    PetscPrintf(PETSC_COMM_WORLD, "TorqueGeneratorMax and GearboxRatio: %le %le \n", fsi_wt[ibi].Torque_generator_max,fsi_wt[ibi].GearBoxRatio );

                    wtm[ibi].Tipspeedratio=fsi_wt[ibi].TSR_max;
                    fsi_wt[ibi].angvel_fixed=0.;
                    fsi_wt[ibi].K_proportional=0.; 
                    fsi_wt[ibi].K_integral=0.; 
                    fsi_wt[ibi].K_derivative=0.; 
                    fsi_wt[ibi].Kratio_torque=0.;
                    fsi_wt[ibi].WindSpeed_rated=0.; 
                    fsi_wt[ibi].WindSpeed_cutin=0.;
                    fsi_wt[ibi].WindSpeed_cutout=0.;
                    fsi_wt[ibi].angvel_axis_err_relax=0.;

                }
                else{
                    fscanf(fd, "%le %le %*s %*s",&(wtm[ibi].indf_axis), &(wtm[ibi].Tipspeedratio)); 
                    fscanf(fd, "%le %*s",&(wtm[ibi].CT)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].CP_max),&(fsi_wt[ibi].TSR_max), &(fsi_wt[ibi].angvel_fixed)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].Torque_generator_max), &(fsi_wt[ibi].GeneratorSpeed_desired), &(fsi_wt[ibi].GearBoxRatio)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].K_proportional), &(fsi_wt[ibi].K_integral), &(fsi_wt[ibi].K_derivative)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].Ki_IPC), &(fsi_wt[ibi].angvel_axis_err_relax), &(fsi_wt[ibi].Kratio_torque));
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].WindSpeed_rated), &(fsi_wt[ibi].WindSpeed_cutin), &(fsi_wt[ibi].WindSpeed_cutout));
                }

                fscanf(fd, "%i %*s",&(fsi_wt[ibi].YawCntrl)); 
                PetscPrintf(PETSC_COMM_WORLD, "Yaw Control: %i \n",fsi_wt[ibi].YawCntrl);
                if(fsi_wt[ibi].YawCntrl==0){
                    PetscPrintf(PETSC_COMM_WORLD, "Yaw Control set to fixed (No Control) \n");
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].nx_tb), &(fsi_wt[ibi].ny_tb), &(fsi_wt[ibi].nz_tb)); 
                    PetscPrintf(PETSC_COMM_WORLD, "Yaw control - Yaw1(nx_tb): %lf Yaw2(ny_tb):%lf Yaw3(nz_tb):%lf \n",fsi_wt[ibi].nx_tb, fsi_wt[ibi].ny_tb, fsi_wt[ibi].nz_tb); 
                } else if (fsi_wt[ibi].YawCntrl==1){
                    PetscPrintf(PETSC_COMM_WORLD, "Yaw Control set to Prescribed \n");
                    PetscPrintf(PETSC_COMM_WORLD, "Input file yaw_prescrb.inp will be read\n");
                    fscanf(fd, "%i %*s",&(fsi_wt[ibi].yawtistrt)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].nx_tb), &(fsi_wt[ibi].ny_tb), &(fsi_wt[ibi].nz_tb)); 
                } else {
                    PetscPrintf(PETSC_COMM_WORLD, "No Yaw control have been implemented yet\n");
                    PetscPrintf(PETSC_COMM_WORLD, "Fixing turbine yaw angle \n");
                    fsi_wt[ibi].YawCntrl=0;
                    fscanf(fd, "%i %*s",&(fsi_wt[ibi].yawtistrt)); 
                    fscanf(fd, "%le %le %le %*s %*s %*s",&(fsi_wt[ibi].nx_tb), &(fsi_wt[ibi].ny_tb), &(fsi_wt[ibi].nz_tb)); 
                }
 
                PetscPrintf(PETSC_COMM_WORLD, "BroadCast \n");

                MPI_Bcast(&(wtm[ibi].num_blade), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].num_foiltype), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                PetscPrintf(PETSC_COMM_WORLD, "Allocate \n");

                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(wtm[ibi].pitch));
                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(wtm[ibi].pitch_old));
                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(wtm[ibi].pitch_IPC));
                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(fsi_wt[ibi].Moment_bladebending));

                PetscPrintf(PETSC_COMM_WORLD, "OldPitch %i %i\n",ibi, wtm[ibi].num_blade);

                wtm[ibi].pitch_old[0] = 0.0;
                wtm[ibi].pitch_old[1] = 0.0;
                wtm[ibi].pitch_old[2] = 0.0;

                PetscPrintf(PETSC_COMM_WORLD, "NewPitch \n");
                fscanf(fd, "%le %le %le %le %*s %*s %*s %*s", &(wtm[ibi].pitch[0]), &(wtm[ibi].pitch[1]), &(wtm[ibi].pitch[2]), &(wtm[ibi].pitch_min));

                PetscPrintf(PETSC_COMM_WORLD, "Pitch set to: %le %le %le %le \n",wtm[ibi].pitch[0],wtm[ibi].pitch[1],wtm[ibi].pitch[2],wtm[ibi].pitch_min);
                PetscPrintf(PETSC_COMM_WORLD, "Done reading Turbine.inp for turbine %i\n", ibi);

                 double rr=sqrt(pow(fsi_wt[ibi].nx_tb,2)+pow(fsi_wt[ibi].ny_tb,2)+pow(fsi_wt[ibi].nz_tb,2));

                fsi_wt[ibi].nx_tb=fsi_wt[ibi].nx_tb/rr; 
                fsi_wt[ibi].ny_tb=fsi_wt[ibi].ny_tb/rr; 
                fsi_wt[ibi].nz_tb=fsi_wt[ibi].nz_tb/rr;

                fsi_wt[ibi].nx_tbo = 0.;
                fsi_wt[ibi].ny_tbo = 0.;
                fsi_wt[ibi].nz_tbo = 1.;
                
                MPI_Bcast(&(fsi_wt[ibi].YawCntrl), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].yawtistrt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(fsi_wt[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                PetscPrintf(PETSC_COMM_WORLD, "The rotating center %f %f %f \n", (fsi_wt[ibi].nx_tb), (fsi_wt[ibi].ny_tb), (fsi_wt[ibi].nz_tb));

                fsi_wt[ibi].x_c=fsi_wt[ibi].x_c/cl;
                fsi_wt[ibi].y_c=fsi_wt[ibi].y_c/cl;
                fsi_wt[ibi].z_c=fsi_wt[ibi].z_c/cl;

                MPI_Bcast(&(fsi_wt[ibi].TPCntrl), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_Gen), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GeneratorSyncSpeed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GeneratorTransSpeed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(fsi_wt[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].indf_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].Tipspeedratio), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].J_rotation), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].r_rotor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].CP_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].TSR_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].angvel_fixed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].Torque_generator_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GeneratorSpeed_desired), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GearBoxRatio), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_proportional), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_integral), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_derivative), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].Ki_IPC), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].Kratio_torque), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].angvel_axis_err_relax), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].WindSpeed_rated), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].WindSpeed_cutin), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].WindSpeed_cutout), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch[0]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch[1]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch[2]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch_min), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].CT), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].r_nearhubinflowcorrection), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].r_nacelle), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                PetscPrintf(PETSC_COMM_WORLD, "Number of blade for %d th turbine  %d \n", ibi, (wtm[ibi].num_blade));
                PetscPrintf(PETSC_COMM_WORLD, "Number of airfoil type for %d th turbine  %d \n", ibi, (wtm[ibi].num_foiltype));
                PetscPrintf(PETSC_COMM_WORLD, "The rotating center for turbine %d is %f %f %f \n", ibi, (fsi_wt[ibi].x_c), (fsi_wt[ibi].y_c), (fsi_wt[ibi].z_c));
                PetscPrintf(PETSC_COMM_WORLD, "Induction factor for %d th turbine  %f \n", ibi, (wtm[ibi].indf_axis));
                PetscPrintf(PETSC_COMM_WORLD, "Tipspeedratio for %d th turbine  %f \n", ibi, (wtm[ibi].Tipspeedratio));
                PetscPrintf(PETSC_COMM_WORLD, "Rotational inertial for %d th turbine  %f \n", ibi, (fsi_wt[ibi].J_rotation));
                PetscPrintf(PETSC_COMM_WORLD, "Radius of rotor for %d th turbine  %f \n", ibi, (fsi_wt[ibi].r_rotor));
                PetscPrintf(PETSC_COMM_WORLD, "Max Cp for %d th turbine  %f \n", ibi, (fsi_wt[ibi].CP_max));
                PetscPrintf(PETSC_COMM_WORLD, "Optimum TSR for %d th turbine  %f \n", ibi, (fsi_wt[ibi].TSR_max));
                PetscPrintf(PETSC_COMM_WORLD, "Fixed angvel for %d th turbine  %f \n", ibi, (fsi_wt[ibi].angvel_fixed));
                PetscPrintf(PETSC_COMM_WORLD, "Maximum Generator Torque for %d th turbine  %f \n", ibi, (fsi_wt[ibi].Torque_generator_max));
                PetscPrintf(PETSC_COMM_WORLD, "Desired Generator Speed for %d th turbine  %f \n", ibi, (fsi_wt[ibi].GeneratorSpeed_desired));
                PetscPrintf(PETSC_COMM_WORLD, "Gear Box Ratio for %d th turbine  %f \n", ibi, (fsi_wt[ibi].GearBoxRatio));
                PetscPrintf(PETSC_COMM_WORLD, "K_proportional for %d th turbine  %f \n", ibi, (fsi_wt[ibi].K_proportional));
                PetscPrintf(PETSC_COMM_WORLD, "K_integral for %d th turbine  %f \n", ibi, (fsi_wt[ibi].K_integral));
                PetscPrintf(PETSC_COMM_WORLD, "K_derivative for %d th turbine  %f \n", ibi, (fsi_wt[ibi].K_derivative));
                PetscPrintf(PETSC_COMM_WORLD, "Ki_IPC for %d th turbine  %f \n", ibi, (fsi_wt[ibi].Ki_IPC));
                PetscPrintf(PETSC_COMM_WORLD, "Kratio_torque for %d th turbine  %f \n", ibi, (fsi_wt[ibi].Kratio_torque));
                PetscPrintf(PETSC_COMM_WORLD, "angvel_axis_err_relax for %d th turbine  %f \n", ibi, (fsi_wt[ibi].angvel_axis_err_relax));
                PetscPrintf(PETSC_COMM_WORLD, "Rated wind speed for %d th turbine  %f \n", ibi, (fsi_wt[ibi].WindSpeed_rated));
                PetscPrintf(PETSC_COMM_WORLD, "Cut-in wind speed for %d th turbine  %f \n", ibi, (fsi_wt[ibi].WindSpeed_cutin));
                PetscPrintf(PETSC_COMM_WORLD, "Cut-out wind speed for %d th turbine  %f \n", ibi, (fsi_wt[ibi].WindSpeed_cutout));
                PetscPrintf(PETSC_COMM_WORLD, "initial pitch1 for %d th turbine  %f \n", ibi, (wtm[ibi].pitch[0]));
                PetscPrintf(PETSC_COMM_WORLD, "initial pitch2 for %d th turbine  %f \n", ibi, (wtm[ibi].pitch[1]));
                PetscPrintf(PETSC_COMM_WORLD, "initial pitch3 for %d th turbine  %f \n", ibi, (wtm[ibi].pitch[2]));
                PetscPrintf(PETSC_COMM_WORLD, "initial pitch_min for %d th turbine  %f \n", ibi, (wtm[ibi].pitch_min));
                PetscPrintf(PETSC_COMM_WORLD, "CT for %d th turbine  %f \n", ibi, (wtm[ibi].CT));
                PetscPrintf(PETSC_COMM_WORLD, "r_nearhubinflowcorrection for %d th turbine  %f \n", ibi, (wtm[ibi].r_nearhubinflowcorrection));
                PetscPrintf(PETSC_COMM_WORLD, "r_nacelle for %d th turbine  %f \n", ibi, (wtm[ibi].r_nacelle));
            }

            for (ibi=0;ibi<NumberOfTurbines;ibi++) {
                fsi_wt[ibi].x_c0=fsi_wt[ibi].x_c;
                fsi_wt[ibi].y_c0=fsi_wt[ibi].y_c;
                fsi_wt[ibi].z_c0=fsi_wt[ibi].z_c;
            }

            fclose(fd);
        }
        else {

            for (ibi=0;ibi<NumberOfTurbines;ibi++) {
                MPI_Bcast(&(wtm[ibi].num_blade), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].num_foiltype), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(wtm[ibi].pitch));
                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(wtm[ibi].pitch_old));
                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(wtm[ibi].pitch_IPC));
                PetscMalloc(wtm[ibi].num_blade*sizeof(PetscReal), &(fsi_wt[ibi].Moment_bladebending));

                wtm[ibi].pitch_old[0] = 0.0;
                wtm[ibi].pitch_old[1] = 0.0;
                wtm[ibi].pitch_old[2] = 0.0;

                MPI_Bcast(&(fsi_wt[ibi].YawCntrl), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].yawtistrt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(fsi_wt[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                fsi_wt[ibi].nx_tbo = 0.;
                fsi_wt[ibi].ny_tbo = 0.;
                fsi_wt[ibi].nz_tbo = 1.;

                MPI_Bcast(&(fsi_wt[ibi].TPCntrl), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_Gen), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GeneratorSyncSpeed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GeneratorTransSpeed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(fsi_wt[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].indf_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].Tipspeedratio), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].J_rotation), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].r_rotor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].CP_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].TSR_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].angvel_fixed), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].Torque_generator_max), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GeneratorSpeed_desired), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].GearBoxRatio), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_proportional), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_integral), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].K_derivative), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].Ki_IPC), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].Kratio_torque), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].angvel_axis_err_relax), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].WindSpeed_rated), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].WindSpeed_cutin), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_wt[ibi].WindSpeed_cutout), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(wtm[ibi].pitch[0]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch[1]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch[2]), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].pitch_min), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].CT), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].r_nearhubinflowcorrection), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(wtm[ibi].r_nacelle), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
            }

            for (ibi=0;ibi<NumberOfTurbines;ibi++) {
                fsi_wt[ibi].x_c0=fsi_wt[ibi].x_c;
                fsi_wt[ibi].y_c0=fsi_wt[ibi].y_c;
                fsi_wt[ibi].z_c0=fsi_wt[ibi].z_c;
            }
        }

        if (rotor_model == 5) {

            for (ibi=0;ibi<NumberOfTurbines;ibi++) {

                ibm_acl2ref[ibi].num_blade = wtm[ibi].num_blade;

                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(ibm_acl2ref[ibi].pitch));
                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(ibm_acl2ref[ibi].pitch_old));
                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(ibm_acl2ref[ibi].pitch_IPC));
                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(fsi_acl2ref[ibi].Moment_bladebending));

                ibm_acl2ref[ibi].pitch_old[0] = 0.0;
                ibm_acl2ref[ibi].pitch_old[1] = 0.0;
                ibm_acl2ref[ibi].pitch_old[2] = 0.0;

                ibm_acl2ref[ibi].num_foiltype = wtm[ibi].num_foiltype;

                fsi_acl2ref[ibi].x_c = fsi_wt[ibi].x_c;
                fsi_acl2ref[ibi].y_c = fsi_wt[ibi].y_c;
                fsi_acl2ref[ibi].z_c = fsi_wt[ibi].z_c;

                ibm_acl2ref[ibi].indf_axis = wtm[ibi].indf_axis;
                ibm_acl2ref[ibi].Tipspeedratio = wtm[ibi].Tipspeedratio;
                fsi_acl2ref[ibi].J_rotation = fsi_wt[ibi].J_rotation;
                fsi_acl2ref[ibi].r_rotor = fsi_wt[ibi].r_rotor;
                fsi_acl2ref[ibi].CP_max = fsi_wt[ibi].CP_max;
                fsi_acl2ref[ibi].TSR_max = fsi_wt[ibi].TSR_max;
                fsi_acl2ref[ibi].angvel_fixed = fsi_wt[ibi].angvel_fixed;
                fsi_acl2ref[ibi].Torque_generator_max = fsi_wt[ibi].Torque_generator_max;
                fsi_acl2ref[ibi].GeneratorSpeed_desired = fsi_wt[ibi].GeneratorSpeed_desired;
                fsi_acl2ref[ibi].GearBoxRatio = fsi_wt[ibi].GearBoxRatio;
                fsi_acl2ref[ibi].K_proportional = fsi_wt[ibi].K_proportional;
                fsi_acl2ref[ibi].K_integral = fsi_wt[ibi].K_integral;
                fsi_acl2ref[ibi].K_derivative = fsi_wt[ibi].K_derivative;
                fsi_acl2ref[ibi].Ki_IPC = fsi_wt[ibi].Ki_IPC;
                fsi_acl2ref[ibi].Kratio_torque = fsi_wt[ibi].Kratio_torque;
                fsi_acl2ref[ibi].angvel_axis_err_relax = fsi_wt[ibi].angvel_axis_err_relax;
                fsi_acl2ref[ibi].WindSpeed_rated = fsi_wt[ibi].WindSpeed_rated;
                fsi_acl2ref[ibi].WindSpeed_cutin = fsi_wt[ibi].WindSpeed_cutin;
                fsi_acl2ref[ibi].WindSpeed_cutout = fsi_wt[ibi].WindSpeed_cutout;
                ibm_acl2ref[ibi].pitch[0] = wtm[ibi].pitch[0];
                ibm_acl2ref[ibi].pitch[1] = wtm[ibi].pitch[1];
                ibm_acl2ref[ibi].pitch[2] = wtm[ibi].pitch[2];
                ibm_acl2ref[ibi].pitch_min = wtm[ibi].pitch_min;
                ibm_acl2ref[ibi].CT = wtm[ibi].CT;
                ibm_acl2ref[ibi].r_nearhubinflowcorrection = wtm[ibi].r_nearhubinflowcorrection;
                ibm_acl2ref[ibi].r_nacelle = wtm[ibi].r_nacelle;

                fsi_acl2ref[ibi].TPCntrl             = fsi_wt[ibi].TPCntrl;
                fsi_acl2ref[ibi].K_Gen               = fsi_wt[ibi].K_Gen;
                fsi_acl2ref[ibi].GeneratorSyncSpeed  = fsi_wt[ibi].GeneratorSyncSpeed;
                fsi_acl2ref[ibi].GeneratorTransSpeed = fsi_wt[ibi].GeneratorTransSpeed;

                fsi_acl2ref[ibi].x_c0=fsi_wt[ibi].x_c0;
                fsi_acl2ref[ibi].y_c0=fsi_wt[ibi].y_c0;
                fsi_acl2ref[ibi].z_c0=fsi_wt[ibi].z_c0;

                fsi_acl2ref[ibi].nx_tb=fsi_wt[ibi].nx_tb;
                fsi_acl2ref[ibi].ny_tb=fsi_wt[ibi].ny_tb;
                fsi_acl2ref[ibi].nz_tb=fsi_wt[ibi].nz_tb;

                fsi_acl2ref[ibi].nx_tbo=0.0;
                fsi_acl2ref[ibi].ny_tbo=0.0;
                fsi_acl2ref[ibi].nz_tbo=1.0;

                PetscPrintf(PETSC_COMM_WORLD, "Reference line file\n");

                PetscPrintf(PETSC_COMM_WORLD, "Number of blade for %d th turbine  %d \n", ibi, (ibm_acl2ref[ibi].num_blade));
                PetscPrintf(PETSC_COMM_WORLD, "Number of airfoil type for %d th turbine  %d \n", ibi, (ibm_acl2ref[ibi].num_foiltype));

                PetscPrintf(PETSC_COMM_WORLD, "The rotating center for %d th turbne %f %f %f \n", ibi, (fsi_acl2ref[ibi].x_c), (fsi_acl2ref[ibi].y_c), (fsi_acl2ref[ibi].z_c));

                PetscPrintf(PETSC_COMM_WORLD, "Induction factor for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].indf_axis));
                PetscPrintf(PETSC_COMM_WORLD, "Tipspeedratio for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].Tipspeedratio));
                PetscPrintf(PETSC_COMM_WORLD, "Rotational inertial for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].J_rotation));
                PetscPrintf(PETSC_COMM_WORLD, "Radius of rotor for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].r_rotor));
                PetscPrintf(PETSC_COMM_WORLD, "Max Cp for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].CP_max));
                PetscPrintf(PETSC_COMM_WORLD, "Optimum TSR for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].TSR_max));
                PetscPrintf(PETSC_COMM_WORLD, "Fixed angvel for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].angvel_fixed));
                PetscPrintf(PETSC_COMM_WORLD, "Maximum Generator Torque for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Torque_generator_max));
                PetscPrintf(PETSC_COMM_WORLD, "Desired Generator Speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].GeneratorSpeed_desired));
                PetscPrintf(PETSC_COMM_WORLD, "Gear Box Ratio for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].GearBoxRatio));

                PetscPrintf(PETSC_COMM_WORLD, "K_proportional for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].K_proportional));
                PetscPrintf(PETSC_COMM_WORLD, "K_integral for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].K_integral));
                PetscPrintf(PETSC_COMM_WORLD, "K_derivative for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].K_derivative));
                PetscPrintf(PETSC_COMM_WORLD, "Ki_IPC for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Ki_IPC));
                PetscPrintf(PETSC_COMM_WORLD, "Kratio_torque for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Kratio_torque));
                PetscPrintf(PETSC_COMM_WORLD, "angvel_axis_err_relax for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].angvel_axis_err_relax));

                PetscPrintf(PETSC_COMM_WORLD, "Rated wind speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].WindSpeed_rated));
                PetscPrintf(PETSC_COMM_WORLD, "Cut-in wind speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].WindSpeed_cutin));
                PetscPrintf(PETSC_COMM_WORLD, "Cut-out wind speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].WindSpeed_cutout));

                PetscPrintf(PETSC_COMM_WORLD, "pitch1 for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[0]));
                PetscPrintf(PETSC_COMM_WORLD, "pitch2 for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[1]));
                PetscPrintf(PETSC_COMM_WORLD, "pitch3 for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[2]));
                PetscPrintf(PETSC_COMM_WORLD, "pitch_min for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch_min));

                PetscPrintf(PETSC_COMM_WORLD, "CT for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].CT));
                PetscPrintf(PETSC_COMM_WORLD, "r_nearhubinflowcorrection for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].r_nearhubinflowcorrection));
                PetscPrintf(PETSC_COMM_WORLD, "r_nacelle for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].r_nacelle));

            }
        }

        if (rotor_model == 6) {

            for (ibi=0;ibi<NumberOfTurbines;ibi++) {

                ibm_acl2ref[ibi].num_blade = wtm[ibi].num_blade;
                ibm_acl2ref[ibi].num_foiltype = wtm[ibi].num_foiltype;

                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(ibm_acl2ref[ibi].pitch));
                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(ibm_acl2ref[ibi].pitch_old));
                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(ibm_acl2ref[ibi].pitch_IPC));
                PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(PetscReal), &(fsi_acl2ref[ibi].Moment_bladebending));

                ibm_acl2ref[ibi].pitch_old[0] = 0.0;
                ibm_acl2ref[ibi].pitch_old[1] = 0.0;
                ibm_acl2ref[ibi].pitch_old[2] = 0.0;

                fsi_acl2ref[ibi].x_c = fsi_wt[ibi].x_c;
                fsi_acl2ref[ibi].y_c = fsi_wt[ibi].y_c;
                fsi_acl2ref[ibi].z_c = fsi_wt[ibi].z_c;

                ibm_acl2ref[ibi].indf_axis = wtm[ibi].indf_axis;
                ibm_acl2ref[ibi].Tipspeedratio = wtm[ibi].Tipspeedratio;
                fsi_acl2ref[ibi].J_rotation = fsi_wt[ibi].J_rotation;
                fsi_acl2ref[ibi].r_rotor = fsi_wt[ibi].r_rotor;
                fsi_acl2ref[ibi].CP_max = fsi_wt[ibi].CP_max;
                fsi_acl2ref[ibi].TSR_max = fsi_wt[ibi].TSR_max;
                fsi_acl2ref[ibi].angvel_fixed = fsi_wt[ibi].angvel_fixed;
                fsi_acl2ref[ibi].Torque_generator_max = fsi_wt[ibi].Torque_generator_max;
                fsi_acl2ref[ibi].GeneratorSpeed_desired = fsi_wt[ibi].GeneratorSpeed_desired;
                fsi_acl2ref[ibi].GearBoxRatio = fsi_wt[ibi].GearBoxRatio;
                fsi_acl2ref[ibi].K_proportional = fsi_wt[ibi].K_proportional;
                fsi_acl2ref[ibi].K_integral = fsi_wt[ibi].K_integral;
                fsi_acl2ref[ibi].K_derivative = fsi_wt[ibi].K_derivative;
                fsi_acl2ref[ibi].Ki_IPC = fsi_wt[ibi].Ki_IPC;
                fsi_acl2ref[ibi].Kratio_torque = fsi_wt[ibi].Kratio_torque;
                fsi_acl2ref[ibi].angvel_axis_err_relax = fsi_wt[ibi].angvel_axis_err_relax;
                fsi_acl2ref[ibi].WindSpeed_rated = fsi_wt[ibi].WindSpeed_rated;
                fsi_acl2ref[ibi].WindSpeed_cutin = fsi_wt[ibi].WindSpeed_cutin;
                fsi_acl2ref[ibi].WindSpeed_cutout = fsi_wt[ibi].WindSpeed_cutout;


                ibm_acl2ref[ibi].pitch[0] = wtm[ibi].pitch[0];
                ibm_acl2ref[ibi].pitch[1] = wtm[ibi].pitch[1];
                ibm_acl2ref[ibi].pitch[2] = wtm[ibi].pitch[2];
                ibm_acl2ref[ibi].pitch_min = wtm[ibi].pitch_min;
                ibm_acl2ref[ibi].CT = wtm[ibi].CT;
                ibm_acl2ref[ibi].r_nearhubinflowcorrection = wtm[ibi].r_nearhubinflowcorrection;
                ibm_acl2ref[ibi].r_nacelle = wtm[ibi].r_nacelle;

                fsi_acl2ref[ibi].x_c0=fsi_wt[ibi].x_c0;
                fsi_acl2ref[ibi].y_c0=fsi_wt[ibi].y_c0;
                fsi_acl2ref[ibi].z_c0=fsi_wt[ibi].z_c0;

                fsi_acl2ref[ibi].nx_tb=fsi_wt[ibi].nx_tb;
                fsi_acl2ref[ibi].ny_tb=fsi_wt[ibi].ny_tb;
                fsi_acl2ref[ibi].nz_tb=fsi_wt[ibi].nz_tb;

                fsi_acl2ref[ibi].nx_tbo=0.0;
                fsi_acl2ref[ibi].ny_tbo=0.0;
                fsi_acl2ref[ibi].nz_tbo=1.0;

                fsi_acl2ref[ibi].TPCntrl             = fsi_wt[ibi].TPCntrl;
                fsi_acl2ref[ibi].K_Gen               = fsi_wt[ibi].K_Gen;
                fsi_acl2ref[ibi].GeneratorSyncSpeed  = fsi_wt[ibi].GeneratorSyncSpeed;
                fsi_acl2ref[ibi].GeneratorTransSpeed = fsi_wt[ibi].GeneratorTransSpeed;

                PetscPrintf(PETSC_COMM_WORLD, "Reference line file\n");

                PetscPrintf(PETSC_COMM_WORLD, "Number of blade for %d th turbine  %d \n", ibi, (ibm_acl2ref[ibi].num_blade));
                PetscPrintf(PETSC_COMM_WORLD, "Number of airfoil type for %d th turbine  %d \n", ibi, (ibm_acl2ref[ibi].num_foiltype));

                PetscPrintf(PETSC_COMM_WORLD, "The rotating center for %d th turbne %f %f %f \n", ibi, (fsi_acl2ref[ibi].x_c), (fsi_acl2ref[ibi].y_c), (fsi_acl2ref[ibi].z_c));
                PetscPrintf(PETSC_COMM_WORLD, "Induction factor for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].indf_axis));
                PetscPrintf(PETSC_COMM_WORLD, "Tipspeedratio for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].Tipspeedratio));
                PetscPrintf(PETSC_COMM_WORLD, "Rotational inertial for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].J_rotation));
                PetscPrintf(PETSC_COMM_WORLD, "Radius of rotor for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].r_rotor));
                PetscPrintf(PETSC_COMM_WORLD, "Max Cp for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].CP_max));
                PetscPrintf(PETSC_COMM_WORLD, "Optimum TSR for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].TSR_max));
                PetscPrintf(PETSC_COMM_WORLD, "Fixed angvel for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].angvel_fixed));
                PetscPrintf(PETSC_COMM_WORLD, "Maximum Generator Torque for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Torque_generator_max));
                PetscPrintf(PETSC_COMM_WORLD, "Desired Generator Speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].GeneratorSpeed_desired));
                PetscPrintf(PETSC_COMM_WORLD, "Gear Box Ratio for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].GearBoxRatio));

                PetscPrintf(PETSC_COMM_WORLD, "K_proportional for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].K_proportional));
                PetscPrintf(PETSC_COMM_WORLD, "K_integral for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].K_integral));
                PetscPrintf(PETSC_COMM_WORLD, "K_derivative for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].K_derivative));
                PetscPrintf(PETSC_COMM_WORLD, "Ki_IPC for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Ki_IPC));
                PetscPrintf(PETSC_COMM_WORLD, "Kratio_torque for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].Kratio_torque));
                PetscPrintf(PETSC_COMM_WORLD, "angvel_axis_err_relax for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].angvel_axis_err_relax));

                PetscPrintf(PETSC_COMM_WORLD, "Rated wind speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].WindSpeed_rated));
                PetscPrintf(PETSC_COMM_WORLD, "Cut-in wind speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].WindSpeed_cutin));
                PetscPrintf(PETSC_COMM_WORLD, "Cut-out wind speed for %d th turbine  %f \n", ibi, (fsi_acl2ref[ibi].WindSpeed_cutout));

                PetscPrintf(PETSC_COMM_WORLD, "pitch1 for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[0]));
                PetscPrintf(PETSC_COMM_WORLD, "pitch2 for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[1]));
                PetscPrintf(PETSC_COMM_WORLD, "pitch3 for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch[2]));
                PetscPrintf(PETSC_COMM_WORLD, "pitch_min for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].pitch_min));
                PetscPrintf(PETSC_COMM_WORLD, "CT for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].CT));
                PetscPrintf(PETSC_COMM_WORLD, "r_nearhubinflowcorrection for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].r_nearhubinflowcorrection));
                PetscPrintf(PETSC_COMM_WORLD, "r_nacelle for %d th turbine  %f \n", ibi, (ibm_acl2ref[ibi].r_nacelle));
            }
        }

        PetscBarrier(PETSC_NULL);

        if (rotor_model==2 || rotor_model==3 || rotor_model==5 || rotor_model==6) {
            for (i=0;i<NumberOfTurbines;i++) {
                PetscMalloc(wtm[i].num_foiltype*sizeof(FOIL), &wtm[i].acl);
            }
        }

        if (rotor_model==5) {
            for (i=0;i<NumberOfTurbines;i++) {
                PetscMalloc(ibm_acl2ref[i].num_foiltype*sizeof(FOIL), &ibm_acl2ref[i].acl);
            }
        }

        ti = 0;
        if (rstart_flg) ti = tistart;

//              for (i=0;i<NumberOfTurbines;i++) {

        PetscPrintf(PETSC_COMM_WORLD, "Turbines read!\n");
    
        PetscPrintf(PETSC_COMM_WORLD, "Rotor Model = %i\n", rotor_model);
    
        if (rotor_model == 1) {
            double reflength = reflength_wt;
            char fname[80];
            sprintf(fname,"acddata");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
            Pre_process(&(user[0]), wtm, NumberOfTurbines); 
        }

        if (rotor_model == 2) {
            double reflength = reflength_wt;
            char fname[80];
            sprintf(fname,"acl2data");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
            Pre_process(&(user[0]), wtm, NumberOfTurbines); 

            sprintf(fname,"Urefdata");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
            Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
            airfoil_ACL(wtm,  fsi_wt, NumberOfTurbines);
            Uref_ACL(user, wtm, ibm_ACD, fsi_wt, NumberOfTurbines);
            TurbineTorqueControlInitialization(fsi_wt, wtm);
        }

        if (rotor_model == 3) {
            double reflength = reflength_wt;
            for (i=0;i<NumberOfTurbines;i++) 
                ACL_read(&wtm[i], i, &fsi_wt[i], reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
            Pre_process(&(user[0]), wtm, NumberOfTurbines); 

            char fname[80];
            sprintf(fname,"Urefdata");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
            Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);

            airfoil_ACL(wtm,  fsi_wt, NumberOfTurbines);
            Uref_ACL(user, wtm, ibm_ACD, fsi_wt, NumberOfTurbines);
            TurbineTorqueControlInitialization(fsi_wt, wtm);
        }

        if (rotor_model == 4) {
            double reflength = reflength_wt;
            char fname[80];
            sprintf(fname,"acddata");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&wtm[i], i, &fsi_wt[i], 0, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
            Pre_process(&(user[0]), wtm, NumberOfTurbines); 

            sprintf(fname,"Urefdata");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
            Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);
        }

        if (rotor_model == 5 ) {
            double reflength = reflength_wt;
            //ACD_read(&wtm[i], i, &fsi_wt[i], 0);
            char fname[80], ascfname[80];
            
            PetscPrintf(PETSC_COMM_WORLD, "Begin reading actuator surface!\n");
            sprintf(ascfname,"acsdata");
            
            for (i=0;i<NumberOfTurbines;i++)
                surface_read(1, &wtm[i], i, &fsi_wt[i], ascfname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Actuator surface pre-processing!\n");
            Pre_process(&(user[0]), wtm, NumberOfTurbines);  

            PetscPrintf(PETSC_COMM_WORLD, "\nBegin reading Uref!\n");
            sprintf(fname,"Urefdata");
            for (i=0;i<NumberOfTurbines;i++)
                disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
            Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);

            PetscPrintf(PETSC_COMM_WORLD, "\nBegin reading actuator lines!\n");
            for (i=0;i<NumberOfTurbines;i++)
                ACL_read(&ibm_acl2ref[i], i, &fsi_acl2ref[i], reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Actuator line pre-processing!\n");
            Pre_process(&(user[0]), ibm_acl2ref, NumberOfTurbines); // 20140807

            PetscPrintf(PETSC_COMM_WORLD, "\nBegin adding colors to  ACS to matach ACL !\n");
            for (i=0;i<NumberOfTurbines;i++)
                matchLineColorToSurface(i, &wtm[i], &ibm_acl2ref[i]);
                        // ascdata write need to be here since matchLineColorToSurface
                        // updates the color information.
            
            PetscPrintf(PETSC_COMM_WORLD, "Writing out actuator surface file!\n");
            for (i=0;i<NumberOfTurbines;i++)
                surface_write(&wtm[i], i, ascfname);

            PetscPrintf(PETSC_COMM_WORLD, "Matching  ACS to ACL!\n");
            calc_s2l(wtm, ibm_acl2ref, fsi_wt, fsi_acl2ref, NumberOfTurbines);
    
            airfoil_ACL(wtm,  fsi_wt, NumberOfTurbines);
            airfoil_ACL(ibm_acl2ref,  fsi_acl2ref, NumberOfTurbines);

            TurbineTorqueControlInitialization(fsi_acl2ref, ibm_acl2ref);

            for (i=0;i<NumberOfTurbines;i++){
                if(fsi_wt[i].YawCntrl == 1){
                    Yaw_Control(ti, i,  fsi_wt );
                    fsi_acl2ref[i].nx_tb=fsi_wt[i].nx_tb;
                    fsi_acl2ref[i].ny_tb=fsi_wt[i].ny_tb;
                    fsi_acl2ref[i].nz_tb=fsi_wt[i].nz_tb;
                }
               //rotorYaw(&wtm[ibi], &fsi_wt[ibi]);
               //rotorYaw(&ibm_acl2ref[ibi], &fsi_acl2ref[ibi]);
            }
        }

        if (rotor_model == 6) {
            double reflength = reflength_wt;
            for (i=0;i<NumberOfTurbines;i++) 
                ACL_read(&wtm[i], i, &fsi_wt[i], reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "rotor model pre-processing!\n");
            Pre_process(&(user[0]), wtm, NumberOfTurbines); 

            char fname[80];
            sprintf(fname,"Urefdata");
            for (i=0;i<NumberOfTurbines;i++) 
                disk_read_ucd(&ibm_ACD[i], i, &fsi_wt[i], 1, fname, reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Uref disk pre-processing!\n");
            Pre_process(&(user[0]), ibm_ACD, NumberOfTurbines);

            airfoil_ACL(wtm,  fsi_wt, NumberOfTurbines);


            PetscPrintf(PETSC_COMM_WORLD, "Read uref line!\n");
            for (i=0;i<NumberOfTurbines;i++) 
                ACL_read(&ibm_acl2ref[i], i, &fsi_acl2ref[i], reflength);
            
            PetscPrintf(PETSC_COMM_WORLD, "Uref line pre-processing!\n");
            Pre_process(&(user[0]), ibm_acl2ref, NumberOfTurbines); // 20140807

            TurbineTorqueControlInitialization(fsi_wt, wtm);

            for (i=0;i<NumberOfTurbines;i++){
                if(fsi_wt[i].YawCntrl == 1){
                    Yaw_Control(ti, i,  fsi_wt );
                    fsi_acl2ref[i].nx_tb=fsi_wt[i].nx_tb;
                    fsi_acl2ref[i].ny_tb=fsi_wt[i].ny_tb;
                    fsi_acl2ref[i].nz_tb=fsi_wt[i].nz_tb;
                }
            }
        }

            //Hossein added from turbine structure
            if (turbinestructuremodel) {

			switch(rotor_model) 
			{ 
				case 5: 
					for (ibi=0;ibi<NumberOfTurbines;ibi++) {
                         			PetscMalloc(ibm_acl2ref[ibi].num_blade*sizeof(Beam), &(ibm_acl2ref[ibi].bladestructure));  
					}
					init_bladestructure(ibm_acl2ref);

					dt_turbinestructure = user[0].dt*reflength_wt/refvel_wt;
					PetscOptionsGetReal(PETSC_NULL, "-dt_turbinestructure", &dt_turbinestructure, PETSC_NULL); 

					break;
				case 6: 
					for (ibi=0;ibi<NumberOfTurbines;ibi++) {
                         			PetscMalloc(wtm[ibi].num_blade*sizeof(Beam), &(wtm[ibi].bladestructure));  
					}
					init_bladestructure(wtm);

					dt_turbinestructure = user[0].dt*reflength_wt/refvel_wt;

					PetscOptionsGetReal(PETSC_NULL, "-dt_turbinestructure", &dt_turbinestructure, PETSC_NULL); 

                                    if (restart_turbinestructure){
                	                PetscPrintf(PETSC_COMM_WORLD, "Displacing blade elements on restart\n");
						 transferbladedisplacement2flowblade(wtm, fsi_wt);
                	                PetscPrintf(PETSC_COMM_WORLD, "Displacing blade elements on restart DONE\n");
                                        }        

					break;
			}
		}
    }

//Meric
//Canopy
    if(cnpy) {
      CanopyInitialization(&user[0], NumberOfCnpy, cnpyHeightDim, cnpyCd);
    }


    if (nacelle_model) {

        PetscReal cl = 1.;
        PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);

        if (!my_rank) {
            FILE *fd;
            char str[256];
            sprintf(str, "%s/Nacelle.inp", path);
            fd = fopen(str, "r");
            if(!fd) PetscPrintf(PETSC_COMM_WORLD, "cannot open %s !\n", str),exit(0);

            char string[256];
            fgets(string, 256, fd);

            for (ibi=0;ibi<NumberOfNacelle;ibi++) {
                PetscPrintf(PETSC_COMM_WORLD, "\n\nBegin reading of the Nacelle.inp file for Nacelle %i!/n", ibi);

                fscanf(fd, "%le %le %le %le %le %le %le %d %le %le %le %le %le %le", &(fsi_nacelle[ibi].nx_tb), &(fsi_nacelle[ibi].ny_tb), &(fsi_nacelle[ibi].nz_tb), &(fsi_nacelle[ibi].x_c), &(fsi_nacelle[ibi].y_c), &(fsi_nacelle[ibi].z_c), &(fsi_nacelle[ibi].angvel_axis), &(fsi_nacelle[ibi].rotate_alongaxis), &(ibm_nacelle[ibi].friction_factor), &(ibm_nacelle[ibi].dh), &(fsi_nacelle[ibi].xnacelle_upstreamend), &(fsi_nacelle[ibi].ynacelle_upstreamend), &(fsi_nacelle[ibi].znacelle_upstreamend), &(ibm_nacelle[ibi].r_nacelle));

                double rr=sqrt(pow(fsi_nacelle[ibi].nx_tb,2)+pow(fsi_nacelle[ibi].ny_tb,2)+pow(fsi_nacelle[ibi].nz_tb,2));

                fsi_nacelle[ibi].nx_tb=fsi_nacelle[ibi].nx_tb/rr; 
                fsi_nacelle[ibi].ny_tb=fsi_nacelle[ibi].ny_tb/rr; 
                fsi_nacelle[ibi].nz_tb=fsi_nacelle[ibi].nz_tb/rr;
                
                MPI_Bcast(&(fsi_nacelle[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                PetscPrintf(PETSC_COMM_WORLD, "The directions of nacelle %f %f %f \n", (fsi_nacelle[ibi].nx_tb), (fsi_nacelle[ibi].ny_tb), (fsi_nacelle[ibi].nz_tb));

                fsi_nacelle[ibi].x_c=fsi_nacelle[ibi].x_c/cl;
                fsi_nacelle[ibi].y_c=fsi_nacelle[ibi].y_c/cl;
                fsi_nacelle[ibi].z_c=fsi_nacelle[ibi].z_c/cl;

                MPI_Bcast(&(fsi_nacelle[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].angvel_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].rotate_alongaxis), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(ibm_nacelle[ibi].friction_factor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(ibm_nacelle[ibi].dh), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);


                MPI_Bcast(&(fsi_nacelle[ibi].xnacelle_upstreamend), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].ynacelle_upstreamend), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].znacelle_upstreamend), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                PetscPrintf(PETSC_COMM_WORLD, "The upstream end coordinates of nacelle %f %f %f \n", (fsi_nacelle[ibi].xnacelle_upstreamend), (fsi_nacelle[ibi].ynacelle_upstreamend), (fsi_nacelle[ibi].znacelle_upstreamend));

                MPI_Bcast(&(ibm_nacelle[ibi].r_nacelle), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                PetscPrintf(PETSC_COMM_WORLD, "The Locations for %d th  nacelle %f %f %f \n", ibi, (fsi_nacelle[ibi].x_c), (fsi_nacelle[ibi].y_c),(fsi_nacelle[ibi].z_c));
                PetscPrintf(PETSC_COMM_WORLD, "The angular velocity for %d th nacelle body %f \n", ibi, (fsi_nacelle[ibi].angvel_axis));
                PetscPrintf(PETSC_COMM_WORLD, "Rotate the nacelle along the axis? %d \n",(fsi_nacelle[ibi].rotate_alongaxis) );
                PetscPrintf(PETSC_COMM_WORLD, "The friction factor of %d the nacelle %f \n", ibi, (ibm_nacelle[ibi].friction_factor));
                PetscPrintf(PETSC_COMM_WORLD, "The wall-normal thickness of %d th nacelle mesh %f \n", ibi, (ibm_nacelle[ibi].dh));
                PetscPrintf(PETSC_COMM_WORLD, "The radius of %d th nacelle %f \n", ibi, (ibm_nacelle[ibi].r_nacelle));

            }

            for (ibi=0;ibi<NumberOfNacelle;ibi++) {
                fsi_nacelle[ibi].x_c0=fsi_nacelle[ibi].x_c;
                fsi_nacelle[ibi].y_c0=fsi_nacelle[ibi].y_c;
                fsi_nacelle[ibi].z_c0=fsi_nacelle[ibi].z_c;
            }

            fclose(fd);

            if (cf_nacelle_fromfile) {
                for (ibi=0;ibi<NumberOfNacelle;ibi++) {
                    FILE *fd;
                    char str[256];
                    sprintf(str, "%s/cf_nacelle%03d", path, ibi);
                    fd = fopen(str, "r");
                    if(!fd) PetscPrintf(PETSC_COMM_WORLD, "cannot open %s !\n", str),exit(0);

                    fscanf(fd, "%d \n", &(ibm_nacelle[ibi].num_cf));
                    MPI_Bcast(&(ibm_nacelle[ibi].num_cf), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                    int i;
                    for (i=0;i<ibm_nacelle[ibi].num_cf;i++) {
                        fscanf(fd, "%le %le \n", &(ibm_nacelle[ibi].r_in[i]), &(ibm_nacelle[ibi].cf_in[i]));

                        ibm_nacelle[ibi].r_in[i] = ibm_nacelle[ibi].r_in[i]/reflength_nacelle;
                    }

                    MPI_Bcast(&(ibm_nacelle[ibi].r_in[0]), ibm_nacelle[ibi].num_cf, MPIU_REAL, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(&(ibm_nacelle[ibi].cf_in[0]), ibm_nacelle[ibi].num_cf, MPIU_REAL, 0, PETSC_COMM_WORLD);
                }
            }
        }
        else {

            for (ibi=0;ibi<NumberOfNacelle;ibi++) {
                MPI_Bcast(&(fsi_nacelle[ibi].nx_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].ny_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].nz_tb), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(fsi_nacelle[ibi].x_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].y_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].z_c), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].angvel_axis), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].rotate_alongaxis), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(ibm_nacelle[ibi].friction_factor), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(ibm_nacelle[ibi].dh), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

                MPI_Bcast(&(fsi_nacelle[ibi].xnacelle_upstreamend), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].ynacelle_upstreamend), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(fsi_nacelle[ibi].znacelle_upstreamend), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(&(ibm_nacelle[ibi].r_nacelle), 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
            }

            for (ibi=0;ibi<NumberOfNacelle;ibi++) {
                fsi_nacelle[ibi].x_c0=fsi_nacelle[ibi].x_c;
                fsi_nacelle[ibi].y_c0=fsi_nacelle[ibi].y_c;
                fsi_nacelle[ibi].z_c0=fsi_nacelle[ibi].z_c;
            }

            if (cf_nacelle_fromfile) {
                for (ibi=0;ibi<NumberOfNacelle;ibi++) {
                    MPI_Bcast(&(ibm_nacelle[ibi].num_cf), 1, MPI_INT, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(&(ibm_nacelle[ibi].r_in[0]), ibm_nacelle[ibi].num_cf, MPIU_REAL, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(&(ibm_nacelle[ibi].cf_in[0]), ibm_nacelle[ibi].num_cf, MPIU_REAL, 0, PETSC_COMM_WORLD);
                }
            }
        }

        PetscPrintf(PETSC_COMM_WORLD, "\n\nBegin reading nacelle000_ UCD file\n");

        double reflength = reflength_nacelle;

        int ipt;
        int NumLoc=(int)NumberOfNacelle/(int)NumNacellePerLoc;
        for (ibi=0;ibi<NumLoc;ibi++) 
        for (ipt=0;ipt<(int)NumNacellePerLoc;ipt++) { 
            int iname=ibi*(int)NumNacellePerLoc+ipt; 
            char fname[80];
            sprintf(fname,"nacelle%3.3d_", ipt);
            disk_read_ucd(&ibm_nacelle[iname], iname, &fsi_nacelle[iname], 0, fname, reflength);    
        }

        PetscPrintf(PETSC_COMM_WORLD, "Nacelle preprocess!\n");
        Pre_process(&(user[0]), ibm_nacelle, NumberOfNacelle);

        PetscPrintf(PETSC_COMM_WORLD, "IP points!\n");
        Coordinates_IP(ibm_nacelle, NumberOfNacelle); 
        PetscPrintf(PETSC_COMM_WORLD, "IP points Pre-process!\n");
        Pre_process_IP(&(user[0]), ibm_nacelle, NumberOfNacelle); 

        //PetscPrintf(PETSC_COMM_WORLD, "Color IB points for zero force!\n");
        //ColorIB_tmp(ibm_nacelle, NumberOfNacelle); //Hossein added from 3.4 code to fix turbine rotation restart issue (temporary)

        ti = 0;
        if (rstart_flg) ti = tistart;

       PetscBarrier(PETSC_NULL);
        PetscPrintf(PETSC_COMM_WORLD, "NacelleYaw Test\n");
        PetscPrintf(PETSC_COMM_WORLD, "Num of Nacelle %d\n", NumberOfNacelle);
        for (ibi=0;ibi<NumberOfNacelle;ibi++) { 
            fsi_nacelle[ibi].nx_tbo = 0.;//orientation assumes nacelle000_ is normal to x-axis
            fsi_nacelle[ibi].ny_tbo = 0.;
            fsi_nacelle[ibi].nz_tbo = 1.;
        }

        PetscBarrier(PETSC_NULL);
        PetscBarrier(PETSC_NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Begin Yaw of Nacelle\n");
        PetscBarrier(PETSC_NULL);

        nacelleYaw_IB(&(user[0]), ibm_nacelle, fsi_nacelle, NumberOfNacelle);

        PetscBarrier(PETSC_NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Finished Yaw of Nacelle\n");

        PetscBarrier(PETSC_NULL);

    }
 
    // rstart not working now
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    if (rstart_flg) {
        ti = tistart; tistart++;

        for (bi=0; bi<block_number; bi++) {
            Ucont_P_Binary_Input(&(user[bi]));
            //!!!!!!!!!!!!!!temp change
            //VecSet(user[bi].Nvert_o,0.);

            DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat);
            DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat);

            DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
            DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

            DAGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
            DAGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);

            DAGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
            DAGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);

            Contra2Cart(&(user[bi]));

            if (rstart_fsi) {
                for (ibi=0;ibi<NumberOfBodies;ibi++) {

                    if(!rotatefsi) FSI_DATA_Input(&fsi[ibi],ibi);

                    if (movefsi) {
                        Elmt_Move_FSI_TRANS(&(user[bi]), &fsi[ibi], &ibm[ibi]);  //Hossein   
                        for (i=0;i<6;i++){
                            fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
                            fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
                        }
                        for (i=0; i<ibm[ibi].n_v; i++) {
                            ibm[ibi].uold[i].x = fsi[ibi].S_real[1];
                            ibm[ibi].uold[i].y = fsi[ibi].S_real[3];
                            ibm[ibi].uold[i].z = fsi[ibi].S_real[5];
                        }
                        for (i=0; i<ibm[ibi].n_v; i++) {
                            ibm[ibi].urm1[i].x = fsi[ibi].S_realm1[1];
                            ibm[ibi].urm1[i].y = fsi[ibi].S_realm1[3];
                            ibm[ibi].urm1[i].z = fsi[ibi].S_realm1[5];
                        }
                    }
                    if (rotatefsi|| MHV) {
                        fsi[ibi].x_c = x_r; //seokkoo
                        fsi[ibi].y_c = y_r;
                        fsi[ibi].z_c = z_r;
                        if(ibi==0) Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
                    }
                    else {
                        for (i=0; i<ibm[ibi].n_v; i++) {
                            ibm[ibi].u[i].x = 0;
                            ibm[ibi].u[i].y = 0;
                            ibm[ibi].u[i].z = 0;
                            ibm[ibi].uold[i] = ibm[ibi].u[i];
                            ibm[ibi].urm1[i] = ibm[ibi].u[i];
                        }
                    }

                    // if read ti, then will start for ti+1
                    for (i=0;i<6;i++){
                        fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
                        fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];          
                    }

                    fsi[ibi].F_x_real=fsi[ibi].F_x;
                    fsi[ibi].F_y_real=fsi[ibi].F_y;
                    fsi[ibi].F_z_real=fsi[ibi].F_z;
             
                    fsi[ibi].M_x_rm3=fsi[ibi].M_x;
                    fsi[ibi].M_y_rm3=fsi[ibi].M_y;
                    fsi[ibi].M_z_rm3=fsi[ibi].M_z;
                    
                    fsi[ibi].M_x_rm2=fsi[ibi].M_x;
                    fsi[ibi].M_y_rm2=fsi[ibi].M_y;
                    fsi[ibi].M_z_rm2=fsi[ibi].M_z;

                    fsi[ibi].M_x_real=fsi[ibi].M_x;
                    fsi[ibi].M_y_real=fsi[ibi].M_y;
                    fsi[ibi].M_z_real=fsi[ibi].M_z;
                }//ibi
            }// if rstart fsi
        } 
    }// bi

// do the search once if elmt is not moving!
    if (immersed) {
        for (level = usermg.mglevels-1; level>=usermg.mglevels-1; level--) {
            user = usermg.mgctx[level].user;
            for (bi=0; bi<block_number; bi++) {
                for (ibi=0;ibi<NumberOfBodies;ibi++) {
                    PetscPrintf(PETSC_COMM_WORLD, "IBM%d Advanced Initial Search from main.c \n", ibi);
                    ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);// 

                    if(sediment){
                        PetscReal ts, te, cputime;     // xiaolei
                        PetscGetTime(&ts);  // xiaolei

                        if(ibi == 0) Connectivity_ib(&(user[bi]), &ibm[ibi]);  // xiaolei add SEDI //Ali Corrected
                        PetscGetTime(&te);  // xiaolei
                        PetscPrintf(PETSC_COMM_WORLD, "Time for creating IB connectivity: %le\n", te-ts);
                    }
                }
                PetscBarrier(PETSC_NULL);

                PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP\n");
                ibm_interpolation_advanced(&user[bi]);
                PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP_ends\n");//}
            }
        }
    }

    // Copy Ucont to Ucont_o for the finest level
    for (bi=0; bi<block_number; bi++) {
        //VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o));
        ti = 0;
        if (rstart_flg) ti = tistart;

        if(ti==tistart && ti==0 && levelset) {
            Levelset_Function_IC(&user[bi]);
            DAGlobalToLocalBegin(user[bi].da, user[bi].Levelset, INSERT_VALUES, user[bi].lLevelset);
            DAGlobalToLocalEnd(user[bi].da, user[bi].Levelset, INSERT_VALUES, user[bi].lLevelset);
            VecCopy(user[bi].Levelset, user[bi].Levelset_o);
        }

        if(ti==tistart) Calc_Inlet_Area(&user[bi]);

        if (ti==0) {
            VecSet(user[bi].Ucont,0.);
            VecSet(user[bi].lUcont,0.);
            VecSet(user[bi].Ucont_o,0.);
            VecSet(user[bi].lUcont_o,0.);
            VecSet(user[bi].Ucat,0.);
            VecSet(user[bi].lUcat,0.);
            VecSet(user[bi].P,0.);
            VecSet(user[bi].lP,0.);
        
            if(initialzero) PetscPrintf(PETSC_COMM_WORLD, "\nInitial Guess is Zero !\n");
            else SetInitialGuessToOne(&(user[bi]));
    
            Contra2Cart(&(user[bi]));
            DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
            DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
        }
      
        VecCopy(user[bi].Ucont, user[bi].Ucont_o);
        //VecCopy(user[bi].Ucont, user[bi].Ucont_rm2);    // allocate at init.c
        VecCopy(user[bi].Ucont, user[bi].Ucont_rm1);
        
        VecCopy(user[bi].Ucat, user[bi].Ucat_o);
        VecCopy(user[bi].P, user[bi].P_o);

        DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
        DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);

        DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
        DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
    }

    PetscInt tisteps = 1000000;
    PetscOptionsGetInt(PETSC_NULL, "-totalsteps", &tisteps, &flg);

    if (tistart==0) tisteps ++;

    /*  put the time accuracy coefficient 1 for the 1st real-time step */
    /*   COEF_TIME_ACCURACY=1.; */

    /*   PreLoadBegin(PETSC_TRUE, "Load"); */

    /* ==================================================================================             */
    /*   pysical time Step Loop */
    for (ti = tistart; ti<tistart + tisteps; ti++) {
        PetscPrintf(PETSC_COMM_WORLD, "Time %d\n", ti);

        if (inletprofile==3) {
            if (MHV && (fsi[1].S_ang_n[0]<0.8*max_angle || fsi[2].S_ang_n[0]>-0.8*max_angle)) 
                angle=angle+1;
            else
                angle=0.;

            fluxin(&(usermg.mgctx[usermg.mglevels-1].user[0]));
        }
        /* ==================================================================================             */
        /*     Strong-Coupling (SC) Loop */
        DoSCLoop= PETSC_TRUE ; itr_sc = 0;
        while (DoSCLoop) {
            itr_sc++;
            PetscPrintf(PETSC_COMM_WORLD, "SC LOOP itr # %d\n", itr_sc);

            /*     Structral Solver! */
            if (immersed){
                //PetscPrintf(PETSC_COMM_WORLD, "IBM 0 Axis:%i at XYAngle:%f and Coord. Angles Rotate: %f %f %f, Center of Rot: %f %f %f\n", fsi[2].rot_dir,fsi[2].XYangle ,fsi[2].S_ang_n[0] , fsi[2].S_ang_n[2] , fsi[0].S_ang_n[4] , fsi[2].x_c, fsi[2].y_c, fsi[2].z_c);
                Struc_Solver(&usermg, ibm, fsi, itr_sc,tistart, &DoSCLoop);
            }else
                DoSCLoop = PETSC_FALSE;
      
            /*     Flow Solver! */
            if(levelset) Calc_Inlet_Area(&(usermg.mgctx[usermg.mglevels-1].user[0]));
            Flow_Solver(&usermg, ibm, fsi, itr_sc, wtm, fsi_wt, ibm_ACD, ibm_acl2ref, fsi_acl2ref, ibm_nacelle, fsi_nacelle); //Hossein added itr_sc

            if(rotatefsi || movefsi) for (ibi=0;ibi<NumberOfBodies;ibi++) ibm_surface_out_with_pressure(&ibm[ibi], ibi);

        }// End of while SC loop
        /* ==================================================================================             */

/*      put the time accuracy coefficient back to 1.5 
        after the 1st real-time step */
/*      COEF_TIME_ACCURACY=1.5; */
    

        /* ==================================================================================             */
        /*     Save the old values (at ti) for later */

        level = usermg.mglevels-1;
        user = usermg.mgctx[level].user;
        for (bi=0; bi<block_number; bi++) {

            if (immersed) {
                VecCopy(user[bi].Nvert, user[bi].Nvert_o);
                DAGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
                DAGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
            }

            //VecCopy(user[bi].Ucont_rm1, user[bi].Ucont_rm2);
            if(levelset) VecCopy(user[bi].Levelset, user[bi].Levelset_o);

            VecCopy(user[bi].Ucont_o, user[bi].Ucont_rm1);
            VecCopy(user[bi].Ucont, user[bi].Ucont_o);
            VecCopy(user[bi].P, user[bi].P_o);

            DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
            DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);

            DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
            DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);

            //seokkoo
            DAGlobalToLocalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
            DAGlobalToLocalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, user[bi].lUcat_old);
        }

        if (immersed && (movefsi || rotatefsi || cop || fish || MHV || sediment)){// && ti<tistart+3) {

            for (ibi=0;ibi<NumberOfBodies;ibi++) {      

                for(i=0;i<ibm[ibi].n_elmt; i++){
                    ibm[ibi].cent_z_old[i] = ibm[ibi].cent_z[i];
                    ibm[ibi].cent_y_old[i] = ibm[ibi].cent_y[i];
                }
         
                for (i=0; i<ibm[ibi].n_v; i++) {
                    ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
                    ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
                    ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];

                    ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;
                    ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
                    ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;

                    ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
                    ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
                    ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
                }

                for (i=0;i<6;i++){
                    fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
                    fsi[ibi].S_real[i]=fsi[ibi].S_new[i];

                    fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
                    fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
                }

                fsi[ibi].F_x_real=fsi[ibi].F_x;
                fsi[ibi].F_y_real=fsi[ibi].F_y;
                fsi[ibi].F_z_real=fsi[ibi].F_z;

                fsi[ibi].M_x_rm3=fsi[ibi].M_x_rm2;
                fsi[ibi].M_y_rm3=fsi[ibi].M_y_rm2;
                fsi[ibi].M_z_rm3=fsi[ibi].M_z_rm2;

                fsi[ibi].M_x_rm2=fsi[ibi].M_x_real;
                fsi[ibi].M_y_rm2=fsi[ibi].M_y_real;
                fsi[ibi].M_z_rm2=fsi[ibi].M_z_real;

                fsi[ibi].M_x_real=fsi[ibi].M_x;
                fsi[ibi].M_y_real=fsi[ibi].M_y;
                fsi[ibi].M_z_real=fsi[ibi].M_z;
            }
        }
/* ==================================================================================             */

////////////////////////////////////---------------------------

    } // ti (physical time) loop
/* ==================================================================================             */
  PetscPrintf(PETSC_COMM_WORLD, "\n\n ******* Finished computation ti=%d ******* \n\n", ti);

  MG_Finalize(&usermg);
  PetscFinalize();

/* ==================================================================================             */
    return(0);

}
