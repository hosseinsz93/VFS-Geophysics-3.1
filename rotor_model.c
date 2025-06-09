/*Modeling the rotor blades */

#include "variables.h"
#include "PrintSequential.h"
#include <algorithm>
#include <limits>
#include <math.h>
#include <vector>
using namespace std;

extern	double dfunc_2h(double r);
extern	double dfunc_h(double r);
extern	double dfunc_2hs1(double r);
extern	double dfunc_2hs2(double r);
extern	double dfunc_2hs3(double r);
extern	double dfunc_2hs4(double r);
extern	double dfunc_2hs5(double r);
extern	double dfunc_2hs6(double r);
extern	double dfunc_4hs1(double r);

extern  double dfunc_nhs2(double r, double n); // with plateau, n half width

extern double dfunc_4h2peak(double r);

extern double dfunc_nhs1(double r, double n);

extern	double dfunc_4h(double r);
extern	double dfunc_4htail(double r);
extern  double dfunc_s3h(double r);
extern  double dfunc_s4h(double r);
extern  double dfunc_sc4h(double r);
extern double dfunc_nh(double r, double n);
extern double dfunc_exp(double r, double n);

extern  double dfunc_4uniform(double r);
extern  double dfunc_6uniform(double r);
extern	PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects); 

extern PetscErrorCode Intp_Nut_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects);
extern Cmpnts ArbitraryRotate(Cmpnts p,double theta,Cmpnts r);
//double  coef_cr1 = 2.2, coef_cr2 = 1.0, coef_cr3 = 4.0;

void color_shells(int n_v, int n_elmt, IBMNodes *ibm);

void to_lower(std::string & s)
{
std::transform(s.begin(), s.end(), s.begin(),
    [](unsigned char c){ return std::tolower(c); });
}

int Itpwidth=(int)halfwidth_dfunc+2;

// in the actuator surface model, the surface meshes for the blades must be the first three.
// If fileType is 0 (the default), read a ucd file. If it's any other value, the xpatch file format is used.
PetscErrorCode surface_read(int fileType, IBMNodes *ibm, int ibi, FSInfo *fsi, char fname[80], double reflength)
{

    int	rank;
    int	n_v , n_elmt ;
    PetscReal	*x_bp , *y_bp , *z_bp ;
    int	*nv1 , *nv2 , *nv3 ;
    PetscReal	*nf_x, *nf_y, *nf_z;
    int	i,ii;
    int	n1e, n2e, n3e;
    PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
    PetscReal     dr;
    PetscReal     *dA ;//area
    PetscReal	*nt_x, *nt_y, *nt_z;
    PetscReal	*ns_x, *ns_y, *ns_z;

    char   ss[20];
    char string[128];

    if (fileType)
        PetscPrintf(PETSC_COMM_WORLD, "\n\nReading Actuator Surface Xpatch file\n");
    else
        PetscPrintf(PETSC_COMM_WORLD, "\n\nReading Actuator Surface ucd file\n");

    PetscPrintf(PETSC_COMM_WORLD, "xc=%lf, yc=%lf, zc=%lf for %i th turbine\n", fsi->x_c, fsi->y_c, fsi->z_c, ibi);
    PetscPrintf(PETSC_COMM_WORLD, "nx=%lf, ny=%lf, nz=%lf for %i th turbine\n", fsi->nx_tb, fsi->ny_tb, fsi->nz_tb, ibi);

  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  	if(!rank) { // root processor read in the data
        FILE *fd;
        char filen[160];  

        if (TheSameObject)
            sprintf(filen,"%s/%s%3.3d" , path, fname, 0);
        else
            sprintf(filen,"%s/%s%3.3d" , path, fname, ibi);

        if (fileType)
          PetscPrintf(PETSC_COMM_SELF, "Reading xpatch file name: %s\n", filen);
        else
          PetscPrintf(PETSC_COMM_SELF, "Reading ucd file name: %s\n", filen);
 
        fd = fopen(filen, "r"); 
        if (!fd) printf("Cannot open %s !!", filen),exit(0);
        else printf("Opened %s !\n", filen);
        n_v =0;

        if (fd) {

            fgets(string, 128, fd);
            fgets(string, 128, fd);
            fgets(string, 128, fd);
            if (fileType) {
                fgets(string, 128, fd);
                fscanf(fd, "%i",&n_v);
            }
            else
                fscanf(fd, "%i %i %*i %*i %*i", &n_v, &n_elmt);
            PetscPrintf(PETSC_COMM_SELF, "number of nodes %d\n",n_v);

            ibm->n_v = n_v;

            MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

            x_bp = ibm->x_bp;	// seokkoo
            y_bp = ibm->y_bp;
            z_bp = ibm->z_bp;

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

            for (i=0; i<n_v; i++) {
                if (fileType)
                    fscanf(fd, "%le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);
                else
                    fscanf(fd, "%*d %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);

                x_bp[i]=x_bp[i]/reflength;
                y_bp[i]=y_bp[i]/reflength;
                z_bp[i]=z_bp[i]/reflength;

                ibm->x_bp_i[i] = x_bp[i];
                ibm->y_bp_i[i] = y_bp[i];
                ibm->z_bp_i[i] = z_bp[i];

                x_bp[i] += fsi->x_c;
                y_bp[i] += fsi->y_c;
                z_bp[i] += fsi->z_c;

                ibm->x_bp0[i] = x_bp[i];
                ibm->y_bp0[i] = y_bp[i];
                ibm->z_bp0[i] = z_bp[i];

                ibm->x_bp[i] = x_bp[i];
                ibm->y_bp[i] = y_bp[i];
                ibm->z_bp[i] = z_bp[i];

                ibm->x_bp_o[i] = x_bp[i];
                ibm->y_bp_o[i] = y_bp[i];
                ibm->z_bp_o[i] = z_bp[i];

                ibm->tmprt[i] = 0.;

                ibm->u[i].x = 0.;
                ibm->u[i].y = 0.;
                ibm->u[i].z = 0.;

                ibm->uold[i].x = 0.;
                ibm->uold[i].y = 0.;
                ibm->uold[i].z = 0.;

                ibm->urm1[i].x = 0.;
                ibm->urm1[i].y = 0.;
                ibm->urm1[i].z = 0.;
            }

            i=0;
            PetscPrintf(PETSC_COMM_WORLD, "xyz_bp[%i] = %lf %lf %lf\n", i, x_bp[i], y_bp[i], z_bp[i]);
            PetscPrintf(PETSC_COMM_WORLD, "Rotating (Yawing) of Actuator Surface\n");

//Temporarily turned off for debugging by KFlora
//rotorYaw(ibm, fsi);
            PetscPrintf(PETSC_COMM_WORLD, "Yawed xyz_bp[%i] = %lf %lf %lf\n", i, x_bp[i], y_bp[i], z_bp[i]);

            MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
            MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            char string1[128];
            if (fileType) {
                fscanf(fd, "%*i");
                PetscPrintf(PETSC_COMM_SELF, "number of -- %d\n", ii);
                fscanf(fd, "%*s");
                fscanf(fd, "%i %*i", &n_elmt);
            }
            ibm->n_elmt=n_elmt;

            PetscPrintf(PETSC_COMM_SELF, "number of elemts %d\n",n_elmt);
            MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

            PetscMalloc(n_elmt*sizeof(int), &nv1);
            PetscMalloc(n_elmt*sizeof(int), &nv2);
            PetscMalloc(n_elmt*sizeof(int), &nv3);
  
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
  
            PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
  
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));


            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_x));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_y));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_z));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->disp_theta));
            //Hossein added from turbine structure
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->utheta_struct));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->ux_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uy_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uz_struct));

            // rotor model
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOAAOA_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFt_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFa_mean));


            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));


            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_z));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));


            PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));


            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->color));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->s2l));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_min));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_max));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->random_color));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Nut_lagr));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->dh_IP));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->P_IPlagr));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->rotationdirection));

            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection_mean));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->circulationSpecified));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->frictionfactor));

            for (i=0; i<n_elmt; i++) {
                if (fileType)
                    fscanf(fd, "%i %i %i %*i %i %*i\n", &nv1[i], &nv2[i], &nv3[i], &(ibm->color[i]));
                else {
                    fscanf(fd, "%*d %*d %s %d %d %d\n", &string1, &nv1[i], &nv2[i], &nv3[i]);
                    std::string cell(string1);
                    to_lower(cell);
                    if (cell != "tri") {
                        PetscPrintf(PETSC_COMM_SELF, "ucd reader doesn't support cell type \"%s\".\n", string1);
                        exit(1);
                    }
                }
                nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
            }
            ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

            if (fileType) {
                int color_translate = ibm->color[0]-1;
                for (i=0; i<n_elmt; i++)
                    ibm->color[i]-=color_translate;
            }
            else {
            // Color is handled differently when using ucd files.
            // Each connected region is found and assigned a
            // sequential color number.
                color_shells(n_v, n_elmt, ibm);
            }

            i=0;
            PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

            fclose(fd);
        }

        for (i=0; i<n_elmt; i++) {

            n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
            dx12 = x_bp[n2e] - x_bp[n1e];
            dy12 = y_bp[n2e] - y_bp[n1e];
            dz12 = z_bp[n2e] - z_bp[n1e];

            dx13 = x_bp[n3e] - x_bp[n1e];
            dy13 = y_bp[n3e] - y_bp[n1e];
            dz13 = z_bp[n3e] - z_bp[n1e];

            nf_x[i] = dy12 * dz13 - dz12 * dy13;
            nf_y[i] = -dx12 * dz13 + dz12 * dx13;
            nf_z[i] = dx12 * dy13 - dy12 * dx13;

            dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

            nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;

            if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
            (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
                ns_x[i] = 1.;     
                ns_y[i] = 0.;     
                ns_z[i] = 0. ;

                nt_x[i] = 0.;
                nt_y[i] = 1.;
                nt_z[i] = 0.;
            } else {
                ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
                ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
                ns_z[i] = 0. ;

                nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
                nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
                nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
            }

            dA[i] = dr/2.; 

            // Calc the center of the element
            ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
            ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
            ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
        }

        ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

        ibm->dA = dA;
        ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        int ti=0;

        // Moved surface write after ?? was called since the color
        // values are corrected by that routine.
        // surface_write(ibm, ibi, fname);
    }
    else if (rank) {
        MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
        ibm->n_v = n_v;

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

        x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
        y_bp = ibm->y_bp;
        z_bp = ibm->z_bp;

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

        for (i=0; i<n_v; i++) {
            ibm->tmprt[i] = 0.;

            ibm->u[i].x = 0.;
            ibm->u[i].y = 0.;
            ibm->u[i].z = 0.;

            ibm->uold[i].x = 0.;
            ibm->uold[i].y = 0.;
            ibm->uold[i].z = 0.;
  
            ibm->urm1[i].x = 0.;
            ibm->urm1[i].y = 0.;
            ibm->urm1[i].z = 0.;      
        }

        MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
        MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
        MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
        n_elmt = ibm->n_elmt;

        PetscMalloc(n_elmt*sizeof(int), &nv1);
        PetscMalloc(n_elmt*sizeof(int), &nv2);
        PetscMalloc(n_elmt*sizeof(int), &nv3);

        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

        PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

        ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
        ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

        ibm->dA = dA;
        ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_x));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_y));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_z));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->disp_theta));
        //Hossein added from turbine structure
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->utheta_struct));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->ux_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uy_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uz_struct));

        // rotor model

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOAAOA_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFt_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFa_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->color));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->s2l));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_min));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_max));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->random_color));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Nut_lagr));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->dh_IP));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->P_IPlagr));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->rotationdirection));

        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection_mean));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->circulationSpecified));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->frictionfactor));

        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    }
	PetscPrintf(PETSC_COMM_WORLD, "Finished Reading %s!\n\n", fname);

    for (i=0; i<ibm->n_elmt; i++) {
        ibm->Urelmag_mean[i]=0.0;
        ibm->Fr_mean[i]=0.0;
        ibm->Ft_mean[i]=0.0;
        ibm->Fa_mean[i]=0.0;
        ibm->Ur_mean[i]=0.0;
        ibm->Ut_mean[i]=0.0;
        ibm->Ua_mean[i]=0.0;
        ibm->AOA_mean[i]=0.0;

        ibm->FFa_mean[i]=0.0;
        ibm->FFt_mean[i]=0.0;
        ibm->AOAAOA_mean[i]=0.0;

        ibm->Urel_mean[i].x=0.0;
        ibm->Urel_mean[i].y=0.0;
        ibm->Urel_mean[i].z=0.0;

        ibm->Uinduced_mean[i].x=0.0;
        ibm->Uinduced_mean[i].y=0.0;
        ibm->Uinduced_mean[i].z=0.0;

        ibm->circulation_mean[i].x=0.0;
        ibm->circulation_mean[i].y=0.0;
        ibm->circulation_mean[i].z=0.0;

        ibm->liftdirection_mean[i].x=0.0;
        ibm->liftdirection_mean[i].y=0.0;
        ibm->liftdirection_mean[i].z=0.0;

        ibm->disp_theta[i]=0;
    }

    for (i=0; i<ibm->n_v; i++) {
        ibm->disp_x[i]=0;
        ibm->disp_y[i]=0;
        ibm->disp_z[i]=0;
    }

    count_AL=0.0;
    return(0);
}


// Write the actuator surface files.
PetscErrorCode surface_write(IBMNodes *ibm, int ibi, char fname[80])
{
    int i;
    FILE *f;
    char filen[256];
    sprintf(filen, "%s/%s_%2.2d_nf.dat",path,fname,ibi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z, color\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"TRIANGLES\", N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-13]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);
    for (i=0; i<ibm->n_v; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->x_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_v; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->y_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_v; i++) {	
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->z_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->nf_x[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->nf_y[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->nf_z[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->nt_x[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->nt_y[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->nt_z[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->ns_x[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->ns_y[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e ", ibm->ns_z[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%d ", ibm->color[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    for (i=0; i<ibm->n_elmt; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
    }

    fclose(f);
}


// Find each connected region and assign a unique color to each region.
// WRO 2022-08-23

struct vertsToShells
{
    vertsToShells() : shells(), used(false)
    {
        shells.reserve(8);
    }

    std::vector<int> shells;
    bool used;
};


void csetFromVset(IBMNodes *ibm, std::vector< vertsToShells > & verts,
    std::vector < bool > & cells, std::vector < int > & vset,
    std::vector < int > & cset, int currColor);
void vsetFromCset(IBMNodes *ibm, std::vector< vertsToShells > & verts,
    std::vector < bool > & cells, std::vector < int > & vset,
    std::vector < int > & cset);


void color_shells(int n_v, int n_elmt, IBMNodes *ibm)
{
    // Build a vector of vertices and the shells they are attached to.
    int * nv1 = ibm->nv1;
    int * nv2 = ibm->nv2;
    int * nv3 = ibm->nv3;
    std::vector< vertsToShells > verts(n_v);
    for (int ie=0; ie<n_elmt; ++ie) {
        verts[nv1[ie]].shells.push_back(ie);
        verts[nv2[ie]].shells.push_back(ie);
        verts[nv3[ie]].shells.push_back(ie);
    }

  // vector to store information about whither a shell has been used or not.
    std::vector < bool > cells(n_elmt, false);

    // vectors to store current vertex and cell sets.
    std::vector < int > vset;
    vset.reserve(n_v);
    std::vector < int > cset;
    cset.reserve(n_elmt);

    int currColor = -1;

    // Start looking for another separate shell region.
    while (true)
    {
        vset.resize(0);
        for (auto & vr : verts)
        if (!vr.used)
        {
            vr.used = true;
            vset.push_back(&vr - &*verts.begin());
            break;
        }

        // Are we done?
        if (vset.empty())
            return;

        ++currColor;

        while (true)
        {
            csetFromVset(ibm, verts, cells, vset, cset, currColor);

            // Are we done with this region?
            if (cset.empty())
            break;

            vsetFromCset(ibm, verts, cells, vset, cset);

            // Are we done with this region?
            if (vset.empty())
                break;
        }
    }
}


void csetFromVset(IBMNodes *ibm, std::vector< vertsToShells > & verts,
    std::vector < bool > & cells, std::vector < int > & vset,
    std::vector < int > & cset, int currColor)
{
    cset.resize(0);

    for (auto & ivs : vset) {
        auto & v = verts[ivs];
        for (auto & ic : v.shells)
        {
            if (!cells[ic])
            {
                cells[ic] = true;
                ibm->color[ic] = currColor;
                cset.push_back(ic);
            }
        }
    }
}


void vsetFromCset(IBMNodes *ibm, std::vector< vertsToShells > & verts,
    std::vector < bool > & cells, std::vector < int > & vset,
    std::vector < int > & cset)
{
    vset.resize(0);

    for (auto & ics : cset)
    {
        if (!verts[ibm->nv1[ics]].used)
        {
            verts[ibm->nv1[ics]].used = true;
            vset.push_back(ibm->nv1[ics]);
        }

        if (!verts[ibm->nv2[ics]].used)
        {
            verts[ibm->nv2[ics]].used = true;
            vset.push_back(ibm->nv2[ics]);
        }

        if (!verts[ibm->nv3[ics]].used)
        {
            verts[ibm->nv3[ics]].used = true;
            vset.push_back(ibm->nv3[ics]);
        }
    }
}


// The colors between the surface and line files probably don't match.
// They don't in the case I'm trying now. This routine matches the
// surface color is the line color by using geometry matching.
void matchLineColorToSurface(int ibi, IBMNodes * srfibm, IBMNodes * linibm)
{
    int colorIdx[linibm->num_blade];
    int nv_blade = linibm->n_v / linibm->num_blade;
    int n_elmt_blade = linibm->n_elmt / linibm->num_blade;

    for (int lnb=0; lnb<linibm->num_blade; ++lnb)
    {
        // Find the last vertex for this blade.
        int li = (lnb+1) * nv_blade - 1;
        Cmpnts lc(linibm->x_bp_i[li], linibm->y_bp_i[li], linibm->z_bp_i[li]);

        // Find the closest vertex to above coordinates in the surface file
        int simin;
        double siDistMin = std::numeric_limits<double>::max();
        for (int si=0; si<srfibm->n_v; ++si)
        {
            Cmpnts sc(srfibm->x_bp_i[si], srfibm->y_bp_i[si], srfibm->z_bp_i[si]);
            double dist = lc.dist(sc);
            if (dist < siDistMin)
            {
                simin = si;
                siDistMin = dist;
            }
        }

        // Find the corresponding surface color by search the surface shell file
        // that contains the above found vertex. That shell points to the color.
        int srfColor, linColor;
        for (int si=0; si<srfibm->n_elmt; ++si)
        {
            if (simin == srfibm->nv1[si] ||
                simin == srfibm->nv2[si] ||
                simin == srfibm->nv3[si])
            {
                srfColor = srfibm->color[si];
                break;
            }
        }

        linColor = linibm->color[(lnb+1)*n_elmt_blade - 1];
        colorIdx[srfColor] = linColor;
    }

    // Update the surface colors
    for (int si=0; si<srfibm->n_elmt; ++si)
        srfibm->color[si] = colorIdx[srfibm->color[si]];
}


PetscErrorCode disk_read_ucd(IBMNodes *ibm, int ibi, FSInfo *fsi, int OneDmZ, char fname[80], double reflength)
{
    int	rank;
    int	n_v , n_elmt ;
    PetscReal	*x_bp , *y_bp , *z_bp ;
    int	*nv1 , *nv2 , *nv3 ;
    PetscReal	*nf_x, *nf_y, *nf_z;
    int	i,ii;
    int	n1e, n2e, n3e;
    PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
    PetscReal     dr;
    PetscReal     *dA ;//area
    PetscReal	*nt_x, *nt_y, *nt_z;
    PetscReal	*ns_x, *ns_y, *ns_z;

    char   ss[20];
    char string[128];

        //PetscPrintf(PETSC_COMM_WORLD, "xc=%le, yc=%le, zc=%le for %i th turbine\n", fsi->x_c, fsi->y_c, fsi->z_c, ibi);
        //PetscPrintf(PETSC_COMM_WORLD, "nx=%le, ny=%le, nz=%le for %i th turbine\n", fsi->nx_tb, fsi->ny_tb, fsi->nz_tb, ibi);
    PetscPrintf(PETSC_COMM_WORLD, "Disk Read UCD: %s with reference length: %f\n", fname,reflength);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(!rank) { // root processor read in the data
        FILE *fd;
        char filen[160];  
        sprintf(filen,"%s/%s" , path, fname);
//      PetscPrintf(PETSC_COMM_SELF, "READ %s\n", filen);

        fd = fopen(filen, "r"); 
        if (!fd) printf("Cannot open %s !!", filen),exit(0);
        else printf("Opened %s !\n", filen);
        n_v =0;

        if (fd) {
            fgets(string, 128, fd);
            fgets(string, 128, fd);
            fgets(string, 128, fd);

            fscanf(fd, "%i %i %*i %*i %*i", &n_v, &n_elmt);
            PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d\n",n_v, n_elmt);
  
            ibm->n_v = n_v;
            ibm->n_elmt = n_elmt;      
      
            MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
            MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
  
            x_bp = ibm->x_bp;	// seokkoo
            y_bp = ibm->y_bp;
            z_bp = ibm->z_bp;

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
  
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

            PetscMalloc(n_v*sizeof(Node2Cell), &ibm->n2c);
            PetscMalloc(n_elmt*sizeof(Cell2Cell), &ibm->c2c);
            PetscMalloc(n_elmt*sizeof(Cell2Cells), &ibm->c2cs);

            PetscMalloc(n_v*sizeof(int), &ibm->node_boundary); // xiaolei add
            PetscMalloc(n_elmt*sizeof(int), &ibm->elmt_boundary); // xiaolei add

            PetscMalloc(3*n_elmt*sizeof(PetscReal), &(ibm->error_force)); // xiaolei add

            // fsi_fem
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u_pre));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_pre));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_pre));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_pre));

            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->af));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->af_old));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->af_pre));

            //  Error Vector
            PetscMalloc(9*n_v*sizeof(PetscReal), &(ibm->error));

            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_z);

            // xiaolei add
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_x_pre);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_y_pre);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_z_pre);

            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumP_Tx);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumP_Ty);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumP_Tz);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->AbsP_tau);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumM_Tx);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumM_Ty);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumM_Tz);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->AbsM_tau);

            for (i=0; i<n_v; i++) {
                fscanf(fd, "%*i %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);

                x_bp[i]=x_bp[i]/reflength;//*13.5/48.;
                y_bp[i]=y_bp[i]/reflength;//*13.5/48.;
                z_bp[i]=z_bp[i]/reflength;//*13.5/48.;

                ibm->x_bp_i[i] = x_bp[i];
                ibm->y_bp_i[i] = y_bp[i];
                ibm->z_bp_i[i] = z_bp[i];

                x_bp[i] += fsi->x_c;
                y_bp[i] += fsi->y_c;
                z_bp[i] += fsi->z_c;

                ibm->x_bp0[i] = x_bp[i];
                ibm->y_bp0[i] = y_bp[i];
                ibm->z_bp0[i] = z_bp[i];

                if (OneDmZ) {
                    double R=fsi->r_rotor/reflength_wt;
                    double rr = loc_refvel*2.0*R;
                    x_bp[i]-=rr*fsi->nx_tb; 
                    y_bp[i]-=rr*fsi->ny_tb; 
                    z_bp[i]-=rr*fsi->nz_tb; 
                }

                ibm->x_bp[i] = x_bp[i];
                ibm->y_bp[i] = y_bp[i];
                ibm->z_bp[i] = z_bp[i];

                ibm->x_bp_o[i] = x_bp[i];
                ibm->y_bp_o[i] = y_bp[i];
                ibm->z_bp_o[i] = z_bp[i];

                ibm->tmprt[i] = 0.;

                ibm->u[i].x = 0.;
                ibm->u[i].y = 0.;
                ibm->u[i].z = 0.;

                ibm->uold[i].x = 0.;
                ibm->uold[i].y = 0.;
                ibm->uold[i].z = 0.;

                ibm->urm1[i].x = 0.;
                ibm->urm1[i].y = 0.;
                ibm->urm1[i].z = 0.;
            }
            i=0;
            PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %lf %lf %lf\n", x_bp[i], y_bp[i], z_bp[i]);

            MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
            MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);


            PetscMalloc(n_elmt*sizeof(int), &nv1);
            PetscMalloc(n_elmt*sizeof(int), &nv2);
            PetscMalloc(n_elmt*sizeof(int), &nv3);
  
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
  
            PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
  
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_x));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_y));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_z));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->disp_theta));
            //Hossein added from turbine structure
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->utheta_struct));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->ux_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uy_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uz_struct));

        // rotor model
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOAAOA_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFt_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFa_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_z));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));


            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_min));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_max));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->random_color));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Nut_lagr));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->dh_IP));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->P_IPlagr));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->rotationdirection));

            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection_mean));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->circulationSpecified));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->frictionfactor));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->color));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->s2l));

            for (i=0; i<n_elmt; i++) {
                fscanf(fd, "%*i %*i %s %i %i %i\n", &ss, nv1+i, nv2+i, nv3+i);
                nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
            }
            ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

            i=0;
            PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

            fclose(fd);
        }

        for (i=0; i<n_elmt; i++) {
      
            n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
            dx12 = x_bp[n2e] - x_bp[n1e];
            dy12 = y_bp[n2e] - y_bp[n1e];
            dz12 = z_bp[n2e] - z_bp[n1e];
  
            dx13 = x_bp[n3e] - x_bp[n1e];
            dy13 = y_bp[n3e] - y_bp[n1e];
            dz13 = z_bp[n3e] - z_bp[n1e];
  
            nf_x[i] = dy12 * dz13 - dz12 * dy13;
            nf_y[i] = -dx12 * dz13 + dz12 * dx13;
            nf_z[i] = dx12 * dy13 - dy12 * dx13;
  
            dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
  
            nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
            if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
                (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
                ns_x[i] = 1.;     
                ns_y[i] = 0.;     
                ns_z[i] = 0. ;

                nt_x[i] = 0.;
                nt_y[i] = 1.;
                nt_z[i] = 0.;
            } else {
                ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
                ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
                ns_z[i] = 0. ;

                nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
                nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
                nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
            }

            dA[i] = dr/2.; 
      
            // Calc the center of the element
            ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
            ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
            ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
        }

        ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

        ibm->dA = dA;
        ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        int ti=0;
        FILE *f;
        sprintf(filen, "%s/%s_%2.2d_nf.dat",path,fname,ibi);
        f = fopen(filen, "w");
        PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
        PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"TRIANGLES\", N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
        for (i=0; i<n_v; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
        }
        for (i=0; i<n_v; i++) {
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
        }
        for (i=0; i<n_v; i++) {	
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]);
        }
        for (i=0; i<n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
        }

        fclose(f);

    }
  	else if (rank) {
        MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
        ibm->n_v = n_v;

        MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
        n_elmt = ibm->n_elmt;

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

        x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
        y_bp = ibm->y_bp;
        z_bp = ibm->z_bp;

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));


        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

        PetscMalloc(n_v*sizeof(Node2Cell), &ibm->n2c);
        PetscMalloc(n_elmt*sizeof(Cell2Cell), &ibm->c2c);
        PetscMalloc(n_elmt*sizeof(Cell2Cells), &ibm->c2cs);

        PetscMalloc(n_v*sizeof(int), &ibm->node_boundary); // xiaolei add
        PetscMalloc(n_elmt*sizeof(int), &ibm->elmt_boundary); // xiaolei add

        PetscMalloc(3*n_elmt*sizeof(PetscReal), &(ibm->error_force)); // xiaolei add

    // fsi_fem
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u_pre));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_pre));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_pre));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_pre));

        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->af));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->af_old));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->af_pre));

        // Error Vector
        PetscMalloc(9*n_v*sizeof(PetscReal), &(ibm->error));

        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_z);

        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_x_pre);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_y_pre);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->Traction_z_pre);

		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumP_Tx);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumP_Ty);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumP_Tz);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->AbsP_tau);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumM_Tx);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumM_Ty);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->SumM_Tz);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->AbsM_tau);

        for (i=0; i<n_v; i++) {
            ibm->tmprt[i] = 0.;
            ibm->u[i].x = 0.;
            ibm->u[i].y = 0.;
            ibm->u[i].z = 0.;

            ibm->uold[i].x = 0.;
            ibm->uold[i].y = 0.;
            ibm->uold[i].z = 0.;
  
            ibm->urm1[i].x = 0.;
            ibm->urm1[i].y = 0.;
            ibm->urm1[i].z = 0.;      
        }

        MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
        MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
        MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        PetscMalloc(n_elmt*sizeof(int), &nv1);
        PetscMalloc(n_elmt*sizeof(int), &nv2);
        PetscMalloc(n_elmt*sizeof(int), &nv3);

        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

        PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

        ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
        ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

        ibm->dA = dA;
        ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_x));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_y));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_z));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->disp_theta));
        //Hossein added from turbine structure
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->utheta_struct));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->ux_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uy_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uz_struct));

        // rotor model

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOAAOA_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFt_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFa_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_min));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_max));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->random_color));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Nut_lagr));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->dh_IP));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->P_IPlagr));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->rotationdirection));

        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection_mean));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->circulationSpecified));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->frictionfactor));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->color));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->s2l));

        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    }
    PetscPrintf(PETSC_COMM_WORLD, "Finished Reading %s !\n\n", fname);

    for (i=0; i<ibm->n_elmt; i++) {
        ibm->Urelmag_mean[i]=0.0;
        ibm->Fr_mean[i]=0.0;
        ibm->Ft_mean[i]=0.0;
        ibm->Fa_mean[i]=0.0;
        ibm->Ur_mean[i]=0.0;
        ibm->Ut_mean[i]=0.0;
        ibm->Ua_mean[i]=0.0;
        ibm->AOA_mean[i]=0.0;
        ibm->color[i]=0;

        ibm->FFa_mean[i]=0.0;
        ibm->FFt_mean[i]=0.0;
        ibm->AOAAOA_mean[i]=0.0;

        ibm->Urel_mean[i].x=0.0;
        ibm->Urel_mean[i].y=0.0;
        ibm->Urel_mean[i].z=0.0;

        ibm->Uinduced_mean[i].x=0.0;
        ibm->Uinduced_mean[i].y=0.0;
        ibm->Uinduced_mean[i].z=0.0;

        ibm->circulation_mean[i].x=0.0;
        ibm->circulation_mean[i].y=0.0;
        ibm->circulation_mean[i].z=0.0;

        ibm->liftdirection_mean[i].x=0.0;
        ibm->liftdirection_mean[i].y=0.0;
        ibm->liftdirection_mean[i].z=0.0;

        ibm->disp_theta[i]=0;
    }

    for (i=0; i<ibm->n_v; i++) {
        ibm->disp_x[i]=0;
        ibm->disp_y[i]=0;
        ibm->disp_z[i]=0;
    }

    for (i=0; i<ibm->n_elmt; i++) {
        ibm->Traction_x[i] = 0;
        ibm->Traction_y[i] = 0;
        ibm->Traction_z[i] = 0;

        ibm->Traction_x_pre[i] = 0;
        ibm->Traction_y_pre[i] = 0;
        ibm->Traction_z_pre[i] = 0;
    }

    count_AL=0.0;

    return(0);
}


PetscErrorCode ACL_read(IBMNodes *ibm, int ibi, FSInfo *fsi, double reflength)
{
    int	rank;
    int	n_v , n_elmt ;
    PetscReal	*x_bp , *y_bp , *z_bp ;
    int	*nv1 , *nv2 , *nv3 ;
    PetscReal	*nf_x, *nf_y, *nf_z;
    int	i,ii;
    int	n1e, n2e, n3e;
    PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
    PetscReal     dr;
    //Added 4/1/06 iman
    PetscReal     *dA ;//area
    PetscReal	*nt_x, *nt_y, *nt_z;
    PetscReal	*ns_x, *ns_y, *ns_z;

    int nv_blade, nelmt_blade;

    char   ss[20];
    //double xt;
    char string[128];

    PetscPrintf(PETSC_COMM_WORLD, "ACL Read \n");

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(!rank) { // root processor read in the data
        FILE *fd;
        PetscPrintf(PETSC_COMM_SELF, "READ acldata\n");
        char filen[80];  

        if (TheSameObject)
            sprintf(filen,"%s/acldata%3.3d" , path, 0);
        else 
            sprintf(filen,"%s/acldata%3.3d" , path, ibi);

        fd = fopen(filen, "r"); 
        if (!fd) 
            printf("Cannot open %s !!", filen),exit(0);
        else 
            printf("Opened %s !\n", filen);

        n_v =0;

        if (fd) {

            fscanf(fd, "%i ",&n_v);
            n_elmt = n_v - 1;
            PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements & blades %d %d %d \n",n_v, n_elmt, ibm->num_blade);
            nv_blade=n_v;
            nelmt_blade=n_elmt; 

            ibm->n_v = n_v * ibm->num_blade;
            ibm->n_elmt = n_elmt * ibm->num_blade;      

            n_v=ibm->n_v;
            n_elmt=ibm->n_elmt;      

            MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

            x_bp = ibm->x_bp; // seokkoo
            y_bp = ibm->y_bp;
            z_bp = ibm->z_bp;

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));
  
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
  
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
            PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

            int nb,j;
            for (nb=0; nb<ibm->num_blade; nb++) {
                if (nb != 0) fscanf(fd, "%*i ");
                for (j=0; j<nv_blade; j++) {
                    i = nb*nv_blade + j;

                    fscanf(fd, "%le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
                
                    x_bp[i]=x_bp[i]/reflength;//*13.5/48.;
                    y_bp[i]=y_bp[i]/reflength;//*13.5/48.;
                    z_bp[i]=z_bp[i]/reflength;//*13.5/48.;

                    ibm->x_bp_i[i] = x_bp[i];
                    ibm->y_bp_i[i] = y_bp[i];
                    ibm->z_bp_i[i] = z_bp[i];

                    x_bp[i] += fsi->x_c;
                    y_bp[i] += fsi->y_c;
                    z_bp[i] += fsi->z_c;

                    ibm->x_bp0[i] = x_bp[i];
                    ibm->y_bp0[i] = y_bp[i];
                    ibm->z_bp0[i] = z_bp[i];


                    ibm->x_bp[i] = x_bp[i];
                    ibm->y_bp[i] = y_bp[i];
                    ibm->z_bp[i] = z_bp[i];

                    ibm->x_bp_o[i] = x_bp[i];
                    ibm->y_bp_o[i] = y_bp[i];
                    ibm->z_bp_o[i] = z_bp[i];


                    ibm->tmprt[i] = 0.;

                    ibm->u[i].x = 0.;
                    ibm->u[i].y = 0.;
                    ibm->u[i].z = 0.;

                    ibm->uold[i].x = 0.;
                    ibm->uold[i].y = 0.;
                    ibm->uold[i].z = 0.;

                    ibm->urm1[i].x = 0.;
                    ibm->urm1[i].y = 0.;
                    ibm->urm1[i].z = 0.;
                }
            }

            i=0;
            PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %lf %lf %lf\n", x_bp[i], y_bp[i], z_bp[i]);

            MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
            MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
            MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

            MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

            PetscMalloc(n_elmt*sizeof(int), &nv1);
            PetscMalloc(n_elmt*sizeof(int), &nv2);
            PetscMalloc(n_elmt*sizeof(int), &nv3);
      
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
  
            // Added 4/1/06 iman
            PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_x));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_y));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_z));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->disp_theta));
            //Hossein added from turbine structure
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->utheta_struct));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->ux_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uy_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uz_struct));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));


            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOAAOA_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFt_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFa_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_z));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
            PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);

            PetscMalloc(n_elmt*sizeof(int), &ibm->color);

            PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_min));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_min));

            PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_max));
            PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_max));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->random_color));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Nut_lagr));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->dh_IP));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->P_IPlagr));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_x));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_y));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_z));

            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->rotationdirection));

            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation_mean));
            PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection_mean));

            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->circulationSpecified));
            PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->frictionfactor));

            // colored as blade index 
            for (nb=0; nb<ibm->num_blade; nb++) {
                for (j=0; j<nelmt_blade; j++) {
                    i = nb*nelmt_blade + j;
                    ibm->color[i] = nb+1;
                }
            }

            for (nb=0; nb<ibm->num_blade; nb++) {
                for (j=0; j<nelmt_blade; j++) {
                    i = nb*nelmt_blade + j;
                    ii = nb*nv_blade + j;
                    nv1[i] = ii; nv2[i] = ii + 1; nv3[i] = ii + 1;
                }
            }

            ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

            i=0;
            PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);

            fclose(fd);
        }
 
        int nb, j;
        // for ACL, the nf_x, nf_y, nf_z denote the direction of the actuator line. and dA denotes the length of each element.
        for (nb=0; nb<ibm->num_blade; nb++) {
  
            for (j=0; j<nelmt_blade; j++) {
                i = nb*nelmt_blade + j;

                n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; 
                dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
                dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
                dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];
  
                nf_x[i] = dx12;
                nf_y[i] = dy12;
                nf_z[i] = dz12;
  
                dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
  
                nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
 
                ns_x[i] = 0.;     
                ns_y[i] = 0.;     
                ns_z[i] = 0. ;

                nt_x[i] = 0.;
                nt_y[i] = 0.;
                nt_z[i] = 0.;
 
                dA[i] = dr;

                // Calc the center of the element
                ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e])/2.;
                ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e])/2.;
                ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e])/2.;

            }
      
        }
     
        ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

        //Added 4/1/06 iman
        ibm->dA = dA;
        ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
        // Added 4/1/06 iman
        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        // Added 4/2/06 iman
        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        // Added 6/4/06 iman
        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        int ti=0;
        FILE *f;
        sprintf(filen, "%s/line%3.3d_%2.2d_nf.dat",path,ti,ibi);
        f = fopen(filen, "w");
        for (i=0; i<ibm->n_elmt; i++) {
            PetscFPrintf(PETSC_COMM_WORLD, f, "%e %le %le \n", ibm->cent_x[i], ibm->cent_y[i], ibm->cent_z[i]);
        }
        fclose(f);
        // kevin debugging??
        PetscPrintf(PETSC_COMM_WORLD, "Printed Line File %s\n", filen);
    }
    else if (rank) {
        MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
        ibm->n_v = n_v;

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

        x_bp = ibm->x_bp;	// added by seokkoo 03.04.2009
        y_bp = ibm->y_bp;
        z_bp = ibm->z_bp;

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_i));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_i));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_i));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
        PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->tmprt));

        for (i=0; i<n_v; i++) {
            ibm->tmprt[i] = 0.;
            ibm->u[i].x = 0.;
            ibm->u[i].y = 0.;
            ibm->u[i].z = 0.;

            ibm->uold[i].x = 0.;
            ibm->uold[i].y = 0.;
            ibm->uold[i].z = 0.;

            ibm->urm1[i].x = 0.;
            ibm->urm1[i].y = 0.;
            ibm->urm1[i].z = 0.;      
        }

        MPI_Bcast(ibm->x_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp_i, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
        MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
        ibm->n_elmt = n_elmt;

        PetscMalloc(n_elmt*sizeof(int), &nv1);
        PetscMalloc(n_elmt*sizeof(int), &nv2);
        PetscMalloc(n_elmt*sizeof(int), &nv3);

        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

        //Added 4/1/06 iman
        PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

        //Added 4/2/06 iman
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

        ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
        ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

        // Added 4/2/06 iman
        ibm->dA = dA;
        ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
        ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

        // Added 6/4/06
        //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_x));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_y));
        PetscMalloc(n_v*sizeof(PetscReal), &(ibm->disp_z));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->disp_theta));
        //Hossein added from turbine structure
            PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->utheta_struct));
            PetscMalloc(n_v*sizeof(PetscReal), &(ibm->ux_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uy_struct));
      		PetscMalloc(n_v*sizeof(PetscReal), &(ibm->uz_struct));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fr_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ft_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Fa_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOAAOA_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFt_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->FFa_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ur_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ut_mean));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ua_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->AOA_mean));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Urelmag));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dF_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_x));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_y));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->U_lagr_z));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Ftmprt_lagr));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Tmprt_lagr));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_min));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->i_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->j_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->k_max));

        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhx));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhy));
        PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->dhz));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_attack));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->angle_twist));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->chord_blade));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CD));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->CL));

        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
        PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);

        PetscMalloc(n_elmt*sizeof(int), &ibm->color);

        PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_min));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_min));

        PetscMalloc(n_elmt*sizeof(int), &(ibm->iIP_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->jIP_max));
        PetscMalloc(n_elmt*sizeof(int), &(ibm->kIP_max));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->random_color));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->centIP_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->U_IPlagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Nut_lagr));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->dh_IP));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->P_IPlagr));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->Shear_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->ShearDesired_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_x));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_y));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->UU_lagr_z));

        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->rotationdirection));

        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Urel_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->Uinduced_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->circulation_mean));
        PetscMalloc(ibm->n_elmt*sizeof(Cmpnts), &(ibm->liftdirection_mean));

        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->circulationSpecified));
        PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &(ibm->frictionfactor));

        MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        //Added 4/2/06 iman
        MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

        MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
        MPI_Bcast(ibm->color, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    }

    for (i=0; i<ibm->n_elmt; i++) {
        ibm->Urelmag_mean[i]=0.0;
        ibm->Fr_mean[i]=0.0;
        ibm->Ft_mean[i]=0.0;
        ibm->Fa_mean[i]=0.0;
        ibm->Ur_mean[i]=0.0;
        ibm->Ut_mean[i]=0.0;
        ibm->Ua_mean[i]=0.0;
        ibm->AOA_mean[i]=0.0;

        ibm->FFa_mean[i]=0.0;
        ibm->FFt_mean[i]=0.0;
        ibm->AOAAOA_mean[i]=0.0;

        ibm->Urel_mean[i].x=0.0;
        ibm->Urel_mean[i].y=0.0;
        ibm->Urel_mean[i].z=0.0;

        ibm->Uinduced_mean[i].x=0.0;
        ibm->Uinduced_mean[i].y=0.0;
        ibm->Uinduced_mean[i].z=0.0;

        ibm->circulation_mean[i].x=0.0;
        ibm->circulation_mean[i].y=0.0;
        ibm->circulation_mean[i].z=0.0;

        ibm->liftdirection_mean[i].x=0.0;
        ibm->liftdirection_mean[i].y=0.0;
        ibm->liftdirection_mean[i].z=0.0;

        ibm->disp_theta[i]=0;
    }

    for (i=0; i<ibm->n_v; i++) {
        ibm->disp_x[i]=0;
        ibm->disp_y[i]=0;
        ibm->disp_z[i]=0;
    }

    count_AL=0.0;
    PetscPrintf(PETSC_COMM_WORLD, "Finished Reading acldata file !\n\n");
    return(0);
}


/***
 * indentify the nearest grid points for ib_delta method. if it is not fsi, the this function will not be called 
 */
// Call in main.c, after reading the grid for turbines
PetscErrorCode Pre_process(UserCtx *user, IBMNodes *ibm, int NumberOfObjects)
{

    DA              da = user->da, fda = user->fda;
    DALocalInfo  	info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    Cmpnts		***coor, ***csi, ***eta, ***zet;

    PetscReal 	***aj;

    Vec		Coor;

    DAGetLocalInfo(da, &info);
    mx = info.mx; my = info.my; mz = info.mz;
    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;

    //PetscPrintf(PETSC_COMM_WORLD, "Debug mx = %i, my= %i, mz= %i, lxs = %i, lxe= %i, lys = %i, lye = %i,  lzs = %i. lze= %i\n", mx, my, mz, lxs, lxe, lys, lye, lzs, lze);

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);

    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);

    int n1e, n2e, n3e;

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        std::vector<double> dhx_ (ibm[ibi].n_elmt);
        std::vector<double> dhy_ (ibm[ibi].n_elmt);
        std::vector<double> dhz_ (ibm[ibi].n_elmt);
        std::vector<double> sum_dhx_ (ibm[ibi].n_elmt);
        std::vector<double> sum_dhy_ (ibm[ibi].n_elmt);
        std::vector<double> sum_dhz_ (ibm[ibi].n_elmt);

        std::vector<int> iclose (ibm[ibi].n_elmt);
        std::vector<int> jclose (ibm[ibi].n_elmt);
        std::vector<int> kclose (ibm[ibi].n_elmt);
        std::vector<int> sum_iclose (ibm[ibi].n_elmt);
        std::vector<int> sum_jclose (ibm[ibi].n_elmt);
        std::vector<int> sum_kclose (ibm[ibi].n_elmt);

        std::vector<int> count (ibm[ibi].n_elmt);
        std::vector<int> sum_count (ibm[ibi].n_elmt);

        std::vector<double> dminl2e (ibm[ibi].n_elmt);
        std::vector<double> dminl2e_glb (ibm[ibi].n_elmt);

        std::vector<int> ic (ibm[ibi].n_elmt);
        std::vector<int> jc (ibm[ibi].n_elmt);
        std::vector<int> kc (ibm[ibi].n_elmt);
        
        //Initialize vectors to zero
        for (l=0; l<ibm[ibi].n_elmt; l++) {

            dhx_[l] = 0.0;
            dhy_[l] = 0.0;
            dhz_[l] = 0.0;

            sum_dhx_[l] = 0.0;
            sum_dhy_[l] = 0.0;
            sum_dhz_[l] = 0.0;

            iclose[l]=0;
            jclose[l]=0;
            kclose[l]=0;

            sum_iclose[l]=0;
            sum_jclose[l]=0;
            sum_kclose[l]=0;

            count[l]=0;
            sum_count[l]=0;
        }

        //Find closest grid node for each IB element
        for (l=0; l<ibm[ibi].n_elmt; l++) {
            int imark=0;
            double dmin=1.e6;
            int icc,jcc,kcc; //interim i,j,k coordinates of closest grid node

            for (i=lxs; i<lxe; i=i+20)
            for (j=lys; j<lye; j=j+20) 
            for (k=lzs; k<lze; k=k+20){ 

                double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
                double d1=sqrt(r1*r1+r2*r2+r3*r3);
                if (d1<dmin) {
                    dmin=d1;
                    icc=i; jcc=j; kcc=k;
                }
            }

            dminl2e[l] = 1.e13;
            for (i=max(icc-21,lxs); i<min(icc+21,lxe); i++)
            for (j=max(jcc-21,lys); j<min(jcc+21,lye); j++) 
            for (k=max(kcc-21,lzs); k<min(kcc+21,lze); k++){ 
                double r1=ibm[ibi].cent_x[l]-coor[k][j][i].x, r2=ibm[ibi].cent_y[l]-coor[k][j][i].y, r3=ibm[ibi].cent_z[l]-coor[k][j][i].z; 
                double d1=sqrt(r1*r1+r2*r2+r3*r3);
                if (d1<dminl2e[l]) {
                    dminl2e[l]=d1;
                    ic[l]=i; jc[l]=j; kc[l]=k;
                }
            }

            //PetscPrintf(PETSC_COMM_WORLD, "Debug - for l=%i, coordinates xyz = %f,%f,%f\n", l, coor[kc[l]][jc[l]][ic[l]].x,coor[kc[l]][jc[l]][ic[l]].y,coor[kc[l]][jc[l]][ic[l]].z);
            //PetscPrintf(PETSC_COMM_WORLD, "ibm centroid xyz = %f,%f,%f and dmin1 = %f,  dmin2=%f\n", ibm[ibi].cent_x[l], ibm[ibi].cent_y[l],ibm[ibi].cent_z[l], dmin, dminl2e[l]);

        }

        double dmin_global;
        MPI_Allreduce (&dminl2e[0], &dminl2e_glb[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);

        for (l=0; l<ibm[ibi].n_elmt; l++) {
            double diff=dminl2e[l]-dminl2e_glb[l];
            if (diff>1.e-6) {
                count[l]=0;	
                ic[l]=0; jc[l]=0; kc[l]=0;
                iclose[l]=0; jclose[l]=0; kclose[l]=0;

                dhx_[l]=0.0;
                dhy_[l]=0.0;
                dhz_[l]=0.0;
            } else {
                count[l]=1;	
                iclose[l]=ic[l]; jclose[l]=jc[l]; kclose[l]=kc[l];
                i=ic[l]; j=jc[l]; k=kc[l];

                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                dhx_[l]=1.0/aj[k][j][i]/area;

                area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                dhy_[l]=1.0/aj[k][j][i]/area;

                area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                dhz_[l]=1.0/aj[k][j][i]/area;	
            }
        }

        MPI_Allreduce (&dhx_[0], &sum_dhx_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce (&dhy_[0], &sum_dhy_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce (&dhz_[0], &sum_dhz_[0], ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        
        MPI_Allreduce (&iclose[0], &sum_iclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce (&jclose[0], &sum_jclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce (&kclose[0], &sum_kclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        
        MPI_Allreduce (&count[0], &sum_count[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            if (sum_count[l]<1) sum_count[l]=1;
            ibm[ibi].dhx[l]=sum_dhx_[l]/(sum_count[l]);
            ibm[ibi].dhy[l]=sum_dhy_[l]/(sum_count[l]);
            ibm[ibi].dhz[l]=sum_dhz_[l]/(sum_count[l]);

            iclose[l]=sum_iclose[l]/(sum_count[l]);
            jclose[l]=sum_jclose[l]/(sum_count[l]);
            kclose[l]=sum_kclose[l]/(sum_count[l]);
        }

        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            int ic=iclose[l];
            int jc=jclose[l];
            int kc=kclose[l];

            //Debug Kflora
            int ifCasek;
            ifCasek=0;

            int ic1=ic-Itpwidth, ic2=ic+Itpwidth; 
        
            if (ic1>lxe||ic2<lxs) {
                ibm[ibi].i_min[l]=lxs;
                ibm[ibi].i_max[l]=lxs;		
            } else {
                ibm[ibi].i_min[l]=PetscMax(ic1, lxs);
                ibm[ibi].i_max[l]=PetscMin(ic2, lxe);
            }

            int jc1=jc-Itpwidth, jc2=jc+Itpwidth; 

            if (jc1>lye||jc2<lys) {
                ibm[ibi].j_min[l]=lys;
                ibm[ibi].j_max[l]=lys;
            } else {
                ibm[ibi].j_min[l]=PetscMax(jc1, lys);
                ibm[ibi].j_max[l]=PetscMin(jc2, lye);
            }

            int kc1=kc-Itpwidth, kc2=kc+Itpwidth; 

            if (kc1>lze||kc2<lzs) {
                ibm[ibi].k_min[l]=lzs;
                ibm[ibi].k_max[l]=lzs;
                //ibm[ibi].k_max[l]=lze;

                ifCasek=1;
            } else {
                ibm[ibi].k_min[l]=PetscMax(kc1, lzs);
                ibm[ibi].k_max[l]=PetscMin(kc2, lze);
            }
            //PetscPrintf(PETSC_COMM_WORLD, "Debug for element %i: Case %i  => imin and imax = %i, %i\n\n", l,ifCasei,  ibm[ibi].i_min[l],ibm[ibi].i_max[l]);
            //PetscPrintf(PETSC_COMM_WORLD, "Debug for element %i: Case %i  => jmin and jmax = %i, %i\n\n", l, ifCasej,  ibm[ibi].j_min[l],ibm[ibi].j_max[l]);
            //PetscPrintf(PETSC_COMM_WORLD, "Debug ic1 = %i, lxe= %i, ic2 = %i, lxs = %i\n", ic1, lxe, ic2, lxs);
            //PetscPrintf(PETSC_COMM_WORLD, "Debug jc1 = %i, lye= %i, jc2 = %i, lys = %i\n", jc1, lye, jc2, lys);

            //PetscPrintf(PETSC_COMM_WORLD, "Debug kc1 = %i, lze= %i, kc2 = %i, lzs = %i\n", kc1, lze, kc2, lzs);
            //PetscPrintf(PETSC_COMM_WORLD, "Debug for element %i: Case %i  => kmin and kmax = %i, %i\n\n", l, ifCasek, ibm[ibi].k_min[l],ibm[ibi].k_max[l]);

        }
    }

    DAVecRestoreArray(fda, Coor, &coor);

    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);

    return(0);

}


/* Read information of ACL, like chord, twist angle, lift and drag coefficients */
// call in the main.c for actuator line simulation only
PetscErrorCode airfoil_ACL(IBMNodes *ibm,  FSInfo *fsi, int NumberOfObjects)
{

    int	rank;
    int	n_CD, n_CL;
    int	ifoil, i, ibi;
    int 	iturbine;

    char string[128];

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if(!rank) { // root processor read in the data
        FILE *fd;
        char filen[80]; 

        for(ibi=0; ibi<NumberOfObjects; ibi++) {
            PetscPrintf(PETSC_COMM_WORLD, "Begin Foil Data %i of %i Objects\n", ibi, NumberOfObjects);
            for(ifoil=0; ifoil<ibm[ibi].num_foiltype; ifoil++) {
                PetscPrintf(PETSC_COMM_WORLD, "Begin Reading airfoil data for Foil %i of %i for ACL\n", ifoil+1, ibm[ibi].num_foiltype);

                if (TheSameObject) 
                    sprintf(filen,"%s/FOIL%3.3d_%2.2d" , path, 0, ifoil);
                else 
                    sprintf(filen,"%s/FOIL%3.3d_%2.2d" , path, ibi, ifoil);

                fd = fopen(filen, "r"); 
                if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open airfoil data file"), exit(0);
                if (fd) {
                    fgets(string, 128, fd);
                    fgets(string, 128, fd);

                    fscanf(fd, "%i ",&(ibm[ibi].acl[ifoil].num_AC));

                    MPI_Bcast(&(ibm[ibi].acl[ifoil].num_AC), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                    PetscMalloc((ibm[ibi].acl[ifoil].num_AC)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].r_AC));  
                    PetscMalloc((ibm[ibi].acl[ifoil].num_AC)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].angle_twistInp));  
                    PetscMalloc((ibm[ibi].acl[ifoil].num_AC)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].chord_bladeInp)); 

                    for(i=0; i<ibm[ibi].acl[ifoil].num_AC; i++) {
                        fscanf(fd, "%le %le %le", &(ibm[ibi].acl[ifoil].r_AC[i]), &(ibm[ibi].acl[ifoil].chord_bladeInp[i]), &(ibm[ibi].acl[ifoil].angle_twistInp[i]));

                        ibm[ibi].acl[ifoil].r_AC[i] = ibm[ibi].acl[ifoil].r_AC[i]/ reflength_wt;
                        ibm[ibi].acl[ifoil].chord_bladeInp[i] = ibm[ibi].acl[ifoil].chord_bladeInp[i] / reflength_wt;

                    } 

                    ibm[ibi].acl[ifoil].r_beg = ibm[ibi].acl[ifoil].r_AC[0];
                    ibm[ibi].acl[ifoil].r_end = ibm[ibi].acl[ifoil].r_AC[ibm[ibi].acl[ifoil].num_AC-1];

                    MPI_Bcast(ibm[ibi].acl[ifoil].r_AC, ibm[ibi].acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                    MPI_Bcast(ibm[ibi].acl[ifoil].angle_twistInp, ibm[ibi].acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                    MPI_Bcast(ibm[ibi].acl[ifoil].chord_bladeInp, ibm[ibi].acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                    MPI_Bcast(&(ibm[ibi].acl[ifoil].r_beg), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                    MPI_Bcast(&(ibm[ibi].acl[ifoil].r_end), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 

                }

                fclose(fd);
                PetscPrintf(PETSC_COMM_WORLD, "Finished Reading airfoil data for Foil %i for ACL\n\n", ifoil);
                
                PetscPrintf(PETSC_COMM_WORLD, "Begin Reading Aerodynamic data  %i\n", ifoil);

                if (TheSameObject) 
                    sprintf(filen,"%s/Aerodym%3.3d_%2.2d" , path, 0, ifoil);
                else
                    sprintf(filen,"%s/Aerodym%3.3d_%2.2d" , path, ibi, ifoil);

                fd = fopen(filen, "r"); 
                if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open Aerodym file"), exit(0);
                if (fd) {
                    fgets(string, 128, fd);
                    fgets(string, 128, fd);

                    fscanf(fd, "%i ",&(ibm[ibi].acl[ifoil].num_Aerodym));

                    MPI_Bcast(&(ibm[ibi].acl[ifoil].num_Aerodym), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                    PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].ang_Aerodym));
                    PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].CLInp));
                    PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].CDInp));
                    PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].CMInp));

                    for(i=0; i<ibm[ibi].acl[ifoil].num_Aerodym; i++) {
                        fscanf(fd, "%le %le %le %le", &(ibm[ibi].acl[ifoil].ang_Aerodym[i]), &(ibm[ibi].acl[ifoil].CLInp[i]), &(ibm[ibi].acl[ifoil].CDInp[i]), &(ibm[ibi].acl[ifoil].CMInp[i]));
                    }

                    MPI_Bcast(ibm[ibi].acl[ifoil].ang_Aerodym, ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(ibm[ibi].acl[ifoil].CLInp,  ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(ibm[ibi].acl[ifoil].CDInp,  ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(ibm[ibi].acl[ifoil].CMInp,  ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);

                }

                fclose(fd);
                PetscPrintf(PETSC_COMM_WORLD, "Finished Reading Aerodynamic data for Foil %i for ACL\n\n", ifoil);

            }
        }

        if (specifycirculation) {
            for (iturbine=0;iturbine<NumberOfTurbines;iturbine++) {
                sprintf(filen,"%s/circulation%2.2d" , path, iturbine);
                fd = fopen(filen, "r");
                if (!fd) PetscPrintf(PETSC_COMM_SELF, "Cannot open circulation file"), exit(0);
                if (fd) {
                    fgets(string, 128, fd);
                    fgets(string, 128, fd);

                    fscanf(fd, "%i ",&(ibm[iturbine].num_circulationInp));

                    MPI_Bcast(&(ibm[iturbine].num_circulationInp), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                    PetscMalloc((ibm[iturbine].num_circulationInp)*sizeof(PetscReal), &(ibm[iturbine].r_circulationInp));
                    PetscMalloc((ibm[iturbine].num_circulationInp)*sizeof(PetscReal), &(ibm[iturbine].circulationInp));

                    for(i=0; i<ibm[iturbine].num_circulationInp; i++) {
                        fscanf(fd, "%le %le", &(ibm[iturbine].r_circulationInp[i]), &(ibm[iturbine].circulationInp[i]));
                    }

                    MPI_Bcast(ibm[iturbine].r_circulationInp, ibm[iturbine].num_circulationInp, MPIU_REAL, 0, PETSC_COMM_WORLD);
                    MPI_Bcast(ibm[iturbine].circulationInp, ibm[iturbine].num_circulationInp, MPIU_REAL, 0, PETSC_COMM_WORLD);
                }
                fclose(fd);
            }
        }

    } else if(rank) {

        for(ibi=0; ibi<NumberOfObjects; ibi++) {
            for(ifoil=0; ifoil<ibm[ibi].num_foiltype; ifoil++) {

                // angle of twist , chord length of blade
                MPI_Bcast(&(ibm[ibi].acl[ifoil].num_AC), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                PetscMalloc((ibm[ibi].acl[ifoil].num_AC)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].r_AC));  
                PetscMalloc((ibm[ibi].acl[ifoil].num_AC)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].angle_twistInp));  
                PetscMalloc((ibm[ibi].acl[ifoil].num_AC)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].chord_bladeInp)); 
        
                MPI_Bcast(ibm[ibi].acl[ifoil].r_AC, ibm[ibi].acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                MPI_Bcast(ibm[ibi].acl[ifoil].angle_twistInp, ibm[ibi].acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                MPI_Bcast(ibm[ibi].acl[ifoil].chord_bladeInp, ibm[ibi].acl[ifoil].num_AC, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                MPI_Bcast(&(ibm[ibi].acl[ifoil].r_beg), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 
                MPI_Bcast(&(ibm[ibi].acl[ifoil].r_end), 1, MPIU_REAL, 0, PETSC_COMM_WORLD); 

                // life drag moment coefficients
                MPI_Bcast(&(ibm[ibi].acl[ifoil].num_Aerodym), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].ang_Aerodym));
                PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].CLInp));
                PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].CDInp));
                PetscMalloc((ibm[ibi].acl[ifoil].num_Aerodym)*sizeof(PetscReal), &(ibm[ibi].acl[ifoil].CMInp));

                MPI_Bcast(ibm[ibi].acl[ifoil].ang_Aerodym, ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(ibm[ibi].acl[ifoil].CLInp,  ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(ibm[ibi].acl[ifoil].CDInp,  ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(ibm[ibi].acl[ifoil].CMInp,  ibm[ibi].acl[ifoil].num_Aerodym, MPIU_REAL, 0, PETSC_COMM_WORLD);
            }
        }

        if (specifycirculation) {
            for (iturbine=0;iturbine<NumberOfTurbines;iturbine++) {
                MPI_Bcast(&(ibm[iturbine].num_circulationInp), 1, MPI_INT, 0, PETSC_COMM_WORLD);

                PetscMalloc((ibm[iturbine].num_circulationInp)*sizeof(PetscReal), &(ibm[iturbine].r_circulationInp));
                PetscMalloc((ibm[iturbine].num_circulationInp)*sizeof(PetscReal), &(ibm[iturbine].circulationInp));

                MPI_Bcast(ibm[iturbine].r_circulationInp, ibm[iturbine].num_circulationInp, MPIU_REAL, 0, PETSC_COMM_WORLD);
                MPI_Bcast(ibm[iturbine].circulationInp, ibm[iturbine].num_circulationInp, MPIU_REAL, 0, PETSC_COMM_WORLD);
            }
        }

    }

    int j;
    double r, fac1, fac2;
    for (ibi=0;ibi<NumberOfTurbines;ibi++) {
        for (i=0;i<ibm[ibi].n_elmt;i++) {
            r = sqrt( pow((ibm[ibi].cent_x[i]-fsi[ibi].x_c),2.0) + pow((ibm[ibi].cent_y[i]-fsi[ibi].y_c),2.0) + pow((ibm[ibi].cent_z[i]-fsi[ibi].z_c),2.0));	

            //if (!rank) printf("Turbine %d #### elmt %d cent_x %le x_c %le \n", ibi, i,  ibm[ibi].cent_x[i], fsi[ibi].x_c);
            //if (!rank) printf("Turbine %d #### elmt %d cent_y %le y_c %le \n", ibi, i,  ibm[ibi].cent_y[i], fsi[ibi].y_c);
            //if (!rank) printf("Turbine %d #### elmt %d cent_z %le z_c %le \n", ibi, i,  ibm[ibi].cent_z[i], fsi[ibi].z_c);
            for(ifoil=0; ifoil<ibm[ibi].num_foiltype; ifoil++) {

                if( r>=ibm[ibi].acl[ifoil].r_beg && r<=ibm[ibi].acl[ifoil].r_end ) {
                    for(j=0; j<ibm[ibi].acl[ifoil].num_AC-1; j++) {
                        if (r>=ibm[ibi].acl[ifoil].r_AC[j] && r<=ibm[ibi].acl[ifoil].r_AC[j+1]) {
                            fac1 = (ibm[ibi].acl[ifoil].r_AC[j+1]-r)/(ibm[ibi].acl[ifoil].r_AC[j+1]-ibm[ibi].acl[ifoil].r_AC[j]);
                            fac2 = (-ibm[ibi].acl[ifoil].r_AC[j]+r)/(ibm[ibi].acl[ifoil].r_AC[j+1]-ibm[ibi].acl[ifoil].r_AC[j]);
                            ibm[ibi].angle_twist[i]	= fac1*ibm[ibi].acl[ifoil].angle_twistInp[j] + fac2*ibm[ibi].acl[ifoil].angle_twistInp[j+1];
                            ibm[ibi].chord_blade[i]	= fac1*ibm[ibi].acl[ifoil].chord_bladeInp[j] + fac2*ibm[ibi].acl[ifoil].chord_bladeInp[j+1];
                        }
                    }
                }
            }

            //if (!rank) printf("Turbine %d #### %d th elmt the chord %le and twist angle %le \n", ibi, i,  ibm[ibi].chord_blade[i], ibm[ibi].angle_twist[i]);
        }

    }

    if (specifycirculation) {
        for (ibi=0;ibi<NumberOfTurbines;ibi++) {
            for (i=0;i<ibm[ibi].n_elmt;i++) {
                r = sqrt( pow((ibm[ibi].cent_x[i]-fsi[ibi].x_c),2.0) + pow((ibm[ibi].cent_y[i]-fsi[ibi].y_c),2.0) + pow((ibm[ibi].cent_z[i]-fsi[ibi].z_c),2.0));	
                for(j=0; j<ibm[ibi].num_circulationInp-1; j++) {
                    if (r>=ibm[ibi].r_circulationInp[j] && r<=ibm[ibi].r_circulationInp[j+1]) {
                        fac1 = (ibm[ibi].r_circulationInp[j+1]-r)/(ibm[ibi].r_circulationInp[j+1]-ibm[ibi].r_circulationInp[j]);
                        fac2 = (-ibm[ibi].r_circulationInp[j]+r)/(ibm[ibi].r_circulationInp[j+1]-ibm[ibi].r_circulationInp[j]);
                        ibm[ibi].circulationSpecified[i] = fac1*ibm[ibi].circulationInp[j] + fac2*ibm[ibi].circulationInp[j+1];
                    }
                }
            }
        }
    }
    PetscPrintf(PETSC_COMM_WORLD, "Finished Reading the airfoil data for ACL\n");

    return(0);

}


// calculate the force on the actuator disk, called in the solvers.c beforce solving the momentum equation
//
// NonUniform_ADModel: assume the velocity on each element has the same direction as the disk-averaged velocity
PetscErrorCode Calc_F_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt	l, ibi;
    double	pi = 3.141592653589793, a = 0.25;
    double 	C_T;  // = 4.0 / 3.0;
    double 	U_ref, A_sum;
    PetscReal	sign;
    double 	Uref_x, Uref_y, Uref_z;


    Calc_U_lagr(user, ibm, fsi, NumberOfObjects);


    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        if (rotor_model==4) {
            C_T=ibm[ibi].CT;
        } else {
            double indf_ax=ibm[ibi].indf_axis;
            C_T = 4.0 * indf_ax * (1-indf_ax);
            C_T = C_T / ( (1.0 - indf_ax)* (1.0 - indf_ax) );
        }

        if (rotor_model==4) U_ref=ibm[ibi].U_ref;
        else {
            U_ref = 0.0;
            A_sum = 0.0;

            double Uref_sum=0.0;

            for (l=0; l<ibm[ibi].n_elmt; l++) {
                double Uref = ibm[ibi].U_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].U_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].U_lagr_z[l]*fsi[ibi].nz_tb;
                Uref_sum += Uref*ibm[ibi].dA[l];
                A_sum += ibm[ibi].dA[l];
            }

            U_ref = Uref_sum/(A_sum+1.e-19);
        }

        PetscPrintf(PETSC_COMM_WORLD, "**** The U_ref at %i th body: %le\n", ibi, U_ref);

        double UUref_sum=0.0;
        A_sum = 0.0;
        if (NonUniform_ADModel) { 
            for (l=0; l<ibm[ibi].n_elmt; l++) {
                double Uref = ibm[ibi].U_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].U_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].U_lagr_z[l]*fsi[ibi].nz_tb;
                UUref_sum += Uref*Uref*ibm[ibi].dA[l];
                A_sum += ibm[ibi].dA[l];
            }
        }

        double UUref = UUref_sum/(A_sum+1.e-19);

        PetscPrintf(PETSC_COMM_WORLD, "**** The UU_ref at %i th body: %le\n", ibi, UUref);

        //U_ref = 1.0; //RTD
        //sign = 1.0;
        for (l=0; l<ibm[ibi].n_elmt; l++) {
            if (NonUniform_ADModel) {
                double Uref = ibm[ibi].U_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].U_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].U_lagr_z[l]*fsi[ibi].nz_tb;
                double UU = Uref*Uref;			

                sign = U_ref / fabs(U_ref+1.e-9);
                sign *= (UU/UUref);
            } else sign = U_ref / fabs(U_ref+1.e-9);

            ibm[ibi].F_lagr_x[l] = -0.5 * C_T * (U_ref * U_ref) * sign * fsi[ibi].nx_tb;
            ibm[ibi].F_lagr_y[l] = -0.5 * C_T * (U_ref * U_ref) * sign * fsi[ibi].ny_tb;
            ibm[ibi].F_lagr_z[l] = -0.5 * C_T * (U_ref * U_ref) * sign * fsi[ibi].nz_tb;  
        }

    }

    return(0);

}


// compute the actual shear on the nacelle surface for the last time step 
PetscErrorCode Comput_actualshear_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt	l, ibi;
    double		pi = 3.141592653589793, a = 0.25;
    double 		C_T;  // = 4.0 / 3.0;
    double 		U_ref, A_sum;
    PetscReal	sign;
    double 		Uref_x, Uref_y, Uref_z;
    int 		nv1, nv2, nv3;
    double 		r1, r2, r3, dh, rx, ry, rz;

    double 		nfx, nfy, nfz;

    Cmpnts 		Ut, Un, Ubt, Ubn, Ub;

    Cmpnts 		Ut_IP, Un_IP;

    PetscInt 	deltafunc_save = deltafunc;
    //PetscBarrier(PETSC_NULL);
    //PetscPrintf(PETSC_COMM_WORLD, "Debug3-1\n");
    Calc_U_lagr(user, ibm, fsi, NumberOfObjects);
    //PetscBarrier(PETSC_NULL);
    //PetscPrintf(PETSC_COMM_WORLD, "Debug3-2\n");
    Calc_U_IPlagr(user, ibm, fsi, NumberOfObjects);
    //PetscBarrier(PETSC_NULL);
    //PetscPrintf(PETSC_COMM_WORLD, "Debug3-3\n");
    Intp_Nut_lagr(user, ibm, fsi, NumberOfObjects);
    //PetscBarrier(PETSC_NULL);
    //PetscPrintf(PETSC_COMM_WORLD, "Debug3-4\n");


  	for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double drag_friction_desired=0.0;
        double drag_friction_actual=0.0;
        double rhs_sum=0.0;

        // Calculate the desired shear, and the estimated shear on the surface 
        double num=0.0;
        double sum_dh=0.0;

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            if (ibm[ibi].color[l]==1000) {

                ibm[ibi].Shear_lagr_x[l] = 0.0;
                ibm[ibi].Shear_lagr_y[l] = 0.0;
                ibm[ibi].Shear_lagr_z[l] = 0.0;

            } else {

                nv1=ibm[ibi].nv1[l];
                nv2=ibm[ibi].nv2[l];
                nv3=ibm[ibi].nv3[l];

                dh = ibm[ibi].dh;

                sum_dh += dh;

                nfx=ibm[ibi].nf_x[l];
                nfy=ibm[ibi].nf_y[l];
                nfz=ibm[ibi].nf_z[l];

                //velocity boundary conditions at surface 
                Ub.x = (ibm[ibi].u[nv1].x +  ibm[ibi].u[nv2].x +  ibm[ibi].u[nv3].x)/3.0; 
                Ub.y = (ibm[ibi].u[nv1].y +  ibm[ibi].u[nv2].y +  ibm[ibi].u[nv3].y)/3.0; 
                Ub.z = (ibm[ibi].u[nv1].z +  ibm[ibi].u[nv2].z +  ibm[ibi].u[nv3].z)/3.0; 
                double UUn = Ub.x*nfx +  Ub.y*nfy +  Ub.z*nfz; 

                Ubn.x = UUn*nfx; 
                Ubn.y = UUn*nfy; 
                Ubn.z = UUn*nfz; 

                Ubt.x = Ub.x - Ubn.x;
                Ubt.y = Ub.y - Ubn.y;
                Ubt.z = Ub.z - Ubn.z;

                //interpolated velocity at surface 
                UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 
                Un.x = UUn*nfx; 
                Un.y = UUn*nfy; 
                Un.z = UUn*nfz; 

                Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
                Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
                Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;

                //interpolated velocity at IP points
                UUn = ibm[ibi].U_IPlagr_x[l]*nfx +  ibm[ibi].U_IPlagr_y[l]*nfy +  ibm[ibi].U_IPlagr_z[l]*nfz; 

                Un_IP.x = UUn*nfx; 
                Un_IP.y = UUn*nfy; 
                Un_IP.z = UUn*nfz; 

                Ut_IP.x = ibm[ibi].U_IPlagr_x[l] - Un.x;
                Ut_IP.y = ibm[ibi].U_IPlagr_y[l] - Un.y;
                Ut_IP.z = ibm[ibi].U_IPlagr_z[l] - Un.z;

                // actual shear for previous time step
                double nu=1./user->ren;
        
                double Urefx = Ut.x-Ut_IP.x;
                double Urefy = Ut.y-Ut_IP.y;
                double Urefz = Ut.z-Ut_IP.z;

                ibm[ibi].Shear_lagr_x[l] = (nu+ibm[ibi].Nut_lagr[l])*Urefx/(ibm[ibi].dh_IP[l]+1.e-19);
                ibm[ibi].Shear_lagr_y[l] = (nu+ibm[ibi].Nut_lagr[l])*Urefy/(ibm[ibi].dh_IP[l]+1.e-19);
                ibm[ibi].Shear_lagr_z[l] = (nu+ibm[ibi].Nut_lagr[l])*Urefz/(ibm[ibi].dh_IP[l]+1.e-19);

//              PetscPrintf(PETSC_COMM_WORLD, "IPlagr_z_Shear = %le\n", ibm[ibi].U_IPlagr_z[l]);
            } 

        }
    }

    return(0);
}


// compute the nut on the nacelle surface for the last time step 

PetscErrorCode Comput_nut_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt	l, ibi;
    double		pi = 3.141592653589793, a = 0.25;
    double 		C_T;  // = 4.0 / 3.0;
    double 		U_ref, A_sum;
    PetscReal	sign;
    double 		Uref_x, Uref_y, Uref_z;
    int 		nv1, nv2, nv3;
    double 		r1, r2, r3, dh, rx, ry, rz;

    double 		nfx, nfy, nfz;

    Cmpnts 		Ut, Un, Ubt, Ubn, Ub;

    Cmpnts 		Ut_IP, Un_IP;

    PetscInt 	deltafunc_save = deltafunc;

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double drag_friction_desired=0.0;
        double drag_friction_actual=0.0;
        double rhs_sum=0.0;

        // Calculate the desired shear, and the estimated shear on the surface 
        double num=0.0;
        double sum_dh=0.0;

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            if (ibm[ibi].color[l]==1000) {

                ibm[ibi].Nut_lagr[l] = 0.0;

            } else {

                nv1=ibm[ibi].nv1[l];
                nv2=ibm[ibi].nv2[l];
                nv3=ibm[ibi].nv3[l];

                //dh = pow(ibm[ibi].dhx[l]*ibm[ibi].dhy[l]*ibm[ibi].dhz[l],1/3); 
                //dh = pow(r1*r2*r3,1.0/3.0);

                dh = ibm[ibi].dh;

                sum_dh += dh;

                nfx=ibm[ibi].nf_x[l];
                nfy=ibm[ibi].nf_y[l];
                nfz=ibm[ibi].nf_z[l];

                //interpolated velocity at surface 
                double UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 
                Un.x = UUn*nfx; 
                Un.y = UUn*nfy; 
                Un.z = UUn*nfz; 

                Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
                Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
                Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;

                //interpolated velocity at IP points
                UUn = ibm[ibi].U_IPlagr_x[l]*nfx +  ibm[ibi].U_IPlagr_y[l]*nfy +  ibm[ibi].U_IPlagr_z[l]*nfz; 

                Un_IP.x = UUn*nfx; 
                Un_IP.y = UUn*nfy; 
                Un_IP.z = UUn*nfz; 

                Ut_IP.x = ibm[ibi].U_IPlagr_x[l] - Un.x;
                Ut_IP.y = ibm[ibi].U_IPlagr_y[l] - Un.y;
                Ut_IP.z = ibm[ibi].U_IPlagr_z[l] - Un.z;

                // actual shear for previous time step
                double nu=1./user->ren;
        
                double Urefx = Ut.x-Ut_IP.x;
                double Urefy = Ut.y-Ut_IP.y;
                double Urefz = Ut.z-Ut_IP.z;


                double ShearDesired_mag = sqrt(pow(ibm[ibi].ShearDesired_lagr_x[l],2)+pow(ibm[ibi].ShearDesired_lagr_y[l],2)+pow(ibm[ibi].ShearDesired_lagr_z[l],2));
                double ShearActual_mag = sqrt(pow(ibm[ibi].Shear_lagr_x[l],2)+pow(ibm[ibi].Shear_lagr_y[l],2)+pow(ibm[ibi].Shear_lagr_z[l],2));
                double dUdn_mag = sqrt(pow(Urefx/(ibm[ibi].dh_IP[l]+1.e-19),2)+pow(Urefy/(ibm[ibi].dh_IP[l]+1.e-19),2)+pow(Urefz/(ibm[ibi].dh_IP[l]+1.e-19),2));
                ibm[ibi].Nut_lagr[l] = fmax(ShearDesired_mag/(dUdn_mag+1.e-19)-nu,0.0);

            } 

        }
    }

    return(0);
}


// compute the model coefficient on the nacelle surface for the current time step 
PetscErrorCode Comput_modelcoef_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt    l, ibi;
    double      pi = 3.141592653589793, a = 0.25;
    double      C_T;  // = 4.0 / 3.0;
    double      U_ref, A_sum;
    PetscReal   sign;
    double      Uref_x, Uref_y, Uref_z;
    int         nv1, nv2, nv3;
    double      r1, r2, r3, dh, rx, ry, rz;

    double      nfx, nfy, nfz;

    Cmpnts      Ut, Un, Ubt, Ubn, Ub;

    Cmpnts      Ut_IP, Un_IP;

    PetscInt    deltafunc_save = deltafunc;


    // compute correction of model coefficient    

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double drag_friction_desired=0.0;
        double drag_friction_actual=0.0;
        double UU_sum=0.0;

        double A_sum=0.0; 

        double num=0.0;
        double sum_dh=0.0;

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            if (ibm[ibi].color[l]!=1000) {

                // actual friction drag for previous time step
                drag_friction_actual += (ibm[ibi].Shear_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].Shear_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].Shear_lagr_z[l]*fsi[ibi].nz_tb)*ibm[ibi].dA[l];

                // desired friction drag for previous time step
                drag_friction_desired += (ibm[ibi].ShearDesired_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].ShearDesired_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].ShearDesired_lagr_z[l]*fsi[ibi].nz_tb)*ibm[ibi].dA[l];
                UU_sum += (ibm[ibi].UU_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].UU_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].UU_lagr_z[l]*fsi[ibi].nz_tb)*ibm[ibi].dA[l];

                A_sum += ibm[ibi].dA[l];

                double ShearDesired_mag = sqrt(pow(ibm[ibi].ShearDesired_lagr_x[l],2)+pow(ibm[ibi].ShearDesired_lagr_y[l],2)+pow(ibm[ibi].ShearDesired_lagr_z[l],2));
                double ShearActual_mag = sqrt(pow(ibm[ibi].Shear_lagr_x[l],2)+pow(ibm[ibi].Shear_lagr_y[l],2)+pow(ibm[ibi].Shear_lagr_z[l],2));
                double UU_mag = sqrt(pow(ibm[ibi].UU_lagr_x[l],2)+pow(ibm[ibi].UU_lagr_y[l],2)+pow(ibm[ibi].UU_lagr_z[l],2));

                if (fabs(UU_mag)>1.e-9) {
                    ibm[ibi].frictionfactor[l] += (ShearDesired_mag-ShearActual_mag)/(UU_mag+1.e-9);
                } else {
                    ibm[ibi].frictionfactor[l] += 0.0;
                }
            }
        }

        PetscPrintf(PETSC_COMM_WORLD, "1Compute Nacelle Coeff  ibi[%i]. friction factor = %f\n",ibi, ibm[ibi].friction_factor );

        if (fabs(UU_sum)>1.e-9) {
            ibm[ibi].friction_factor += (drag_friction_desired-drag_friction_actual)/(UU_sum+1.e-9);
        } else {
            ibm[ibi].friction_factor += 0.0;
        }
        
        PetscPrintf(PETSC_COMM_WORLD, "2Compute Nacelle Coeff  ibi[%i]. friction factor = %f\n",ibi, ibm[ibi].friction_factor );


        ibm[ibi].friction_factor = fmax(ibm[ibi].friction_factor, 0);
        ibm[ibi].friction_factor = fmin(ibm[ibi].friction_factor, 5);

        double Rex = ibm[ibi].U_ref*ibm[ibi].r_nacelle*user->ren/reflength_wt+1.0;

        double frictionfactor=0.370*pow(log10(Rex),-2.584); 
        PetscPrintf(PETSC_COMM_WORLD, "3Computed Nacelle Coeff  friction factor2 = %f\n",frictionfactor);

        double ustar = sqrt(0.5*frictionfactor*pow(ibm[ibi].U_ref,2));

        PetscPrintf(PETSC_COMM_WORLD, "Ustar = %le, ibm.dh=%f and ren=%le\n", ustar, ibm[ibi].dh, user->ren);
        
        double Uib = ustar*(log(ibm[ibi].dh*ustar*user->ren)/0.4+5.1);

        double delta = 0.382*ibm[ibi].r_nacelle/pow(Rex,0.2);

        double Udelta = ustar*(log(delta*ustar*user->ren)/0.4+5.1);

        PetscPrintf(PETSC_COMM_WORLD, "Uib = %le\n", Uib);

        ibm[ibi].coeff_SG = pow(Udelta/Uib,2);

        ibm[ibi].coeff_SG = fmax(ibm[ibi].coeff_SG, 1);
        ibm[ibi].coeff_SG = fmin(ibm[ibi].coeff_SG, 10);

        //if (fsi[ibi].rotate_alongaxis) ibm[ibi].friction_factor = 0;
    }

    return(0);
    }


// compute the desired shear on the nacelle surface for the current time step 
    PetscErrorCode Comput_desiredshear_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
    {

    PetscInt	l, ibi;
    double		pi = 3.141592653589793, a = 0.25;
    double 		C_T;  // = 4.0 / 3.0;
    double 		U_ref, A_sum;
    PetscReal	sign;
    double 		Uref_x, Uref_y, Uref_z;
    int 		nv1, nv2, nv3;
    double 		r1, r2, r3, dh, rx, ry, rz;

    double 		nfx, nfy, nfz;

    Cmpnts 		Ut, Un, Ubt, Ubn, Ub;

    Cmpnts 		Ut_IP, Un_IP;

    PetscInt 	deltafunc_save = deltafunc;

    // calculate the desired shear for the current time step	
    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double drag_friction_desired=0.0;
        double drag_friction_actual=0.0;
        double rhs_sum=0.0;

        //ibm[ibi].friction_factor=0.0;

        // Calculate the desired shear, and the estimated shear on the surface 
        double num=0.0;
        double sum_dh=0.0;

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            if (ibm[ibi].color[l]==1000) {

                ibm[ibi].ShearDesired_lagr_x[l] = 0.0;
                ibm[ibi].ShearDesired_lagr_y[l] = 0.0;
                ibm[ibi].ShearDesired_lagr_z[l] = 0.0;

            } else {

                nv1=ibm[ibi].nv1[l];
                nv2=ibm[ibi].nv2[l];
                nv3=ibm[ibi].nv3[l];

                dh = ibm[ibi].dh;

                sum_dh += dh;

                nfx=ibm[ibi].nf_x[l];
                nfy=ibm[ibi].nf_y[l];
                nfz=ibm[ibi].nf_z[l];

                //velocity boundary conditions at surface 
                Ub.x = (ibm[ibi].u[nv1].x +  ibm[ibi].u[nv2].x +  ibm[ibi].u[nv3].x)/3.0; 
                Ub.y = (ibm[ibi].u[nv1].y +  ibm[ibi].u[nv2].y +  ibm[ibi].u[nv3].y)/3.0; 
                Ub.z = (ibm[ibi].u[nv1].z +  ibm[ibi].u[nv2].z +  ibm[ibi].u[nv3].z)/3.0; 
                double UUn = Ub.x*nfx +  Ub.y*nfy +  Ub.z*nfz; 

                Ubn.x = UUn*nfx; 
                Ubn.y = UUn*nfy; 
                Ubn.z = UUn*nfz; 

                Ubt.x = Ub.x - Ubn.x;
                Ubt.y = Ub.y - Ubn.y;
                Ubt.z = Ub.z - Ubn.z;

                //interpolated velocity at surface 
                UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 
                Un.x = UUn*nfx; 
                Un.y = UUn*nfy; 
                Un.z = UUn*nfz; 

                Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
                Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
                Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;

                //interpolated velocity at IP points
                
                double UIPx = ibm[ibi].U_IPlagr_x[l];
                double UIPy = ibm[ibi].U_IPlagr_y[l];
                double UIPz = ibm[ibi].U_IPlagr_z[l];
                
                UUn = UIPx*nfx + UIPy*nfy + UIPz*nfz; 

                Un_IP.x = UUn*nfx; 
                Un_IP.y = UUn*nfy; 
                Un_IP.z = UUn*nfz; 

                Ut_IP.x = UIPx - Un.x;
                Ut_IP.y = UIPy - Un.y;
                Ut_IP.z = UIPz - Un.z;

                double umag_t = sqrt(pow(Ut_IP.x,2)+pow(Ut_IP.y,2)+pow(Ut_IP.z,2))+1.e-9;

                double tx=Ut_IP.x/umag_t;
                double ty=Ut_IP.y/umag_t;
                double tz=Ut_IP.z/umag_t;

                Ut_IP.x = ibm[ibi].U_ref*tx;
                Ut_IP.y = ibm[ibi].U_ref*ty;
                Ut_IP.z = ibm[ibi].U_ref*tz;

                // friction factor
                double lx = ibm[ibi].cent_x[l]-fsi[ibi].xnacelle_upstreamend;
                double ly = ibm[ibi].cent_y[l]-fsi[ibi].ynacelle_upstreamend;
                double lz = ibm[ibi].cent_z[l]-fsi[ibi].znacelle_upstreamend;

                double L = fabs(lx*fsi[ibi].nx_tb+ly*fsi[ibi].ny_tb+lz*fsi[ibi].nz_tb);

                double frictionfactor;

                if (cf_nacelle_fromfile) {
                    int i, ifind;
                    double fac1, fac2;
                    for (i=1;i<ibm[ibi].num_cf;i++) {
                        if (L>ibm[ibi].r_in[i-1] && L<=ibm[ibi].r_in[i]) {
                            ifind = i;
                        }
                    }

                    fac1 = (ibm[ibi].r_in[ifind]-L)/(ibm[ibi].r_in[ifind]-ibm[ibi].r_in[ifind-1]);
                    fac2 = (-ibm[ibi].r_in[ifind-1]+L)/(ibm[ibi].r_in[ifind]-ibm[ibi].r_in[ifind-1]);

                    frictionfactor = ibm[ibi].cf_in[ifind-1]*fac1 + ibm[ibi].cf_in[ifind]*fac2;

                } else {
                    double Rex = ibm[ibi].U_ref*L*user->ren+1.0;

                    frictionfactor=0.370*pow(log10(Rex),-2.584); 
                    if (l==0) PetscPrintf(PETSC_COMM_WORLD, "Desired shear friction factor = %f\n",frictionfactor);
                }


                // desired shear for current time step 
                double Urefx = Ubt.x-Ut_IP.x;
                double Urefy = Ubt.y-Ut_IP.y;
                double Urefz = Ubt.z-Ut_IP.z;

                double Umag = sqrt(Urefx*Urefx+Urefy*Urefy+Urefz*Urefz)+1.e-9;

                ibm[ibi].ShearDesired_lagr_x[l] = 0.5 * frictionfactor * Urefx * Umag;
                ibm[ibi].ShearDesired_lagr_y[l] = 0.5 * frictionfactor * Urefy * Umag;
                ibm[ibi].ShearDesired_lagr_z[l] = 0.5 * frictionfactor * Urefz * Umag;  

                // without coefficient
                Urefx = Ubt.x-Ut.x;
                Urefy = Ubt.y-Ut.y;
                Urefz = Ubt.z-Ut.z;

                Umag = sqrt(Urefx*Urefx+Urefy*Urefy+Urefz*Urefz);

                ibm[ibi].UU_lagr_x[l] = 0.5 * Urefx * Umag;
                ibm[ibi].UU_lagr_y[l] = 0.5 * Urefy * Umag;
                ibm[ibi].UU_lagr_z[l] = 0.5 * Urefz * Umag;  

            }
        }

    }


    return(0);
}


// still working on 
// actuator nacelle model
PetscErrorCode Calc_F_lagr_nacelle1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt	l, ibi;
    double		pi = 3.141592653589793, a = 0.25;
    double 		C_T;  // = 4.0 / 3.0;
    double 		U_ref, A_sum;
    PetscReal	sign;
    double 		Uref_x, Uref_y, Uref_z;
    int 		nv1, nv2, nv3;
    double 		r1, r2, r3, dh, rx, ry, rz;

    double 		nfx, nfy, nfz;

    Cmpnts 		Ut, Un, Ubt, Ubn, Ub;

    Cmpnts 		Ut_IP, Un_IP;

    PetscInt 	deltafunc_save = deltafunc;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double sumforce_pressure=0.0, sumforce_friction=0.0;
        double sumforce1_pressure=0.0, sumforce1_friction=0.0;
        double sumforce2_friction=0.0;

        double sumforce_x=0.0, sumforce_y=0.0, sumforce_z=0.0;
        double sumforce1_x=0.0, sumforce1_y=0.0, sumforce1_z=0.0;


        double sum_dh=0.0;
        for (l=0; l<ibm[ibi].n_elmt; l++) {

            nv1=ibm[ibi].nv1[l];
            nv2=ibm[ibi].nv2[l];
            nv3=ibm[ibi].nv3[l];

            dh = ibm[ibi].dh;

            sum_dh += dh;

            nfx=ibm[ibi].nf_x[l];
            nfy=ibm[ibi].nf_y[l];
            nfz=ibm[ibi].nf_z[l];

            //velocity boundary conditions at surface 
            Ub.x = (ibm[ibi].u[nv1].x +  ibm[ibi].u[nv2].x +  ibm[ibi].u[nv3].x)/3.0; 
            Ub.y = (ibm[ibi].u[nv1].y +  ibm[ibi].u[nv2].y +  ibm[ibi].u[nv3].y)/3.0; 
            Ub.z = (ibm[ibi].u[nv1].z +  ibm[ibi].u[nv2].z +  ibm[ibi].u[nv3].z)/3.0; 
            double UUn = Ub.x*nfx +  Ub.y*nfy +  Ub.z*nfz; 

            Ubn.x = UUn*nfx; 
            Ubn.y = UUn*nfy; 
            Ubn.z = UUn*nfz; 

            Ubt.x = Ub.x - Ubn.x;
            Ubt.y = Ub.y - Ubn.y;
            Ubt.z = Ub.z - Ubn.z;

            double fx_n, fy_n, fz_n, fx_t, fy_t, fz_t;

            double Urefx, Urefy, Urefz, Umag;


            if (nacelle_model == 1) {
                //interpolated velocity at surface 
                double UUn;
                UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 
                Un.x = UUn*nfx; 
                Un.y = UUn*nfy; 
                Un.z = UUn*nfz; 

                Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
                Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
                Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;

                // normal force 
                Urefx = Ubn.x-Un.x;
                Urefy = Ubn.y-Un.y;
                Urefz = Ubn.z-Un.z;

                fx_n =  dh*Urefx/user->dt;
                fy_n =  dh*Urefy/user->dt;
                fz_n =  dh*Urefz/user->dt;  

                // tangential force 
                Urefx = Ubt.x-Ut.x;
                Urefy = Ubt.y-Ut.y;
                Urefz = Ubt.z-Ut.z;

                fx_t =  dh*Urefx/user->dt;
                fy_t =  dh*Urefy/user->dt;
                fz_t =  dh*Urefz/user->dt;  
            } else 
            {
                //interpolated velocity at surface 
                double UUn;
                UUn = ibm[ibi].U_lagr_x[l]*nfx +  ibm[ibi].U_lagr_y[l]*nfy +  ibm[ibi].U_lagr_z[l]*nfz; 
                Un.x = UUn*nfx; 
                Un.y = UUn*nfy; 
                Un.z = UUn*nfz; 

                Ut.x = ibm[ibi].U_lagr_x[l] - Un.x;
                Ut.y = ibm[ibi].U_lagr_y[l] - Un.y;
                Ut.z = ibm[ibi].U_lagr_z[l] - Un.z;


                // normal force
                Urefx = Ubn.x-Un.x;
                Urefy = Ubn.y-Un.y;
                Urefz = Ubn.z-Un.z;

                fx_n =  dh*Urefx/user->dt;
                fy_n =  dh*Urefy/user->dt;
                fz_n =  dh*Urefz/user->dt;  

                // tangential force 
                if (nacelle_model==2) {
                    fx_t=0;
                    fy_t=0;
                    fz_t=0;
                } else if (nacelle_model==3) {
                    Urefx = Ubt.x-Ut.x;
                    Urefy = Ubt.y-Ut.y;
                    Urefz = Ubt.z-Ut.z;

                    Umag = sqrt(Urefx*Urefx+Urefy*Urefy+Urefz*Urefz);

                    fx_t = ibm[ibi].ShearDesired_lagr_x[l];
                    fy_t = ibm[ibi].ShearDesired_lagr_y[l];
                    fz_t = ibm[ibi].ShearDesired_lagr_z[l];  
                } else if (nacelle_model==4) {
                    Urefx = Ubt.x-Ut.x;
                    Urefy = Ubt.y-Ut.y;
                    Urefz = Ubt.z-Ut.z;

                    Umag = sqrt(Urefx*Urefx+Urefy*Urefy+Urefz*Urefz);

                    fx_t = ibm[ibi].ShearDesired_lagr_x[l];
                    fy_t = ibm[ibi].ShearDesired_lagr_y[l];
                    fz_t = ibm[ibi].ShearDesired_lagr_z[l];  

                    fx_t += 0.5 * ibm[ibi].frictionfactor[l] * Urefx * Umag;
                    fy_t += 0.5 * ibm[ibi].frictionfactor[l] * Urefy * Umag;
                    fz_t += 0.5 * ibm[ibi].frictionfactor[l] * Urefz * Umag;  
                } else if (nacelle_model==5) {
                    Urefx = Ubt.x-Ut.x;
                    Urefy = Ubt.y-Ut.y;
                    Urefz = Ubt.z-Ut.z;

                    Umag = sqrt(Urefx*Urefx+Urefy*Urefy+Urefz*Urefz);

                    fx_t = ibm[ibi].ShearDesired_lagr_x[l];
                    fy_t = ibm[ibi].ShearDesired_lagr_y[l];
                    fz_t = ibm[ibi].ShearDesired_lagr_z[l];  

                    if (l==0) PetscPrintf(PETSC_COMM_WORLD, "Nacelle friction factor = %f\n",ibm[ibi].friction_factor );

                    fx_t += 0.5 * ibm[ibi].friction_factor * Urefx * Umag;
                    fy_t += 0.5 * ibm[ibi].friction_factor * Urefy * Umag;
                    fz_t += 0.5 * ibm[ibi].friction_factor * Urefz * Umag;  
                }

                if (ibm[ibi].color[l]==1000) {
                    fx_n=0.0; fy_n=0.0; fz_n=0.0;
                    fx_t=0.0; fy_t=0.0; fz_t=0.0;
                    ibm[ibi].Shear_lagr_x[l]=0.0;
                    ibm[ibi].Shear_lagr_y[l]=0.0;
                    ibm[ibi].Shear_lagr_z[l]=0.0;
                    ibm[ibi].P_IPlagr[l]=0.0;
                    ibm[ibi].Nut_lagr[l]=0.0;
                } 	

                ibm[ibi].F_lagr_x[l] = fx_n+fx_t;
                ibm[ibi].F_lagr_y[l] = fy_n+fy_t;
                ibm[ibi].F_lagr_z[l] = fz_n+fz_t;  

                sumforce_pressure += (fx_n*fsi[ibi].nx_tb + fy_n*fsi[ibi].ny_tb + fz_n*fsi[ibi].nz_tb)*ibm[ibi].dA[l];
                sumforce_friction += (fx_t*fsi[ibi].nx_tb + fy_t*fsi[ibi].ny_tb + fz_t*fsi[ibi].nz_tb)*ibm[ibi].dA[l];

                sumforce1_pressure += (ibm[ibi].P_IPlagr[l]*nfx*fsi[ibi].nx_tb + ibm[ibi].P_IPlagr[l]*nfy*fsi[ibi].ny_tb + ibm[ibi].P_IPlagr[l]*nfz*fsi[ibi].nz_tb)*ibm[ibi].dA[l];
                sumforce1_friction += (ibm[ibi].Shear_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].Shear_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].Shear_lagr_z[l]*fsi[ibi].nz_tb)*ibm[ibi].dA[l];
                sumforce2_friction += (ibm[ibi].ShearDesired_lagr_x[l]*fsi[ibi].nx_tb + ibm[ibi].ShearDesired_lagr_y[l]*fsi[ibi].ny_tb + ibm[ibi].ShearDesired_lagr_z[l]*fsi[ibi].nz_tb)*ibm[ibi].dA[l];


                sumforce_x += (fx_n + fx_t)*ibm[ibi].dA[l];
                sumforce_y += (fy_n + fy_t)*ibm[ibi].dA[l];
                sumforce_z += (fz_n + fz_t)*ibm[ibi].dA[l];

                sumforce1_x += (ibm[ibi].P_IPlagr[l]*nfx + ibm[ibi].Shear_lagr_x[l])*ibm[ibi].dA[l];
                sumforce1_y += (ibm[ibi].P_IPlagr[l]*nfy + ibm[ibi].Shear_lagr_y[l])*ibm[ibi].dA[l];
                sumforce1_z += (ibm[ibi].P_IPlagr[l]*nfz + ibm[ibi].Shear_lagr_z[l])*ibm[ibi].dA[l];
            }

            if (!rank) {
                FILE *f;
                char filen[80];
                sprintf(filen, "nacelleforcecoefficients%2.2d_%2.2d", nacelle_model, ibi);
                if (ti==1) {
                    f = fopen(filen, "w");
                    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=time coeff_SG Uref sumaxialforce_pressure_applied sumaxialforce_friction_applied sumaxialforce_pressure_actual sumaxialforce_friction_actual sumaxialforce_friction_desired sumforce_x_applied sumforce_y_applied sumforce_z_applied sumforce_x_actual sumforce_y_actual sumforce_z_actual dh\n");
                } else f = fopen(filen, "a");

                PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n",ti, ibm[ibi].coeff_SG, ibm[ibi].U_ref, sumforce_pressure, sumforce_friction, sumforce1_pressure, sumforce1_friction, sumforce2_friction, sumforce_x, sumforce_y, sumforce_z, sumforce1_x, sumforce1_y, sumforce1_z, sum_dh/(double)ibm[ibi].n_elmt );
                fclose(f);

            }
        }
    }
    return(0);
}


// calculate force at largrangian points using actuator line model
//
// called in solvers.c before solving the momentum equation
PetscErrorCode Export_ForceOnBlade(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt    l, ibi, j;
    double      pi = 3.141592653589793;
    PetscReal   A_sum, U_ref, rr, rx, ry, rz, tmp, u_tangent;
    Cmpnts      U_b, n_relV, n_L, n_blade, n_rot;	
    int         nv1, nv2, nv3, ifoil;
    PetscReal   fac1, fac2, r;

    int istall;

    int rank=0;

    PetscPrintf(PETSC_COMM_WORLD, "Export ForceOnBlade ti %i  tiout %i \n", ti, tiout);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (ti>tistart+1) {

        //PetscPrintf(PETSC_COMM_WORLD, "AL: export\n");
        count_AL += 1.0;
        FILE *f;
        char filen[80];
        
        double fac = 1.0 / count_AL;

        for (ibi=0; ibi<NumberOfObjects; ibi++) {

            if (ti == (ti/tiout) * tiout) {
                PetscPrintf(PETSC_COMM_WORLD, "Export ForceAOA 1 \n");
                if (!rank) {
                    sprintf(filen, "ForceAOA_%2.2d", ibi);
                    f = fopen(filen, "w");

                    PetscFPrintf(PETSC_COMM_WORLD, f, "r Ft Fa AOA C Urel FFt FFa AOAAOA Urel.x Urel.y Urel.z Uinduced.x Uinduced.y Uinduced.z circulation.x circulation.y circulation.z liftdirection.x liftdirection.y liftdirection.z CL CD\n");
                }
            }

            for (l=0; l<ibm[ibi].n_elmt; l++) {
                nv1 = ibm[ibi].nv1[l];
                nv2 = ibm[ibi].nv2[l];
                nv3 = ibm[ibi].nv3[l];

                rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
                ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
                rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;

                double r = sqrt(rx*rx+ry*ry+rz*rz);
                n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;

                double ux=0.5*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x); 
                double uy=0.5*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y); 
                double uz=0.5*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z); 
                
                double Ublade=sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));

                n_rot.x=ux/(Ublade+1.e-20);
                n_rot.y=uy/(Ublade+1.e-20);
                n_rot.z=uz/(Ublade+1.e-20);

                ibm[ibi].Ft_mean[l] += ibm[ibi].F_lagr_x[l]*n_rot.x+ibm[ibi].F_lagr_y[l]*n_rot.y+ibm[ibi].F_lagr_z[l]*n_rot.z;
                ibm[ibi].Fa_mean[l] += ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb;
                
                ibm[ibi].FFt_mean[l] += pow(ibm[ibi].F_lagr_x[l]*n_rot.x+ibm[ibi].F_lagr_y[l]*n_rot.y+ibm[ibi].F_lagr_z[l]*n_rot.z,2);
                ibm[ibi].FFa_mean[l] += pow(ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb,2);

                ibm[ibi].Urelmag_mean[l] += ibm[ibi].Urelmag[l];
                ibm[ibi].AOA_mean[l] += ibm[ibi].angle_attack[l];
                ibm[ibi].AOAAOA_mean[l] += pow(ibm[ibi].angle_attack[l],2);

                ibm[ibi].Urel_mean[l].x+=ibm[ibi].Urel[l].x;
                ibm[ibi].Urel_mean[l].y+=ibm[ibi].Urel[l].y;
                ibm[ibi].Urel_mean[l].z+=ibm[ibi].Urel[l].z;

                ibm[ibi].Uinduced_mean[l].x+=ibm[ibi].Uinduced[l].x;
                ibm[ibi].Uinduced_mean[l].y+=ibm[ibi].Uinduced[l].y;
                ibm[ibi].Uinduced_mean[l].z+=ibm[ibi].Uinduced[l].z;

                ibm[ibi].circulation_mean[l].x+=ibm[ibi].circulation[l].x;
                ibm[ibi].circulation_mean[l].y+=ibm[ibi].circulation[l].y;
                ibm[ibi].circulation_mean[l].z+=ibm[ibi].circulation[l].z;

                ibm[ibi].liftdirection_mean[l].x+=ibm[ibi].liftdirection[l].x;
                ibm[ibi].liftdirection_mean[l].y+=ibm[ibi].liftdirection[l].y;
                ibm[ibi].liftdirection_mean[l].z+=ibm[ibi].liftdirection[l].z;

                double Ft = -ibm[ibi].Ft_mean[l]*fac;
                double Fa = -ibm[ibi].Fa_mean[l]*fac;
                double AOA = ibm[ibi].AOA_mean[l]*fac;

                double FFt = ibm[ibi].FFt_mean[l]*fac-ibm[ibi].Ft_mean[l]*ibm[ibi].Ft_mean[l]*fac*fac;
                double FFa = ibm[ibi].FFa_mean[l]*fac-ibm[ibi].Fa_mean[l]*ibm[ibi].Fa_mean[l]*fac*fac;
                double AOAAOA = ibm[ibi].AOAAOA_mean[l]*fac-ibm[ibi].AOA_mean[l]*ibm[ibi].AOA_mean[l]*fac*fac;

                
                if (ti == (ti/tiout) * tiout) {
                    if (!rank) PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", r, Ft, Fa, AOA, ibm[ibi].chord_blade[l], ibm[ibi].Urelmag_mean[l]*fac, FFt, FFa, AOAAOA, ibm[ibi].Urel_mean[l].x*fac, ibm[ibi].Urel_mean[l].y*fac,  ibm[ibi].Urel_mean[l].z*fac,  ibm[ibi].Uinduced_mean[l].x*fac, ibm[ibi].Uinduced_mean[l].y*fac, ibm[ibi].Uinduced_mean[l].z*fac, ibm[ibi].circulation_mean[l].x*fac, ibm[ibi].circulation_mean[l].y*fac, ibm[ibi].circulation_mean[l].z*fac,  ibm[ibi].liftdirection_mean[l].x*fac, ibm[ibi].liftdirection_mean[l].y*fac,  ibm[ibi].liftdirection_mean[l].z*fac, ibm[ibi].CL[l], ibm[ibi].CD[l]);
                }
            }
            
            if (ti == (ti/tiout) * tiout) {
                if (!rank) fclose(f);
            }
        }
    }

    if (ti==tistart+1 || ti==tistart || ti == (ti) ) {
        PetscPrintf(PETSC_COMM_WORLD, "Export ForceAOA 2 \n");
        FILE *f;
        char filen[80];

        for (ibi=0; ibi<NumberOfObjects; ibi++) {
            if (!rank) {
                sprintf(filen, "./ForceAOA/ForceAOA_%06d_%2.2d", ti, ibi);
                f = fopen(filen, "w");

                PetscFPrintf(PETSC_COMM_WORLD, f, "Variables = r Ft Fa AOA C Urel Urel.x Urel.y Urel.z Uinduced.x Uinduced.y Uinduced.z circulation.x circulation.y circulation.z liftdirection.x liftdirection.y liftdirection.z Fx Fy Fz CL CD\n");
            }

            for (l=0; l<ibm[ibi].n_elmt; l++) {
                nv1 = ibm[ibi].nv1[l];
                nv2 = ibm[ibi].nv2[l];
                nv3 = ibm[ibi].nv3[l];


                rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
                ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
                rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;

                double r = sqrt(rx*rx+ry*ry+rz*rz);
                n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;

                double ux=0.5*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x); 
                double uy=0.5*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y); 
                double uz=0.5*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z); 

                double Ublade=sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));

                n_rot.x=ux/(Ublade+1.e-20);
                n_rot.y=uy/(Ublade+1.e-20);
                n_rot.z=uz/(Ublade+1.e-20);

                
                double Ft = -(ibm[ibi].F_lagr_x[l]*n_rot.x+ibm[ibi].F_lagr_y[l]*n_rot.y+ibm[ibi].F_lagr_z[l]*n_rot.z);
                double Fa = -(ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb);
                
                double AOA = ibm[ibi].angle_attack[l];

                
                if (!rank) PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", r, Ft, Fa, AOA, ibm[ibi].chord_blade[l], ibm[ibi].Urelmag[l], ibm[ibi].Urel[l].x, ibm[ibi].Urel[l].y,  ibm[ibi].Urel[l].z,  ibm[ibi].Uinduced[l].x, ibm[ibi].Uinduced[l].y, ibm[ibi].Uinduced[l].z, ibm[ibi].circulation[l].x, ibm[ibi].circulation[l].y, ibm[ibi].circulation[l].z,  ibm[ibi].liftdirection[l].x, ibm[ibi].liftdirection[l].y,  ibm[ibi].liftdirection[l].z, ibm[ibi].F_lagr_x[l], ibm[ibi].F_lagr_y[l], ibm[ibi].F_lagr_z[l], ibm[ibi].CL[l], ibm[ibi].CD[l]);
            }

            if (!rank) fclose(f);

        }

    }

    return(0);
}


// calculate force at largrangian points using actuator line model
//
// called in solvers.c before solving the momentum equation
PetscErrorCode Calc_F_lagr_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    PetscInt      	l, ibi, j;
    double		pi = 3.141592653589793;
    PetscReal	A_sum, U_ref, rr, rx, ry, rz, tmp, u_tangent;
    Cmpnts	      	U_b, n_relV, n_L, n_blade, n_rot;	
    int		nv1, nv2, nv3, ifoil;
    PetscReal 	fac1, fac2, r;

    int istall;

    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


    double sign_L;

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        int l_nearhub;
        for (l=1; l<ibm[ibi].n_elmt; l++) {

            rx = ibm[ibi].cent_x[l-1]-fsi[ibi].x_c;
            ry = ibm[ibi].cent_y[l-1]-fsi[ibi].y_c;
            rz = ibm[ibi].cent_z[l-1]-fsi[ibi].z_c;

            double r1 = sqrt(rx*rx+ry*ry+rz*rz)+1.e-20;

            rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
            ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
            rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;

            double r2 = sqrt(rx*rx+ry*ry+rz*rz)+1.e-20;

            //PetscPrintf(PETSC_COMM_WORLD, " %f %f %i %i\n", r1-ibm[ibi].r_nearhubinflowcorrection/reflength_wt, r2-ibm[ibi].r_nearhubinflowcorrection/reflength_wt, l, (r1-ibm[ibi].r_nearhubinflowcorrection/reflength_wt < 1.e-20 && r2-ibm[ibi].r_nearhubinflowcorrection/reflength_wt > 1.e-20) );
            if (r1-ibm[ibi].r_nearhubinflowcorrection/reflength_wt < 1.e-20 || r2-ibm[ibi].r_nearhubinflowcorrection/reflength_wt > 1.e-20){
                l_nearhub=l;
                // PetscPrintf(PETSC_COMM_WORLD, "LNEARHUB %i\n",l);
            }

            if (r1-ibm[ibi].r_nearhubinflowcorrection/reflength_wt < 1.e-20 && r2-ibm[ibi].r_nearhubinflowcorrection/reflength_wt > 1.e-20) l_nearhub=l;

        }

        double dmax_AOA=10000;
        int icount=0; 

        Cmpnts coor_refblade[ibm[ibi].n_elmt];

        Cmpnts p, q, nr, na, nt; 

        int i;	
        for (i=0; i<ibm[ibi].n_elmt; i++) {
        
            p.x=ibm[ibi].cent_x[i]-fsi[ibi].x_c;
            p.y=ibm[ibi].cent_y[i]-fsi[ibi].y_c;
            p.z=ibm[ibi].cent_z[i]-fsi[ibi].z_c;

            na.x=fsi[ibi].nx_tb;	
            na.y=fsi[ibi].ny_tb;	
            na.z=fsi[ibi].nz_tb;	

            double isign = (1.e-19+fsi[ibi].angvel_axis)/(fabs(fsi[ibi].angvel_axis)+1.e-19);

            if (isign>0.0) isign = 1.0;
            else isign = -1.0;

            double theta = isign * refangle_AL*M_PI/180.0;

            if (rotor_model==5) theta=-theta;

            q=ArbitraryRotate(p,theta,na);

            coor_refblade[i].x=q.x+fsi[ibi].x_c;
            coor_refblade[i].y=q.y+fsi[ibi].y_c;
            coor_refblade[i].z=q.z+fsi[ibi].z_c;

        }
        
        double AOA_old[ibm[ibi].n_elmt];

        for (l=0; l<ibm[ibi].n_elmt; l++) {
            ibm[ibi].Uinduced[l].x=0.0;
            ibm[ibi].Uinduced[l].y=0.0;
            ibm[ibi].Uinduced[l].z=0.0;

            AOA_old[l]=1000;
        }

        do {

            icount++; 
            dmax_AOA=0;
            for (l=0; l<ibm[ibi].n_elmt; l++) {
                dmax_AOA+=fabs(AOA_old[l]-ibm[ibi].angle_attack[l])/ibm[ibi].n_elmt;
            }

            PetscBarrier(PETSC_NULL);
            for (l=0; l<ibm[ibi].n_elmt; l++) {
                AOA_old[l]=ibm[ibi].angle_attack[l];
            }

            int iblade;
            int NelmtPerBlade = (int)((double)(ibm[ibi].n_elmt/ibm[ibi].num_blade)+0.1);
            //PetscPrintf(PETSC_COMM_WORLD, "NelmtPerBlade %d %d %d\n", NelmtPerBlade, ibm[ibi].n_elmt, ibm[ibi].num_blade);
            for (iblade=0; iblade<ibm[ibi].num_blade; iblade++){ 
                for (l=iblade*NelmtPerBlade; l<(iblade+1)*NelmtPerBlade; l++) {

                    nv1 = ibm[ibi].nv1[l];
                    nv2 = ibm[ibi].nv2[l];
                    nv3 = ibm[ibi].nv3[l];

                    rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
                    ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
                    rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;

                    double r = sqrt(rx*rx+ry*ry+rz*rz)+1.e-20;
                    n_blade.x = rx/r; n_blade.y = ry/r; n_blade.z = rz/r;

                    double ux, uy, uz;

                    if (rotor_model == 3 || rotor_model == 5 || rotor_model == 6) {
                        ux=0.5*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x); 
                        uy=0.5*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y); 
                        uz=0.5*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z); 
                    } else if (rotor_model == 2) {
                        double fac=1.0/3.0;
                        ux=fac*(ibm[ibi].u[nv1].x+ibm[ibi].u[nv2].x+ibm[ibi].u[nv3].x); 
                        uy=fac*(ibm[ibi].u[nv1].y+ibm[ibi].u[nv2].y+ibm[ibi].u[nv3].y); 
                        uz=fac*(ibm[ibi].u[nv1].z+ibm[ibi].u[nv2].z+ibm[ibi].u[nv3].z); 
                    }

                    double Ublade=sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));
                    //PetscPrintf(PETSC_COMM_WORLD, "Ublade %f \n", Ublade);
                    n_rot.x=ux/(Ublade+1.e-20);
                    n_rot.y=uy/(Ublade+1.e-20);
                    n_rot.z=uz/(Ublade+1.e-20);

                    ibm[ibi].rotationdirection[l].x=n_rot.x;
                    ibm[ibi].rotationdirection[l].y=n_rot.y;
                    ibm[ibi].rotationdirection[l].z=n_rot.z;

                    ibm[ibi].Urel[l].x=ibm[ibi].U_lagr_x[l]-ux;
                    ibm[ibi].Urel[l].y=ibm[ibi].U_lagr_y[l]-uy;
                    ibm[ibi].Urel[l].z=ibm[ibi].U_lagr_z[l]-uz;

                    PetscBarrier(PETSC_NULL);
                    if (r<ibm[ibi].r_nearhubinflowcorrection/reflength_wt) {
                        ibm[ibi].Urel[l].x=ibm[ibi].U_lagr_x[l_nearhub]-ux;
                        ibm[ibi].Urel[l].y=ibm[ibi].U_lagr_y[l_nearhub]-uy;
                        ibm[ibi].Urel[l].z=ibm[ibi].U_lagr_z[l_nearhub]-uz;
                    } else {
                        ibm[ibi].Urel[l].x=ibm[ibi].U_lagr_x[l]-ux;
                        ibm[ibi].Urel[l].y=ibm[ibi].U_lagr_y[l]-uy;
                        ibm[ibi].Urel[l].z=ibm[ibi].U_lagr_z[l]-uz;
                    }

                    ibm[ibi].Urel[l].x=ibm[ibi].Urel[l].x-ibm[ibi].Uinduced[l].x;
                    ibm[ibi].Urel[l].y=ibm[ibi].Urel[l].y-ibm[ibi].Uinduced[l].y;
                    ibm[ibi].Urel[l].z=ibm[ibi].Urel[l].z-ibm[ibi].Uinduced[l].z;

                    double Ur=ibm[ibi].Urel[l].x*n_blade.x+ibm[ibi].Urel[l].y*n_blade.y+ibm[ibi].Urel[l].z*n_blade.z;

                    ibm[ibi].Urel[l].x-=Ur*n_blade.x;
                    ibm[ibi].Urel[l].y-=Ur*n_blade.y;
                    ibm[ibi].Urel[l].z-=Ur*n_blade.z;

                    U_ref=sqrt(pow(ibm[ibi].Urel[l].x,2)+pow(ibm[ibi].Urel[l].y,2)+pow(ibm[ibi].Urel[l].z,2));

                    n_relV.x = ibm[ibi].Urel[l].x/(U_ref+1.e-20); n_relV.y = ibm[ibi].Urel[l].y/(U_ref+1.e-20); n_relV.z = ibm[ibi].Urel[l].z/(U_ref+1.e-20); 
                    tmp = -n_relV.x*n_rot.x-n_relV.y*n_rot.y-n_relV.z*n_rot.z;
                    tmp /= (1+1e-9);

                    ibm[ibi].angle_attack[l] = acos(tmp) * 180.0 / pi - ibm[ibi].angle_twist[l] - ibm[ibi].pitch[iblade];

                    if (turbinestructuremodel && torsion_turbinestructure) ibm[ibi].angle_attack[l] += - ibm[ibi].disp_theta[l]*180/pi; // not sure about the sign

                    ibm[ibi].liftdirection[l].x = n_relV.y*n_blade.z - n_relV.z*n_blade.y; 
                    ibm[ibi].liftdirection[l].y = n_relV.z*n_blade.x - n_relV.x*n_blade.z; 
                    ibm[ibi].liftdirection[l].z = n_relV.x*n_blade.y - n_relV.y*n_blade.x;
        
                    tmp = fsi[ibi].nx_tb*ibm[ibi].liftdirection[l].x + fsi[ibi].ny_tb*ibm[ibi].liftdirection[l].y + fsi[ibi].nz_tb*ibm[ibi].liftdirection[l].z;

                    if (tmp < 0.0) ibm[ibi].liftdirection[l].x = -ibm[ibi].liftdirection[l].x, ibm[ibi].liftdirection[l].y = -ibm[ibi].liftdirection[l].y, ibm[ibi].liftdirection[l].z = -ibm[ibi].liftdirection[l].z;   

                    for(ifoil=0; ifoil<ibm[ibi].num_foiltype; ifoil++) {
                
                        if( r>=ibm[ibi].acl[ifoil].r_beg && r<=ibm[ibi].acl[ifoil].r_end ) {

                            int inrange=0;

                            for(j=0; j<ibm[ibi].acl[ifoil].num_Aerodym-1; j++) {
                                if (ibm[ibi].angle_attack[l]>=ibm[ibi].acl[ifoil].ang_Aerodym[j] && ibm[ibi].angle_attack[l]<=ibm[ibi].acl[ifoil].ang_Aerodym[j+1]) {
                                    fac1 = (ibm[ibi].acl[ifoil].ang_Aerodym[j+1]-ibm[ibi].angle_attack[l])/(ibm[ibi].acl[ifoil].ang_Aerodym[j+1]-ibm[ibi].acl[ifoil].ang_Aerodym[j]);
                                    fac2 = (-ibm[ibi].acl[ifoil].ang_Aerodym[j]+ibm[ibi].angle_attack[l])/(ibm[ibi].acl[ifoil].ang_Aerodym[j+1]-ibm[ibi].acl[ifoil].ang_Aerodym[j]);
                                    ibm[ibi].CD[l] = fac1*ibm[ibi].acl[ifoil].CDInp[j] + fac2*ibm[ibi].acl[ifoil].CDInp[j+1];
                                    inrange=1;
                                }
                            }

                            if (!inrange) {
                                if (ibm[ibi].angle_attack[l]<ibm[ibi].acl[ifoil].ang_Aerodym[0] && ibm[ibi].angle_attack[l]>-45.0) ibm[ibi].CD[l] = (-1.2 - ibm[ibi].acl[ifoil].CDInp[0]) * (ibm[ibi].angle_attack[l] - ibm[ibi].acl[ifoil].ang_Aerodym[0]) / (-45.0-ibm[ibi].acl[ifoil].ang_Aerodym[0]) + ibm[ibi].acl[ifoil].CDInp[0];
                                if ( ibm[ibi].angle_attack[l]<=-45.0 && ibm[ibi].angle_attack[l]>=-90.0 ) ibm[ibi].CD[l] = -1.2;
                                if ( ibm[ibi].angle_attack[l] < -90.0 ) ibm[ibi].CD[l] = -1.2;

                                if (ibm[ibi].angle_attack[l]>ibm[ibi].acl[ifoil].ang_Aerodym[ibm[ibi].acl[ifoil].num_Aerodym-1] && ibm[ibi].angle_attack[l]<45.0) ibm[ibi].CD[l] = (1.2 - ibm[ibi].acl[ifoil].CDInp[ibm[ibi].acl[ifoil].num_Aerodym-1]) * (ibm[ibi].angle_attack[l] - ibm[ibi].acl[ifoil].ang_Aerodym[ibm[ibi].acl[ifoil].num_Aerodym-1]) / (45.0-ibm[ibi].acl[ifoil].ang_Aerodym[ibm[ibi].acl[ifoil].num_Aerodym-1]) + ibm[ibi].acl[ifoil].CDInp[ibm[ibi].acl[ifoil].num_Aerodym-1];
                                if ( ibm[ibi].angle_attack[l]>=45.0 && ibm[ibi].angle_attack[l]<=90.0 ) ibm[ibi].CD[l] = 1.2;
                            }

                            inrange=0;
                            for(j=0; j<ibm[ibi].acl[ifoil].num_Aerodym-1; j++) {
                                if (ibm[ibi].angle_attack[l]>=ibm[ibi].acl[ifoil].ang_Aerodym[j] && ibm[ibi].angle_attack[l]<=ibm[ibi].acl[ifoil].ang_Aerodym[j+1]) {
                                    fac1 = (ibm[ibi].acl[ifoil].ang_Aerodym[j+1]-ibm[ibi].angle_attack[l])/(ibm[ibi].acl[ifoil].ang_Aerodym[j+1]-ibm[ibi].acl[ifoil].ang_Aerodym[j]);
                                    fac2 = (-ibm[ibi].acl[ifoil].ang_Aerodym[j]+ibm[ibi].angle_attack[l])/(ibm[ibi].acl[ifoil].ang_Aerodym[j+1]-ibm[ibi].acl[ifoil].ang_Aerodym[j]);
                                    ibm[ibi].CL[l] = fac1*ibm[ibi].acl[ifoil].CLInp[j] + fac2*ibm[ibi].acl[ifoil].CLInp[j+1];
                                    inrange=1;
                                }
                            }

                            if (!inrange) {

                                if (ibm[ibi].angle_attack[l]<ibm[ibi].acl[ifoil].ang_Aerodym[0] && ibm[ibi].angle_attack[l]>-45.0) ibm[ibi].CL[l] = (-1.05 - ibm[ibi].acl[ifoil].CLInp[0] ) * (ibm[ibi].angle_attack[l] - ibm[ibi].acl[ifoil].ang_Aerodym[0] ) / (-45.0-ibm[ibi].acl[ifoil].ang_Aerodym[0] ) + ibm[ibi].acl[ifoil].CLInp[0] ;
                                if ( ibm[ibi].angle_attack[l]<=-45.0 && ibm[ibi].angle_attack[l]>=-90.0 ) ibm[ibi].CL[l] = 1.05*sin(2.0*ibm[ibi].angle_attack[l]*pi/180.0);
                                if ( ibm[ibi].angle_attack[l] < -90.0 ) ibm[ibi].CL[l] = 0.0;
                                if (ibm[ibi].angle_attack[l]>ibm[ibi].acl[ifoil].ang_Aerodym[ibm[ibi].acl[ifoil].num_Aerodym-1] && ibm[ibi].angle_attack[l]<45.0) ibm[ibi].CL[l] = (1.05 - ibm[ibi].acl[ifoil].CLInp[ibm[ibi].acl[ifoil].num_Aerodym-1] ) * (ibm[ibi].angle_attack[l] - ibm[ibi].acl[ifoil].ang_Aerodym[ibm[ibi].acl[ifoil].num_Aerodym-1] ) / (45.0-ibm[ibi].acl[ifoil].ang_Aerodym[ibm[ibi].acl[ifoil].num_Aerodym-1] ) + ibm[ibi].acl[ifoil].CLInp[ibm[ibi].acl[ifoil].num_Aerodym-1] ;
                                if ( ibm[ibi].angle_attack[l]>=45.0 && ibm[ibi].angle_attack[l]<=90.0 ) ibm[ibi].CL[l] = 1.05*sin(2.0*ibm[ibi].angle_attack[l]*pi/180.0);
                            }
                        }
                    }

                    //PetscPrintf(PETSC_COMM_WORLD, "AL: 3D_CH\n");
                    // add 3D and rotational effects  // Chariaropoulos PK and Hansen MOL, JFE 2000:122:330-6
                    if (correction3D_CH) {
        
                        double  coef_cr1 = c1_CH, coef_cr2 = c2_CH, coef_cr3 = c3_CH;
                        double angle_twist = ibm[ibi].angle_twist[l]*pi/180;
                        double angle_AOA = ibm[ibi].angle_attack[l]*pi/180;
                        double angle_inflow = angle_twist + angle_AOA;
                        ibm[ibi].CD[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * pow(cos(angle_twist), coef_cr3) * (ibm[ibi].CD[l] - CD_0);
                        ibm[ibi].CL[l] += coef_cr1 * pow( ibm[ibi].chord_blade[l]/r, coef_cr2 ) * pow(cos(angle_twist), coef_cr3) * (2.0*pi*(ibm[ibi].angle_attack[l]-AOA_0)/180.0 - ibm[ibi].CL[l]);
                    } 

                    // add 3D and rotational effects  // Du, Zhaohui, and Michael S. Selig. "A 3-D stall-delay model for horizontal axis wind turbine
                    // performance prediction." AIAA Paper 21 (1998).

                    if (correction3D_DS) {
                        double R=fsi[ibi].r_rotor/reflength_wt;
                        double tipspeed = fsi[ibi].angvel_axis*R;
                        double Lambda = tipspeed/sqrt(pow(ibm[ibi].U_ref,2)+pow(tipspeed,2)); 
                        double f1=1.6*(ibm[ibi].chord_blade[l]/r);
                        double f2=pow(ibm[ibi].chord_blade[l]/r,d_DS*R/(fabs(Lambda)*r+1.e-19));
                        double fL=0.5*((f1*a_DS-f2)/(0.1267*b_DS+f2)-1)/pi;
        
                        f2=pow(ibm[ibi].chord_blade[l]/r,0.5*d_DS*R/(fabs(Lambda)*r+1.e-19));
                        double fD=0.5*((f1*a_DS-f2)/(0.1267*b_DS+f2)-1)/pi;
        
                        ibm[ibi].CD[l] -= fD *( ibm[ibi].CD[l] - CD_0);
                        ibm[ibi].CL[l] += fL *( 2.0*pi*(ibm[ibi].angle_attack[l]-AOA_0)/180.0 - ibm[ibi].CL[l]);
                    } 

                    // Modification to 2D force
                    double F1=1.0;
                    double F1_Fa=1.0;
                    // Shen Model 
                    if (Shen_AL) {
                        double R=fsi[ibi].r_rotor/reflength_wt;
                        double pi = 3.1415926;
                        double phi1 = (ibm[ibi].angle_attack[l] + ibm[ibi].angle_twist[l] + ibm[ibi].pitch[iblade])*pi/180.0;
                        if (turbinestructuremodel && torsion_turbinestructure) phi1 += ibm[ibi].disp_theta[l]; 
                        double Coef_a = a_shen; //0.125;
                        double Coef_b = b_shen; //21;
                        double Coef_c = c_shen; //0.1;
                        double gg = -Coef_a*(ibm[ibi].num_blade*fabs(ibm[ibi].Tipspeedratio)-Coef_b);
                        double g1 = exp(gg) + Coef_c; 
                        //g1=g1_AL;
                        F1 = 2.0*acos(exp(-g1*ibm[ibi].num_blade*(max(R-r,0.0))/2.0/r/fabs(sin(phi1))))/pi;
                        F1_Fa = F1*tipcorrectionratio_Fa;
                    }
                    
                    // Prandtl model 
                    if (Prandtl_AL) {

                        double R=fsi[ibi].r_rotor/reflength_wt;
                        double Phi = (ibm[ibi].angle_attack[l] + ibm[ibi].angle_twist[l] + ibm[ibi].pitch[iblade])*pi/180.0;
                        if (turbinestructuremodel && torsion_turbinestructure) Phi += ibm[ibi].disp_theta[l]; 
                        double f = 0.5*(double)ibm[ibi].num_blade*(R-r)/(r*sin(Phi)+1.e-19);
                        F1 = 2.0*acos(exp(-f))/M_PI;
                        F1_Fa = F1;
                    }

                    // Shen1_AL is calibrated for 2.5 MW clipper turbine 
                    if (Shen1_AL) {
                        double R=fsi[ibi].r_rotor/reflength_wt;
                        double pi = 3.1415926;
                        double phi1 = (ibm[ibi].angle_attack[l] + ibm[ibi].angle_twist[l] + ibm[ibi].pitch[iblade])*pi/180.0;
                        if (turbinestructuremodel && torsion_turbinestructure) phi1 += ibm[ibi].disp_theta[l]; 
                        double Lambda = fabs(ibm[ibi].Tipspeedratio);
                        double g1;

                        if (Lambda<7.0636) {
                            g1 = -1.7636*Lambda+11.6112;
                        } else if (Lambda>7.0636 && Lambda<9.8146) {
                            g1 = -0.5805*Lambda+3.2542;
                        } else {
                            g1 = -0.5077*Lambda+2.5397;
                        }
                    
                        g1*=correction_ALShen1;
                        double _g1 = exp(g1);
                        F1 = 2.0*acos(exp(-_g1*ibm[ibi].num_blade*(max(R-r,0.0))/2.0/r/fabs(sin(phi1))))/pi;
                        F1_Fa = F1*tipcorrectionratio_Fa;
                    }

                    double factor;

                    if (rotor_model==2) {
                        factor=1.0;
                    } else if (rotor_model == 3 || rotor_model == 5 || rotor_model == 6) {
                        factor=ibm[ibi].chord_blade[l];
                    }

                    ibm[ibi].Urelmag[l] = U_ref;

                    if (specifycirculation==1) {
                        ibm[ibi].F_lagr_x[l] = -U_ref * ibm[ibi].circulationSpecified[l] * ibm[ibi].liftdirection[l].x;
                        ibm[ibi].F_lagr_y[l] = -U_ref * ibm[ibi].circulationSpecified[l] * ibm[ibi].liftdirection[l].y;
                        ibm[ibi].F_lagr_z[l] = -U_ref * ibm[ibi].circulationSpecified[l] * ibm[ibi].liftdirection[l].z;
                    } 
                    else if (specifycirculation==2) {
                        ibm[ibi].F_lagr_x[l] = -U_ref * ibm[ibi].circulationSpecified[l] * ibm[ibi].liftdirection[l].x;
                        ibm[ibi].F_lagr_y[l] = -U_ref * ibm[ibi].circulationSpecified[l] * ibm[ibi].liftdirection[l].y;
                        ibm[ibi].F_lagr_z[l] = -U_ref * ibm[ibi].circulationSpecified[l] * ibm[ibi].liftdirection[l].z;

                        // force from C_D
                        ibm[ibi].F_lagr_x[l] += -0.5 * U_ref * U_ref * (ibm[ibi].CD[l]) * n_relV.x * factor * F1;
                        ibm[ibi].F_lagr_y[l] += -0.5 * U_ref * U_ref * (ibm[ibi].CD[l]) * n_relV.y * factor * F1;
                        ibm[ibi].F_lagr_z[l] += -0.5 * U_ref * U_ref * (ibm[ibi].CD[l]) * n_relV.z * factor * F1;
                    }
                    else{

                        // force from C_D
                        ibm[ibi].F_lagr_x[l] = -0.5 * U_ref * U_ref * (ibm[ibi].CD[l]) * n_relV.x * factor;
                        ibm[ibi].F_lagr_y[l] = -0.5 * U_ref * U_ref * (ibm[ibi].CD[l]) * n_relV.y * factor;
                        ibm[ibi].F_lagr_z[l] = -0.5 * U_ref * U_ref * (ibm[ibi].CD[l]) * n_relV.z * factor; 
            
                        // add force from C_L
                        ibm[ibi].F_lagr_x[l] += -0.5 * U_ref * U_ref * (ibm[ibi].CL[l]) * ibm[ibi].liftdirection[l].x * factor;
                        ibm[ibi].F_lagr_y[l] += -0.5 * U_ref * U_ref * (ibm[ibi].CL[l]) * ibm[ibi].liftdirection[l].y * factor;
                        ibm[ibi].F_lagr_z[l] += -0.5 * U_ref * U_ref * (ibm[ibi].CL[l]) * ibm[ibi].liftdirection[l].z * factor;

                        double Fa = (ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb)*F1_Fa;
                        double Ft = (ibm[ibi].F_lagr_x[l]*n_rot.x+ibm[ibi].F_lagr_y[l]*n_rot.y+ibm[ibi].F_lagr_z[l]*n_rot.z)*F1;

                        ibm[ibi].F_lagr_x[l] = Fa*fsi[ibi].nx_tb + Ft*n_rot.x;
                        ibm[ibi].F_lagr_y[l] = Fa*fsi[ibi].ny_tb + Ft*n_rot.y;
                        ibm[ibi].F_lagr_z[l] = Fa*fsi[ibi].nz_tb + Ft*n_rot.z;

                        // PetscPrintf(PETSC_COMM_WORLD,"Turb_Force %le %le %le %i %i\n",ibm[ibi].F_lagr_z[l],U_ref,factor,ibi,l);

                        if (correction3D_CL) {
                            if (r<0.75*fsi[ibi].r_rotor/reflength_wt) {
                                double angle_twist = ibm[ibi].angle_twist[l]*pi/180;
                                double angle_AOA = ibm[ibi].angle_attack[l]*pi/180;
                                double angle_inflow = angle_twist + angle_AOA;

                                double nx=fsi[ibi].nx_tb, ny=fsi[ibi].ny_tb, nz=fsi[ibi].nz_tb;
                                double F_correc = -0.5*c0_CL * pow( ibm[ibi].chord_blade[l]/r, 1) * pow(cos(angle_inflow), 2) * U_ref * U_ref * factor;
                
                                ibm[ibi].F_lagr_x[l]+=F_correc*nx;	
                                ibm[ibi].F_lagr_y[l]+=F_correc*ny;	
                                ibm[ibi].F_lagr_z[l]+=F_correc*nz;	

                                //PetscPrintf(PETSC_COMM_WORLD,"Turb_Force_Correction %le %le %le %i %i\n",ibm[ibi].F_lagr_z[l],U_ref,factor,ibi,l);
                            }
                        }
                    }
                }

                if (smoothforce_AL) {
                    int i,j,k;
                    int Elmt_blade=ibm[ibi].n_elmt/ibm[ibi].num_blade;
                    double tmp_forcex[ibm[ibi].n_elmt];
                    double tmp_forcey[ibm[ibi].n_elmt];
                    double tmp_forcez[ibm[ibi].n_elmt];
        
                    int ii=0;
                    do {
                        ii++;
                        for (l=0; l<ibm[ibi].n_elmt; l++) {
                            tmp_forcex[l]=ibm[ibi].F_lagr_x[l];
                            tmp_forcey[l]=ibm[ibi].F_lagr_y[l];
                            tmp_forcez[l]=ibm[ibi].F_lagr_z[l];
                        }

                        for (k=0; k<ibm[ibi].num_blade; k++) {
                            int j1=k*Elmt_blade+2;
                            int j2=(k+1)*Elmt_blade-2;
                            for (j=j1; j<j2; j++) {
                                double sum_forcex=0.0;
                                double sum_forcey=0.0;
                                double sum_forcez=0.0;
                                for (i=j-2; i<j+2; i++) {			
                                    double weight=(double)(i-j);
                                    sum_forcex+=tmp_forcex[i]*dfunc_4h(weight);
                                    sum_forcey+=tmp_forcey[i]*dfunc_4h(weight);
                                    sum_forcez+=tmp_forcez[i]*dfunc_4h(weight);
                                }
                                ibm[ibi].F_lagr_x[j]=sum_forcex;
                                ibm[ibi].F_lagr_y[j]=sum_forcey;
                                ibm[ibi].F_lagr_z[j]=sum_forcez;
                            }	
                        }

                        for (l=0; l<ibm[ibi].n_elmt; l++) {
                            tmp_forcex[l]=ibm[ibi].F_lagr_x[l];
                            tmp_forcey[l]=ibm[ibi].F_lagr_y[l];
                            tmp_forcez[l]=ibm[ibi].F_lagr_z[l];
                        }

                        for (k=0; k<ibm[ibi].num_blade; k++) {
                            int j1=k*Elmt_blade+1;
                            int j2=(k+1)*Elmt_blade-1;
                            for (j=j1; j<j2; j++) {
                                double sum_forcex=0.0;
                                double sum_forcey=0.0;
                                double sum_forcez=0.0;
                                for (i=j-1; i<j+1; i++) {			
                                    double weight=(double)(i-j);
                                    sum_forcex+=tmp_forcex[i]*dfunc_2h(weight);
                                    sum_forcey+=tmp_forcey[i]*dfunc_2h(weight);
                                    sum_forcez+=tmp_forcez[i]*dfunc_2h(weight);
                                }
                                ibm[ibi].F_lagr_x[j]=sum_forcex;
                                ibm[ibi].F_lagr_y[j]=sum_forcey;
                                ibm[ibi].F_lagr_z[j]=sum_forcez;
                            }	
                        }

                        for (k=0; k<ibm[ibi].num_blade; k++) {
                            int jj=k*Elmt_blade; 
                            int j1=k*Elmt_blade+1;
                            int j2=k*Elmt_blade+2;
                            ibm[ibi].F_lagr_x[jj]=2*ibm[ibi].F_lagr_x[j1]-ibm[ibi].F_lagr_x[j2];
                            ibm[ibi].F_lagr_y[jj]=2*ibm[ibi].F_lagr_y[j1]-ibm[ibi].F_lagr_y[j2];
                            ibm[ibi].F_lagr_z[jj]=2*ibm[ibi].F_lagr_z[j1]-ibm[ibi].F_lagr_z[j2];

                            jj=(k+1)*Elmt_blade-1; 
                            j1=(k+1)*Elmt_blade-2;
                            j2=(k+1)*Elmt_blade-3;
                            ibm[ibi].F_lagr_x[jj]=2*ibm[ibi].F_lagr_x[j1]-ibm[ibi].F_lagr_x[j2];
                            ibm[ibi].F_lagr_y[jj]=2*ibm[ibi].F_lagr_y[j1]-ibm[ibi].F_lagr_y[j2];
                            ibm[ibi].F_lagr_z[jj]=2*ibm[ibi].F_lagr_z[j1]-ibm[ibi].F_lagr_z[j2];
                        }	
                    } while(ii<=smoothforce_AL);
                }

                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    r = sqrt( pow((ibm[ibi].cent_x[l]-fsi[ibi].x_c),2) + pow((ibm[ibi].cent_y[l]-fsi[ibi].y_c),2) + pow((ibm[ibi].cent_z[l]-fsi[ibi].z_c),2));
                    if (r<ibm[ibi].r_nacelle/reflength_wt) {
                        ibm[ibi].F_lagr_x[l]=0.0;	
                        ibm[ibi].F_lagr_y[l]=0.0;	
                        ibm[ibi].F_lagr_z[l]=0.0;	
                        ibm[ibi].angle_attack[l]=0.0;
                    }
                }

                // circulation on the blade
                for (l=0; l<ibm[ibi].n_elmt; l++) {

                    r = sqrt( pow((ibm[ibi].cent_x[l]-fsi[ibi].x_c),2) + pow((ibm[ibi].cent_y[l]-fsi[ibi].y_c),2) + pow((ibm[ibi].cent_z[l]-fsi[ibi].z_c),2));
                    Cmpnts Ureldirection, circulationdirection;

                    double UU = sqrt(pow(ibm[ibi].Urel[l].x+ibm[ibi].Uinduced[l].x,2) + pow(ibm[ibi].Urel[l].y+ibm[ibi].Uinduced[l].y,2) + pow(ibm[ibi].Urel[l].z+ibm[ibi].Uinduced[l].z,2))+1.e-19;
                    Ureldirection.x=(ibm[ibi].Urel[l].x+ibm[ibi].Uinduced[l].x)/UU;
                    Ureldirection.y=(ibm[ibi].Urel[l].y+ibm[ibi].Uinduced[l].y)/UU;
                    Ureldirection.z=(ibm[ibi].Urel[l].z+ibm[ibi].Uinduced[l].z)/UU;
                    
                    circulationdirection.x=ibm[ibi].liftdirection[l].y*Ureldirection.z-ibm[ibi].liftdirection[l].z*Ureldirection.y;
                    circulationdirection.y=ibm[ibi].liftdirection[l].z*Ureldirection.x-ibm[ibi].liftdirection[l].x*Ureldirection.z;
                    circulationdirection.z=ibm[ibi].liftdirection[l].x*Ureldirection.y-ibm[ibi].liftdirection[l].y*Ureldirection.x;

                    double Taumag= -(ibm[ibi].F_lagr_x[l]*ibm[ibi].liftdirection[l].x + ibm[ibi].F_lagr_y[l]*ibm[ibi].liftdirection[l].y + ibm[ibi].F_lagr_z[l]*ibm[ibi].liftdirection[l].z)/UU;
                
                    ibm[ibi].circulation[l].x = Taumag*circulationdirection.x;
                    ibm[ibi].circulation[l].y = Taumag*circulationdirection.y;
                    ibm[ibi].circulation[l].z = Taumag*circulationdirection.z;
                }

                // induced velocity
                int ll;

                if (rotor_model==6) {
                    for (l=0; l<ibm[ibi].n_elmt; l++) {

                        ibm[ibi].Uinduced[l].x=0.0;
                        ibm[ibi].Uinduced[l].y=0.0;
                        ibm[ibi].Uinduced[l].z=0.0;

                        for (ll=0; ll<ibm[ibi].n_elmt; ll++) {
                            Cmpnts r;
                            r.x=coor_refblade[l].x-ibm[ibi].cent_x[ll];
                            r.y=coor_refblade[l].y-ibm[ibi].cent_y[ll];
                            r.z=coor_refblade[l].z-ibm[ibi].cent_z[ll];

                            double rr=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);

                            Cmpnts Tau;
                            Tau.x=ibm[ibi].circulation[ll].x*ibm[ibi].dA[ll];
                            Tau.y=ibm[ibi].circulation[ll].y*ibm[ibi].dA[ll];
                            Tau.z=ibm[ibi].circulation[ll].z*ibm[ibi].dA[ll];

                            double fac=1.0/(4*M_PI*rr*rr*rr+1.e-19);

                            if (rr>(ibm[ibi].r_nacelle/reflength_wt)) {
                                ibm[ibi].Uinduced[l].x+=(Tau.y*r.z-Tau.z*r.y)*fac;
                                ibm[ibi].Uinduced[l].y+=(Tau.z*r.x-Tau.x*r.z)*fac;
                                ibm[ibi].Uinduced[l].z+=(Tau.x*r.y-Tau.y*r.x)*fac;
                            }
                        }
                    }
                }

                if (rotor_model==5) {
                    for (l=0; l<ibm[ibi].n_elmt; l++) {

                        ibm[ibi].Uinduced[l].x=0.0;
                        ibm[ibi].Uinduced[l].y=0.0;
                        ibm[ibi].Uinduced[l].z=0.0;

                        for (ll=0; ll<ibm[ibi].n_elmt; ll++) {

                            Cmpnts r;
                            r.x=-coor_refblade[ll].x+ibm[ibi].cent_x[l];
                            r.y=-coor_refblade[ll].y+ibm[ibi].cent_y[l];
                            r.z=-coor_refblade[ll].z+ibm[ibi].cent_z[l];

                            double rr=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);

                            Cmpnts Tau;
                            Tau.x=ibm[ibi].circulation[ll].x*ibm[ibi].dA[ll];
                            Tau.y=ibm[ibi].circulation[ll].y*ibm[ibi].dA[ll];
                            Tau.z=ibm[ibi].circulation[ll].z*ibm[ibi].dA[ll];

                            double fac=1.0/(4*M_PI*rr*rr*rr+1.e-19);

                            if (rr>(ibm[ibi].r_nacelle/reflength_wt)) {
                                ibm[ibi].Uinduced[l].x+=(Tau.y*r.z-Tau.z*r.y)*fac;
                                ibm[ibi].Uinduced[l].y+=(Tau.z*r.x-Tau.x*r.z)*fac;
                                ibm[ibi].Uinduced[l].z+=(Tau.x*r.y-Tau.y*r.x)*fac;
                            }
                        }
                    }
                }
            }
        } while(fabs(dmax_AOA)>0.01 && icount<maxiteraction_rotormodel);

        // compute the bending moment on each blade
        
        int iblade;
        
        int NelmtPerBlade = (int)((double)(ibm[ibi].n_elmt/ibm[ibi].num_blade)+0.1);
        for (iblade=0; iblade<ibm[ibi].num_blade; iblade++) {

            fsi[ibi].Moment_bladebending[iblade] = 0;
            for (l=iblade*NelmtPerBlade; l<(iblade+1)*NelmtPerBlade; l++) {

                nv1 = ibm[ibi].nv1[l];
                nv2 = ibm[ibi].nv2[l];
                nv3 = ibm[ibi].nv3[l];

                rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;
                ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;
                rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;

                double r = sqrt(rx*rx+ry*ry+rz*rz)+1.e-20;

                double Fa = -(ibm[ibi].F_lagr_x[l]*fsi[ibi].nx_tb+ibm[ibi].F_lagr_y[l]*fsi[ibi].ny_tb+ibm[ibi].F_lagr_z[l]*fsi[ibi].nz_tb)*ibm[ibi].dA[l];

                fsi[ibi].Moment_bladebending[iblade] += Fa*r; 
            }
        }
    }

    return(0);
}


/* Interpolate the velocity at the Lagrangian points */
// subroutine for Calc_F_lagr_ACL and Calc_F_lagr
PetscErrorCode Calc_U_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    DA              da = user->da, fda = user->fda;
    DALocalInfo     info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    Cmpnts	***ucat, ***coor, ***csi, ***eta, ***zet;

    PetscReal 	***aj;

    Vec		Coor;

    PetscReal	dfunc;

    PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

    double r1, r2, r3;

    DAGetLocalInfo(da, &info);
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


    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);
    DAVecGetArray(fda, user->lUcat, &ucat);
    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);

    for (ibi=0; ibi<NumberOfObjects; ibi++) {
        for (l=0; l<ibm[ibi].n_elmt; l++) {
            ibm[ibi].U_lagr_x[l] = 0.0;
            ibm[ibi].U_lagr_y[l] = 0.0;
            ibm[ibi].U_lagr_z[l] = 0.0;
        }
    }

    clock_t start, end;
    double elapsed;

    start = clock();
    double ni[3], nj[3], nk[3];

    for (ibi=0; ibi<NumberOfObjects; ibi++) {
        for (l=0; l<ibm[ibi].n_elmt; l++) {
            //PetscPrintf(PETSC_COMM_WORLD, "Debug in Calc_U_lagr for element %i with kmin and kmax = %i, %i\n\n", l, ibm[ibi].k_min[l],ibm[ibi].k_max[l]);

            for (k = ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
            for (j = ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
            for (i = ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++){
                //PetscPrintf(PETSC_COMM_WORLD, "Debug in Calc_U_lagr inner loop with delta function\n");

                xc = 0.125 *
                    (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                    coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
                    coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
                    coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
                yc = 0.125 *
                    (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                    coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
                    coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
                    coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
                zc = 0.125 *
                    (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                    coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
                    coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
                    coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);

                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                double dhx_=1.0/aj[k][j][i]/area;

                area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                double dhy_=1.0/aj[k][j][i]/area;

                area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                double dhz_=1.0/aj[k][j][i]/area;	

                double rx= (xc - ibm[ibi].cent_x[l]), ry = (yc - ibm[ibi].cent_y[l]), rz = (zc - ibm[ibi].cent_z[l]);

                Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                r1=(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
                r2=(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
                r3=(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 

                if (deltafunc == 0) {
                    dfunc = dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);
                } else if (deltafunc == 6) {		
                    double n = halfwidth_dfunc;
                    dfunc = dfunc_nh(r1,n) * dfunc_nh(r2,n) * dfunc_nh(r3,n);
                } else if (deltafunc == 7) {		
                    double n = halfwidth_dfunc;
                    dfunc = dfunc_exp(r1,n) * dfunc_exp(r2,n) * dfunc_exp(r3,n);
                } else if (deltafunc == 8) { 
                    dfunc = dfunc_s4h(r1) * dfunc_s4h(r2) * dfunc_s4h(r3);
                } else { 
                    dfunc = dfunc_sc4h(r1) * dfunc_sc4h(r2) * dfunc_sc4h(r3);//Same as in Yang Wind Energy 2018
                }

                ibm[ibi].U_lagr_x[l] += ucat[k][j][i].x * dfunc;
                ibm[ibi].U_lagr_y[l] += ucat[k][j][i].y * dfunc;
                ibm[ibi].U_lagr_z[l] += ucat[k][j][i].z * dfunc;            
            }
        }
    }

    for (ibi=0; ibi<NumberOfObjects; ibi++){
        std::vector<Cmpnts> u_local (ibm[ibi].n_elmt); 
        std::vector<Cmpnts> u_sum (ibm[ibi].n_elmt); 

        for (i=0; i<ibm[ibi].n_elmt; i++ ) {
            u_local[i].x = ibm[ibi].U_lagr_x[i];
            u_local[i].y = ibm[ibi].U_lagr_y[i];
            u_local[i].z = ibm[ibi].U_lagr_z[i];
        }

        int count_vector=3*ibm[ibi].n_elmt;
        MPI_Allreduce( &(u_local[0]), &(u_sum[0]), count_vector, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

        for (i=0; i<ibm[ibi].n_elmt; i++ ) {
            ibm[ibi].U_lagr_x[i] = u_sum[i].x;
            ibm[ibi].U_lagr_y[i] = u_sum[i].y;
            ibm[ibi].U_lagr_z[i] = u_sum[i].z;
        }
    }

    if (ii_periodicWT || jj_periodicWT || kk_periodicWT) {

        double xc_min = 1.0e6;        
        double yc_min = 1.0e6;
        double zc_min = 1.0e6;

        for (ibi=0; ibi<NumberOfObjects; ibi++) {
            xc_min = PetscMin(xc_min, fsi[ibi].x_c);
            yc_min = PetscMin(yc_min, fsi[ibi].y_c);
            zc_min = PetscMin(zc_min, fsi[ibi].z_c);
        }

        int IndWT[Nz_WT][Ny_WT][Nx_WT];

        double fac_x=1.0/Sx_WT, fac_y=1.0/Sy_WT, fac_z=1.0/Sz_WT;
        for (ibi=0; ibi<NumberOfObjects; ibi++){
            int ii=(int)((fsi[ibi].x_c-xc_min+1.e-9)*fac_x);
            int jj=(int)((fsi[ibi].y_c-yc_min+1.e-9)*fac_y);
            int kk=(int)((fsi[ibi].z_c-zc_min+1.e-9)*fac_z);

            PetscPrintf(PETSC_COMM_WORLD, "ibi ii jj kk %d %d %d %d \n", ibi, ii, jj, kk);
            IndWT[kk][jj][ii]=ibi;
        }

        int i, j, k;

        for (k=0;k<Nz_WT;k++)
        for (j=0;j<Ny_WT;j++)
        for (i=0;i<Nx_WT;i++){

            if (ii_periodicWT && i==0) {
                ibi=IndWT[k][j][i];
                int ibi1=IndWT[k][j][Nx_WT-1];
                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    ibm[ibi].U_lagr_x[l]=ibm[ibi].U_lagr_x[l]+ibm[ibi1].U_lagr_x[l];
                    ibm[ibi].U_lagr_y[l]=ibm[ibi].U_lagr_y[l]+ibm[ibi1].U_lagr_y[l];
                    ibm[ibi].U_lagr_z[l]=ibm[ibi].U_lagr_z[l]+ibm[ibi1].U_lagr_z[l];
                }
            }

            if (jj_periodicWT && j==0) {
                ibi=IndWT[k][j][i];
                int ibi1=IndWT[k][Ny_WT-1][i];
                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    ibm[ibi].U_lagr_x[l]=ibm[ibi].U_lagr_x[l]+ibm[ibi1].U_lagr_x[l];
                    ibm[ibi].U_lagr_y[l]=ibm[ibi].U_lagr_y[l]+ibm[ibi1].U_lagr_y[l];
                    ibm[ibi].U_lagr_z[l]=ibm[ibi].U_lagr_z[l]+ibm[ibi1].U_lagr_z[l];
                }
            }

            if (kk_periodicWT && k==0){
                ibi=IndWT[k][j][i];
                int ibi1=IndWT[Nz_WT-1][j][i];
                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    ibm[ibi].U_lagr_x[l]=ibm[ibi].U_lagr_x[l]+ibm[ibi1].U_lagr_x[l];
                    ibm[ibi].U_lagr_y[l]=ibm[ibi].U_lagr_y[l]+ibm[ibi1].U_lagr_y[l];
                    ibm[ibi].U_lagr_z[l]=ibm[ibi].U_lagr_z[l]+ibm[ibi1].U_lagr_z[l];
                }
            }
        }

        for (k=0;k<Nz_WT;k++)
        for (j=0;j<Ny_WT;j++)
        for (i=0;i<Nx_WT;i++){
            if (ii_periodicWT && i==Nx_WT-1){
                ibi=IndWT[k][j][i];
                int ibi1=IndWT[k][j][0];
                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    ibm[ibi].U_lagr_x[l]=ibm[ibi1].U_lagr_x[l];
                    ibm[ibi].U_lagr_y[l]=ibm[ibi1].U_lagr_y[l];
                    ibm[ibi].U_lagr_z[l]=ibm[ibi1].U_lagr_z[l];
                }
            }

            if (jj_periodicWT && j==Ny_WT-1) {
                ibi=IndWT[k][j][i];
                int ibi1=IndWT[k][0][i];
                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    ibm[ibi].U_lagr_x[l]=ibm[ibi1].U_lagr_x[l];
                    ibm[ibi].U_lagr_y[l]=ibm[ibi1].U_lagr_y[l];
                    ibm[ibi].U_lagr_z[l]=ibm[ibi1].U_lagr_z[l];
                }
            }

            if (kk_periodicWT && k==Nz_WT-1) {
                ibi=IndWT[k][j][i];
                int ibi1=IndWT[0][j][i];
                for (l=0; l<ibm[ibi].n_elmt; l++) {
                    ibm[ibi].U_lagr_x[l]=ibm[ibi1].U_lagr_x[l];
                    ibm[ibi].U_lagr_y[l]=ibm[ibi1].U_lagr_y[l];
                    ibm[ibi].U_lagr_z[l]=ibm[ibi1].U_lagr_z[l];
                }
            }
        }
    }

    if (MoveFrame){
        for (ibi=0; ibi<NumberOfObjects; ibi++)
        for (l=0; l<ibm[ibi].n_elmt; l++){
            ibm[ibi].U_lagr_x[l] += u_frame;
            ibm[ibi].U_lagr_y[l] += v_frame;
            ibm[ibi].U_lagr_z[l] += w_frame;
        }
    }

    DAVecRestoreArray(fda, Coor, &coor);
    DAVecRestoreArray(fda, user->lUcat, &ucat);
    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);

    return(0);
}


// distribute the force on the turbines to background grid.
//
// called after Calc_F_lagr_ACL or Calc_F_lagr before solving the momentum equatiion
PetscErrorCode Calc_F_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df)
{

    //PetscPrintf(PETSC_COMM_WORLD, "Debug1 in Calc_F_eul  in rotor_model.c\n");
    int my_rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);

    DA              da = user->da, fda = user->fda;
    DALocalInfo     info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    Cmpnts		***lf_eul, ***coor, ***csi, ***eta, ***zet;

    PetscReal 	***aj, ***nvert;

    Vec		Coor;

    PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

    double		dfunc;

    double r1, r2, r3;

    DAGetLocalInfo(da, &info);
    mx = info.mx; my = info.my; mz = info.mz;
    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;

    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);
    DAVecGetArray(fda, user->lF_eul, &lf_eul);
    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);

    DAVecGetArray(da, user->lNvert, &nvert);

    //PetscPrintf(PETSC_COMM_WORLD, "Debug2 in Calc_F_eul  in rotor_model.c\n");

    double tempEulx, tempEuly, tempEulz;
    tempEulx = 0.0;
    tempEuly = 0.0;
    tempEulz = 0.0;
    int	tempElement;
    char my_msg[256];
    MPI_Comm  comm;
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    comm = MPI_COMM_WORLD;
    double ni[3], nj[3], nk[3];
    for (ibi=0; ibi<NumberOfObjects; ibi++) {
        //PetscPrintf(PETSC_COMM_WORLD, "Debug3 in Calc_F_eul in rotor_model.c with ibm %i having %i elements\n", ibi, ibm[ibi].n_elmt);
        
        int cnt = 0;
        for (l=0; l<ibm[ibi].n_elmt; l++) {
            if (l==0) {
                //PetscPrintf(PETSC_COMM_SELF, "Debug4-SELF rank %i cnt %i for k, j, i min/max: %i, %i, %i, %i,%i, %i\n",my_rank, cnt, ibm[ibi].k_min[l], 
                //	ibm[ibi].k_max[l],ibm[ibi].j_min[l], ibm[ibi].j_max[l],ibm[ibi].i_min[l], ibm[ibi].i_max[l]);
                //++cnt;
                for(int p=0;p<worldSize;p++){
                    if (my_rank==p){
                        ++cnt;
                        //printf("rank %d, cnt %d  for element %d\n", my_rank, cnt, l);
                    }
    //              MPI_Barrier(comm);
                }
                MPI_Barrier(comm);

            }
            //if(ibm[ibi].k_min[l] != ibm[ibi].k_max[l] && ibm[ibi].j_min[l] != ibm[ibi].j_max[l] && ibm[ibi].i_min[l] != ibm[ibi].i_max[l]) {
                //printf( "Debug4 For rank %i,  k, j, i min/max: %i, %i,   %i, %i,   %i, %i\n",my_rank, ibm[ibi].k_min[l], ibm[ibi].k_max[l],ibm[ibi].j_min[l], ibm[ibi].j_max[l],ibm[ibi].i_min[l], ibm[ibi].i_max[l]);
            //}
            for (k=ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
            for (j=ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
            for (i=ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) {

                double xi = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25;	
                double yi = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25;	
                double zi = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; 

                double xj = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25;	
                double yj = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25;	
                double zj = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25; 

                double xk = (coor[k  ][j  ][i].x + coor[k  ][j-1][i].x + coor[k  ][j][i-1].x + coor[k][j-1][i-1].x) * 0.25;	
                double yk = (coor[k  ][j  ][i].y + coor[k  ][j-1][i].y + coor[k  ][j][i-1].y + coor[k][j-1][i-1].y) * 0.25;	
                double zk = (coor[k  ][j  ][i].z + coor[k  ][j-1][i].z + coor[k  ][j][i-1].z + coor[k][j-1][i-1].z) * 0.25; 
        
                xc=coor[k][j][i].x;
                yc=coor[k][j][i].y;
                zc=coor[k][j][i].z;

                double rxi= (xi - ibm[ibi].cent_x[l]), ryi = (yi - ibm[ibi].cent_y[l]), rzi = (zi - ibm[ibi].cent_z[l]);
                double rxj= (xj - ibm[ibi].cent_x[l]), ryj = (yj - ibm[ibi].cent_y[l]), rzj = (zj - ibm[ibi].cent_z[l]);
                double rxk= (xk - ibm[ibi].cent_x[l]), ryk = (yk - ibm[ibi].cent_y[l]), rzk = (zk - ibm[ibi].cent_z[l]);

                Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
        
                double dhx_, dhy_, dhz_;
                //PetscPrintf(PETSC_COMM_WORLD, "Debug6 in Calc_F_eul  in rotor_model.c\n");

                if (forcewidthfixed) {
                    dhx_=dhi_fixed;
                    dhy_=dhj_fixed;
                    dhz_=dhk_fixed;
                } 
                else if (forcedistrwidth_surfacescale) {
                    if (rotor_model == 5) dhx_=sqrt(ibm[ibi].dA[l]);
                    else dhx_=ibm[ibi].dA[l];

                    dhy_=dhx_;
                    dhz_=dhx_;
                }
                else {
                    double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                    dhx_=1.0/aj[k][j][i]/area;

                    area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                    dhy_=1.0/aj[k][j][i]/area;

                    area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                    dhz_=1.0/aj[k][j][i]/area;	
                }

                double r1i=fabs(rxi*ni[0]+ryi*ni[1]+rzi*ni[2])/dhx_; 
                double r2i=fabs(rxi*nj[0]+ryi*nj[1]+rzi*nj[2])/dhy_; 
                double r3i=fabs(rxi*nk[0]+ryi*nk[1]+rzi*nk[2])/dhz_; 

                double r1j=fabs(rxj*ni[0]+ryj*ni[1]+rzj*ni[2])/dhx_; 
                double r2j=fabs(rxj*nj[0]+ryj*nj[1]+rzj*nj[2])/dhy_; 
                double r3j=fabs(rxj*nk[0]+ryj*nk[1]+rzj*nk[2])/dhz_; 

                double r1k=fabs(rxk*ni[0]+ryk*ni[1]+rzk*ni[2])/dhx_; 
                double r2k=fabs(rxk*nj[0]+ryk*nj[1]+rzk*nj[2])/dhy_; 
                double r3k=fabs(rxk*nk[0]+ryk*nk[1]+rzk*nk[2])/dhz_; 

                vol_eul = 1.0/(dhx_*dhy_*dhz_);
        
                double dfunci, dfuncj, dfunck;
                if (df == 0) {
                    dfunci = vol_eul * dfunc_2h(r1i) * dfunc_2h(r2i) * dfunc_2h(r3i);
                    dfuncj = vol_eul * dfunc_2h(r1j) * dfunc_2h(r2j) * dfunc_2h(r3j);
                    dfunck = vol_eul * dfunc_2h(r1k) * dfunc_2h(r2k) * dfunc_2h(r3k);
                } else if (df == 6) {		
                    double n = halfwidth_dfunc;
                    dfunci = vol_eul * dfunc_nh(r1i,n) * dfunc_nh(r2i,n) * dfunc_nh(r3i,n);
                    dfuncj = vol_eul * dfunc_nh(r1j,n) * dfunc_nh(r2j,n) * dfunc_nh(r3j,n);
                    dfunck = vol_eul * dfunc_nh(r1k,n) * dfunc_nh(r2k,n) * dfunc_nh(r3k,n);

                } else if (df == 7) {		
                    double n = halfwidth_dfunc;
                    dfunci = vol_eul * dfunc_exp(r1i,n) * dfunc_exp(r2i,n) * dfunc_exp(r3i,n);
                    dfuncj = vol_eul * dfunc_exp(r1j,n) * dfunc_exp(r2j,n) * dfunc_exp(r3j,n);
                    dfunck = vol_eul * dfunc_exp(r1k,n) * dfunc_exp(r2k,n) * dfunc_exp(r3k,n);
                } else if (df == 8) { 
                    dfunci = vol_eul * dfunc_s4h(r1i) * dfunc_s4h(r2i) * dfunc_s4h(r3i);
                    dfuncj = vol_eul * dfunc_s4h(r1j) * dfunc_s4h(r2j) * dfunc_s4h(r3j);
                    dfunck = vol_eul * dfunc_s4h(r1k) * dfunc_s4h(r2k) * dfunc_s4h(r3k);
                } else { 
                    dfunci = vol_eul * dfunc_sc4h(r1i) * dfunc_sc4h(r2i) * dfunc_sc4h(r3i);
                    dfuncj = vol_eul * dfunc_sc4h(r1j) * dfunc_sc4h(r2j) * dfunc_sc4h(r3j);
                    dfunck = vol_eul * dfunc_sc4h(r1k) * dfunc_sc4h(r2k) * dfunc_sc4h(r3k);
                }

                if(l==4046 || l==4832 || l==435){
                    tempElement = l;
                    tempEulx =lf_eul[k][j][i].x;
                    tempEuly =lf_eul[k][j][i].y;
                    tempEulz=lf_eul[k][j][i].z;
                }
                 lf_eul[k][j][i].x += ibm[ibi].F_lagr_x[l] * dfunci * ibm[ibi].dA[l] * csi[k][j][i].x +
                    ibm[ibi].F_lagr_y[l] * dfunci * ibm[ibi].dA[l] * csi[k][j][i].y +
                    ibm[ibi].F_lagr_z[l] * dfunci * ibm[ibi].dA[l] * csi[k][j][i].z;

                 lf_eul[k][j][i].y += ibm[ibi].F_lagr_x[l] * dfuncj * ibm[ibi].dA[l] *  eta[k][j][i].x +
                    ibm[ibi].F_lagr_y[l] * dfuncj * ibm[ibi].dA[l] *  eta[k][j][i].y +
                    ibm[ibi].F_lagr_z[l] * dfuncj * ibm[ibi].dA[l] *  eta[k][j][i].z;

                 lf_eul[k][j][i].z += ibm[ibi].F_lagr_x[l] * dfunck * ibm[ibi].dA[l] *  zet[k][j][i].x +
                    ibm[ibi].F_lagr_y[l] * dfunck * ibm[ibi].dA[l] *  zet[k][j][i].y +
                    ibm[ibi].F_lagr_z[l] * dfunck * ibm[ibi].dA[l] *  zet[k][j][i].z;

                if(( k==44 && j==3 && i == 17 && l==4046) || (k==44 && j==6 && i==17 && l==435)  ||(k==44 && j==7 && i==17 && l==4832)){
                    /*PetscPrintf(PETSC_COMM_SELF,"\n*****Debug on rank %i where the LagrangianF for element %i in x,y,z direction = %le, %le, %le  *****\n", 
                        my_rank, l, ibm[ibi].F_lagr_x[l], ibm[ibi].F_lagr_y[l], ibm[ibi].F_lagr_z[l]);

                    PetscPrintf(PETSC_COMM_SELF,"rank =%d,lfx=%le, ilfy=%le, lfz=%le, dfuncj=%le, dA=%le, etax=%le, etay=%le, etaz=%le\n", 
                        my_rank,ibm[ibi].F_lagr_x[l], ibm[ibi].F_lagr_y[l], ibm[ibi].F_lagr_z[l], dfuncj, ibm[ibi].dA[l], eta[k][j][i].x, 
                        eta[k][j][i].y, eta[k][j][i].z);

                    PetscPrintf(PETSC_COMM_SELF,"r1j = %le, r2j=%le, r3j=%le, vol_eul=%le\n", r1j, r2j, r3j, vol_eul);

                    PetscPrintf(PETSC_COMM_SELF,"Rank %i PreEularianF for element %i on grid pt at       k,j,i (%i, %i, %i) in x,y,z = %le, %le, %le\n", 
                        my_rank, l,k,j,i, tempEulx, tempEuly, tempEulz);

                    PetscPrintf(PETSC_COMM_SELF,"Rank %i Modified EularianF for element %i on grid pt at k,j,i (%i, %i, %i) in x,y,z = %le, %le, %le\n", 
                        my_rank, l,k,j,i, lf_eul[k][j][i].x , lf_eul[k][j][i].y , lf_eul[k][j][i].z  );*/
                }
            }
        }
    }

    for (k=lzs;k<lze;k++) 
    for (j=lys;j<lye;j++) 
    for (i=lxs;i<lxe;i++) {
    //	PetscPrintf(PETSC_COMM_WORLD, "Debug7 in Calc_F_eul  in rotor_model.c\n");

        int ii, jj, kk;
        double _nvert;
        _nvert = 0.0;

        for (kk=max(k-1,lzs);kk<min(k+2,lze);kk++) 
        for (jj=max(j-1,lys);jj<min(j+2,lye);jj++) 
        for (ii=max(i-1,lxs);ii<min(i+2,lxe);ii++) {
            _nvert += nvert[kk][jj][ii];
            //PetscPrintf(PETSC_COMM_WORLD, "Debug8 with _nvert[%i][%i][%i]=%f in Calc_F_eul  in rotor_model.c\n", kk,jj,ii,_nvert);
        }

        //If solid node, set velocity to zero
        if ( _nvert > 2.9 ) { 
            lf_eul[k][j][i].x=0.0;
            lf_eul[k][j][i].y=0.0;
            lf_eul[k][j][i].z=0.0;
    //		PetscPrintf(PETSC_COMM_WORLD, "Debug9 in Calc_F_eul  in rotor_model.c\n");
        }
    }

    //Set velocity at boundaries
    for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
    for (i=xs; i<xe; i++) {
        int a=i, b=j, c=k;
        if(i==1) {
            if (i_periodicIB) a=mx-1;
            else a=0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug10A lf_eul[k][j][i].z =%f at a,b,c (%i,%i,%i)\n",lf_eul[k][j][i].z,a,b,c );

            lf_eul[k][j][i].x += lf_eul[c][b][a].x;
            lf_eul[k][j][i].y += lf_eul[c][b][a].y;
            lf_eul[k][j][i].z += lf_eul[c][b][a].z;
            lf_eul[c][b][a].x = 0.0;
            lf_eul[c][b][a].y = 0.0;
            lf_eul[c][b][a].z = 0.0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug10B lf_eul[k][j][i].z =%f at a,b,c (%i,%i,%i)\n",lf_eul[k][j][i].z,a,b,c );

        }

        if(i==mx-2) {
            if (i_periodicIB) a=0;
            else a=mx-1;
            lf_eul[k][j][i].x += lf_eul[c][b][a].x;
            lf_eul[k][j][i].y += lf_eul[c][b][a].y;
            lf_eul[k][j][i].z += lf_eul[c][b][a].z;
            lf_eul[c][b][a].x = 0.0;
            lf_eul[c][b][a].y = 0.0;
            lf_eul[c][b][a].z = 0.0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug11 in Calc_F_eul  in rotor_model.c\n");
        }

        if(j==1) {
            if (j_periodicIB) b=my-1;
            else b=0;
            lf_eul[k][j][i].x += lf_eul[c][b][a].x;
            lf_eul[k][j][i].y += lf_eul[c][b][a].y;
            lf_eul[k][j][i].z += lf_eul[c][b][a].z;
            lf_eul[c][b][a].x = 0.0;
            lf_eul[c][b][a].y = 0.0;
            lf_eul[c][b][a].z = 0.0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug12 in Calc_F_eul  in rotor_model.c\n");
        }

        if(j==my-2) {
            if (j_periodicIB) b=0;
            else b=my-1;
            lf_eul[k][j][i].x += lf_eul[c][b][a].x;
            lf_eul[k][j][i].y += lf_eul[c][b][a].y;
            lf_eul[k][j][i].z += lf_eul[c][b][a].z;
            lf_eul[c][b][a].x = 0.0;
            lf_eul[c][b][a].y = 0.0;
            lf_eul[c][b][a].z = 0.0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug13 in Calc_F_eul  in rotor_model.c\n");
        }

        if(k==1) {
            if (k_periodicIB) c=mz-1;
            else c=0;
            lf_eul[k][j][i].x += lf_eul[c][b][a].x;
            lf_eul[k][j][i].y += lf_eul[c][b][a].y;
            lf_eul[k][j][i].z += lf_eul[c][b][a].z;
            lf_eul[c][b][a].x = 0.0;
            lf_eul[c][b][a].y = 0.0;
            lf_eul[c][b][a].z = 0.0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug14 in Calc_F_eul  in rotor_model.c\n");
        }

        if(k==mz-2) {
            if (k_periodicIB) c=0;
            else c=mz-1;
            lf_eul[k][j][i].x += lf_eul[c][b][a].x;
            lf_eul[k][j][i].y += lf_eul[c][b][a].y;
            lf_eul[k][j][i].z += lf_eul[c][b][a].z;
            lf_eul[c][b][a].x = 0.0;
            lf_eul[c][b][a].y = 0.0;
            lf_eul[c][b][a].z = 0.0;
            //PetscPrintf(PETSC_COMM_WORLD, "Debug15 in Calc_F_eul  in rotor_model.c\n");
        }
    }

    DAVecRestoreArray(fda, Coor, &coor);
    DAVecRestoreArray(fda, user->lF_eul, &lf_eul);
    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);

    DAVecRestoreArray(da, user->lNvert, &nvert);
    DALocalToLocalBegin(fda, user->lF_eul, INSERT_VALUES, user->lF_eul);
    DALocalToLocalEnd(fda, user->lF_eul, INSERT_VALUES, user->lF_eul);

    DALocalToGlobal(fda, user->lF_eul, INSERT_VALUES, user->F_eul);

    return(0);
}


// exporting the forces on the turbines actuator disk model
PetscErrorCode Calc_forces_rotor(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi, char fname[80], int NumberOfObjects)
{
    PetscInt      	l, ibi;
    double        	pi = 3.141592653589793;
    PetscReal	A_Sum, F_Sum, P_Sum, U_Sum, rx, ry, rz, M_Sum;

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double C_T;	
        if (rotor_model==1) {	
        double indf_ax=ibm[ibi].indf_axis;

        C_T = 4.0 * indf_ax * (1-indf_ax);

        C_T = C_T / ( (1.0 - indf_ax)* (1.0 - indf_ax) );
        }

        if (rotor_model==4) {
            C_T=ibm[ibi].CT;
        }

        A_Sum = 0.0; P_Sum = 0.0; U_Sum = 0.0, F_Sum=0.0, M_Sum=0.0;
        double nx=fsi[ibi].nx_tb, ny=fsi[ibi].ny_tb, nz=fsi[ibi].nz_tb;

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            A_Sum += ibm[ibi].dA[l] ;

            double F_axis=ibm[ibi].F_lagr_x[l]*nx+ibm[ibi].F_lagr_y[l]*ny+ibm[ibi].F_lagr_z[l]*nz;
            double U_axis=ibm[ibi].U_lagr_x[l]*nx+ibm[ibi].U_lagr_y[l]*ny+ibm[ibi].U_lagr_z[l]*nz;

            F_Sum += F_axis*ibm[ibi].dA[l] ;
            P_Sum += F_axis*ibm[ibi].dA[l]*U_axis;

            U_Sum += U_axis*ibm[ibi].dA[l] ;

            rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;	
            ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;	
            rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;	

            double M_x= ry*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] - rz*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l];  
            double M_y= rz*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] - rx*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l];  
            double M_z= rx*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] - ry*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l];  

            M_Sum+=M_x*nx+M_y*ny+M_z*nz;
        }

        U_Sum=U_Sum/A_Sum;
    }

    return(0);
}

// exporting the forces on the turbines, actuator line model
PetscErrorCode Calc_forces_ACL(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int bi)
{
    PetscInt      l, ibi;
    double        pi = 3.141592653589793;
    PetscReal A_Sum = 0.0, P_Sum = 0.0, U_Sum = 0.0, F_Sum=0.0, M_Sum=0.0, Fx_Sum=0.0, Fy_Sum=0.0, Fz_Sum=0.0;
    PetscReal F_z, P, r, rx, ry, rz; 

    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


    for (ibi=0; ibi<NumberOfTurbines; ibi++) {
        A_Sum = 0.0; P_Sum = 0.0; U_Sum = 0.0; F_Sum=0.0; M_Sum=0.0; Fx_Sum=0.0; Fy_Sum=0.0; Fz_Sum=0.0;

        double nx=fsi[ibi].nx_tb, ny=fsi[ibi].ny_tb, nz=fsi[ibi].nz_tb;

        for (l=0; l<ibm[ibi].n_elmt; l++) {

            A_Sum += ibm[ibi].dA[l] ;

            double F_axis=ibm[ibi].F_lagr_x[l]*nx+ibm[ibi].F_lagr_y[l]*ny+ibm[ibi].F_lagr_z[l]*nz;
            double U_axis=ibm[ibi].U_lagr_x[l]*nx+ibm[ibi].U_lagr_y[l]*ny+ibm[ibi].U_lagr_z[l]*nz;

            F_Sum += F_axis*ibm[ibi].dA[l] ;
            P_Sum += F_axis*ibm[ibi].dA[l]*U_axis;

            U_Sum += U_axis*ibm[ibi].dA[l] ;

            Fx_Sum += ibm[ibi].F_lagr_x[l]*ibm[ibi].dA[l] ;
            Fy_Sum += ibm[ibi].F_lagr_y[l]*ibm[ibi].dA[l] ;
            Fz_Sum += ibm[ibi].F_lagr_z[l]*ibm[ibi].dA[l] ;

            rx = ibm[ibi].cent_x[l]-fsi[ibi].x_c;	
            ry = ibm[ibi].cent_y[l]-fsi[ibi].y_c;	
            rz = ibm[ibi].cent_z[l]-fsi[ibi].z_c;	

            double M_x= ry*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l] - rz*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l];  
            double M_y= rz*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l] - rx*ibm[ibi].F_lagr_z[l] * ibm[ibi].dA[l];  
            double M_z= rx*ibm[ibi].F_lagr_y[l] * ibm[ibi].dA[l] - ry*ibm[ibi].F_lagr_x[l] * ibm[ibi].dA[l];  

            M_Sum+=M_x*nx+M_y*ny+M_z*nz;
        }

        U_Sum=U_Sum/A_Sum;

        fsi[ibi].Torque_aero=-M_Sum;    
        fsi[ibi].Force_axis=-F_Sum;   

        double P_moment=M_Sum*fsi[ibi].angvel_axis;
        double R0=fsi[0].r_rotor/reflength_wt;
        double rho_cfd=1.0; // no levelset 	
        double T_wt=fsi[0].r_rotor/refvel_wt; // the 0 turbine cannot be in the other turbines' wakes.  
        double T_cfd=R0/refvel_cfd; // the 0 turbine cannot be in the other turbines' wakes.  
        double ratio_rho=rho_air/rho_cfd;
        double ratio_L=reflength_wt;	
        double ratio_T=T_wt/T_cfd;	
        double ratio_V=refvel_wt/refvel_cfd;

        double Time_real=ti*user->dt*ratio_T;	
        double ang_real=fsi[ibi].ang_axis;
        double angvel_real=fsi[ibi].angvel_axis/ratio_T;
        double Force_real=ratio_rho*fsi[ibi].Force_axis*pow(ratio_V,2)*pow(ratio_L,2);
        double Torquea_real=ratio_rho*fsi[ibi].Torque_aero*pow(ratio_V,2)*pow(ratio_L,3);
        double Torquec_real=fsi[ibi].Torque_generator;
        double Torquec_cfd=fsi[ibi].Torque_generator/(ratio_rho*pow(ratio_V,2)*pow(ratio_L,3));
        double Uref_real=ibm[ibi].U_ref*ratio_V;
        double Ud_real=U_Sum*ratio_V;

        double ratio_torque = ratio_rho*pow(ratio_V,2)*pow(ratio_L,3);
        //PetscPrintf(PETSC_COMM_WORLD, "ti=%i, dt=%f, ratio_T=%f, R0=%f, T_wt=%f, ratio_T=%f\n", ti, user->dt,R0, T_wt, ratio_T );

        if (!rank) {
            FILE *f;
            char filen[80];
            sprintf(filen, "Turbine_AL%2.2d_%2.2d",rotor_model,ibi);
            if (ti==1) {
                f = fopen(filen, "w");
                PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"time\", \"angle\", \"angvel_axis\", \"Force_axis\", \"Torque_fluid\", \"Torque_generator\", \"Uref\", \"Ud\" ,\"TSR\", \"pitch1\", \"pitch2\", \"pitch3\", \"pitch1_IPC\", \"pitch2_IPC\", \"pitch3_IPC\", \"Moment1_bladebending\", \"Moment2_bladebending\", \"Moment3_bladebending\", \"angvel_axis_err\" \n");
            } else f = fopen(filen, "a");

            PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n",ti*user->dt, fsi[ibi].ang_axis, fsi[ibi].angvel_axis,fsi[ibi].Force_axis, fsi[ibi].Torque_aero, Torquec_cfd, ibm[ibi].U_ref, U_Sum, ibm[ibi].Tipspeedratio, ibm[ibi].pitch[0], ibm[ibi].pitch[1], ibm[ibi].pitch[2], ibm[ibi].pitch_IPC[0], ibm[ibi].pitch_IPC[1], ibm[ibi].pitch_IPC[2], fsi[ibi].Moment_bladebending[0], fsi[ibi].Moment_bladebending[1], fsi[ibi].Moment_bladebending[2], fsi[ibi].angvel_axis_err*ratio_T);
            fclose(f);

            sprintf(filen, "Turbine_AL_real%2.2d_%2.2d",rotor_model,ibi);
            if (ti==1) {
                f = fopen(filen, "w");
                PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=\"time (s)\", \"angle (rad)\", \"angvel_axis (s<sup>-1</sup>)\", \"Force_axis (N)\", \"Torque_fluid (Nm)\", \"Torque_generator (Nm)\", \"Uref (m/s)\", \"Ud (m/s)\", \"TSR \", \"pitch1\", \"pitch2\", \"pitch3\", \"pitch1_IPC\", \"pitch2_IPC\", \"pitch3_IPC\", \"Moment1_bladebending (Nm)\", \"Moment2_bladebending (Nm)\", \"Moment3_bladebending (Nm)\", \"angvel_axis_err\" \n");
            } else f = fopen(filen, "a");

            PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n",Time_real, ang_real, angvel_real, Force_real, Torquea_real, Torquec_real, Uref_real, Ud_real, ibm[ibi].Tipspeedratio, ibm[ibi].pitch[0], ibm[ibi].pitch[1], ibm[ibi].pitch[2], ibm[ibi].pitch_IPC[0], ibm[ibi].pitch_IPC[1], ibm[ibi].pitch_IPC[2], fsi[ibi].Moment_bladebending[0]*ratio_torque, fsi[ibi].Moment_bladebending[1]*ratio_torque, fsi[ibi].Moment_bladebending[2]*ratio_torque, fsi[ibi].angvel_axis_err);
            fclose(f);
        }
    }
    return(0);
}


// calculating the reference velocity for actuator line model
PetscErrorCode Uref_ACL(UserCtx *user, IBMNodes *ibm_ACL, IBMNodes *ibm_ACD, FSInfo *fsi_wt, int NumberOfObjects)
{
    PetscInt      l, ibi;
    PetscReal     A_Sum, U_Sum;

    PetscPrintf(PETSC_COMM_WORLD, "Start Uref_ACL for Number of Objects =%i and fsi_wt=%f\n", NumberOfObjects, fsi_wt);

    Calc_U_lagr(user, ibm_ACD, fsi_wt, NumberOfObjects);

    for (ibi=0; ibi<NumberOfObjects; ibi++) {

        double nx=fsi_wt[ibi].nx_tb, ny=fsi_wt[ibi].ny_tb, nz=fsi_wt[ibi].nz_tb;
        PetscPrintf(PETSC_COMM_WORLD, "nx ny nz: %le %le %le \n", nx, ny ,nz);

            U_Sum = 0.0; A_Sum = 0.0; 
            for (l=0; l<ibm_ACD[ibi].n_elmt; l++) {
                double U_axis=ibm_ACD[ibi].U_lagr_x[l]*nx+ibm_ACD[ibi].U_lagr_y[l]*ny+ibm_ACD[ibi].U_lagr_z[l]*nz;
                U_Sum += U_axis*ibm_ACD[ibi].dA[l] ;
                        A_Sum += ibm_ACD[ibi].dA[l] ;
            }

            U_Sum /= A_Sum;

        ibm_ACL[ibi].U_ref = U_Sum; // / (1.0 - indf_a);

    }

    PetscPrintf(PETSC_COMM_WORLD, "Finish Uref_ACL with U_Sum=%f and A_Sum=%f\n", U_Sum, A_Sum);
    return(0);
}


/* ==================================================================================             */
PetscErrorCode Export_lagrdata(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension)
{

    int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
    PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
    int i,j,k;
    int n1e, n2e, n3e;
    PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
    PetscReal rx,ry,rz;
    int n1;

    if (dimension == 2) {
        if ( ti==tistart || ti % tiout == 0) {
            int rank=0;
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            int i;

            if (!rank) {
                FILE *f;
                char filen[80];
                sprintf(filen, "%s/%s%06d_%03d_nf.dat",path,fname,ti,ibi);
                f = fopen(filen, "w");
                int n_v=ibm->n_v; 	
                int n_elmt=ibm->n_elmt; 	

                PetscFPrintf(PETSC_COMM_WORLD, f, "TITLE = \"Actuator surface mesh\" \n ");
                PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x, y, z, ub_x, ub_y, ub_z, tmprtb, n_x,n_y,n_z, F_lagr_x, F_lagr_y, F_lagr_z, s2l, color, U_lagr_x, U_lagr_y, U_lagr_z, Tmprt_lagr, Ftmprt_lagr\n");
                PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"TRIANGLES\", N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-7]=NODAL,[8-20]=CELLCENTERED)\n", n_v, n_elmt);
                PetscFPrintf(PETSC_COMM_WORLD, f, "STRANDID=0.1 SOLUTIONTIME=%le \n", ((double)ti)*dt);

                for (i=0; i<n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->x_bp[i]);
                }
                for (i=0; i<n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->y_bp[i]);
                }
                for (i=0; i<n_v; i++) {	
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->z_bp[i]);
                }

                for (i=0; i<n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->u[i].x);
                }
                for (i=0; i<n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->u[i].y);
                }
                for (i=0; i<n_v; i++) {	
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->u[i].z);
                }
                for (i=0; i<n_v; i++) {	
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->tmprt[i]);
                }
                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nf_x[i]);
                }
                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nf_y[i]);
                }
                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nf_z[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->F_lagr_x[i]);
                }
                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->F_lagr_y[i]);
                }
                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->F_lagr_z[i]);
                }
                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", ibm->s2l[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", ibm->color[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->U_lagr_x[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->U_lagr_y[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->U_lagr_z[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Tmprt_lagr[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Ftmprt_lagr[i]);
                }

                for (i=0; i<n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
                }

                fclose(f);
            }
        }
    }

    if (dimension == 1) {
        if ( ti==tistart || ti % tiout == 0 /*ti == (ti/tiout) * tiout */) {
            int rank=0;
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

            if (!rank) {
                FILE *f;
                char filen[80];  
                sprintf(filen, "%s/%s%06d_%03d_nf.dat", path,fname, ti,ibi);
                f = fopen(filen, "w");

                PetscFPrintf(PETSC_COMM_WORLD, f, "TITLE = \"Actuator line mesh\" \n ");
                PetscFPrintf(PETSC_COMM_WORLD, f, "VARIABLES = \"X\", \"Y\", \"Z\", \"ub_x\", \"ub_y\", \"ub_z\", \"tmprtb\", \"color\", \"U_lagr_x\", \"U_lagr_y\", \"U_lagr_z\", \"F_lagr_x\", \"F_lagr_y\", \"F_lagr_z\", \"Tmprt_lagr\", \"Ftmprt_lagr\"\n");
                PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"P_1\", DATAPACKING=BLOCK, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG, VARLOCATION=([1-7]=NODAL,[8-16]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);
                PetscFPrintf(PETSC_COMM_WORLD, f, "STRANDID=0.1 SOLUTIONTIME=%le \n", ((double)ti)*dt);

                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->x_bp[i]);
                }
                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->y_bp[i]);
                }
                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->z_bp[i]);
                }
                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->u[i].x);
                }
                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->u[i].y);
                }
                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->u[i].z);
                }
                for (i=0; i<ibm->n_v; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->tmprt[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%d \n", ibm->color[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->U_lagr_x[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->U_lagr_y[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->U_lagr_z[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->F_lagr_x[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->F_lagr_y[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->F_lagr_z[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->Tmprt_lagr[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%le \n", ibm->Ftmprt_lagr[i]);
                }
                for (i=0; i<ibm->n_elmt; i++) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d \n", ibm->nv1[i]+1, ibm->nv2[i]+1);
                }
                fclose(f);
            }
        } 
    }

    return(0);
}

/**
 * rotate the blade along the ref line for pitch 
 */

PetscErrorCode bladepitch_Rot(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects) 
{

    int ibi, elmt_s, elmt_l, l;
    int n1e_s, n2e_s, n3e_s;
    int n1e_l, n2e_l, n3e_l;
    double nfx_l, nfy_l, nfz_l;
    Cmpnts p, q, nr, na, nt; 

    int ne[3];

    double xp, yp, zp;
    double xr, yr, zr;

    for (ibi=0; ibi<NumberOfObjects; ibi++) {
        int indicator[ibm_surface[ibi].n_v];

        for (int nl=0; nl < ibm_surface[ibi].n_v; nl++)
            indicator[nl] = 0;

        for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) {
            elmt_l = ibm_surface[ibi].s2l[elmt_s];

            int NelmtPerBlade = (int)((double)(ibm_line[ibi].n_elmt/ibm_line[ibi].num_blade)+0.1);

            int iblade = (double)(elmt_l/NelmtPerBlade);

            n1e_l = ibm_line[ibi].nv1[elmt_l]; 
            n2e_l = ibm_line[ibi].nv2[elmt_l]; 

            double x1_l, y1_l, z1_l, x2_l, y2_l, z2_l; 

            x1_l = ibm_line[ibi].x_bp[n1e_l];
            y1_l = ibm_line[ibi].y_bp[n1e_l];
            z1_l = ibm_line[ibi].z_bp[n1e_l];

            x2_l = ibm_line[ibi].x_bp[n2e_l];
            y2_l = ibm_line[ibi].y_bp[n2e_l];
            z2_l = ibm_line[ibi].z_bp[n2e_l];

            double dx12 = x2_l - x1_l;
            double dy12 = y2_l - y1_l;
            double dz12 = z2_l - z1_l;
  
            double dr = sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)+1.e-19;
      
            nfx_l = dx12/dr; 
            nfy_l = dy12/dr; 
            nfz_l = dz12/dr;

            //PetscPrintf(PETSC_COMM_WORLD, "B \n");
            
            // surface points  
            n1e_s = ibm_surface[ibi].nv1[elmt_s]; 
            n2e_s = ibm_surface[ibi].nv2[elmt_s]; 
            n3e_s = ibm_surface[ibi].nv3[elmt_s];

            ne[0] = n1e_s;
            ne[1] = n2e_s;
            ne[2] = n3e_s;

            int  k;

            for (k=0; k<3; k++) {

                //PetscPrintf(PETSC_COMM_WORLD, "C \n");

                if (indicator[ne[k]] == 0) {
                    indicator[ne[k]] = 1;
                    //PetscPrintf(PETSC_COMM_WORLD, "xr=%.3f yr=%.3f zr=%.3f \n", xr, yr, zr);
                    xr = ibm_surface[ibi].x_bp[ne[k]]; 
                    yr = ibm_surface[ibi].y_bp[ne[k]]; 
                    zr = ibm_surface[ibi].z_bp[ne[k]]; 

                    dx12 = xr - x1_l;
                    dy12 = yr - y1_l;
                    dz12 = zr - z1_l;
      
                    double ds = dx12*nfx_l + dy12*nfy_l + dz12*nfz_l;

                    xp = x1_l + ds*nfx_l; 
                    yp = y1_l + ds*nfy_l; 
                    zp = z1_l + ds*nfz_l; 

                    // rotate 
                    p.x = xr - xp;
                    p.y = yr - yp;
                    p.z = zr - zp;

                    na.x = nfx_l;	
                    na.y = nfy_l;	
                    na.z = nfz_l;	

                    double theta=(ibm_surface[ibi].pitch[iblade]-ibm_surface[ibi].pitch_old[iblade])*M_PI/180.0;

                    //PetscPrintf(PETSC_COMM_WORLD, "x1=%.3f y1=%.3f z1=%.3f \n", x1_l, y1_l, z1_l);
                    //PetscPrintf(PETSC_COMM_WORLD, "x2=%.3f y2=%.3f z2=%.3f \n", x2_l, y2_l, z2_l);
                    //PetscPrintf(PETSC_COMM_WORLD, "xp=%.3f yp=%.3f zp=%.3f \n", xp, yp, zp);
                    //PetscPrintf(PETSC_COMM_WORLD, "theta=%.3f \n", ibm_surface[ibi].pitch[k]);
                    //PetscPrintf(PETSC_COMM_WORLD, "theta_old=%.3f \n", ibm_surface[ibi].pitch_old[k]);
                    //PetscPrintf(PETSC_COMM_WORLD, "nax=%.3f nay=%.3f naz=%.3f \n", na.x, na.y, na.z);
                    //PetscPrintf(PETSC_COMM_WORLD, "px=%.3f py=%.3f pz=%.3f \n", p.x, p.y, p.z);

                    q=ArbitraryRotate(p,theta,na);

                    //PetscPrintf(PETSC_COMM_WORLD, "qx=%.3f qy=%.3f qz=%.3f \n", q.x, q.y, q.z);
                    xr = q.x+xp;
                    yr = q.y+yp;
                    zr = q.z+zp;

                    //PetscPrintf(PETSC_COMM_WORLD, "xr=%.3f yr=%.3f zr=%.3f \n", xr, yr, zr);

                    //PetscPrintf(PETSC_COMM_WORLD, "\n");
                    //PetscPrintf(PETSC_COMM_WORLD, "\n");

                    ibm_surface[ibi].x_bp[ne[k]] = xr;
                    ibm_surface[ibi].y_bp[ne[k]] = yr;
                    ibm_surface[ibi].z_bp[ne[k]] = zr;
                }
                //PetscPrintf(PETSC_COMM_WORLD, "D \n");
            }
        }
    }

    for (ibi=0; ibi<NumberOfObjects; ibi++) {
        for (int k=0; k<3; k++)
            ibm_surface[ibi].pitch_old[k] = ibm_surface[ibi].pitch[k];
    }

    return(0);
}


/* ==================================================================================             */
PetscErrorCode rotor_Rot(FSInfo *FSinfo, IBMNodes *ibm, PetscReal dt, int ibi, char fname[80], int dimension)
{

    int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
    PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
    int i,j,k;
    int n1e, n2e, n3e;
    PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
    PetscReal rx,ry,rz;
    int n1;

    double rot_angle;

    Cmpnts p, q, nr, na, nt; 

    if (ti==tistart && rstart_turbinerotation) 
    {
        for (i=0; i<n_v; i++) 
        {
            p.x=ibm->x_bp[i]-FSinfo->x_c;
            p.y=ibm->y_bp[i]-FSinfo->y_c;
            p.z=ibm->z_bp[i]-FSinfo->z_c;

            na.x=FSinfo->nx_tb;	
            na.y=FSinfo->ny_tb;	
            na.z=FSinfo->nz_tb;	
    PetscPrintf(PETSC_COMM_WORLD, "rotor_Rot1 %le %le %le\n", p.x,p.y,p.z);

            double theta=FSinfo->ang_axis;
    PetscPrintf(PETSC_COMM_WORLD, "rotor_Rot2 %le %le %le %le\n", na.x,na.y,na.z, theta);

            q=ArbitraryRotate(p,theta,na);

            ibm->x_bp[i]=q.x+FSinfo->x_c;
            ibm->y_bp[i]=q.y+FSinfo->y_c;
            ibm->z_bp[i]=q.z+FSinfo->z_c;
        }
    }


    PetscPrintf(PETSC_COMM_WORLD, "angvel of turbine %le \n", FSinfo->angvel_axis);

    for (i=0; i<n_v; i++) 
    {
        p.x=ibm->x_bp[i]-FSinfo->x_c;
        p.y=ibm->y_bp[i]-FSinfo->y_c;
        p.z=ibm->z_bp[i]-FSinfo->z_c;

        na.x=FSinfo->nx_tb;	
        na.y=FSinfo->ny_tb;	
        na.z=FSinfo->nz_tb;	

        double theta=FSinfo->angvel_axis*dt;

        q=ArbitraryRotate(p,theta,na);

        ibm->x_bp[i]=q.x+FSinfo->x_c;
        ibm->y_bp[i]=q.y+FSinfo->y_c;
        ibm->z_bp[i]=q.z+FSinfo->z_c;

        double rx = ibm->x_bp[i]-FSinfo->x_c;
        double ry = ibm->y_bp[i]-FSinfo->y_c;
        double rz = ibm->z_bp[i]-FSinfo->z_c;

        double rr = sqrt(rx*rx+ry*ry+rz*rz)+1.e-19;
        nr.x = rx/rr; 
        nr.y = ry/rr; 
        nr.z = rz/rr;

        nt.x=na.y*nr.z-na.z*nr.y;
        nt.y=na.z*nr.x-na.x*nr.z;
        nt.z=na.x*nr.y-na.y*nr.x;
        
        double Ut=FSinfo->angvel_axis*rr;

        ibm->u[i].x = Ut*nt.x;
        ibm->u[i].y = Ut*nt.y;
        ibm->u[i].z = Ut*nt.z;

    }

    if (turbinestructuremodel) 
    {
        ibm->x_bp[i]+=ibm->disp_x[i];
        ibm->y_bp[i]+=ibm->disp_y[i];
        ibm->z_bp[i]+=ibm->disp_z[i];
    }

    if (dimension == 2) 
    {
     
        for (i=0; i<n_elmt; i++) 
        {

            n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
            dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
            dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
            dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 

            dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
            dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
            dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 

            ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
            ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
            ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

            dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
            ibm->nf_z[i]*ibm->nf_z[i]);

            ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;

            if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
            (((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) 
            {
                ibm->ns_x[i] = 1.;
                ibm->ns_y[i] = 0.;
                ibm->ns_z[i] = 0 ;

                // nt = ns x nf
                ibm->nt_x[i] = 0.;
                ibm->nt_y[i] = 1.;
                ibm->nt_z[i] = 0.;
            } else {
                ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
                ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
                ibm->ns_z[i] = 0 ;

                // nt = ns x nf
                ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
                ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
                ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
            }

            ibm->dA[i] = dr/2.;

            ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e]+ibm->x_bp[n3e])/3.;
            ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e]+ibm->y_bp[n3e])/3.;
            ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e]+ibm->z_bp[n3e])/3.;
        }
    }

    // for ACL, the nf_x, nf_y, nf_z denote the direction of the actuator line. and dA denotes the length of each element.
    int nb;

    int n_elmt_1 = (ibm->n_elmt)/ibm->num_blade;

    //      PetscPrintf(PETSC_COMM_WORLD, "rotormodel %d \n", rotor_model);
    if (dimension == 1) {
        for (nb=0; nb<ibm->num_blade; nb++) {

            for (j=0; j<n_elmt_1; j++) 
            {
                i = nb * n_elmt_1 + j;
      
                n1e = ibm->nv1[i]; n2e = ibm->nv2[i]; 
                dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
                dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
                dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];
      
                ibm->nf_x[i] = dx12;
                ibm->nf_y[i] = dy12;
                ibm->nf_z[i] = dz12;
      
                dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + ibm->nf_z[i]*ibm->nf_z[i]);
      
                ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
          
                ibm->ns_x[i] = 0.;     
                ibm->ns_y[i] = 0.;     
                ibm->ns_z[i] = 0. ;

                ibm->nt_x[i] = 0.;
                ibm->nt_y[i] = 0.;
                ibm->nt_z[i] = 0.;
         
                ibm->dA[i] = dr;

                ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e])/2.;
                ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e])/2.;
                ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e])/2.;
            }
        }
    }
    return(0);
}


PetscErrorCode refAL_Rot(FSInfo *FSinfo, IBMNodes *ibm, IBMNodes *ibm_ref, int ibi)
{

    int n_v = ibm->n_v, n_elmt = ibm->n_elmt;
    PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
    int i,j,k;
    int n1e, n2e, n3e;
    PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
    PetscReal rx,ry,rz;
    int n1;

    double rot_angle;

    Cmpnts p, q, nr, na, nt; 


    for (i=0; i<n_v; i++) 
    {
        p.x=ibm->x_bp[i]-FSinfo->x_c;
        p.y=ibm->y_bp[i]-FSinfo->y_c;
        p.z=ibm->z_bp[i]-FSinfo->z_c;
        
        na.x=FSinfo->nx_tb;	
        na.y=FSinfo->ny_tb;	
        na.z=FSinfo->nz_tb;	
        
        double isign = (1.e-19+FSinfo->angvel_axis)/(fabs(FSinfo->angvel_axis)+1.e-19);
        
        if (isign>0.0) isign = 1.0;
        else isign = -1.0;

        double theta = isign * refangle_AL*M_PI/180.0;
        
        //PetscPrintf(PETSC_COMM_WORLD, "theta %le refangle_AL %le isign %le\n", theta, refangle_AL, isign);

        q=ArbitraryRotate(p,theta,na);
        
        ibm_ref->x_bp[i]=q.x+FSinfo->x_c;
        ibm_ref->y_bp[i]=q.y+FSinfo->y_c;
        ibm_ref->z_bp[i]=q.z+FSinfo->z_c;
        
    }
    

    int nb;

    int n_elmt_1 = (ibm->n_elmt)/ibm->num_blade;

    for (nb=0; nb<ibm->num_blade; nb++) 
    {

        for (j=0; j<n_elmt_1; j++) 
        {
            i = nb * n_elmt_1 + j;

            n1e = ibm_ref->nv1[i]; n2e = ibm_ref->nv2[i]; 
            dx12 = ibm_ref->x_bp[n2e] - ibm_ref->x_bp[n1e];
            dy12 = ibm_ref->y_bp[n2e] - ibm_ref->y_bp[n1e];
            dz12 = ibm_ref->z_bp[n2e] - ibm_ref->z_bp[n1e];

            ibm_ref->nf_x[i] = dx12;
            ibm_ref->nf_y[i] = dy12;
            ibm_ref->nf_z[i] = dz12;

            dr = sqrt(ibm_ref->nf_x[i]*ibm_ref->nf_x[i] + ibm_ref->nf_y[i]*ibm_ref->nf_y[i] + ibm_ref->nf_z[i]*ibm_ref->nf_z[i]);

            ibm_ref->nf_x[i] /=dr; ibm_ref->nf_y[i]/=dr; ibm_ref->nf_z[i]/=dr;
      
            ibm_ref->ns_x[i] = 0.;     
            ibm_ref->ns_y[i] = 0.;     
            ibm_ref->ns_z[i] = 0. ;

            ibm_ref->nt_x[i] = 0.;
            ibm_ref->nt_y[i] = 0.;
            ibm_ref->nt_z[i] = 0.;
     
            ibm_ref->dA[i] = dr;

            ibm_ref->cent_x[i]= (ibm_ref->x_bp[n1e]+ibm_ref->x_bp[n2e])/2.;
            ibm_ref->cent_y[i]= (ibm_ref->y_bp[n1e]+ibm_ref->y_bp[n2e])/2.;
            ibm_ref->cent_z[i]= (ibm_ref->z_bp[n1e]+ibm_ref->z_bp[n2e])/2.;
        }
    }
    
    return(0);
}


PetscErrorCode calc_s2l(IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects) 
{

    int ibi, elmt_s, elmt_l;
    double rs, rl, dd_min, dd;

    double xs, ys, zs, xl, yl, zl;	
    int colors, colorl;
    double nx, ny, nz;

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) 
        {
            dd_min = 100000;
            xs = ibm_surface[ibi].cent_x[elmt_s] - fsi_surface[ibi].x_c;
            ys = ibm_surface[ibi].cent_y[elmt_s] - fsi_surface[ibi].y_c;
            zs = ibm_surface[ibi].cent_z[elmt_s] - fsi_surface[ibi].z_c;
            colors = ibm_surface[ibi].color[elmt_s];
            // kevin debugging???			
            int temp_elmt_l;
            temp_elmt_l = 99;
            int temp_elmt_colorl;
            temp_elmt_colorl = 99;
            for (elmt_l=0; elmt_l<ibm_line[ibi].n_elmt; elmt_l++) 
            {
                xl = ibm_line[ibi].cent_x[elmt_l] - fsi_line[ibi].x_c;
                yl = ibm_line[ibi].cent_y[elmt_l] - fsi_line[ibi].y_c;
                zl = ibm_line[ibi].cent_z[elmt_l] - fsi_line[ibi].z_c;

                rl = sqrt(xl*xl+yl*yl+zl*zl);

                nx = xl/rl; ny = yl/rl; nz = zl/rl;

                rs = fabs(xs*nx+ys*ny+zs*nz);

                colorl = ibm_line[ibi].color[elmt_l];

                //PetscPrintf(PETSC_COMM_WORLD, "rs=%le, rl=%le, colors=%d, colorl=%d for %d th line segment\n", rs, rl, colors, colorl, elmt_l);
                double icolor = double(colors-colorl);
                dd=fabs(rs-rl);

                //Temporary edit by KFlora to check color matching
                if (dd-dd_min<1.e-9 && colors==colorl) 
                {
                
                // if (dd-dd_min<1.e-9 ) {
                    temp_elmt_l = elmt_l;
                    temp_elmt_colorl = colorl;
                    dd_min=dd;
                    ibm_surface[ibi].s2l[elmt_s] = elmt_l;

                    //	if (dd-dd_min<1.e-9 && colors==colorl) {
                    //		dd_min=dd;
                    //		ibm_surface[ibi].s2l[elmt_s] = elmt_l;
                }
            }
            // kevin debugging??
        //	PetscPrintf(PETSC_COMM_WORLD, "For elmt_s = %i with color= %i, closest elmt_l = %i with color = %i at distance = %f\n", elmt_s, colors, temp_elmt_l, temp_elmt_colorl, dd_min);
        }
    }
    	return(0);
}


PetscErrorCode UlagrProjection_s2l(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects) 
{

    int ibi, elmt_s, elmt_l, l;

    double count, Ux, Uy, Uz;
    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (elmt_l=0; elmt_l<ibm_line[ibi].n_elmt; elmt_l++) 
        {
            Ux=0.0; Uy=0.0; Uz=0.0; count=0.0;
            for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) 
            {
                if (elmt_l==ibm_surface[ibi].s2l[elmt_s]) 
                {
                    Ux+=ibm_surface[ibi].U_lagr_x[elmt_s];
                    Uy+=ibm_surface[ibi].U_lagr_y[elmt_s];
                    Uz+=ibm_surface[ibi].U_lagr_z[elmt_s];
                    count+=1.0;
                } 
            }

            if (count>0.5) 
            {
                ibm_line[ibi].U_lagr_x[elmt_l]=Ux/count;
                ibm_line[ibi].U_lagr_y[elmt_l]=Uy/count;
                ibm_line[ibi].U_lagr_z[elmt_l]=Uz/count;
            }
        }
    }

        return(0);
}


//Projection of turbine forces from the actuator line to the actuator surface
PetscErrorCode ForceProjection_l2s(UserCtx *user, IBMNodes *ibm_surface, IBMNodes *ibm_line, FSInfo *fsi_surface, FSInfo *fsi_line, int NumberOfObjects) 
{

    int ibi, elmt_s, elmt_l, l;

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {

        for (elmt_s=0; elmt_s<ibm_surface[ibi].n_elmt; elmt_s++) 
        {
            elmt_l=ibm_surface[ibi].s2l[elmt_s];
                        ibm_surface[ibi].F_lagr_x[elmt_s] = ibm_line[ibi].F_lagr_x[elmt_l] /ibm_line[ibi].chord_blade[elmt_l];
                        ibm_surface[ibi].F_lagr_y[elmt_s] = ibm_line[ibi].F_lagr_y[elmt_l] /ibm_line[ibi].chord_blade[elmt_l];
                        ibm_surface[ibi].F_lagr_z[elmt_s] = ibm_line[ibi].F_lagr_z[elmt_l] /ibm_line[ibi].chord_blade[elmt_l];
        }
    }
    return(0);
}



double dfunc_4htail(double r)
{
    if (fabs(r) <= 1.0) 
    {
        return 1-r*r;
    } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) 
    {
        return 2-3*fabs(r)+r*r;
    } else 
    {
        return 0.0;
    }
}

double dfunc_4h(double r)
{
    if (fabs(r) <= 1.0) 
    {
        return (3.0-2.0*fabs(r)+sqrt(1.0+4.0*fabs(r)-4.0*r*r))/8.0;
    } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) 
    {
        return (5.0-2.0*fabs(r)-sqrt(-7.0+12.0*fabs(r)-4.0*r*r))/8.0;
    } else 
    {
        return 0.0;
    }
}

double dfunc_s3h(double r)
{
    if (fabs(r) <= 1.0) 
    {
        return 17.0/48.0 + sqrt(3.0)*3.14159265/108.0 + fabs(r)/4.0 - r*r/4.0 + (1.0-2.0*fabs(r))*sqrt(-12.0*r*r+12.0*fabs(r)+1.0)/16.0 - sqrt(3.0)*asin(sqrt(3.0)*(2.0*fabs(r)-1.0)/2.0)/12.0;
    } else if (fabs(r) >= 1.0 && fabs(r) <= 2.0) 
    {
        return 55.0/48.0 - sqrt(3.0)*3.14159265/108.0 - 13.0*fabs(r)/12.0 + r*r/4.0 + (2.0*fabs(r)-3.0)*sqrt(-12.0*r*r+36.0*fabs(r)-23.0)/48.0 + sqrt(3.0)*asin(sqrt(3.0)*(2.0*fabs(r)-3.0)/2.0)/36.0;
    } else 
    { 
        return 0.0;
    }
}


double dfunc_2h(double r)
{
    if (fabs(r) < 1.0) 
    {
        return 1.0-fabs(r);
    } else 
    {
        return 0.0;
    }
}


double dfunc_s4h(double r)
{

    if (fabs(r) <= 0.5) 
    {
        return 3.0/8.0+3.14159265/32.0-pow(r,2)/4.0;
    } else if (fabs(r) >= 0.5 && fabs(r) <= 1.5) 
    {
        return 1.0/4.0+(1.0-fabs(r))*sqrt(-2.0+8.0*fabs(r)-4.0*pow(r,2))/8.0-asin(sqrt(2.0)*(fabs(r)-1.0))/8.0;
    } else if (fabs(r) >= 1.5 && fabs(r) <= 2.5) 
    {
        return
        17.0/16.0-3.14159265/64.0-3.0*fabs(r)/4.0+pow(r,2)/8.0+(fabs(r)-2.0)*sqrt(-14.0+16.0*fabs(r)-4.0*pow(r,2))/16.0+asin(sqrt(2.0)*(fabs(r)-2.0))/16.0;
    } else {
        return 0.0;
    }
}



double dfunc_sc4h(double r)
{

    if (fabs(r) <= 1.5) 
    {
        return 0.25*(M_PI+2*sin(0.25*M_PI*(2*fabs(r)+1))-2*sin(0.25*M_PI*(2*fabs(r)-1)))/M_PI;
    } else if (fabs(r) >= 1.5 && fabs(r) <= 2.5) 
    {
        return -0.125*(-5*M_PI+2*M_PI*fabs(r)+4*sin(0.25*M_PI*(2*fabs(r)-1)))/M_PI;
    } else 
    {
        return 0.0;
    }

}


Cmpnts ArbitraryRotate(Cmpnts p,double theta,Cmpnts r) 
{
    Cmpnts q = {0.0,0.0,0.0};
    double costheta,sintheta;

    double rr=sqrt(r.x*r.x+r.y*r.y+r.z*r.z)+1.0e-11;
    r.x=r.x/rr; r.y=r.y/rr; r.z=r.z/rr;

    costheta = cos(theta);
    sintheta = sin(theta);

    q.x += (costheta + (1 - costheta) * r.x * r.x) * p.x;
    q.x += ((1 - costheta) * r.x * r.y - r.z * sintheta) * p.y;
    q.x += ((1 - costheta) * r.x * r.z + r.y * sintheta) * p.z;

    q.y += ((1 - costheta) * r.x * r.y + r.z * sintheta) * p.x;
    q.y += (costheta + (1 - costheta) * r.y * r.y) * p.y;
    q.y += ((1 - costheta) * r.y * r.z - r.x * sintheta) * p.z;

    q.z += ((1 - costheta) * r.x * r.z - r.y * sintheta) * p.x;
    q.z += ((1 - costheta) * r.y * r.z + r.x * sintheta) * p.y;
    q.z += (costheta + (1 - costheta) * r.z * r.z) * p.z;

    return(q);
}


double dfunc_4uniform(double r)
{

    if (fabs(r) < 2.1) {
        return 0.25;
    } else {
        return 0.0;
    }
}


double dfunc_6uniform(double r)
{

    if (fabs(r) < 3.00001) {
        return 0.166666667;
    } else {
        return 0.0;
    }
}


double dfunc_nh(double r, double n)
{
    if (fabs(r) < n) {
        return (n-fabs(r))/pow(n,2);
    } else {
        return 0.0;
    }
}


double dfunc_exp(double r, double n)
{
    return exp(-(r/n)*(r/n))/(pow(n,1)*pow(3.1415926,0.5));
}



double dfunc_2hs1(double r)
{
    if (fabs(r) <=0.5){
        return 3.0/4.0-pow(r,2.0);
    } else if (fabs(r) >= 0.5 && fabs(r) <=1.5) {
        return pow((2*fabs(r) - 3.0),2.0)/8.0;
    } else {
        return 0.0;
    }
}


double dfunc_2hs2(double r)
{
    if (fabs(r) <=1.0){
        return pow(fabs(r),3.0)/2.0 - pow(r,2) + 2.0/3.0;
    } else if (fabs(r) >= 1.0 && fabs(r) <=2.0) {
        return -pow((fabs(r) - 2.0),3.0)/6.0;
    } else {
        return 0.0;
    }
}



double dfunc_2hs3(double r)
{
    if (fabs(r) <=0.5){
        return pow(r,4.0)/4. - (5.*pow(r,2.))/8. + 115./192.;
    } else if (fabs(r) >= 0.5 && fabs(r) <=1.5) {
        return - pow(r,4.)/6. + (5.*pow(fabs(r),3.))/6. - (5.*pow(r,2.))/4. + (5.*fabs(r))/24. + 55./96.;
    } else if (fabs(r) >= 1.5 && fabs(r) <=2.5) {
        return pow((2.*fabs(r) - 5.),4.)/384.;
    } else {
        return 0.0;
    }
}



double dfunc_2hs4(double r)
{
    if (fabs(r) <=1.){
        return - pow(fabs(r),5.)/12. + pow(r,4.)/4. - pow(r,2.)/2. + 11./20.;
    } else if (fabs(r) >= 1. && fabs(r) <=2.) {
        return pow(fabs(r),5.)/24. - (3.*pow(r,4.))/8. + (5.*pow(fabs(r),3.))/4. - (7.*pow(r,2.))/4. + (5.*fabs(r))/8. + 17./40.;
    } else if (fabs(r) >= 2. && fabs(r) <=3.) {
        return -pow((fabs(r) - 3.),5.)/120.;
    } else {
        return 0.0;
    }
}


double dfunc_2hs5(double r)
{
    if (fabs(r) <=0.5){
        return - pow(r,6.)/36. + (7.*pow(r,4.))/48. - (77.*pow(r,2.))/192. + 5887./11520.;
    } else if (fabs(r) >= 0.5 && fabs(r) <=1.5) {
        return pow(r,6.)/48. - (7.*pow(fabs(r),5.))/48. + (21.*pow(r,4.))/64. - (35.*pow(fabs(r),3.))/288. - (91.*pow(r,2.))/256. - (7.*fabs(r))/768. + 7861./15360.;
    } else if (fabs(r) >= 1.5 && fabs(r) <=2.5) {
        return - pow(r,6.)/120. + (7.*pow(fabs(r),5.))/60. - (21.*pow(r,4.))/32. + (133.*pow(fabs(r),3.))/72. - (329.*pow(r,2.))/128. + (1267.*fabs(r))/960. + 1379./7680.;
    } else if (fabs(r) >= 2.5 && fabs(r) <=3.5) {
        return pow((2.*fabs(r) - 7.),6.)/46080.;
    } else {
        return 0.0;
    }
}


double dfunc_2hs6(double r)
{
    if (fabs(r) <=1.){
        return pow(fabs(r),7.)/144. - pow(r,6.)/36. + pow(r,4.)/9. - pow(r,2.)/3. + 151./315.;
    } else if (fabs(r) >= 1. && fabs(r) <=2.) {
        return - pow(fabs(r),7.)/240. + pow(r,6.)/20. - (7.*pow(fabs(r),5.))/30. + pow(r,4.)/2. - (7.*pow(fabs(r),3.))/18. - pow(r,2.)/10. - (7.*fabs(r))/90. + 103./210.;
    } else if (fabs(r) >= 2. && fabs(r) <=3.) {
        return pow(fabs(r),7.)/720. - pow(r,6.)/36. + (7.*pow(fabs(r),5.))/30. - (19.*pow(r,4.))/18. + (49.*pow(fabs(r),3.))/18. - (23.*pow(r,2.))/6. + (217.*fabs(r))/90. - 139./630.;
    } else if (fabs(r) >= 3. && fabs(r) <=4.) {
        return -pow((fabs(r) - 4.),7.)/5040.;
    } else {
        return 0.0;
    }
}


double dfunc_4hs1(double r)
{
    if (fabs(r) <=0.5){
        return 15./64. - pow(r,2)/16.;
    } else if (fabs(r) >= 0.5 && fabs(r) <=3.5) {
        return 1.0/4.0 - fabs(r)/16.0;
    } else if (fabs(r) >= 3.5 && fabs(r) <=4.5) {
        return pow((2.*fabs(r) - 9.),2.)/128.;
    } else {
        return 0.0;
    }
}


double dfunc_nhs2(double r, double n)
{
    if (fabs(r) <= n-1.){
        return 1./(2.*n);
    } else if (fabs(r) >= n-1. && fabs(r) <= n) {
        return fabs(r) - n/2. - (2.*pow(r,2.) + 3.*fabs(r) - 1.)/(4.*n) + 3./4.;
    } else if (fabs(r) >= n && fabs(r) <= n+1.) {
        return ((n - fabs(r) + 1.)*(2.*n - 2.*fabs(r) + 1.))/(4.*n);
    } else {
        return 0.0;
    }

}



double dfunc_nhs1(double r, double n)
{
    if (fabs(r) <= n-0.5){
        return 1./(2.*n);
    } else if (fabs(r) >= n-0.5 && fabs(r) <= n+0.5) {
        return 1./2. - (fabs(r) - 1./2.)/(2.*n);
    } else {
        return 0.0;
    }

}


double dfunc_h(double r)
{
    if (fabs(r) < 0.5) {
        return 1.0;
    } else {
        return 0.0;
    }
}


double dfunc_4h2peak(double r)
{
    if (fabs(r)>2.0) {
        return 0.0;
    } else if (r>=-2.0 && r<=0.0) {
        return 0.5*(1.0-fabs(r+1.0));
    } else if (r>=0.0 && r<=2.0) {
        return 0.5*(1.0-fabs(r-1.0));
    }
}


// get coordinates of the interpolation points in the wall normal direction
PetscErrorCode Coordinates_IP(IBMNodes *ibm, int NumberOfObjects)
{

    PetscInt        i, j, k, l, ibi;

    int n1e, n2e, n3e;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {

        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {

            //double dh = 1.5*fmax(fmax(ibm[ibi].dhx[l], ibm[ibi].dhy[l]), ibm[ibi].dhz[l]);
            double dh = ibm[ibi].dh;

            //double dh = pow(ibm[ibi].dhx[l]*ibm[ibi].dhy[l]*ibm[ibi].dhz[l],1/3); 

            ibm[ibi].dh_IP[l] = dh;

            ibm[ibi].centIP_x[l] = ibm[ibi].cent_x[l] + dh*ibm[ibi].nf_x[l];
            ibm[ibi].centIP_y[l] = ibm[ibi].cent_y[l] + dh*ibm[ibi].nf_y[l];
            ibm[ibi].centIP_z[l] = ibm[ibi].cent_z[l] + dh*ibm[ibi].nf_z[l];

        }

        if (!rank) 
        {
            FILE *f;
            char filen[80];
            sprintf(filen, "IPNodes_%2.2d", ibi);
            f = fopen(filen, "w");
            PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x, y, z\n");

            for (l=0; l<ibm[ibi].n_elmt; l++) PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", ibm[ibi].centIP_x[l], ibm[ibi].centIP_y[l], ibm[ibi].centIP_z[l]);
            fclose(f);
        }

    }

    return(0);
}


PetscErrorCode Pre_process_IP(UserCtx *user, IBMNodes *ibm, int NumberOfObjects)
{
    DA              da = user->da, fda = user->fda;
    DALocalInfo     info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt        lxs, lxe, lys, lye, lzs, lze;

    Cmpnts          ***coor, ***csi, ***eta, ***zet;

    PetscReal       ***aj, ***nvert;

    Vec         Coor;


    DAGetLocalInfo(da, &info);
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
PetscPrintf(PETSC_COMM_WORLD, "PreProcess IP 1\n");

    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);

PetscPrintf(PETSC_COMM_WORLD, "PreProcess IP 2\n");


    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);
PetscPrintf(PETSC_COMM_WORLD, "PreProcess IP 3\n");

    DAVecGetArray(da, user->lNvert_4diffusedinterface, &nvert);

PetscPrintf(PETSC_COMM_WORLD, "PreProcess IP 4\n");


    int n1e, n2e, n3e;

    PetscPrintf(PETSC_COMM_WORLD, "PreProcess IP first loop\n");
    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {

        std::vector<int> iclose (ibm[ibi].n_elmt);
        std::vector<int> jclose (ibm[ibi].n_elmt);
        std::vector<int> kclose (ibm[ibi].n_elmt);
        std::vector<int> sum_iclose (ibm[ibi].n_elmt);
        std::vector<int> sum_jclose (ibm[ibi].n_elmt);
        std::vector<int> sum_kclose (ibm[ibi].n_elmt);

        std::vector<int> count (ibm[ibi].n_elmt);
        std::vector<int> sum_count (ibm[ibi].n_elmt);

        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            iclose[l]=0;
            jclose[l]=0;
            kclose[l]=0;

            sum_iclose[l]=0;
            sum_jclose[l]=0;
            sum_kclose[l]=0;
        }

        //int count[ibm[ibi].n_elmt], sum_count[ibm[ibi].n_elmt];
        PetscPrintf(PETSC_COMM_WORLD, "PreProcess IP second loop\n");
        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {

            int imark=0;
            double dmin=1.e4;
            int ic,jc,kc;
            
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        double r1=ibm[ibi].centIP_x[l]-coor[k][j][i].x, r2=ibm[ibi].centIP_y[l]-coor[k][j][i].y, r3=ibm[ibi].centIP_z[l]-coor[k][j][i].z; 
                        double d1=sqrt(r1*r1+r2*r2+r3*r3);
                        if (d1<dmin)
                        {
                            dmin=d1;
                            ic=i;jc=j;kc=k;
                            count[l]=1;	
                        }
                    }
                }
            }

            i=ic; j=jc; k=kc;
            double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
            double dhx=1.0/aj[k][j][i]/area;

            area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
            double dhy=1.0/aj[k][j][i]/area;

            area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
            double dhz=1.0/aj[k][j][i]/area;	

            double dh_grid = 2.0*sqrt(dhx*dhx+dhy*dhy+dhz*dhz);

            double diff=dmin-dh_grid;

            if (diff>1.e-4) 
            {
                count[l]=0;	
                ic=0; jc=0; kc=0;
                iclose[l]=0; jclose[l]=0; kclose[l]=0;
            } else 
            {
                iclose[l]=ic; jclose[l]=jc; kclose[l]=kc;
                i=ic; j=jc; k=kc;
            }
        }

        PetscBarrier(PETSC_NULL);
        PetscPrintf(PETSC_COMM_WORLD, "MPI_ALLREDUCE\n");

        //PetscBarrier(PETSC_NULL);
        MPI_Allreduce (&iclose[0], &sum_iclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        //PetscBarrier(PETSC_NULL);
        MPI_Allreduce (&jclose[0], &sum_jclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        //PetscBarrier(PETSC_NULL);
        MPI_Allreduce (&kclose[0], &sum_kclose[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        //PetscBarrier(PETSC_NULL);

        MPI_Allreduce (&count[0], &sum_count[0], ibm[ibi].n_elmt, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {          
            iclose[l]=sum_iclose[l]/(sum_count[l]);
            jclose[l]=sum_jclose[l]/(sum_count[l]);
            kclose[l]=sum_kclose[l]/(sum_count[l]);
        }

        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {

            int ic=iclose[l];
            int jc=jclose[l];
            int kc=kclose[l];

            int ic1=ic-Itpwidth, ic2=ic+Itpwidth; 
            if (ic1>lxe||ic2<lxs) 
            {
                ibm[ibi].iIP_min[l]=lxs;
                ibm[ibi].iIP_max[l]=lxs;
            } else 
            {
                ibm[ibi].iIP_min[l]=PetscMax(ic1, lxs);
                ibm[ibi].iIP_max[l]=PetscMin(ic2, lxe);
            }

            int jc1=jc-Itpwidth, jc2=jc+Itpwidth; 
            if (jc1>lye||jc2<lys) 
            {
                ibm[ibi].jIP_min[l]=lys;
                ibm[ibi].jIP_max[l]=lys;
            } else 
            {
                ibm[ibi].jIP_min[l]=PetscMax(jc1, lys);
                ibm[ibi].jIP_max[l]=PetscMin(jc2, lye);
            }
            int kc1=kc-Itpwidth, kc2=kc+Itpwidth; 

            if (kc1>lze||kc2<lzs) 
            {
                ibm[ibi].kIP_min[l]=lzs;
                ibm[ibi].kIP_max[l]=lzs;
                //PetscPrintf(PETSC_COMM_WORLD, "Debug-B Set kmin and kmax wotj lzs = %i and lze = %i\n", lzs, lze);
            } else 
            {
                ibm[ibi].kIP_min[l]=PetscMax(kc1, lzs);
                ibm[ibi].kIP_max[l]=PetscMin(kc2, lze);
            }
        }
    }

    DAVecRestoreArray(fda, Coor, &coor);
    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);
    DAVecRestoreArray(da, user->lNvert_4diffusedinterface, &nvert);

    return(0);
}


/* Interpolate the velocity at the Lagrangian Interpolation points */
PetscErrorCode Calc_U_IPlagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    DA              da = user->da, fda = user->fda;
    DALocalInfo     info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt    lxs, lxe, lys, lye, lzs, lze;

    Cmpnts  ***ucat, ***coor, ***csi, ***eta, ***zet;

    PetscReal   ***aj, ***p;

    Vec     Coor;

    PetscReal   dfunc;

    PetscReal   xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

    double r1, r2, r3;

    DAGetLocalInfo(da, &info);
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

    //PetscPrintf(PETSC_COMM_WORLD, "Calc_U_IPlagr In\n");

    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);
    DAVecGetArray(fda, user->lUcat, &ucat);
    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);

    DAVecGetArray(da, user->lP, &p);

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            ibm[ibi].U_IPlagr_x[l] = 0.0;
            ibm[ibi].U_IPlagr_y[l] = 0.0;
            ibm[ibi].U_IPlagr_z[l] = 0.0;
            ibm[ibi].P_IPlagr[l] = 0.0;
        }
    }

    int my_rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);

    clock_t start, end;
    double elapsed;

    start = clock();
    double ni[3], nj[3], nk[3];

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            for (k = ibm[ibi].kIP_min[l]; k<ibm[ibi].kIP_max[l]; k++) 
            for (j = ibm[ibi].jIP_min[l]; j<ibm[ibi].jIP_max[l]; j++) 
            for (i = ibm[ibi].iIP_min[l]; i<ibm[ibi].iIP_max[l]; i++) 
            {
                xc = 0.125 *(coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
                coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
                coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
                
                yc = 0.125 *(coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
                coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
                coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
                
                zc = 0.125 *(coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
                coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
                coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);

                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                double dhx_=1.0/aj[k][j][i]/area;

                area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                double dhy_=1.0/aj[k][j][i]/area;

                area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                double dhz_=1.0/aj[k][j][i]/area;

                double rx= (xc - ibm[ibi].centIP_x[l]), ry = (yc - ibm[ibi].centIP_y[l]), rz = (zc - ibm[ibi].centIP_z[l]);

                Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                r1=(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
                r2=(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
                r3=(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 

                dfunc =  dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);
                ibm[ibi].U_IPlagr_x[l] += ucat[k][j][i].x * dfunc;
                ibm[ibi].U_IPlagr_y[l] += ucat[k][j][i].y * dfunc;
                ibm[ibi].U_IPlagr_z[l] += ucat[k][j][i].z * dfunc;          

                ibm[ibi].P_IPlagr[l] += p[k][j][i] * dfunc;          
            }
        }
    }
    start = clock();
    for (ibi=0; ibi<NumberOfObjects; ibi++)  
    {
        std::vector<Cmpnts> u_local (ibm[ibi].n_elmt); 
        std::vector<Cmpnts> u_sum (ibm[ibi].n_elmt); 
        std::vector<double> p_local (ibm[ibi].n_elmt); 
        std::vector<double> p_sum (ibm[ibi].n_elmt); 

        for (i=0; i<ibm[ibi].n_elmt; i++ ) 
        {
            u_local[i].x = ibm[ibi].U_IPlagr_x[i];
            u_local[i].y = ibm[ibi].U_IPlagr_y[i];
            u_local[i].z = ibm[ibi].U_IPlagr_z[i];
        }

        int count_vector=3*ibm[ibi].n_elmt;

        //PetscBarrier(PETSC_NULL);
        MPI_Allreduce( &(u_local[0]), &(u_sum[0]), count_vector, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        //MPI_Allreduce( &u_local, &u_sum, count_vector, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        //PetscBarrier(PETSC_NULL);

        for (i=0; i<ibm[ibi].n_elmt; i++ ) 
        {
            MPI_Allreduce( &u_local[i].x, &u_sum[i].x, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            MPI_Allreduce( &u_local[i].y, &u_sum[i].y, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            MPI_Allreduce( &u_local[i].z, &u_sum[i].z, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            ibm[ibi].U_IPlagr_x[i] = u_sum[i].x;
            ibm[ibi].U_IPlagr_y[i] = u_sum[i].y;
            ibm[ibi].U_IPlagr_z[i] = u_sum[i].z;
        }

        //PetscBarrier(PETSC_NULL);

        for (i=0; i<ibm[ibi].n_elmt; i++ ) 
        {
            p_local[i] = ibm[ibi].P_IPlagr[i];
        }

        //PetscBarrier(PETSC_NULL);
        MPI_Allreduce( &(p_local[0]), &(p_sum[0]), ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        //PetscBarrier(PETSC_NULL);

        for (i=0; i<ibm[ibi].n_elmt; i++ ) 
        {
            ibm[ibi].P_IPlagr[i] = p_sum[i];
        }

        //PetscBarrier(PETSC_NULL);
    }

    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;

    //PetscPrintf(PETSC_COMM_WORLD, "Time for U_larg global sum %le \n", elapsed);
    //PetscBarrier(PETSC_NULL);

    if (ii_periodicWT || jj_periodicWT || kk_periodicWT) 
    {
        double xc_min = 1.0e6;        
        double yc_min = 1.0e6;
        double zc_min = 1.0e6;

        for (ibi=0; ibi<NumberOfObjects; ibi++) 
        {
            xc_min = PetscMin(xc_min, fsi[ibi].x_c);
            yc_min = PetscMin(yc_min, fsi[ibi].y_c);
            zc_min = PetscMin(zc_min, fsi[ibi].z_c);
        }

        //PetscPrintf(PETSC_COMM_WORLD, "HERE 1 \n");
        int IndWT[Nz_WT][Ny_WT][Nx_WT];

        double fac_x=1.0/Sx_WT, fac_y=1.0/Sy_WT, fac_z=1.0/Sz_WT;
        for (ibi=0; ibi<NumberOfObjects; ibi++) 
        {
            int ii=(int)((fsi[ibi].x_c-xc_min+1.e-9)*fac_x);
            int jj=(int)((fsi[ibi].y_c-yc_min+1.e-9)*fac_y);
            int kk=(int)((fsi[ibi].z_c-zc_min+1.e-9)*fac_z);
            PetscPrintf(PETSC_COMM_WORLD, "ibi ii jj kk %d %d %d %d \n", ibi, ii, jj, kk);
            IndWT[kk][jj][ii]=ibi;
        }

        int i, j, k;
    }

    if (MoveFrame) 
    {
        for (ibi=0; ibi<NumberOfObjects; ibi++)
        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            ibm[ibi].U_IPlagr_x[l] += u_frame;
            ibm[ibi].U_IPlagr_y[l] += v_frame;
            ibm[ibi].U_IPlagr_z[l] += w_frame;
        }
    }

    DAVecRestoreArray(fda, Coor, &coor);
    DAVecRestoreArray(fda, user->lUcat, &ucat);
    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);
    DAVecRestoreArray(da,  user->lP,  &p);

    return(0);
}

// Hossein temporarily added from 3.4 code to solve the restart rotation issue
// Hard Coding !!!!!
// coloring each IB nodes, color=1000, no force on the IB nodes
PetscErrorCode ColorIB_tmp(IBMNodes *ibm, int NumberOfObjects)
{

  	PetscInt        i, j, k, l, ibi;


  	for (ibi=0; ibi<NumberOfObjects; ibi++) {

	        for (l=0; l<ibm[ibi].n_elmt; l++) {

			//if (fabs(ibm[ibi].cent_z[l]-(-0.03051))<1.e-3 && fabs(fabs(ibm[ibi].nf_z[l])-1.0)<1.e-4) ibm[ibi].color[l]=1000;
			//if (fabs(ibm[ibi].cent_z[l]-(1.5-0.03051))<1.e-3 && fabs(fabs(ibm[ibi].nf_z[l])-1.0)<1.e-4) ibm[ibi].color[l]=1000;
			//if (fabs(ibm[ibi].cent_z[l]-0.0352)<1.e-3 && fabs(fabs(ibm[ibi].nf_z[l])-1.0)<1.e-4) ibm[ibi].color[l]=1000;
//			if (fabs(ibm[ibi].nf_y[l])<1.e-3) ibm[ibi].color[l]=1000;
//			if (fabs(ibm[ibi].cent_z[l]-(1.4695))<1.e-3 && fabs(fabs(ibm[ibi].nf_z[l])-1.0)<1.e-4) ibm[ibi].color[l]=1000;

		}
	}



  	return(0);

}

/* Interpolate the eddy viscosity at the Lagrangian points */
PetscErrorCode Intp_Nut_lagr(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{

    DA              da = user->da, fda = user->fda;
    DALocalInfo     info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    Cmpnts	***ucat, ***coor, ***csi, ***eta, ***zet;

    PetscReal 	***aj;

    Vec		Coor;

    PetscReal	dfunc;

    PetscReal	xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

    double r1, r2, r3;

    DAGetLocalInfo(da, &info);
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

    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);
    DAVecGetArray(fda, user->lUcat, &ucat);
    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);

    PetscReal ***lnu_t;

    if(les) { DAVecGetArray(da, user->lNu_t, &lnu_t); //Hossein-added for rans to fix nacelle model and rans problem
    } else if(rans) {
               DAVecGetArray(da, user->lNu_t, &lnu_t);
    }

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            ibm[ibi].Nut_lagr[l] = 0.0;
        }
    }


    clock_t start, end;
    double elapsed;

    start = clock();
    double ni[3], nj[3], nk[3];

    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (l=0; l<ibm[ibi].n_elmt; l++) {

            for (k = ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
            for (j = ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
            for (i = ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) 
            {

                xc = 0.125 *(coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                    coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
                    coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
                    coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
                    
                yc = 0.125 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                    coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
                    coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
                    coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
                    
                zc = 0.125 *(coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                    coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
                    coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
                    coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);

                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                double dhx_=1.0/aj[k][j][i]/area;

                area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                double dhy_=1.0/aj[k][j][i]/area;

                area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                double dhz_=1.0/aj[k][j][i]/area;	

                double rx= (xc - ibm[ibi].cent_x[l]), ry = (yc - ibm[ibi].cent_y[l]), rz = (zc - ibm[ibi].cent_z[l]);

                Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                r1=(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
                r2=(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
                r3=(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 

                dfunc = dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);

                ibm[ibi].Nut_lagr[l] += lnu_t[k][j][i] * dfunc;

            }

        }
    }

    start = clock();
    for (ibi=0; ibi<NumberOfObjects; ibi++)  
    {

        //PetscBarrier(PETSC_NULL);
        std::vector<double> nut_local (ibm[ibi].n_elmt); 
        std::vector<double> nut_sum (ibm[ibi].n_elmt); 

        for (i=0; i<ibm[ibi].n_elmt; i++ ) 
        {
            nut_local[i] = ibm[ibi].Nut_lagr[i];
        }

        //PetscBarrier(PETSC_NULL);
        MPI_Allreduce( &(nut_local[0]), &(nut_sum[0]), ibm[ibi].n_elmt, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        //PetscBarrier(PETSC_NULL);

        for (i=0; i<ibm[ibi].n_elmt; i++ ) 
        {
            ibm[ibi].Nut_lagr[i] = nut_sum[i];
        }

        //PetscBarrier(PETSC_NULL);
    }

    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;

    //PetscPrintf(PETSC_COMM_WORLD, "Time for U_larg global sum %le \n", elapsed);

    DAVecRestoreArray(fda, Coor, &coor);
    DAVecRestoreArray(fda, user->lUcat, &ucat);
    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);

    if(les) { DAVecRestoreArray(da, user->lNu_t, &lnu_t); //Hossein-added for rans to fix nacelle model and rans problem
    } else if(rans) {
                DAVecRestoreArray(da, user->lNu_t, &lnu_t);
    }

    return(0);
}


// distribute the force on the turbines to background grid.
//
// called after Calc_F_lagr_ACL or Calc_F_lagr before solving the momentum equatiion
PetscErrorCode Calc_Nut_eul(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, PetscInt NumberOfObjects, double dh, int df)
{

    DA              da = user->da, fda = user->fda;
    DALocalInfo     info;
    PetscInt        xs, xe, ys, ye, zs, ze; // Local grid information
    PetscInt        mx, my, mz; // Dimensions in three directions
    PetscInt        i, j, k, l, ibi;
    PetscInt        lxs, lxe, lys, lye, lzs, lze;

    Cmpnts          ***coor, ***csi, ***eta, ***zet;

    PetscReal       ***aj, ***nvert;

    PetscReal       ***lnut_eul;

    Vec         Coor;

    PetscReal   xc, yc, zc, hx, hy, hz, vol_eul, vol_lagr;

    double      dfunc;

    double      r1, r2, r3;

    DAGetLocalInfo(da, &info);
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

    DAGetGhostedCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);
    DAVecGetArray(da, user->lNut_eul, &lnut_eul);
    DAVecGetArray(fda, user->lCsi,  &csi);
    DAVecGetArray(fda, user->lEta,  &eta);
    DAVecGetArray(fda, user->lZet,  &zet);
    DAVecGetArray(da,  user->lAj,  &aj);

    DAVecGetArray(da, user->lNvert, &nvert);

    double ni[3], nj[3], nk[3];
    for (ibi=0; ibi<NumberOfObjects; ibi++) 
    {
        for (l=0; l<ibm[ibi].n_elmt; l++) 
        {
            for (k=ibm[ibi].k_min[l]; k<ibm[ibi].k_max[l]; k++) 
            for (j=ibm[ibi].j_min[l]; j<ibm[ibi].j_max[l]; j++) 
            for (i=ibm[ibi].i_min[l]; i<ibm[ibi].i_max[l]; i++) 
            {
                xc=coor[k][j][i].x;
                yc=coor[k][j][i].y;
                zc=coor[k][j][i].z;
                double rx= (xc - ibm[ibi].cent_x[l]), ry = (yc - ibm[ibi].cent_y[l]), rz = (zc - ibm[ibi].cent_z[l]);

                Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                double dhx_=1.0/aj[k][j][i]/area;

                area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                double dhy_=1.0/aj[k][j][i]/area;

                area = sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );
                double dhz_=1.0/aj[k][j][i]/area;	

                double r1=fabs(rx*ni[0]+ry*ni[1]+rz*ni[2])/dhx_; 
                double r2=fabs(rx*nj[0]+ry*nj[1]+rz*nj[2])/dhy_; 
                double r3=fabs(rx*nk[0]+ry*nk[1]+rz*nk[2])/dhz_; 

                vol_eul = 1;

                double dfunc;
                if (df == 0) 
                {
                    dfunc = vol_eul * dfunc_2h(r1) * dfunc_2h(r2) * dfunc_2h(r3);
                } else if (df == 7) 
                {
                    double n = halfwidth_dfunc;
                    dfunc = vol_eul * dfunc_exp(r1,n) * dfunc_exp(r2,n) * dfunc_exp(r3,n);
                } else 
                { 
                    dfunc = vol_eul * dfunc_sc4h(r1) * dfunc_sc4h(r2) * dfunc_sc4h(r3);
                }

                lnut_eul[k][j][i] += ibm[ibi].Nut_lagr[l] * dfunc;

            }
        }
    }

    for (k=lzs;k<lze;k++) 
    for (j=lys;j<lye;j++) 
    for (i=lxs;i<lxe;i++) 
    {
        int ii, jj, kk;
        double _nvert;
        _nvert = 0.0;

        for (kk=k-1;kk<k+2;kk++) 
        for (jj=j-1;jj<j+2;jj++) 
        for (ii=i-1;ii<i+2;ii++) 
        {
            _nvert += nvert[kk][jj][ii];
        }

        //_nvert = nvert[k][j][i];
        if ( _nvert >2.9 ) 
        { 
                        lnut_eul[k][j][i]=0.0;
        }

        if (i==1) {
            lnut_eul[k][j][i-1]=0.0;
        }

        if (j==1) {
            lnut_eul[k][j-1][i]=0.0;
        }

        if (k==1) {
            lnut_eul[k-1][j][i]=0.0;
        }

        if (i==mx-2) {
            lnut_eul[k][j][i+1]=0.0;
        }

        if (j==my-2) {
            lnut_eul[k][j+1][i]=0.0;
        }

        if (k==mz-2) {
            lnut_eul[k+1][j][i]=0.0;
        }

    }

    DAVecRestoreArray(fda, Coor, &coor);
    DAVecRestoreArray(da, user->lNut_eul, &lnut_eul);
    DAVecRestoreArray(fda, user->lCsi,  &csi);
    DAVecRestoreArray(fda, user->lEta,  &eta);
    DAVecRestoreArray(fda, user->lZet,  &zet);
    DAVecRestoreArray(da,  user->lAj,  &aj);
    DAVecRestoreArray(da, user->lNvert, &nvert);
    DALocalToLocalBegin(da, user->lNut_eul, INSERT_VALUES, user->lNut_eul);
    DALocalToLocalEnd(da, user->lNut_eul, INSERT_VALUES, user->lNut_eul);
    DALocalToGlobal(da, user->lNut_eul, INSERT_VALUES, user->Nut_eul);

    if (les) 
    {

        PetscReal ***lnu_t;
        DAVecGetArray(da, user->lNut_eul, &lnut_eul);
        DAVecGetArray(da, user->lNu_t, &lnu_t);

        for(k=zs; k<ze; k++)
        for(j=ys; j<ye; j++)
        for(i=xs; i<xe; i++) 
        {
            lnu_t[k][j][i] += lnut_eul[k][j][i];
        }

        DAVecRestoreArray(da, user->lNut_eul, &lnut_eul);
        DAVecRestoreArray(da, user->lNu_t, &lnu_t);

        DALocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
        DALocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
    }

    return(0);
}


PetscErrorCode nacelleYaw_IB(UserCtx *user, IBMNodes *ibm, FSInfo *fsi, int NumberOfObjects)
{
    // Compute rotation matrix using Rodrigues equation

    for(int ibi=0; ibi<NumberOfObjects; ibi++)
    {

        double nx=fsi[ibi].nx_tb; 
        double ny=fsi[ibi].ny_tb; 
        double nz=fsi[ibi].nz_tb; 

        double nxo=fsi[ibi].nx_tbo; 
        double nyo=fsi[ibi].ny_tbo; 
        double nzo=fsi[ibi].nz_tbo; 

        if(nx==nxo){
          return(0);
        }

        double rr = sqrt(nx*nx+ny*ny+nz*nz);
        nx = nx/rr;
        ny = ny/rr;
        nz = nz/rr;

        rr = sqrt(nxo*nxo+nyo*nyo+nzo*nzo);
        nxo = nxo/rr;
        nyo = nyo/rr;
        nzo = nzo/rr;

        double vx = nyo*nz-nzo*ny;
        double vy = nzo*nx-nxo*nz;
        double vz = nxo*ny-nyo*nx;

        double snth = sqrt(vx*vx+vy*vy+vz*vz);

        double csth = nxo*nx+nyo*ny+nzo*nz;
        double icsth = 1./(1.+csth);

        double v11 = 0.0; double v12 = -vz; double v13 =  vy;
        double v21 =  vz; double v22 = 0.0; double v23 = -vx;
        double v31 = -vy; double v32 =  vx; double v33 = 0.0;

        double vv11 = v11*v11+v12*v21+v13*v31; 
        double vv12 = v11*v12+v12*v22+v13*v32; 
        double vv13 = v11*v13+v12*v23+v13*v33;
        double vv21 = v21*v11+v22*v21+v23*v31; 
        double vv22 = v21*v12+v22*v22+v23*v32; 
        double vv23 = v21*v13+v22*v23+v23*v33; 
        double vv31 = v31*v11+v32*v21+v33*v31; 
        double vv32 = v31*v12+v32*v22+v33*v32; 
        double vv33 = v31*v13+v32*v23+v33*v33; 

        double R11 = 1.0 + v11 + vv11*icsth;
        double R12 =       v12 + vv12*icsth;
        double R13 =       v13 + vv13*icsth;
        double R21 =       v21 + vv21*icsth;
        double R22 = 1.0 + v22 + vv22*icsth;
        double R23 =       v23 + vv23*icsth;
        double R31 =       v31 + vv31*icsth;
        double R32 =       v32 + vv32*icsth;
        double R33 = 1.0 + v33 + vv33*icsth;

        for (int l=0; l<ibm[ibi].n_elmt; l++) 
        {
            double x_c=fsi[ibi].x_c;
            double y_c=fsi[ibi].y_c;
            double z_c=fsi[ibi].z_c;
            double xcent = R11*(ibm[ibi].cent_x[l]-x_c) + R12*(ibm[ibi].cent_y[l]-y_c) + R13*(ibm[ibi].cent_z[l]-z_c);
            double ycent = R21*(ibm[ibi].cent_x[l]-x_c) + R22*(ibm[ibi].cent_y[l]-y_c) + R23*(ibm[ibi].cent_z[l]-z_c);
            double zcent = R31*(ibm[ibi].cent_x[l]-x_c) + R32*(ibm[ibi].cent_y[l]-y_c) + R33*(ibm[ibi].cent_z[l]-z_c);

            //PetscPrintf(PETSC_COMM_WORLD, "cent_x %le %le %i \n",xcent,ibm[ibi].cent_x[l],l);

            ibm[ibi].cent_x[l] = xcent+x_c;
            ibm[ibi].cent_y[l] = ycent+y_c;
            ibm[ibi].cent_z[l] = zcent+z_c;
            //PetscPrintf(PETSC_COMM_WORLD, "cent_x %le %i \n",xcent,l);
        }

        for (int l=0; l<ibm[ibi].n_elmt; l++) 
        {
            double xnf = R11*ibm[ibi].nf_x[l] + R12*ibm[ibi].nf_y[l] + R13*ibm[ibi].nf_z[l];
            double ynf = R21*ibm[ibi].nf_x[l] + R22*ibm[ibi].nf_y[l] + R23*ibm[ibi].nf_z[l];
            double znf = R31*ibm[ibi].nf_x[l] + R32*ibm[ibi].nf_y[l] + R33*ibm[ibi].nf_z[l];
            ibm[ibi].nf_x[l] = xnf;
            ibm[ibi].nf_y[l] = ynf;
            ibm[ibi].nf_z[l] = znf;
            //     PetscPrintf(PETSC_COMM_WORLD, "xnf_x %le %i \n",xnf,l);
        }
    }

    Pre_process(&(user[0]), ibm, NumberOfObjects);

    Coordinates_IP(ibm, NumberOfObjects);

    //There is a here in this subroutine
    Pre_process_IP(&(user[0]), ibm, NumberOfObjects);

    return(0);
}    


PetscErrorCode rotorYaw(IBMNodes *ibm, FSInfo *fsi)
{
    // Compute rotation matrix using Rodrigues equation

    double nx=fsi->nx_tb; 
    double ny=fsi->ny_tb; 
    double nz=fsi->nz_tb; 

    double nxo=fsi->nx_tbo; 
    double nyo=fsi->ny_tbo; 
    double nzo=fsi->nz_tbo; 

    if(nx==nxo){
        return(0);
    }

    //  PetscPrintf(PETSC_COMM_WORLD, "Rotor Yaw %le %le %le \n",nx,ny,nz);
    //  PetscPrintf(PETSC_COMM_WORLD, "Rotor Yaw Old %le %le %le \n",nxo,nyo,nzo);

    double rr = sqrt(nx*nx+ny*ny+nz*nz);
    nx = nx/rr;
    ny = ny/rr;
    nz = nz/rr;

    rr = sqrt(nxo*nxo+nyo*nyo+nzo*nzo);
    nxo = nxo/rr;
    nyo = nyo/rr;
    nzo = nzo/rr;

    if (nzo==nz) return(0);

    double vx = nyo*nz-nzo*ny;
    double vy = nzo*nx-nxo*nz;
    double vz = nxo*ny-nyo*nx;

    double snth = sqrt(vx*vx+vy*vy+vz*vz);

    double csth = nxo*nx+nyo*ny+nzo*nz;
    double icsth = 1./(1.+csth);

    double v11 = 0.0; double v12 = -vz; double v13 =  vy;
    double v21 =  vz; double v22 = 0.0; double v23 = -vx;
    double v31 = -vy; double v32 =  vx; double v33 = 0.0;

    double vv11 = v11*v11+v12*v21+v13*v31; 
    double vv12 = v11*v12+v12*v22+v13*v32; 
    double vv13 = v11*v13+v12*v23+v13*v33;
    double vv21 = v21*v11+v22*v21+v23*v31; 
    double vv22 = v21*v12+v22*v22+v23*v32; 
    double vv23 = v21*v13+v22*v23+v23*v33; 
    double vv31 = v31*v11+v32*v21+v33*v31; 
    double vv32 = v31*v12+v32*v22+v33*v32; 
    double vv33 = v31*v13+v32*v23+v33*v33; 

    double R11 = 1.0 + v11 + vv11*icsth;
    double R12 =       v12 + vv12*icsth;
    double R13 =       v13 + vv13*icsth;
    double R21 =       v21 + vv21*icsth;
    double R22 = 1.0 + v22 + vv22*icsth;
    double R23 =       v23 + vv23*icsth;
    double R31 =       v31 + vv31*icsth;
    double R32 =       v32 + vv32*icsth;
    double R33 = 1.0 + v33 + vv33*icsth;

    if(ti==tistart && rstart_turbinerotation)
    {
        for (int l=0; l<ibm->n_v; l++) 
        {
            double x_c=fsi->x_c;
            double y_c=fsi->y_c;
            double z_c=fsi->z_c;
            double xbp = R11*(ibm->x_bp0[l]-x_c) + R12*(ibm->y_bp0[l]-y_c) + R13*(ibm->z_bp0[l]-z_c);
            double ybp = R21*(ibm->x_bp0[l]-x_c) + R22*(ibm->y_bp0[l]-y_c) + R23*(ibm->z_bp0[l]-z_c);
            double zbp = R31*(ibm->x_bp0[l]-x_c) + R32*(ibm->y_bp0[l]-y_c) + R33*(ibm->z_bp0[l]-z_c);

            //PetscPrintf(PETSC_COMM_WORLD, "cent_x %le %le %i \n",xcent,ibm[ibi].cent_x[l],l);
            ibm->x_bp[l] = xbp+x_c;
            ibm->y_bp[l] = ybp+y_c;
            ibm->z_bp[l] = zbp+z_c;
            //PetscPrintf(PETSC_COMM_WORLD, "cent_x %le %i \n",xcent,l);
        }
    } else 
    {
        for (int l=0; l<ibm->n_v; l++) 
        {
          double x_c=fsi->x_c;
          double y_c=fsi->y_c;
          double z_c=fsi->z_c;
          double xbp = R11*(ibm->x_bp[l]-x_c) + R12*(ibm->y_bp[l]-y_c) + R13*(ibm->z_bp[l]-z_c);
          double ybp = R21*(ibm->x_bp[l]-x_c) + R22*(ibm->y_bp[l]-y_c) + R23*(ibm->z_bp[l]-z_c);
          double zbp = R31*(ibm->x_bp[l]-x_c) + R32*(ibm->y_bp[l]-y_c) + R33*(ibm->z_bp[l]-z_c);

        //    PetscPrintf(PETSC_COMM_WORLD, "cent_x %le %le %i \n",xcent,ibm[ibi].cent_x[l],l);
          ibm->x_bp[l] = xbp+x_c;
          ibm->y_bp[l] = ybp+y_c;
          ibm->z_bp[l] = zbp+z_c;
        //     PetscPrintf(PETSC_COMM_WORLD, "cent_x %le %i \n",xcent,l);
        }

    }

    fsi->nx_tbo=fsi->nx_tb; 
    fsi->ny_tbo=fsi->ny_tb; 
    fsi->nz_tbo=fsi->nz_tb; 

    return(0);

}
