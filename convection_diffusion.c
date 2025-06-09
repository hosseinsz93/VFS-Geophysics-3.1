#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern PetscInt immersed, NumberOfBodies, ti, tistart, wallfunction, y_direction;
extern PetscErrorCode MyKSPMonitor1(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy);
extern PetscInt inletprofile;
extern PetscInt conv_diff, zero_grad;  
extern PetscInt SuspendedParticles; //  0: dye or contaminant with no settling velocity; 1: suspended sediment and/or any other particles nd materials with settling velocity
extern PetscReal w_s,U_Bulk,d50,deltab,Background_Conc,Inlet_Conc;
extern PetscInt mobile_bed;

double sigma_phi = 1.0;

extern double kappa;

#define PI 3.14159265

void RHS_Conv_Diff(UserCtx *user, Vec ConvDiff_RHS)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	
	PetscReal	***aj;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***Conc, ***Conc_o, ***convdiff_rhs;
	Cmpnts	***ucont, ***ucat;
	Cmpnts	***csi, ***eta, ***zet;
        PetscReal	***nvert, ***lnu_t;
 
	Vec Fp1, Fp2, Fp3;
	Vec Visc1, Visc2, Visc3;
	PetscReal ***fp1, ***fp2, ***fp3;
	PetscReal ***visc1, ***visc2, ***visc3;
	
	Cmpnts	***icsi, ***ieta, ***izet;
	Cmpnts	***jcsi, ***jeta, ***jzet;
	Cmpnts	***kcsi, ***keta, ***kzet;
	PetscReal	***iaj, ***jaj, ***kaj, ***rho, ***mu;


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
	
        VecDuplicate(user->lConc, &Fp1);
	VecDuplicate(user->lConc, &Fp2);
	VecDuplicate(user->lConc, &Fp3);
	VecDuplicate(user->lConc, &Visc1);
	VecDuplicate(user->lConc, &Visc2);
	VecDuplicate(user->lConc, &Visc3);
	
	VecSet(Fp1,0);
	VecSet(Fp2,0);
	VecSet(Fp3,0);
	VecSet(Visc1,0);
	VecSet(Visc2,0);
	VecSet(Visc3,0);
	
	
	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
	}
	
        if(rans || les)	DAVecGetArray(da, user->lNu_t, &lnu_t);
	DAVecGetArray(fda, user->lUcont, &ucont);
	DAVecGetArray(fda, user->lUcat,  &ucat);
	
	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	
	DAVecGetArray(fda, user->lICsi, &icsi);
	DAVecGetArray(fda, user->lIEta, &ieta);
	DAVecGetArray(fda, user->lIZet, &izet);
	DAVecGetArray(fda, user->lJCsi, &jcsi);
	DAVecGetArray(fda, user->lJEta, &jeta);
	DAVecGetArray(fda, user->lJZet, &jzet);
	DAVecGetArray(fda, user->lKCsi, &kcsi);
	DAVecGetArray(fda, user->lKEta, &keta);
	DAVecGetArray(fda, user->lKZet, &kzet);
	
	DAVecGetArray(da, user->lNvert, &nvert);
	
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lIAj, &iaj);
	DAVecGetArray(da, user->lJAj, &jaj);
	DAVecGetArray(da, user->lKAj, &kaj);
	
	DAVecGetArray(da, user->lConc, &Conc);
	DAVecGetArray(da, user->lConc_o, &Conc_o);
	DAVecGetArray(da, ConvDiff_RHS, &convdiff_rhs);
		
	DAVecGetArray(da, Fp1, &fp1);
	DAVecGetArray(da, Fp2, &fp2);
	DAVecGetArray(da, Fp3, &fp3);
	DAVecGetArray(da, Visc1, &visc1);
	DAVecGetArray(da, Visc2, &visc2);
	DAVecGetArray(da, Visc3, &visc3);		
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs-1; i<lxe; i++) {
		double csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
		double eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
		double zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
		double ajc = iaj[k][j][i];
		
		double g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
		double g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
		double g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
		
		double dCdc, dCde, dCdz;
		double nu_t;
		
	        double diff_coef;
		
		dCdc = Conc[k][j][i+1] - Conc[k][j][i];

		if ((nvert[k][j+1][i])> 1.1 || (nvert[k][j+1][i+1])> 1.1) {
			dCde = (Conc[k][j  ][i+1] + Conc[k][j  ][i] - Conc[k][j-1][i+1] - Conc[k][j-1][i]) * 0.5;
		}
		else if  ((nvert[k][j-1][i])> 1.1 || (nvert[k][j-1][i+1])> 1.1) {
			dCde = (Conc[k][j+1][i+1] + Conc[k][j+1][i] - Conc[k][j  ][i+1] - Conc[k][j  ][i]) * 0.5;
		}
		else {
			dCde = (Conc[k][j+1][i+1] + Conc[k][j+1][i] - Conc[k][j-1][i+1] - Conc[k][j-1][i]) * 0.25;
		}	  

		if ((nvert[k+1][j][i])> 1.1 || (nvert[k+1][j][i+1])> 1.1) {
			dCdz = (Conc[k  ][j][i+1] + Conc[k  ][j][i] - Conc[k-1][j][i+1] - Conc[k-1][j][i]) * 0.5;
		}
		else if ((nvert[k-1][j][i])> 1.1 || (nvert[k-1][j][i+1])> 1.1) {
			dCdz = (Conc[k+1][j][i+1] + Conc[k+1][j][i] - Conc[k  ][j][i+1] - Conc[k  ][j][i]) * 0.5;
		}
		else {
			dCdz = (Conc[k+1][j][i+1] + Conc[k+1][j][i] - Conc[k-1][j][i+1] - Conc[k-1][j][i]) * 0.25;
		}
		
		if(rans || les){ nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j][i+1]);
		nu_t = PetscMax(nu_t, 0);}


		if( nvert[k][j][i]+nvert[k][j][i+1]>0.1 || i==0 || i==mx-2 || periodic ) {
			fp1[k][j][i] = -ucont[k][j][i].x * Upwind ( Conc[k][j][i], Conc[k][j][i+1], ucont[k][j][i].x);
		}
		else {
			fp1[k][j][i] = -ucont[k][j][i].x * weno3 ( Conc[k][j][i-1], Conc[k][j][i], Conc[k][j][i+1], Conc[k][j][i+2], ucont[k][j][i].x );
		}
		
	                //PetscPrintf(PETSC_COMM_WORLD, "Concentration %e\n", Conc[k][j][i] );
		double nu = 1./user->ren;
		if(levelset) {
			if(nvert[k][j][i]>0.1 || i==0) nu=mu[k][j][i+1];
			else if(nvert[k][j][i+1]>0.1 || i==mx-2) nu=mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k][j][i+1] );
		}
			diff_coef = nu + sigma_phi * nu_t;
                        if (density_current == 2) diff_coef = 0.00000001;
			//if(rans || les) {visc1[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * (nu + sigma_phi * nu_t);}
                        //else {visc1[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * (nu + sigma_phi * nu_t);}
			if(rans || les) {visc1[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * diff_coef;}
                        else {visc1[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * diff_coef;}
	}
	
	for (k=lzs; k<lze; k++)
	for (j=lys-1; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
		double eta0= jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
		double zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
		double ajc = jaj[k][j][i];
		
		double g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
		double g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
		double g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
		
		double dCdc, dCde, dCdz;
	 	double nu_t;		

	        double diff_coef;
	
	
		if ((nvert[k][j][i+1])> 1.1 || (nvert[k][j+1][i+1])> 1.1) {
			dCdc = (Conc[k][j+1][i  ]+ Conc[k][j][i  ] - Conc[k][j+1][i-1] - Conc[k][j][i-1]) * 0.5;
		}
		else if ((nvert[k][j][i-1])> 1.1 || (nvert[k][j+1][i-1])> 1.1) {
			dCdc = (Conc[k][j+1][i+1] + Conc[k][j][i+1] - Conc[k][j+1][i  ] - Conc[k][j][i  ]) * 0.5;
		}
		else {
			dCdc = (Conc[k][j+1][i+1] + Conc[k][j][i+1] - Conc[k][j+1][i-1] - Conc[k][j][i-1]) * 0.25;
		}

		dCde = Conc[k][j+1][i] - Conc[k][j][i];

		if ((nvert[k+1][j][i])> 1.1 || (nvert[k+1][j+1][i])> 1.1) {
			dCdz = (Conc[k  ][j+1][i] + Conc[k  ][j][i] - Conc[k-1][j+1][i] - Conc[k-1][j][i]) * 0.5;
		}
		else if ((nvert[k-1][j][i])> 1.1 || (nvert[k-1][j+1][i])> 1.1) {
			dCdz = (Conc[k+1][j+1][i] + Conc[k+1][j][i] - Conc[k  ][j+1][i] - Conc[k  ][j][i]) * 0.5;
		}
		else {
			dCdz = (Conc[k+1][j+1][i] + Conc[k+1][j][i] - Conc[k-1][j+1][i] - Conc[k-1][j][i]) * 0.25;
		}
		
                if(rans || les) {nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i]);
		nu_t = PetscMax(nu_t, 0);}
		
		if( nvert[k][j][i]+nvert[k][j+1][i]>0.1 || j==0 || j==my-2 || periodic ) {
			fp2[k][j][i] = -ucont[k][j][i].y * Upwind ( Conc[k][j][i], Conc[k][j+1][i], ucont[k][j][i].y);
		}
		else {
			fp2[k][j][i] = -ucont[k][j][i].y * weno3 ( Conc[k][j-1][i], Conc[k][j][i], Conc[k][j+1][i], Conc[k][j+2][i], ucont[k][j][i].y );
		}
		
		
		double nu = 1./user->ren;
		if(levelset) {
			if(nvert[k][j][i]>0.1 || j==0) nu=mu[k][j+1][i];
			else if(nvert[k][j+1][i]>0.1 || j==my-2) nu=mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k][j+1][i] );
		}
		
			diff_coef = nu + sigma_phi * nu_t;
                        if (density_current == 2) diff_coef = 0.00000001;
			//if(rans || les) {visc2[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * (nu + sigma_phi * nu_t);}
			//else {visc2[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * (nu + sigma_phi * nu_t);}
			if(rans || les) {visc2[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * diff_coef;}
			else {visc2[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * diff_coef;}
	}
	
	for (k=lzs-1; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
		double eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
		double zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
		double ajc = kaj[k][j][i];
		
		double g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
		double g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
		double g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
		
		double dCdc, dCde, dCdz;
		double nu_t;
		
	        double diff_coef;
		
                if ((nvert[k][j][i+1])> 1.1 || (nvert[k+1][j][i+1])> 1.1) {
			dCdc = (Conc[k+1][j][i  ] + Conc[k][j][i  ] - Conc[k+1][j][i-1] - Conc[k][j][i-1]) * 0.5;
		}
		else if ((nvert[k][j][i-1])> 1.1 || (nvert[k+1][j][i-1])> 1.1) {
			dCdc = (Conc[k+1][j][i+1] + Conc[k][j][i+1] - Conc[k+1][j][i  ] - Conc[k][j][i  ]) * 0.5;
		}
		else {
			dCdc = (Conc[k+1][j][i+1] + Conc[k][j][i+1] - Conc[k+1][j][i-1] - Conc[k][j][i-1]) * 0.25;
		}

		if ((nvert[k][j+1][i])> 1.1 || (nvert[k+1][j+1][i])> 1.1) {
			dCde = (Conc[k+1][j  ][i] + Conc[k][j  ][i] - Conc[k+1][j-1][i] - Conc[k][j-1][i]) * 0.5;
		}
		else if ((nvert[k][j-1][i])> 1.1 || (nvert[k+1][j-1][i])> 1.1) {
			dCde = (Conc[k+1][j+1][i] + Conc[k][j+1][i] - Conc[k+1][j  ][i] - Conc[k][j  ][i]) * 0.5;
		}
		else {
			dCde = (Conc[k+1][j+1][i] + Conc[k][j+1][i] - Conc[k+1][j-1][i] - Conc[k][j-1][i]) * 0.25;
		}

		dCdz = Conc[k+1][j][i] - Conc[k][j][i];
		
                if(rans || les) {nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k+1][j][i]);
		nu_t = PetscMax(nu_t, 0);}
		
		if( nvert[k][j][i]+nvert[k+1][j][i]>0.1 || k==0 || k==mz-2 || periodic ) {
			fp3[k][j][i] = -ucont[k][j][i].z * Upwind ( Conc[k][j][i], Conc[k+1][j][i], ucont[k][j][i].z);
		}
		else {
			fp3[k][j][i] = -ucont[k][j][i].z * weno3 ( Conc[k-1][j][i], Conc[k][j][i], Conc[k+1][j][i], Conc[k+2][j][i], ucont[k][j][i].z );
		}
		
		double nu = 1./user->ren;
		if(levelset) {
			if(nvert[k][j][i]>0.1 || k==0) nu=mu[k+1][j][i];
			else if(nvert[k+1][j][i]>0.1 || k==mz-2) nu=mu[k][j][i];
			else nu = 0.5 * ( mu[k][j][i] + mu[k+1][j][i] );
		}
		
			diff_coef = nu + sigma_phi * nu_t;
                        if (density_current ==2) diff_coef = 0.00000001;
                        //if(rans || les) {visc3[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * (nu + sigma_phi * nu_t);}
			//else {visc3[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * (nu + sigma_phi * nu_t);}
                        if(rans || les) {visc3[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * diff_coef;}
			else {visc3[k][j][i] = (g11 * dCdc + g21 * dCde + g31 * dCdz) * ajc * diff_coef;}
	}
	
	DAVecRestoreArray(da, Fp1, &fp1);
	DAVecRestoreArray(da, Fp2, &fp2);
	DAVecRestoreArray(da, Fp3, &fp3);
	DAVecRestoreArray(da, Visc1, &visc1);
	DAVecRestoreArray(da, Visc2, &visc2);
	DAVecRestoreArray(da, Visc3, &visc3);

	DALocalToLocalBegin(da, Fp1, INSERT_VALUES, Fp1);
	DALocalToLocalEnd(da, Fp1, INSERT_VALUES, Fp1);
	
	DALocalToLocalBegin(da, Fp2, INSERT_VALUES, Fp2);
	DALocalToLocalEnd(da, Fp2, INSERT_VALUES, Fp2);
	
	DALocalToLocalBegin(da, Fp3, INSERT_VALUES, Fp3);
	DALocalToLocalEnd(da, Fp3, INSERT_VALUES, Fp3);
	
	DALocalToLocalBegin(da, Visc1, INSERT_VALUES, Visc1);
	DALocalToLocalEnd(da, Visc1, INSERT_VALUES, Visc1);
	
	DALocalToLocalBegin(da, Visc2, INSERT_VALUES, Visc2);
	DALocalToLocalEnd(da, Visc2, INSERT_VALUES, Visc2);
	
	DALocalToLocalBegin(da, Visc3, INSERT_VALUES, Visc3);
	DALocalToLocalEnd(da, Visc3, INSERT_VALUES, Visc3);
	
	DAVecGetArray(da, Fp1, &fp1);
	DAVecGetArray(da, Fp2, &fp2);
	DAVecGetArray(da, Fp3, &fp3);
	DAVecGetArray(da, Visc1, &visc1);
	DAVecGetArray(da, Visc2, &visc2);
	DAVecGetArray(da, Visc3, &visc3);

	if(periodic)
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int a=i, b=j, c=k;

		int flag=0;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
				
		if(flag) {
			fp1[k][j][i] = fp1[c][b][a];
			fp2[k][j][i] = fp2[c][b][a];
			fp3[k][j][i] = fp3[c][b][a];
			visc1[k][j][i] = visc1[c][b][a];
			visc2[k][j][i] = visc2[c][b][a];
			visc3[k][j][i] = visc3[c][b][a];
		}
	}
	


	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if ( i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1 ) { // it was 0.1 for k-omega  not 1.1
			convdiff_rhs[k][j][i] = 0;
			//komega_rhs[k][j][i].x = komega_rhs[k][j][i].y = 0;
			continue;
		}
		
                double ajc = aj[k][j][i];
                 
		if ( nvert[k][j][i] < 0.1 ) {
			double r = 1.;
			
			if(levelset) r = rho[k][j][i];

                  convdiff_rhs[k][j][i] = ( fp1[k][j][i] - fp1[k][j][i-1] + fp2[k][j][i] - fp2[k][j-1][i] + fp3[k][j][i] - fp3[k-1][j][i] ) * ajc;// advection
		  convdiff_rhs[k][j][i] += ( visc1[k][j][i] - visc1[k][j][i-1] + visc2[k][j][i] - visc2[k][j-1][i] + visc3[k][j][i] - visc3[k-1][j][i] ) * ajc / r;// diffusion
		}
	}
	
	if(rans || les) DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
        DAVecRestoreArray(da, Fp1, &fp1);
	DAVecRestoreArray(da, Fp2, &fp2);
	DAVecRestoreArray(da, Fp3, &fp3);
	DAVecRestoreArray(da, Visc1, &visc1);
	DAVecRestoreArray(da, Visc2, &visc2);
	DAVecRestoreArray(da, Visc3, &visc3);
	
	DAVecRestoreArray(fda, user->lUcont, &ucont);
	DAVecRestoreArray(fda, user->lUcat,  &ucat);
	
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	
	DAVecRestoreArray(fda, user->lICsi, &icsi);
	DAVecRestoreArray(fda, user->lIEta, &ieta);
	DAVecRestoreArray(fda, user->lIZet, &izet);
	DAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DAVecRestoreArray(fda, user->lJEta, &jeta);
	DAVecRestoreArray(fda, user->lJZet, &jzet);
	DAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DAVecRestoreArray(fda, user->lKEta, &keta);
	DAVecRestoreArray(fda, user->lKZet, &kzet);
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lIAj, &iaj);
	DAVecRestoreArray(da, user->lJAj, &jaj);
	DAVecRestoreArray(da, user->lKAj, &kaj);
	
	DAVecRestoreArray(da, user->lConc, &Conc);
	DAVecRestoreArray(da, user->lConc_o, &Conc_o);
	DAVecRestoreArray(da, ConvDiff_RHS, &convdiff_rhs);
	
        if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
	}
	
	VecDestroy(Fp1);
	VecDestroy(Fp2);
	VecDestroy(Fp3);
	VecDestroy(Visc1);
	VecDestroy(Visc2);
	VecDestroy(Visc3);
};


void Force_Current(UserCtx *user)
{
	DA		da = user->da;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	double Ri = user->Ri;

	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	PetscReal 	***conc;
	PetscReal	***nvert;

	Cmpnts		***fcurrent, ***csi, ***eta, ***zet;

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

	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(user->da, user->lConc, &conc);
	DAVecGetArray(user->fda, user->lFCurrent, &fcurrent);
  	DAVecGetArray(user->fda, user->lCsi,  &csi);
  	DAVecGetArray(user->fda, user->lEta,  &eta);
  	DAVecGetArray(user->fda, user->lZet,  &zet);


	double force_x, force_y, force_z;
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	// pressure node; cell center

		force_x = 0.0;
		force_y = 0.0;
		force_z = 0.0;

//		if (user->bctype_tmprt[0] || user->bctype_tmprt[1]) force_x = Ri * (tmprt[k][j][i]-tmprt_xzAvg[j]); 
//		if (user->bctype_tmprt[2] || user->bctype_tmprt[3]) force_y = Ri * tmprt[k][j][i]; 
//		if (user->bctype_tmprt[4] || user->bctype_tmprt[5]) force_z = Ri * (tmprt[k][j][i]-tmprt_xzAvg[j]); 
		force_z = Ri * conc[k][j][i]; 

	      	fcurrent[k][j][i].x = force_x *  csi[k][j][i].x + force_y *  csi[k][j][i].y + force_z *  csi[k][j][i].z;
	      	fcurrent[k][j][i].y = force_x *  eta[k][j][i].x + force_y *  eta[k][j][i].y + force_z *  eta[k][j][i].z;
	      	fcurrent[k][j][i].z = force_x *  zet[k][j][i].x + force_y *  zet[k][j][i].y + force_z *  zet[k][j][i].z;

//		if (k==3 && i==3) printf("the force  %d %d %le %le %le %le\n", my-2, j, ftmprt[k][j][i].y, force_x, force_y, force_z);



                if( nvert[k][j][i]>0.1) {
                  	fcurrent[k][j][i].x = 0.0;
                	fcurrent[k][j][i].y = 0.0;
                	fcurrent[k][j][i].z = 0.0;

                }
                else if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || j==my-2) {
                        fcurrent[k][j][i].x = 0.0;
                        fcurrent[k][j][i].y = 0.0;
                        fcurrent[k][j][i].z = 0.0;

                }

//		if (user->bctype_tmprt[1] && i==mx-2) ftmprt[k][j][i].x = 0.0;
//		if (user->bctype_tmprt[3] && j==my-2) ftmprt[k][j][i].y = 0.0;
//		if (user->bctype_tmprt[5] && k==mz-2) ftmprt[k][j][i].z = 0.0;


	}
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(user->da, user->lConc, &conc);
	DAVecRestoreArray(user->fda, user->lFCurrent, &fcurrent);

	DALocalToGlobal(user->fda, user->lFCurrent, INSERT_VALUES, user->FCurrent);

//	TECIOOut_rhs(user,  user->FTmprt);

        DAVecRestoreArray(user->fda, user->lCsi,  &csi);
        DAVecRestoreArray(user->fda, user->lEta,  &eta);
        DAVecRestoreArray(user->fda, user->lZet,  &zet);
	
};

void Conv_Diff_IC(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
        Vec Coor;
	Cmpnts	***coor;

	double nu = 1./user->ren;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
	PetscReal	***nvert;
	PetscReal	***Conc;

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

	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lConc, &Conc);
	
        DAGetGhostedCoordinates(da, &Coor);
        DAVecGetArray(fda, Coor, &coor);
	// I.C. for Conc
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {	// concentration node; cell center
		
		//Conc[k][j][i] = 0.55;
		if(density_current == 1){
                 double xc  = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25 ;
                    if (xc < 1. && nvert[k][j][i]<1.1) {
                         Conc[k][j][i]=1.0;
                         if(j==my-2) Conc[k][j+1][i]=1.0;
                         if(j==1) Conc[k][j-1][i]=1.0;       
                         if(i==mx-2) Conc[k][j][i+1]=1.0;
                         if(i==1) Conc[k][j][i-1]=1.0;       
                         if(k==1) Conc[k-1][j][i]=1.0;       
                              }
                    else   {
                         Conc[k][j][i]=0.0;
		         if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) Conc[k][j][i] = 0.;
                              }
		}
		else if(density_current == 2) {
			if( nvert[k][j][i]<5.1) Conc[k][j][i] = 0.0;
		}
		else if(density_current == 3) {
			if( nvert[k][j][i]<5.1) Conc[k][j][i] = 0.0;
		}
                else {
                      if( nvert[k][j][i]<5.1) Conc[k][j][i] = 0.0;
                     }
	
        	if( nvert[k][j][i]>1.6) {   //Solid nodes common command for all
			Conc[k][j][i] = 0.0;
		}
		else if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) {   //commom command for all
			Conc[k][j][i] = 0.;
		}

        	if( nvert[k][j][i]<= 1.1) {   //fluid node and ib
			Conc[k][j][i] = Background_Conc;}
	}
	
	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lConc, &Conc);
	DAVecRestoreArray(fda, Coor, &coor);
	
	DALocalToGlobal(da, user->lConc, INSERT_VALUES, user->Conc);
	VecCopy(user->Conc, user->Conc_o);
	
};

void SettlingVelocityCorrection (UserCtx *user)
{
	PetscInt      i, j, k;
	DALocalInfo	info ;
	PetscInt	xs, xe, ys, ye, zs, ze;
	PetscInt  	mx,my,mz;	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts        ***ucont;
	PetscReal     ***nvert;

	DA            da = user->da,fda = user->fda;
	info = user->info;
  
	xs = info.xs; xe = info.xs + info.xm;
	ys = info.ys; ye = info.ys + info.ym;
	zs = info.zs; ze = info.zs + info.zm;
	mx = info.mx; my = info.my; mz = info.mz;
  
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	Cmpnts ***ucat, ***icsi, ***jeta, ***kzet;
	
	/*if(!immersed && ti==tistart) {
          DAVecGetArray(da, user->lNvert, &nvert);
          for (k=lzs; k<lze; k++) {
            for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
                if(user->bctype[0]==-1 && i==1) nvert[k][j][i]=1;
                if(user->bctype[1]==-1 && i==mx-2) nvert[k][j][i]=1;
                if(user->bctype[2]==-1 && j==1) nvert[k][j][i]=1;
                if(user->bctype[3]==-1 && j==my-2) nvert[k][j][i]=1;
                if(user->bctype[4]==-1 && k==1) nvert[k][j][i]=1;
                if(user->bctype[5]==-1 && k==mz-2) nvert[k][j][i]=1;

                if(user->bctype[0]==-2 && i==1) nvert[k][j][i]=1;
                if(user->bctype[1]==-2 && i==mx-2) nvert[k][j][i]=1;
                if(user->bctype[2]==-2 && j==1) nvert[k][j][i]=1;
		if(user->bctype[3]==-2 && j==my-2) nvert[k][j][i]=1;
                if(user->bctype[4]==-2 && k==1) nvert[k][j][i]=1;
                if(user->bctype[5]==-2 && k==mz-2) nvert[k][j][i]=1;
              }
	    }
	  }
          DAVecRestoreArray(da, user->lNvert, &nvert);
	  DALocalToLocalBegin(da, user->lNvert, INSERT_VALUES, user->lNvert);
	  DALocalToLocalEnd(da, user->lNvert, INSERT_VALUES, user->lNvert);
	  DALocalToGlobal(da, user->lNvert, INSERT_VALUES, user->Nvert);
	}*/
	
	DAVecGetArray(fda, user->lUcat, &ucat);
	DAVecGetArray(fda, user->lUcont, &ucont);//
	DAVecGetArray(fda, user->lICsi, &icsi);//
	DAVecGetArray(fda, user->lJEta, &jeta);//
	DAVecGetArray(fda, user->lKZet, &kzet);//
	DAVecGetArray(da, user->lNvert, &nvert);//
	
/*	if(periodic)
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	
		int flag=0, a=i, b=j, c=k;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		
		if(flag) {
			ucat[k][j][i] = ucat[c][b][a];
		}
	}*/

	double  ucx, ucy, ucz;
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	  
		//if(immersed) {
			//double f = 1.0;

			//if(immersed==3) f = 0;
			
	//		if ((int)nvert[k][j][i]==1) {
				ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
				if(!y_direction)
				{
					ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
					ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5 - w_s;
				}
				else {
					ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5 - w_s;
					ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
				}
				ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z);
					
				ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
				if(!y_direction)
				{
					ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
					ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5 - w_s;
				}
				else {
					ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5 - w_s;
				    ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
                }

				ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z);
					
				ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
				if (!y_direction)
				{
					ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
					ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5 - w_s;
				}
				else {
					ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5 - w_s;
				    ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
                     }
				ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z);
		//	}
		/*		
			if ((int)(nvert[k][j][i+1])==1) {
				ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
				ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
				ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
				ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z);
			}
						
			if ((int)(nvert[k][j+1][i])==1) {
				ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z);
			}
						
			if ((int)(nvert[k+1][j][i])==1) {
				ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z);
			}
			
			if(!movefsi && !rotatefsi && immersed==3) {
				if ((nvert[k][j][i+1]+nvert[k][j][i])>1.1) ucont[k][j][i].x = 0;
				if ((nvert[k][j+1][i]+nvert[k][j][i])>1.1) ucont[k][j][i].y = 0;
				if ((nvert[k+1][j][i]+nvert[k][j][i])>1.1) ucont[k][j][i].z = 0;
			}
		*/
	  
		// wall func
		
		if ( (user->bctype[0]==-1 && i==1) ||( user->bctype[1]==-1 && i==mx-2) ) {
			if(i==1) ucont[k][j][i-1].x=0;
			else ucont[k][j][i].x=0;
		}
		
		if ( (user->bctype[0]==-2 && i==1) ||( user->bctype[1]==-2 && i==mx-2) ) {
			if(i==1) ucont[k][j][i-1].x=0;
			else ucont[k][j][i].x=0;
		}
		
		if ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 && j==my-2) ) {
			if(j==1) ucont[k][j-1][i].y=0;
                        else ucont[k][j][i].y=0;
		}
		
		if ( (user->bctype[2]==-2 && j==1) || (user->bctype[3]==-2 && j==my-2) ) {
			if(j==1) ucont[k][j-1][i].y=0;
                        else ucont[k][j][i].y=0;
		}
		
	}
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if ( (user->bctype[0]==10) && i==1) ucont[k][j][0].x = 0;
		if ( (user->bctype[1]==10) && i==mx-2) ucont[k][j][mx-2].x = 0;
		if ( (user->bctype[2]==10) && j==1) ucont[k][0][i].y = 0;
		if ( (user->bctype[3]==10) && j==my-2) ucont[k][my-2][i].y = 0;
		if ( (user->bctype[4]==10) && k==1) ucont[0][j][i].z = 0;
		if ( (user->bctype[5]==10) && k==mz-2) ucont[mz-2][j][i].z = 0;
		
		if ( i_periodic && i==0 ) ucont[k][j][0].x = ucont[k][j][mx-2].x;
		if ( i_periodic && i==mx-1 ) ucont[k][j][mx-1].x = ucont[k][j][1].x;
		if ( j_periodic && j==0 ) ucont[k][0][i].y = ucont[k][my-2][i].y;
		if ( j_periodic && j==my-1 ) ucont[k][my-1][i].y = ucont[k][1][i].y;
		if ( k_periodic && k==0 ) ucont[0][j][i].z = ucont[mz-2][j][i].z;
		if ( k_periodic && k==mz-1 ) ucont[mz-1][j][i].z = ucont[1][j][i].z;
		
		if ( ii_periodic && i==0 ) ucont[k][j][0].x = ucont[k][j][-2].x;
		if ( ii_periodic && i==mx-1 ) ucont[k][j][mx-1].x = ucont[k][j][mx+1].x;
		
		if ( jj_periodic && j==0 ) ucont[k][0][i].y = ucont[k][-2][i].y;
		if ( jj_periodic && j==my-1 ) ucont[k][my-1][i].y = ucont[k][my+1][i].y;
		
		if ( kk_periodic && k==0 ) ucont[0][j][i].z = ucont[-2][j][i].z;
		if ( kk_periodic && k==mz-1 ) ucont[mz-1][j][i].z = ucont[mz+1][j][i].z;
	}
	DAVecRestoreArray(fda, user->lUcat, &ucat);//
	DAVecRestoreArray(fda, user->lICsi, &icsi);//
	DAVecRestoreArray(fda, user->lJEta, &jeta);//
	DAVecRestoreArray(fda, user->lKZet, &kzet);//
	DAVecRestoreArray(fda, user->lUcont, &ucont);//
	DAVecRestoreArray(da, user->lNvert, &nvert);//
	
	
	DALocalToLocalBegin(fda, user->lUcont, INSERT_VALUES, user->lUcont);
	DALocalToLocalEnd(fda, user->lUcont, INSERT_VALUES, user->lUcont);
	
	return;
}

void SettlingVelocityCorrectionBack (UserCtx *user)
{
	PetscInt      i, j, k;
	DALocalInfo	info ;
	PetscInt	xs, xe, ys, ye, zs, ze;
	PetscInt  	mx,my,mz;	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	Cmpnts        ***ucont;
	PetscReal     ***nvert;

	DA            da = user->da,fda = user->fda;
	info = user->info;
  
	xs = info.xs; xe = info.xs + info.xm;
	ys = info.ys; ye = info.ys + info.ym;
	zs = info.zs; ze = info.zs + info.zm;
	mx = info.mx; my = info.my; mz = info.mz;
  
	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;
  
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
  
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	Cmpnts ***ucat, ***icsi, ***jeta, ***kzet;
	
	
	DAVecGetArray(fda, user->lUcat, &ucat);
	DAVecGetArray(fda, user->lUcont, &ucont);//
	DAVecGetArray(fda, user->lICsi, &icsi);//
	DAVecGetArray(fda, user->lJEta, &jeta);//
	DAVecGetArray(fda, user->lKZet, &kzet);//
	DAVecGetArray(da, user->lNvert, &nvert);//
	
	double  ucx, ucy, ucz;
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	  
		//if(immersed) {
			//double f = 1.0;

			//if(immersed==3) f = 0;
			
		//	if ((int)nvert[k][j][i]==1) {
				ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
				ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
				ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
				ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z);
					
				ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z);
					
				ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z);
		//	}
				
		/*	if ((int)(nvert[k][j][i+1])==1) {
				ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
				ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
				ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
				ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z);
			}
						
			if ((int)(nvert[k][j+1][i])==1) {
				ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z);
			}
						
			if ((int)(nvert[k+1][j][i])==1) {
				ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
				ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
				ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
				ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z);
			}
			
			if(!movefsi && !rotatefsi && immersed==3) {
				if ((nvert[k][j][i+1]+nvert[k][j][i])>1.1) ucont[k][j][i].x = 0;
				if ((nvert[k][j+1][i]+nvert[k][j][i])>1.1) ucont[k][j][i].y = 0;
				if ((nvert[k+1][j][i]+nvert[k][j][i])>1.1) ucont[k][j][i].z = 0;
			}
		*/
	  
		// wall func
		
		if ( (user->bctype[0]==-1 && i==1) ||( user->bctype[1]==-1 && i==mx-2) ) {
			if(i==1) ucont[k][j][i-1].x=0;
			else ucont[k][j][i].x=0;
		}
		
		if ( (user->bctype[0]==-2 && i==1) ||( user->bctype[1]==-2 && i==mx-2) ) {
			if(i==1) ucont[k][j][i-1].x=0;
			else ucont[k][j][i].x=0;
		}
		
		if ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 && j==my-2) ) {
			if(j==1) ucont[k][j-1][i].y=0;
                        else ucont[k][j][i].y=0;
		}
		
		if ( (user->bctype[2]==-2 && j==1) || (user->bctype[3]==-2 && j==my-2) ) {
			if(j==1) ucont[k][j-1][i].y=0;
                        else ucont[k][j][i].y=0;
		}
		
	}
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if ( (user->bctype[0]==10) && i==1) ucont[k][j][0].x = 0;
		if ( (user->bctype[1]==10) && i==mx-2) ucont[k][j][mx-2].x = 0;
		if ( (user->bctype[2]==10) && j==1) ucont[k][0][i].y = 0;
		if ( (user->bctype[3]==10) && j==my-2) ucont[k][my-2][i].y = 0;
	if ( (user->bctype[4]==10) && k==1) ucont[0][j][i].z = 0;
		if ( (user->bctype[5]==10) && k==mz-2) ucont[mz-2][j][i].z = 0;
		
		if ( i_periodic && i==0 ) ucont[k][j][0].x = ucont[k][j][mx-2].x;
		if ( i_periodic && i==mx-1 ) ucont[k][j][mx-1].x = ucont[k][j][1].x;
		if ( j_periodic && j==0 ) ucont[k][0][i].y = ucont[k][my-2][i].y;
		if ( j_periodic && j==my-1 ) ucont[k][my-1][i].y = ucont[k][1][i].y;
		if ( k_periodic && k==0 ) ucont[0][j][i].z = ucont[mz-2][j][i].z;
		if ( k_periodic && k==mz-1 ) ucont[mz-1][j][i].z = ucont[1][j][i].z;
		
		if ( ii_periodic && i==0 ) ucont[k][j][0].x = ucont[k][j][-2].x;
		if ( ii_periodic && i==mx-1 ) ucont[k][j][mx-1].x = ucont[k][j][mx+1].x;
		
		if ( jj_periodic && j==0 ) ucont[k][0][i].y = ucont[k][-2][i].y;
		if ( jj_periodic && j==my-1 ) ucont[k][my-1][i].y = ucont[k][my+1][i].y;
		
		if ( kk_periodic && k==0 ) ucont[0][j][i].z = ucont[-2][j][i].z;
		if ( kk_periodic && k==mz-1 ) ucont[mz-1][j][i].z = ucont[mz+1][j][i].z;
	}
	DAVecRestoreArray(fda, user->lUcat, &ucat);//
	DAVecRestoreArray(fda, user->lICsi, &icsi);//
	DAVecRestoreArray(fda, user->lJEta, &jeta);//
	DAVecRestoreArray(fda, user->lKZet, &kzet);//
	DAVecRestoreArray(fda, user->lUcont, &ucont);//
	DAVecRestoreArray(da, user->lNvert, &nvert);//
	
	
	DALocalToLocalBegin(fda, user->lUcont, INSERT_VALUES, user->lUcont);
	DALocalToLocalEnd(fda, user->lUcont, INSERT_VALUES, user->lUcont);
	
	return;
}


double Depth_Averaged_Conc( PetscReal ***Conc, PetscReal ***nvert, int i, int j, int k, int lye)
{
/*
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions 

	PetscReal	***Conc;
	PetscReal	***nvert;

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

	DAVecGetArray(da, user->lConc, &Conc);
	DAVecGetArray(da, user->lNvert, &nvert);
  */      
	PetscInt	jj; 
        double c_v;
        double sum = 0.;
               c_v = 0.;
	for (jj=j-1; jj<lye; jj++)
             {
		if ( nvert[k][jj][i]<1.1) {
                         sum += 1.;
			 c_v += Conc[k][jj][i]; 
                         }
             }
         if(sum>0.) c_v = c_v/sum;                    

/*
	DAVecRestoreArray(da, user->lConc, &Conc);
        DALocalToLocalBegin(da, user->lConc, INSERT_VALUES, user->lConc);
        DALocalToLocalEnd(da, user->lConc, INSERT_VALUES, user->lConc);
*/
return c_v;
}


void Conv_Diff_BC(UserCtx *user)
{
	DA		da = user->da, fda = user->fda;
	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k, ibi;
	
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

        //double nu = 1./user->ren;
	PetscReal	***aj, ***rho, ***mu;
	PetscReal	***lnu_t;
        extern PetscInt ti;	
	//Cmpnts2 ***K_Omega;
        Vec Coor;
	Cmpnts	***csi, ***eta, ***zet, ***ucat, ***coor;
	PetscReal	***nvert, ***ustar;
	PetscReal	***Conc;
        
        double nu_t, nu_tt;
        
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

	DAVecGetArray(fda, user->lCsi, &csi);
	DAVecGetArray(fda, user->lEta, &eta);
	DAVecGetArray(fda, user->lZet, &zet);
	DAVecGetArray(fda, user->lUcat, &ucat);

	if(rans || les) DAVecGetArray(da, user->lNu_t, &lnu_t);

        DAGetGhostedCoordinates(da, &Coor);
        DAVecGetArray(fda, Coor, &coor);
	
	DAVecGetArray(da, user->lNvert, &nvert);
	DAVecGetArray(da, user->lAj, &aj);
	DAVecGetArray(da, user->lUstar, &ustar);

	DAVecGetArray(da, user->lConc, &Conc);
	
	if(levelset) {
		DAVecGetArray(da, user->lDensity, &rho);
		DAVecGetArray(da, user->lMu, &mu);
	}
	
	// BC for Concentration_pollution_sespended sediment
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {	// pressure node; cell center

		double ren = user->ren;
  
		// from saved inflow file
		if(inletprofile==1000 && k==1) {
		//	Conc[k-1][j][i] = user->concentration_plane[j][i];
		//	Conc[k-1][j][i] = std::max ( Conc[k-1][j][i], 1. );
			Conc[k-1][j][i] = Inlet_Conc; 
		}
		else if ( user->bctype[4]==5 && k==1 && nvert[k][j][i]<1.3) {	// inflow
			if(!density_current)Conc[k-1][j][i] =Inlet_Conc;  
               // Jobson's dye injection materials input @ free surface
               // Nitrogen injection input @ free surface of OSL

                // double xc  = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25 ;
                 //if(ti==0)
                 //{
                   // if (xc < 2.) {
                     //    if(j==my-2) Conc[k][j+1][i]=1.0;
                               if (density_current == 1) Conc[k][j][i]= Conc[k+1][j][i];
                               else if (density_current == 2){
                                                             if(coor[k][j][i].z < 6.6666666667){Conc[k][j][i]= 1.0;}
                                                             else {Conc[k][j][i]= 0.0;}
                                                             }
                               else if (density_current == 3) Conc[k-1][j][i]= 1.0;
                               
                              //}  else   {Conc[k][j][i]=0.0;}}
                  //else { Conc[k][j][i]=Conc[k+1][j][i];}

                /* Used for the Dye Material Validation inlet @ x/y=3 prifile
                 double zc  = (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25 ;
                 if ( zc < 0.3 ){Conc[k][j][i] = 0.0;}
                 else {Conc[k][j][i] = exp ( (zc - 0.7643) / 0.1194 );}
                */                             
                 
                /* The below used for net deposition with sed. particles injected from free-surfacec at the inlet */              
                 /*if (j==my-2) {
                               Conc[k][j+1][i]=10.;
                               Conc[k][j][i]=10.0;
                              }
                       else   {
                               Conc[k][j][i]=0.0;
                              }*/
                // OSL sediment injection from whole the inlet area 
                      /* if (j==my-2) {
                               Conc[k][j+1][i]=30.*(3./0.01)*1.e-4;
                               Conc[k][j][i]=30.*(3./0.01)*1.e-4;
                              }
                       else   {
                               Conc[k][j][i]=0.0;
                              }*/
                 // Conc[k][j][i] = 1.e-4; inlet sediment concentration into the OSL
                //    Conc[k][j][i] = 1.0; // Nitrogen_inflow though the inlet plane
         	}	
		
		// "slip"  same as "Symetrical B.C."
		if ( user->bctype[0] == 10 && i==1 ) Conc[k][j][i-1] = Conc[k][j][i];
		if ( user->bctype[1] == 10 && i==mx-2 ) Conc[k][j][i+1] = Conc[k][j][i];
		if ( user->bctype[2] == 10 && j==1 ) Conc[k][j-1][i] = Conc[k][j][i];
		if ( SuspendedParticles && user->bctype[3] == 10 && j==my-2 ){
                       double delz = coor[k][j][i].z - coor[k][j-1][i].z;
					   double dely = coor[k][j][i].y - coor[k][j-1][i].y;
					   if (y_direction) delz = dely;
		       double nu  = 1./ren;
                       if (rans || les) {nu_t = PetscMax(0.0 , 0.5*(lnu_t[k][j][i]+lnu_t[k][j-1][i])); nu_tt = nu + sigma_phi * nu_t;}
                       else {nu_t = nu; nu_tt = nu;}
                      
                       Conc[k][j+1][i] = Conc[k][j][i]/(1. + delz * w_s / nu_tt);}
		if ( !SuspendedParticles && user->bctype[3] == 10 && j==my-2 ) Conc[k][j+1][i] = Conc[k][j][i];
		if ( user->bctype[4] == 10 && k==1 ) Conc[k-1][j][i] = Conc[k][j][i];
		if ( user->bctype[5] == 10 && k==mz-2 ) Conc[k+1][j][i] = Conc[k][j][i];
	
                // for density current in a closed reservoir	
		if ( user->bctype[4] == 1 && k==1 ) Conc[k][j][i] = Conc[k+1][j][i];
		//if ( user->bctype[4] == 1 && k==1 ) Conc[k-1][j][i] = Conc[k][j][i];
		if ( user->bctype[5] == 1 && k==mz-2 ) Conc[k+1][j][i] = Conc[k][j][i];



	
		// outflow
		if ( user->bctype[5] == 4 && k==mz-2 ) {
			Conc[k+1][j][i] = Conc[k][j][i];
		}
		
		// couette
		if ( user->bctype[3] == 12 && j==my-2 ) {
			//double dist = distance[k][j][i];
			//K_Omega[k][j][i].y = wall_omega(ren, dist);
			//if(j==my-2) K_Omega[k][j+1][i].x = - K_Omega[k][j][i].x;
		}
		
		// wall
		if ( user->bctype[0] == 1 && i<=1 ) {
			//double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			//double dist = distance[k][j][i];
			//dist = 0.5/aj[k][j][i]/area;
			//K_Omega[k][j][i].y = wall_omega(ren, dist);
			//if(i==1) K_Omega[k][j][i-1].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j][i-1].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
			if(i==1) Conc[k][j][i] = Conc[k][j][i+1];
			if(i==1) Conc[k][j][i-1] = Conc[k][j][i];
		}
		if ( user->bctype[1] == 1 && i>=mx-2 ) {
			//double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			//double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			//dist = 0.5/aj[k][j][i]/area;
			//K_Omega[k][j][i].y = wall_omega(ren, dist);
			//if(i==mx-2) K_Omega[k][j][i+1].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j][i+1].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
			if(i==mx-2) Conc[k][j][i] =  Conc[k][j][i-1];
			if(i==mx-2) Conc[k][j][i+1] =  Conc[k][j][i];
		}
		if ( user->bctype[2] == 1 && j<=1 ) {
			//double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			//double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			//dist = 0.5/aj[k][j][i]/area;
			//K_Omega[k][j][i].y = wall_omega(ren, dist);
			//if(j==1) K_Omega[k][j-1][i].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j-1][i].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
			if(j==1) Conc[k][j][i] =  Conc[k][j+1][i];
			if(j==1) Conc[k][j-1][i] =  Conc[k][j][i];
		}
		if ( user->bctype[3] == 1 && j>=my-2 ) {
			//double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			//double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			//dist = 0.5/aj[k][j][i]/area;
			//K_Omega[k][j][i].y = wall_omega(ren, dist);
			//if(j==my-2) K_Omega[k][j+1][i].x = - K_Omega[k][j][i].x + 1.e-5;
			//K_Omega[k][j+1][i].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
			if(j==my-2) Conc[k][j][i] =  Conc[k][j-1][i];
			if(j==my-2) Conc[k][j+1][i] =  Conc[k][j][i];
		}
				
		// wall function
		if( nvert[k][j][i]<1.1 && ( ( (user->bctype[0]==-1 || user->bctype[0]==-2) && i==1) || ( (user->bctype[1]==-1 || user->bctype[1]==-2) &&  i==mx-2) ) && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) {
		        /*
			double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
			double sb, sc; 
			Cmpnts Ua, Uc;
			Ua.x = Ua.y = Ua.z = 0;

			sb = 0.5/aj[k][j][i]/area;
			
			if(i==1) {
				sc = 2*sb + 0.5/aj[k][j][i+1]/area;
				Uc = ucat[k][j][i+1];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j][i-1]/area;
				Uc = ucat[k][j][i-1];
			}*/
			/*
			double ni[3], nj[3], nk[3];
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
			if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;

			wall_function (1./ren, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
			*/
			//double utau = ustar[k][j][i];
			
			//double Kc = utau*utau/sqrt(0.09);
			//double Ob = utau/sqrt(0.09)/(kappa*sb);
			
			//K_Omega[k][j][i].x = Kc;
			//K_Omega[k][j][i].y = Ob;
                        if(i==1)Conc[k][j][i] = Conc[k][j][i+1];
                        //if(i==1)Conc[k][j][i-1] = Conc[k][j][i];
                        if(i==mx-2)Conc[k][j][i] = Conc[k][j][i-1];
                        //if(i==mx-2)Conc[k][j][i+1] = Conc[k][j][i];
		}
		
		if( nvert[k][j][i]<1.1 && ( ( (user->bctype[2]==-1 || user->bctype[2]==-2) && j==1) || ( (user->bctype[3]==-1 || user->bctype[3]==-2) &&  j==my-2) ) && (i!=0 && i!=mx-1 && k!=0 && k!=mz-1)) {
		        if(mobile_bed && !immersed){
			double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
			double sb, sc; 
			Cmpnts Ua, Uc;
			Ua.x = Ua.y = Ua.z = 0;

			sb = 0.5/aj[k][j][i]/area;
			
			if(j==1) {
				sc = 2*sb + 0.5/aj[k][j+1][i]/area;
				Uc = ucat[k][j+1][i];
			}
			else {
				sc = 2*sb + 0.5/aj[k][j-1][i]/area;
				Uc = ucat[k][j-1][i];
			}
			
			double ni[3], nj[3], nk[3];
			Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
                        double nfx = nj[0], nfy = nj[1], nfz = nj[2];
			if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;

			wall_function (1./ren, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
			
			double utau = ustar[k][j][i];
			//Cmpnts uc;
			//if(j==1) {uc = ucat[k][j+1][i];}
			//else {uc = ucat[k][j][i];}

                         // Ali bed B.C.
                        PetscInt VanRijn=0, Fredsoe=1, non_equilibrium = 1;
                        //double U_bulk_vel = 1.0;//1.87;//0.593;//1.5;
                        double U_bulk_vel = U_Bulk;//1.87;//0.593;//1.5;
                        double ho = 3.;//8.;//.3;
                        double fluid_density = 1000.0;
                        //double D50 = 0.0005;//0.0007;
                        double D50 = d50;//0.0005;//0.0007;
                        //double Deltab = 10.*D50;//0.08;//2.*D50;
                        double Deltab = deltab;//10.*D50;//0.08;//2.*D50;
                        double particle_sg = 2.5;
                        double gi = 9.81;
                        double visk = 1.e-6; 
                        double Angle_repose = 45. * 3.1415 / 180.;
                        double CShields_param, Ctau_0, Custar_0, tau_0, tt, cb, bita, tmp;
                        double betta, alfa, gamma, nu_tt;
                        double sgd_tmp   = ( particle_sg - 1. ) * gi * D50;
                        double Dstar = D50 * pow(((particle_sg-1.) * gi / ( visk * visk )), 0.3333333333 );
                              
                        int shields = 0, soulsby = 1; // soulseby and Whitehouse methods: see J Geoph. Res. vol 115, 2010, Chou & Fringer

                        if(shields){

                        	if( Dstar <= 4. ) CShields_param = 0.24 / Dstar;
                           	if( Dstar > 4. && Dstar <= 10. ) CShields_param = 0.14 / pow( Dstar, 0.64 );
                         	if( Dstar > 10. && Dstar <= 20. ) CShields_param = 0.04 / pow( Dstar, 0.1 );
                        	if( Dstar > 20. && Dstar <= 150. ) CShields_param = 0.013 * pow( Dstar, 0.29 );
                            	if( Dstar > 150. ) CShields_param = 0.055;
                                  } 

                       if (soulsby) { CShields_param = 0.3/ ( 1. + 1.2 * Dstar) + 0.055 * (1. - exp (-0.02 * Dstar)) ; }

                         //double tita_x = atan ( - nfx / nfz );
                         //double tita_y = atan ( - nfy / nfz );
                         //double tita1  = uc.x * sin (tita_x) + uc.y * sin (tita_y); 
                         //double tita2  = sqrt (uc.x * uc.x + uc.y * uc.y);
                         //double tita3  = tita1 / tita2;
                         //if (fabs(tita3) >= 1. || fabs(tita3)<= -1.) { tmp =1.;} else {
                         //double tita   = asin (tita3);
                         //       tmp    = sin ( tita + Angle_repose ) / sin (Angle_repose);}
                         
                        //CShields_param = CShields_param * tmp;
                                       
        
	                
                               if(VanRijn) {
                                        Ctau_0 = CShields_param * sgd_tmp * fluid_density;
	                                    tau_0 =  fluid_density*(utau*U_bulk_vel)*(utau*U_bulk_vel);
		                       	tt = (tau_0 - Ctau_0)/Ctau_0;
                                        tt = PetscMax( tt, 0. );
                                        cb =  0.015 *(D50/Deltab)*pow( tt, 1.5 )/pow( Dstar,0.3);
                                           } 
      
              	               if(Fredsoe) {
                                        Ctau_0 = CShields_param; 
                                        tau_0 = (utau*U_bulk_vel)*(utau*U_bulk_vel); // the inlet vel. increase by 1.5 times
                                        tau_0 = tau_0/(sgd_tmp);
                                        tt = tau_0 - Ctau_0;
                                        tt = PetscMax( tt, 1.e-7 );
                                        bita = PI * 0.51 / 6.;
                                        tt= tt*tt*tt*tt;
                                        bita=bita*bita*bita*bita;
                                        cb = pow((1.+ bita/tt),-0.25)*PI/6.;
                                        cb = PetscMax (0.0, cb);
                                           }
                              if(non_equilibrium){
                                                  double delz = coor[k][j+1][i].z - coor[k][j][i].z;
												  double dely = coor[k][j+1][i].y - coor[k][j][i].y;
												  if(y_direction) delz = dely;
		                                  double nu  = 1./ren;
                                                  if(rans || les) {nu_t = PetscMax(0.0 , 0.5*(lnu_t[k][j][i]+lnu_t[k][j-1][i])); nu_tt = nu + sigma_phi * nu_t;} 
                                                  else {nu_t = nu; nu_tt = nu;}
                                                  
                                                         alfa = 1. - exp ( PetscMax ( 0.0, w_s * ( delz - Deltab / ho ) / nu_tt ) );
                                                         alfa = PetscMax (0.0, alfa);
                                                         betta = (1. - alfa) * ( w_s * delz / nu_tt); 
                                                         gamma = 1. - w_s * delz / nu_tt;
                                                  }
	                //PetscPrintf(PETSC_COMM_WORLD, "tmp_j=1   cb %e %e \n", tmp,cb );
	                //PetscPrintf(PETSC_COMM_WORLD, "u*, tau_0, Ctau_0, Cb %e %e %e %e\n", utau, tau_0, Ctau_0, cb );
	                //PetscPrintf(PETSC_COMM_WORLD, "Cb alfa betta gamma %e %e %e %e \n", cb, alfa, betta, gamma);
		        if(j==1 && !non_equilibrium) Conc[k][j][i] = cb;
		        if(j==1 && non_equilibrium) Conc[k][j][i] = Conc[k][j+1][i] + fabs(betta/gamma) * cb;
	
                        } else {
                         
                        if(j==1)Conc[k][j][i] = Conc[k][j+1][i];
                        //if(j==1)Conc[k][j-1][i] = Conc[k][j][i];
                        if(j==my-2)Conc[k][j][i] = Conc[k][j-1][i];
                        //if(j==my-2)Conc[k][j+1][i] = Conc[k][j][i];
                       }
		}
		
		if ( nvert[k][j][i] > 1.1 ) Conc[k][j][i]= 0.;
	}

	DAVecRestoreArray(da, user->lConc, &Conc);

        DALocalToLocalBegin(da, user->lConc, INSERT_VALUES, user->lConc);
        DALocalToLocalEnd(da, user->lConc, INSERT_VALUES, user->lConc);

        DAVecGetArray(da, user->lConc, &Conc);
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int flag=0, a=i, b=j, c=k;
		
		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
	
		if(flag) Conc[k][j][i] = Conc[c][b][a];
	}
	
        if(immersed){
	for(ibi=0; ibi<NumberOfBodies; ibi++)
	{
		extern IBMNodes *ibm_ptr;
		IBMNodes *ibm = &ibm_ptr[ibi];
		
		IBMListNode *current;
		current = user->ibmlist[ibi].head;
		while (current) {
			IBMInfo *ibminfo = &current->ibm_intp;
			int ni = ibminfo->cell;
			current = current->next;
			double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
			int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
			int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
			int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
			i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
			double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
			//double Kc = (K_Omega[kp1][jp1][ip1].x * sk1 + K_Omega[kp2][jp2][ip2].x * sk2 + K_Omega[kp3][jp3][ip3].x * sk3);
			double C = (Conc[kp1][jp1][ip1] * sk1 + Conc[kp2][jp2][ip2] * sk2 + Conc[kp3][jp3][ip3] * sk3);
			double ren = user->ren;
			if(levelset) ren = rho[k][j][i]/ mu[k][j][i];
				
			double nfx = ibm->nf_x[ni];
			double nfy = ibm->nf_y[ni];
			double nfz = ibm->nf_z[ni];
			Cmpnts uc = ucat[k][j][i];
		            
			if(wallfunction && ti>0) {
				//double utau = ustar[k][j][i];
                         // Ali bed B.C.
 
                        //int rigid = ibm->Rigidity[ni];
                        int mobile = ibm->Mobility[ni];
                        PetscInt Fringer = 1, MikiBC = 0, _MikiBC = 0;
			//if(mobile_bed && !rigid)
			//if(mobile_bed && mobile)
			//if(mobile_bed && ibm->Rigidity[ni]==0)
			if(mobile_bed && (ibm->nf_z[ni]> 0.99 || (sediment && ibm->Rigidity[ni]==0)))
                        {
                        double utau = ustar[k][j][i];
                       
                        PetscInt VanRijn=0, Fredsoe=1, non_equilibrium = 0;
                       
                        //if(LiveBed) non_equilibrium = 1;
                        
                        //double U_bulk_vel = 1.0;//1.87;//1.1;// Max Value0.593; //Ubulk=0.501m/s;
                        double U_bulk_vel = U_Bulk;//1.0;//1.87;//1.1;// Max Value0.593; //Ubulk=0.501m/s;
			double ho = 0.3;//8.0;//0.152;
                        double fluid_density = 1000.0;
                        //double D50 = 0.0005; //0.00025;
                        double D50 = d50;//0.0005; //0.00025;
                        //double Deltab = 10.*D50;//2.*D50;
                        double Deltab = deltab;//10.*D50;//2.*D50;
                        double particle_sg = 2.5;
                        double porosity = 0.45;
                        double gi = 9.81;
                        double visk = 1.e-6; 
                        double Angle_repose = 45. * 3.1415 / 180.;
                        double CShields_param, Ctau_0, Custar_0, tau_0, tt, cb, bita, tmp, coeff;
                        double betta, alfa, gamma;
                        double sgd_tmp   = ( particle_sg - 1. ) * gi * D50;
                        double Dstar = D50 * pow(((particle_sg-1.) * gi / ( visk * visk )), 0.3333333333 );
                		int shields = 0, soulsby = 1, old_method = 0, whitehouse = 1; // soulseby and Whitehouse methods: see J Geoph. Res. vol 115, 2010, Chou & Fringer
                        

                        if(shields){
                              	if( Dstar <= 4. ) CShields_param = 0.24 / Dstar;
	                        if( Dstar > 4. && Dstar <= 10. ) CShields_param = 0.14 / pow( Dstar, 0.64 );
                              	if( Dstar > 10. && Dstar <= 20. ) CShields_param = 0.04 / pow( Dstar, 0.1 );
                          	if( Dstar > 20. && Dstar <= 150. ) CShields_param = 0.013 * pow( Dstar, 0.29 );
                          	if( Dstar > 150. ) CShields_param = 0.055;
                                  } 
                        if (soulsby) { CShields_param = 0.3/ ( 1. + 1.2 * Dstar) + 0.055 * (1. - exp (-0.02 * Dstar)) ; }

                         /*double tita_x = atan ( - nfx / nfz );
                         double tita_y = atan ( - nfy / nfz );
                         double tita1  = uc.x * sin (tita_x) + uc.y * sin (tita_y); 
                         double tita2  = sqrt (uc.x * uc.x + uc.y * uc.y);
                         double tita3  = tita1 / tita2;
                         if (fabs(tita3) >= 1. || fabs(tita3)<= -1.) { tmp =1.;} else {
                         double tita   = asin (tita3);
                                tmp    = sin ( tita + Angle_repose ) / sin (Angle_repose);}
                         CShields_param *= tmp; */
                         
						if(!Fringer){
                                       
                               if(VanRijn) {
                                        Ctau_0 = CShields_param * sgd_tmp * fluid_density;
	                                tau_0 =  fluid_density*(utau*U_bulk_vel)*(utau*U_bulk_vel);
		                       	tt = (tau_0 - Ctau_0)/Ctau_0;
                                        tt = PetscMax( tt, 0. );
                                        cb =  0.015 *(D50/Deltab)*pow( tt, 1.5 )/pow( Dstar,0.3);
                                           } 
      
              	               if(Fredsoe) {
                                        Ctau_0 = CShields_param; 
                                        tau_0 = (utau*U_bulk_vel)*(utau*U_bulk_vel); // the inlet vel. increase by 1.5 times
                                        tau_0 = tau_0/(sgd_tmp);
                                        tt = tau_0 - Ctau_0;
                                        tt = PetscMax( tt, 1.e-7 );
                                        bita = PI * 0.51 / 6.;
                                        tt= tt*tt*tt*tt;
                                        bita=bita*bita*bita*bita;
                                        cb = pow((1.+ bita/tt),-0.25)*PI/6.;
                                        cb = PetscMax (0.0, cb);
                                           }
                               
                               if(non_equilibrium){
                                                  //double delz = coor[k][j+1][i].z - coor[k][j][i].z;
                                                  double delz = sc - sb;
		                                  double nu  = 1./ren;
                                                  if(rans || les) {nu_t = PetscMax(0.0 , 0.5*(lnu_t[k][j][i]+lnu_t[k][j-1][i])); nu_tt = nu + sigma_phi * nu_t;}
                                                  else {nu_t = nu; nu_tt = nu;}
                                                         /*alfa = 1. - exp ( w_s * PetscMax( 0., delz - ( Deltab / ho ) ) / nu_tt);
                                                         alfa = PetscMax ( 0.0, alfa);
                                                         alfa = PetscMin ( 1.0, alfa);
                                                         betta = (1. - alfa) * ( w_s * delz / nu_tt); 
                                                         gamma = 1. - w_s * delz / nu_tt;
                                                             coeff = PetscMin (fabs(betta/gamma), 1.0); */
                                                             coeff = PetscMin (fabs(w_s * delz / nu_tt), 1.0);
                                                  }
		                if(!non_equilibrium) Conc[k][j][i] = cb;
		               /* if(non_equilibrium){ 
                                    if(C > cb) Conc[k][j][i] = C - coeff * cb;
                                    if(C <= cb) Conc[k][j][i] = C + coeff * cb;
                                                   }*/
		                if(non_equilibrium) Conc[k][j][i] = C + coeff * cb;
		                //if(non_equilibrium) Conc[k][j][i] = Conc[k][j+1][i] + (betta/gamma) * cb;
						
                                //const double yplus_min = 0.25;
				//sb = PetscMax( sb, 1./ren * yplus_min / utau );	// prevent the case sb=0
				//double K, Ob;//, Om;
				// Wilcox pp.108-109
				//K = utau*utau/sqrt(0.09);
				//Ob = utau/sqrt(0.09)/(kappa*sb);
				//K_Omega[k][j][i].x = K;
				//K_Omega[k][j][i].y = Ob;

						} // if not-Fringer
						
		  if(Fringer)    { //This will consider the Pick up function and other stuff entraining sediment into the flow. 
                                           //to account for the concentration change effect on the mommentom---density current
						
                            Ctau_0 = CShields_param * sgd_tmp * fluid_density;
	                    tau_0 =  fluid_density*(utau*U_bulk_vel)*(utau*U_bulk_vel);
		            tt = (tau_0 - Ctau_0)/Ctau_0;
                            tt = PetscMax( tt, 0. );
                            double P_k = 0.00033 * pow( tt, 1.5 ) * pow( Dstar,0.3) * sqrt(sgd_tmp);

                            double C_extr;
                            if(C >= Conc[k][j][i]) {C_extr = PetscMax(C * (sc/(sc-sb)) - Conc[k][j][i] * (sb/(sc-sb)), 0.0);}
                            else {C_extr = PetscMax(Conc[k][j][i] * (sc/(sc-sb)) - C * (sb/(sc-sb)), 0.0);}
                            
		            double delz = coor[k][j+1][i].z - coor[k][j][i].z;
					double dely = coor[k][j+1][i].y - coor[k][j][i].y;
					if(y_direction) delz = dely;

                            double nu  = 1./ren;
                            if(rans || les) {
				             nu_t = PetscMax(0.0 , 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i])); 
					     nu_tt = nu + sigma_phi * nu_t;
			                    } else { 
						    nu_t = nu; 
			                            nu_tt = nu;
			                           }
			    double wbar = 0.5 * (ucat[k][j][i].z + ucat[k][j+1][i].z);
					if(y_direction)wbar = 0.5 * (ucat[k][j][i].y + ucat[k][j+1][i].y);
			    //double wbar = 0.0; 
			    double cbar = 0.5 * (Conc[k][j][i] + Conc[k][j+1][i]);
			    double LeoTerm = wbar * cbar - nu_tt * (Conc[k][j+1][i] - Conc[k][j][i])/delz - w_s * cbar;

                            Conc[k][j][i] = Conc[k][j][i] - user->dt * LeoTerm / delz + user->dt * (P_k - w_s * C_extr) / delz;

                            //Streeter Boundary Condition: (1/nu_t).partial(C)/Partial(z) = -P_k
                            //Conc[k][j][i] = C + P_k * (sc-sb) / nu_tt;
                            
                            if (Conc[k][j][i] > porosity) Conc[k][j][i] = porosity;
                            if (Conc[k][j][i] < 0.0) Conc[k][j][i] = 0.0;
                            if (Conc[k][j][i] < C) Conc[k][j][i] = C;
						 }

		         	}  //if mobile-bed
  
				else if(mobile_bed && ibm->Rigidity[ni]) {Conc[k][j][i] = C;}

                                else if(!mobile_bed && zero_grad){Conc[k][j][i] = C;}
                                

                                else if(!mobile_bed && !zero_grad && MikiBC){
                                  
                               // Miki_Hondzo's B.C. just on sediments
			       if(ibm->nf_z[ni]<0.8 || ((ibm->cent_x[ni]<45. && ibm->cent_z[ni]>6.98) || (ibm->cent_y[ni]<55.5 && ibm->cent_x[ni]>75.2 && ibm->cent_z[ni]>6.62)))
                                  {
                                    Conc[k][j][i] = C;
                                  } else {
                                  
                                   double utau = ustar[k][j][i];
                                   double U_bulk_vel = 0.16; //0.032;
                                   double ho = 0.16;
                                   double Re_tau = (utau * U_bulk_vel) * ho * 1000000.;

                                   double c_v = Depth_Averaged_Conc(Conc, nvert, i, j, k, lye);
                                    
                                   double c_s = c_v / ( 1.0 + 0.67 * exp ((Re_tau - 17300.) * (Re_tau - 17300.)/(-2.1 * 100000000.)));
                                   
	                           //PetscPrintf(PETSC_COMM_WORLD, "c_v, c_s Re %e %e %e\n", c_v, c_s, Re_tau);
   
                                   Conc[k][j][i] = c_s;}

                                }
                                else if(!mobile_bed && !zero_grad && _MikiBC){
                                  
                               // Miki_Hondzo's Newly developed B.C. for NO3 (flux based BC) on the sediment bed material
			       if(ibm->nf_z[ni]<0.9 || ibi)
                                  {
                                    Conc[k][j][i] = C;
                                  } else {
                                  
                                   double utau = ustar[k][j][i];
                                   double U_bulk_vel = 0.16667; 
                                   double ho = (coor[k][lye][i].z - coor[k][j][i].z) * 0.16;

                                   if(ho < 0.0) PetscPrintf(PETSC_COMM_WORLD, "\n error ****** Flow depth cannot be  negative **** error\n");

                                   double Re_tau = (utau * U_bulk_vel) * ho * 1000000.;

                                   double c_v = Depth_Averaged_Conc(Conc, nvert, i, j, k, lye);
                                    
	                           //PetscPrintf(PETSC_COMM_WORLD, "c_v, c_s Re %e %e %e\n", c_v, c_s, Re_tau);
		                   double delz = coor[k][j+1][i].z - coor[k][j][i].z;
                                   double nu  = 1./ren;
                                   if(rans || les) {
				      nu_t = PetscMax(0.0 , 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i])); 
				      nu_tt = nu + sigma_phi * nu_t;
			                           } else { 
						    nu_t = nu; 
			                            nu_tt = nu;
			                                  }
                                  double rho_b = 1.87; //=1.63 +- 0.49 g/cm3
                                  double BOM = 0.162;//=0.102 +- 0.12 g/cm2
                                  
                                  double u_star = utau * U_bulk_vel;
                                  
                                  if(ho = 0.0) PetscPrintf(PETSC_COMM_WORLD, "\n error ****** Flow depth is zero **** error\n");
                                  if(c_v = 0.0) PetscPrintf(PETSC_COMM_WORLD, "\n error ****** Averaged No3 concentration is zero **** error\n");

                                  //double Flux_No3 = 1.23 * pow(10.,-6.) * (u_star * c_v) * pow((u_star * ho * 1000000.),-(9./11.)) * pow((BOM/(ho*c_v)),(3./5.)) * pow((rho_b/c_v),(5./7.)) + 0.024; former simulation done with this
                                  //double Flux_No3 = 1.23 * pow(10.,-6.) * (u_star * c_v) * pow((u_star * ho * 1000000.),-(9./11.)) * pow((BOM/(ho*c_v)),(3./5.)) * pow((rho_b/c_v),(5./7.)) - 0.024;
                                  double Flux_No3 = 1.0 * pow(10.,-9.6) * pow((u_star * ho * 1000000.),-(9./11.)) * pow((10000.*BOM/(ho*c_v)),(3./5.)) * pow((1000000.*rho_b/c_v),(5./7.));
                                  
                                  if (nu_tt > 0.) {Conc[k][j][i] = PetscMax (0.0, (C-fabs(Flux_No3 * (sc - sb) / nu_tt)));
                                                   if(Conc[k][j][i] > C) Conc[k][j][i] = C;
                                                  }
                                  if (nu_tt <= 0.) Conc[k][j][i] = C;
                                                                  }

                                }
				  else {Conc[k][j][i] = C;}
 
                                }  //WallFunction
       
		                if(!wallfunction) {
				//const double yplus_min = 0.25;
				//double utau = ustar[k][j][i];
				//K_Omega[k][j][i].x = Kc * sb / sc;
				//sb = PetscMax( sb, 1./ren * yplus_min / utau );	// prevent the case sb=0
				//K_Omega[k][j][i].y = wall_omega(ren, sb);	
				
				//if ( K_Omega[k][j][i].x < 0 ) K_Omega[k][j][i].x = utau*utau/sqrt(0.09);
				Conc[k][j][i] = C;
		        	}
			if(user->bctype[4]==5 && k==1) Conc[k][j][i] = Conc[k-1][j][i];
		};
	}
}
	
	DAVecRestoreArray(fda, user->lCsi, &csi);
	DAVecRestoreArray(fda, user->lEta, &eta);
	DAVecRestoreArray(fda, user->lZet, &zet);
	DAVecRestoreArray(fda, user->lUcat, &ucat);

	if(rans || les) DAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
	DAVecRestoreArray(fda, Coor, &coor);

	DAVecRestoreArray(da, user->lNvert, &nvert);
	DAVecRestoreArray(da, user->lAj, &aj);
	DAVecRestoreArray(da, user->lUstar, &ustar);
	
	DAVecRestoreArray(da, user->lConc, &Conc);
	
	if(levelset) {
		DAVecRestoreArray(da, user->lDensity, &rho);
		DAVecRestoreArray(da, user->lMu, &mu);
	}
	DALocalToLocalBegin(da, user->lConc, INSERT_VALUES, user->lConc);
	DALocalToLocalEnd(da, user->lConc, INSERT_VALUES, user->lConc);
};

PetscErrorCode FormFunction_Conv_Diff(SNES snes, Vec Conc, Vec Rhs, void *ptr)
{
	UserCtx *user = (UserCtx*)ptr;

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
	PetscReal	***nvert;
	//Cmpnts2 ***rhs;
	PetscReal  ***rhs;

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
	

	DAGlobalToLocalBegin(user->da, Conc, INSERT_VALUES, user->lConc);
	DAGlobalToLocalEnd(user->da, Conc, INSERT_VALUES, user->lConc);
        
	Conv_Diff_BC(user);
 
        if(SuspendedParticles)SettlingVelocityCorrection(user);

	RHS_Conv_Diff(user, Rhs);
		
	
	VecAXPY(Rhs, -1/user->dt, Conc);
	VecAXPY(Rhs, 1/user->dt, user->Conc_o);
	
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->da, Rhs, &rhs);
	
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1) { //it was 0.1 for k-omega
			rhs[k][j][i] = 0.;
		}
		
		// couette
		//if ( user->bctype[3] == 12 && j==my-2 ) rhs[k][j][i].y = 0;
	
		//wall_omega
		/*
                if( i<=1 && user->bctype[0]==1 ) rhs[k][j][i].y = 0;
		if( i>=mx-2 && user->bctype[1]==1 ) rhs[k][j][i].y = 0;
		
		if( j<=1 && user->bctype[2]==1 ) rhs[k][j][i].y = 0;
		if( j>=my-2 && user->bctype[3]==1 ) rhs[k][j][i].y = 0;
		
		if( k==1 && user->bctype[4]==1 ) rhs[k][j][i].y = 0;
		if( k==mz-2 && user->bctype[5]==1 ) rhs[k][j][i].y = 0;
		*/

		if( i>=mx-2 && user->bctype[1]==1 ) rhs[k][j][i] = 0.;
		if( j<=1 && user->bctype[2]==1 ) rhs[k][j][i] = 0.;

		if( j<=1 && user->bctype[2]==1 ) rhs[k][j][i] = 0.;
		if( j>=my-2 && user->bctype[3]==1 ) rhs[k][j][i] = 0.;
		
		if( k==1 && user->bctype[4]==1 ) rhs[k][j][i] = 0.;
		if( k==mz-2 && user->bctype[5]==1 ) rhs[k][j][i] = 0.;

		// wall function k, omega
		if ( 	( i==1 && user->bctype[0] == -1 ) || ( i==mx-2 && user->bctype[1] == -1 ) ||
			( j==1 && user->bctype[2] == -1 ) || ( j==my-2 && user->bctype[3] == -1 ) ||
			( k==1 && user->bctype[4] == -1 ) || ( k==mz-2 && user->bctype[5] == -1 ) ||
			( i==1 && user->bctype[0] == -2 ) || ( i==mx-2 && user->bctype[1] == -2 ) ||
                        ( j==1 && user->bctype[2] == -2 ) || ( j==my-2 && user->bctype[3] == -2 ) ||
                        ( k==1 && user->bctype[4] == -2 ) || ( k==mz-2 && user->bctype[5] == -2 )
			) {
			//rhs[k][j][i].x = 0;    if you compute this value using an equation for the i==1 then the Rhs of it should be set to 0
			//rhs[k][j][i].y = 0;
			rhs[k][j][i] = 0.;
          		  }
	}
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->da, Rhs, &rhs);
	
	return(0);
}

int snes_convdiff_created=0;
Vec r_convdiff;
Mat J_convdiff;
SNES snes_convdiff;

void Solve_Conv_Diff(UserCtx *user)
{
	
	KSP ksp;
	PC pc;
	
	
	double norm;
	
	int bi=0;
	double tol=1.e-6;//5.e-5;//1.e-6
	
	if(!snes_convdiff_created) {
		snes_convdiff_created=1;
		
		VecDuplicate(user[bi].Conc, &r_convdiff);
		SNESCreate(PETSC_COMM_WORLD,&snes_convdiff);
		SNESSetFunction(snes_convdiff,r_convdiff,FormFunction_Conv_Diff,(void *)&user[bi]);
		MatCreateSNESMF(snes_convdiff, &J_convdiff);
		SNESSetJacobian(snes_convdiff,J_convdiff,J_convdiff,MatMFFDComputeJacobian,(void *)&user[bi]);
		SNESSetType(snes_convdiff, SNESTR);			//SNESTR,SNESLS	
		SNESSetMaxLinearSolveFailures(snes_convdiff,10000);
		SNESSetMaxNonlinearStepFailures(snes_convdiff,10000);		
		SNESKSPSetUseEW(snes_convdiff, PETSC_TRUE);
		SNESKSPSetParametersEW(snes_convdiff,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		SNESSetTolerances(snes_convdiff,PETSC_DEFAULT,tol,PETSC_DEFAULT,50,50000);
			
		SNESGetKSP(snes_convdiff, &ksp);
		KSPSetType(ksp, KSPGMRES);
		
		KSPGetPC(ksp,&pc);
		PCSetType(pc,PCNONE);
		
		int maxits=50/*20 or 10000*/;	double rtol=tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;
		KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
	}
	
	extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy);
	SNESMonitorSet(snes_convdiff,MySNESMonitor,PETSC_NULL,PETSC_NULL);
	
	PetscPrintf(PETSC_COMM_WORLD, "\nSolving Convection-Diffusion ...\n");
	
	SNESSolve(snes_convdiff, PETSC_NULL, user[bi].Conc);
	
	SNESGetFunctionNorm(snes_convdiff, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "\nConcentration SNES residual norm=%.5e\n\n", norm);

	DALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

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
	
	DAGlobalToLocalBegin(user->da, user->Conc, INSERT_VALUES, user->lConc);
	DAGlobalToLocalEnd(user->da, user->Conc, INSERT_VALUES, user->lConc);

	PetscReal ***conc, ***lconc;
	
	DAVecGetArray(user->da, user->Conc, &conc);
	DAVecGetArray(user->da, user->lConc, &lconc);

	if(periodic)
	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		int flag=0, a=i, b=j, c=k;

		if(i_periodic && i==0) a=mx-2, flag=1;
		else if(i_periodic && i==mx-1) a=1, flag=1;
		
		if(j_periodic && j==0) b=my-2, flag=1;
		else if(j_periodic && j==my-1) b=1, flag=1;
		
		if(k_periodic && k==0) c=mz-2, flag=1;
		else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		
		
		if(flag) conc[k][j][i] = lconc[c][b][a];
	}
	DAVecRestoreArray(user->da, user->Conc, &conc);
	DAVecRestoreArray(user->da, user->lConc, &lconc);

	DAGlobalToLocalBegin(user->da, user->Conc, INSERT_VALUES, user->lConc);
	DAGlobalToLocalEnd(user->da, user->Conc, INSERT_VALUES, user->lConc);

	Conv_Diff_BC(user);

        if(SuspendedParticles) SettlingVelocityCorrectionBack(user);

	DALocalToGlobal(user->da, user->lConc, INSERT_VALUES, user->Conc);
};

