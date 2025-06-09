#include "variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c, L_dim, Nseg;
extern PetscInt  bed_cell_depth, rstart, periodic_morpho, no_ibm_search, Paraview, levelset;
extern char path[256];

PetscErrorCode ibm_surface_out(IBMNodes *ibm, PetscInt ti,
			       PetscInt ibi)
{
  PetscPrintf(PETSC_COMM_WORLD, "Write surface file for ibm%i from ibm_surface_out\n", ibi);
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == ti) {
      FILE *f;
      char filen[80];
      sprintf(filen, "%s/surface%3.3d_%2.2d.dat",path,ti,ibi);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"TRIANGLES\", N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", ibm->n_v, ibm->n_elmt);
      for (i=0; i<ibm->n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }
  return(0);
}

PetscErrorCode ibm_read(IBMNodes *ibm1, IBMNodes *ibm2)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal	t, dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscReal     cl=30.;

  double xt;
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata0", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file")
    n_v =0;
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%le", &xt);
    ibm1->n_v = n_v/2;
    ibm2->n_v = n_v/2;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->urm1));

    MPI_Bcast(&(ibm2->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->urm1));

    for (i=0; i<n_v; i++) {
      fscanf(fd, "%le %le %le %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i], &t, &t, &t);
      if(no_ibm_search){
		x_bp[i] = x_bp[i];// / cl;// 28.
		y_bp[i] = y_bp[i];// / cl;
		z_bp[i] = z_bp[i];// / cl;
      }
	  else{
		x_bp[i] = x_bp[i] / cl;// 28.
		y_bp[i] = y_bp[i] / cl;
		z_bp[i] = z_bp[i] / cl;
	  }
    }
/*     ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */   

    PetscReal temp;
    for (i=0; i<n_v/2; i++) {
      temp = y_bp[i];
      ibm1->x_bp0[i] = z_bp[i];
      ibm1->y_bp0[i] = x_bp[i]-CMy_c;
      ibm1->z_bp0[i] = -temp+CMz_c;

      ibm1->x_bp[i] = z_bp[i];
      ibm1->y_bp[i] = x_bp[i]-CMy_c;
      ibm1->z_bp[i] = -temp+CMz_c;

      ibm1->x_bp_o[i] = z_bp[i];
      ibm1->y_bp_o[i] = x_bp[i]-CMy_c;
      ibm1->z_bp_o[i] = -temp+CMz_c;   

      ibm1->u[i].x = 0.;
      ibm1->u[i].y = 0.;
      ibm1->u[i].z = 0.;
      
      ibm1->uold[i].x = 0.;
      ibm1->uold[i].y = 0.;
      ibm1->uold[i].z = 0.;
      
      ibm1->urm1[i].x = 0.;
      ibm1->urm1[i].y = 0.;
      ibm1->urm1[i].z = 0.;   
    }

    for (i=0; i<n_v/2; i++) {
      temp = y_bp[i+n_v/2];
      ibm2->x_bp0[i] = z_bp[i+n_v/2];
      ibm2->y_bp0[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp0[i] = -temp+CMz_c;

      ibm2->x_bp[i] = z_bp[i+n_v/2];
      ibm2->y_bp[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp[i] = -temp+CMz_c;

      ibm2->x_bp_o[i] = z_bp[i+n_v/2];
      ibm2->y_bp_o[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp_o[i] = -temp+CMz_c;

      ibm2->u[i].x = 0.;
      ibm2->u[i].y = 0.;
      ibm2->u[i].z = 0.;
      
      ibm2->uold[i].x = 0.;
      ibm2->uold[i].y = 0.;
      ibm2->uold[i].z = 0.;
      
      ibm2->urm1[i].x = 0.;
      ibm2->urm1[i].y = 0.;
      ibm2->urm1[i].z = 0.;   
    }

/*     ibm1->x_bp=ibm1->x_bp0;ibm1->y_bp=ibm1->y_bp0;ibm1->z_bp=ibm1->z_bp0; */
/*     ibm2->x_bp=ibm2->x_bp0;ibm2->y_bp=ibm2->y_bp0;ibm2->z_bp=ibm2->z_bp0; */
/*     ibm1->x_bp_o=ibm1->x_bp0;ibm1->y_bp_o=ibm1->y_bp0;ibm1->z_bp_o=ibm1->z_bp0; */
/*     ibm2->x_bp_o=ibm2->x_bp0;ibm2->y_bp_o=ibm2->y_bp0;ibm2->z_bp_o=ibm2->z_bp0; */

    MPI_Bcast(ibm1->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    fscanf(fd, "%i\n", &n_elmt);
    ibm1->n_elmt = n_elmt/2;
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm2->n_elmt = n_elmt/2;
 
    PetscPrintf(PETSC_COMM_WORLD, "elmts , %d \n", n_elmt);

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

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
    
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_z));

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_z));

    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i %i %i\n", nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
    }


    fclose(fd);

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

      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
		(((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) 
      {
		if (nf_z[i]>0) {
			ns_x[i] = 1.;     
			ns_y[i] = 0.;     
			ns_z[i] = 0. ;
			
			nt_x[i] = 0.;
			nt_y[i] = 1.;
			nt_z[i] = 0.;
	} 	else {
			ns_x[i] = -1.;     
			ns_y[i] = 0.;     
			ns_z[i] = 0. ;
			
			nt_x[i] = 0.;
			nt_y[i] = -1.;
			nt_z[i] = 0.;
		}	  
      } else {
			ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
			ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
			ns_z[i] = 0. ;
			
			nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
			nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
			nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
		}
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
    }

    for (i=0; i<n_elmt/2; i++) {
      ibm1->nv1[i]=nv1[i];
      ibm1->nv2[i]=nv2[i];
      ibm1->nv3[i]=nv3[i];

      ibm1->nf_x[i]=  nf_z[i];
      ibm1->nf_y[i]=  nf_x[i];
      ibm1->nf_z[i]= -nf_y[i];

      ibm1->ns_x[i]=  ns_z[i];
      ibm1->ns_y[i]=  ns_x[i];
      ibm1->ns_z[i]= -ns_y[i];

      ibm1->nt_x[i]=  nt_z[i];
      ibm1->nt_y[i]=  nt_x[i];
      ibm1->nt_z[i]= -nt_y[i];
    }

    for (i=0; i<n_elmt/2; i++) {
      ibm2->nv1[i]=nv1[i+n_elmt/2]-n_v/2;
      ibm2->nv2[i]=nv2[i+n_elmt/2]-n_v/2;
      ibm2->nv3[i]=nv3[i+n_elmt/2]-n_v/2;

      ibm2->nf_x[i]=  nf_z[i+n_elmt/2];
      ibm2->nf_y[i]=  nf_x[i+n_elmt/2];
      ibm2->nf_z[i]= -nf_y[i+n_elmt/2];

      ibm2->ns_x[i]=  ns_z[i+n_elmt/2];
      ibm2->ns_y[i]=  ns_x[i+n_elmt/2];
      ibm2->ns_z[i]= -ns_y[i+n_elmt/2];

      ibm2->nt_x[i]=  nt_z[i+n_elmt/2];
      ibm2->nt_y[i]=  nt_x[i+n_elmt/2];
      ibm2->nt_z[i]= -nt_y[i+n_elmt/2];

    }

    MPI_Bcast(ibm1->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm1->n_v = n_v/2;

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->urm1));

    MPI_Bcast(&(ibm2->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->urm1));

    for (i=0; i<n_v/2; i++) {

      ibm1->u[i].x = 0.;
      ibm1->u[i].y = 0.;
      ibm1->u[i].z = 0.;
      
      ibm1->uold[i].x = 0.;
      ibm1->uold[i].y = 0.;
      ibm1->uold[i].z = 0.;
      
      ibm1->urm1[i].x = 0.;
      ibm1->urm1[i].y = 0.;
      ibm1->urm1[i].z = 0.;   

      ibm2->u[i].x = 0.;
      ibm2->u[i].y = 0.;
      ibm2->u[i].z = 0.;
      
      ibm2->uold[i].x = 0.;
      ibm2->uold[i].y = 0.;
      ibm2->uold[i].z = 0.;
      
      ibm2->urm1[i].x = 0.;
      ibm2->urm1[i].y = 0.;
      ibm2->urm1[i].z = 0.;   

    }

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);


    MPI_Bcast(ibm1->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm1->n_elmt = n_elmt/2;
    ibm2->n_elmt = n_elmt/2;
    
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_z));

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_z));

    i=10;
      PetscPrintf(PETSC_COMM_SELF, " ibm1 xbp %d %le %le %le\n", i, ibm1->z_bp0[i], ibm1->z_bp_o[i], ibm1->z_bp[i]);
/*     } */

    MPI_Bcast(ibm1->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

  }

  return(0);
}



PetscErrorCode ibm_read_ucd(IBMNodes *ibm, PetscInt ibi, PetscReal CMx_c=0.0, PetscReal CMy_c=0.0, \
	PetscReal CMz_c=0.0, PetscReal x_r=0.0, PetscReal y_r=0.0, PetscReal z_r=0.0, PetscReal XYangle=0.0)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;

  char   ss[20];
  char string[128];

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
		//PetscPrintf(PETSC_COMM_SELF, "From ucd read, IBM%i has CMx_c=%f, CMy_c=%f, CMz_c=%f, x_r=%f, y_r=%f, z_r=%f\n", ibi, CMx_c, CMy_c, CMz_c, x_r, y_r, z_r);

    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "\nBegin reading ibmdata for IB %i\n", ibi);
    char filen[80];  
    sprintf(filen,"%s/ibmdata%2.2d" , path, ibi);
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file")
    n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      
      fscanf(fd, "%i %i %i %i %i",&n_v,&n_elmt,&ii,&ii,&ii);
      PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d\n",n_v, n_elmt);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;      
      
      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
	    
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));	// added by seokkoo 03.04.2009
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      x_bp = ibm->x_bp;	// seokkoo
      y_bp = ibm->y_bp;
      z_bp = ibm->z_bp;
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
      
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
      PetscReal cl = 1.;
      PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);




     
	  //Adjust for center of mass translation in fsi_rot.dat input file
		PetscPrintf(PETSC_COMM_WORLD, "Before transation2, x_r=%f, y_r=%f, z_r= %f\n", x_r, y_r, z_r);
		x_r += CMx_c;
		y_r += CMy_c;
		z_r += CMz_c;
		PetscPrintf(PETSC_COMM_WORLD, "After transation2, x_r=%f, y_r=%f, z_r= %f\n", x_r, y_r, z_r);
		PetscReal x_temp = 0;
	  
	


    for (i=0; i<n_v; i++) {
		fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	
		if(no_ibm_search){
			x_bp[i] = x_bp[i];///cl*L_dim + CMx_c;//0.25 ;// 24.;	
			y_bp[i] = y_bp[i];///cl*L_dim + CMy_c ;//+ ibi*2.;//5.;//2.;//8.;//6.;//2.   ;// 24.;
			z_bp[i] = z_bp[i];///cl*L_dim + CMz_c ;//+ ibi*1.5;//2.;//8.;//15.;//2.   ;// 24.;
		}
		else {
			x_bp[i] = x_bp[i]/cl*L_dim + CMx_c/cl*L_dim;//0.25 ;// 24.;	
			y_bp[i] = y_bp[i]/cl*L_dim + CMy_c/cl*L_dim;//+ ibi*2.;//5.;//2.;//8.;//6.;//2.   ;// 24.;
			z_bp[i] = z_bp[i]/cl*L_dim + CMz_c/cl*L_dim;//+ ibi*1.5;//2.;//8.;//15.;//2.   ;// 24.;
			//if (i == 10) PetscPrintf(PETSC_COMM_WORLD, "Before Rotation of %f radians, x_bp[10] = %f and y_bp[10] = %f\n", XYangle, x_bp[10], y_bp[10]);
			x_temp = x_r + (x_bp[i]-x_r)*cos(XYangle) - (y_bp[i]-y_r)*sin(XYangle);
			y_bp[i] = y_r + (x_bp[i]-x_r)*sin(XYangle) + (y_bp[i]-y_r)*cos(XYangle);
			x_bp[i] = x_temp;
		//	if (i == 10) PetscPrintf(PETSC_COMM_WORLD, "After Rotation with xr=%f and yr=%f, x_bp[10] = %f and y_bhp[10] = %f\n", x_r, y_r, x_bp[10], y_bp[10]);

		}

        ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

	ibm->x_bp0[i] = x_bp[i];
	ibm->y_bp0[i] = y_bp[i];
	ibm->z_bp0[i] = z_bp[i];

	ibm->x_bp_o[i] = x_bp[i];
	ibm->y_bp_o[i] = y_bp[i];
	ibm->z_bp_o[i] = z_bp[i];

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
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

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

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
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
      
      // Added 6/4/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added
      
      // Added 1/11/2010 ali
      PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->Bvel));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_grad_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_grad_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_grad_z)); //Hossein
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->v_grad_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->v_grad_y));

      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_grad_x)); //Hossein
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_grad_z)); //Hossein

      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_grad_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_grad_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_grad_z)); //Hossein

      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->sbb));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->scc));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_max));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->v_max));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_max)); //Hossein
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_max));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_u));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_v));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_w));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vx));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vy));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vz));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Shvel));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Shvel_v));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->BShS));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->BShS_v));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->BCShS));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_flux_12));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_flux_13));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_flux_23));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->SCont));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Delz));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Dely)); //Hossein
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dVol));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dzbp));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->elmt_depth));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z_old));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_zl));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z_AVE));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y_AVE)); //Hossein
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y_old)); //Hossein
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_yl)); //Hossein

      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->atke_old));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->netflux_old));

      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Dflux));  //Hossein
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Eflux)); //Hossein
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->qflux)); //Hossein

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_l));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->deltz_p_ds));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->deltz_p_us));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->A_us));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->A_ds));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->max_bed_angle));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->C));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_ave));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->k_ave));
      PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->Rigidity));
      PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->Mobility));
      if(periodic_morpho) PetscMalloc(Nseg*sizeof(PetscReal), &(ibm->SedFlux));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->fluid_den));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->fluid_vis));
      
     // end added 
      
     PetscMalloc(n_v*sizeof(Node2Cell), &ibm->n2c); 
     PetscMalloc(n_elmt*sizeof(Cell2Cell), &ibm->c2c); 
     PetscMalloc(n_elmt*sizeof(Cell2Cells), &ibm->c2cs); 
     
      
	//seokkoo begin
	{	//only for rank 0
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->_nv1);
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->_nv2);
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->_nv3);
		
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->_x_bp);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->_y_bp);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->_z_bp);
	}
	
	PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->count);
	PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->local2global_elmt);
	
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
	PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);
	//seokkoo end

      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, &ss, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
	      
		// seokkoo
	      ibm->_nv1[i] = nv1[i];
	      ibm->_nv2[i] = nv2[i];
	      ibm->_nv3[i] = nv3[i];
	      // seokkoo
	      ibm->_x_bp[i] = ibm->x_bp[i];
	      ibm->_y_bp[i] = ibm->y_bp[i];
	      ibm->_z_bp[i] = ibm->z_bp[i];

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
      

      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
		(((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
		ns_x[i] = 1.;     
		ns_y[i] = 0.;     
		ns_z[i] = 0. ;
		
		nt_x[i] = 0.;
		nt_y[i] = 1.;
		nt_z[i] = 0.;
      } 
	  else {
		ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
		ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
		ns_z[i] = 0. ;
		
		nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
		nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
		nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
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

    PetscInt ti=0;
   

  extern PetscInt Paraview;
  FILE *f;
	if(Paraview) {

		// vtk file name
		PetscInt n_cells=3;
		  //PetscInt rank,i;
		PetscPrintf(PETSC_COMM_WORLD, "Write surface file for ibm%i from ibm_io for Paraview1\n", ibi);
		//char filen[80];
		sprintf(filen, "%s/surface%3.3d_%2.2d.vtk",path,ibi,ti);
		f = fopen(filen, "w"); // open file

		PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
		PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
		PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
		PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
	   
		//    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  5523993 double\n");
		PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibm->n_v);
		for (i=0; i<ibm->n_v; i++) {
		  PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->x_bp[i],ibm->y_bp[i],ibm->z_bp[i]);
		}

		PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibm->n_elmt, (n_cells+1)*ibm->n_elmt);
		for (i=0; i<ibm->n_elmt; i++) {
		  PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);
		}
		PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibm->n_elmt);
		for (i=0; i<ibm->n_elmt; i++) {
		  PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
		}

		PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", ibm->n_v);
		PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
		for (i=0; i<ibm->n_v; i++) {
		  PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->u[i].x,ibm->u[i].y,ibm->u[i].z);
		}
		PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n", ibm->n_elmt);
		PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nf float\n");
		for (i=0; i<ibm->n_elmt; i++) {
		  PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nf_x[i],ibm->nf_y[i],ibm->nf_z[i]);
		}

		PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS p double\n");    
		PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
		for (i=0; i<ibm->n_elmt; i++) {
		  PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", ibm->nt_y[i]);
		}
		fclose(f);
	} else {

   //FILE *f;
    //char filen[80];
    sprintf(filen, "%s/surface%3.3d_%2.2d_nf.dat",path,ti,ibi);
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
	  
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
    
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

    for (i=0; i<n_v; i++) {
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

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    //Added 4/1/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

    // Added 1/11/2010 ali
    PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->Bvel));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_grad_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_grad_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_grad_z)); //Hossein

    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->v_grad_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->v_grad_y));

    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_grad_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_grad_z));

    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_grad_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_grad_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_grad_z)); //Hossein

    //PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->z_grad_x));
    //PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->z_grad_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->sbb));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->scc));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->u_max));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->v_max));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_max));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->c_max));
    //PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->z_max));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_u));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_v));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_w));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vx));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vy));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->vz));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Shvel));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Shvel_v));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->BShS));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->BShS_v));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->BCShS));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_flux_12));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_flux_13));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->F_flux_23));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->SCont));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Delz));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Dely)); //Hossein
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dVol));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->dzbp));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->elmt_depth));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z_old));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z_AVE));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_zl));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y_AVE)); //Hossein
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y_old)); //Hossein
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_yl)); //Hossein

    // PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->atke));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->atke_old));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->netflux_old));

    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Dflux));  //Hossein
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Eflux)); //Hossein
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->qflux)); //Hossein

    //PetscMalloc(100000*sizeof(PetscReal), &(ibm->dtime));
    //PetscMalloc(100000*sizeof(PetscReal), &(ibm->time_bedchange));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_l));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->deltz_p_ds));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->deltz_p_us));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->A_us));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->A_ds));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->max_bed_angle));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->C));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->w_ave));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->k_ave));
    PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->Rigidity));
    PetscMalloc(n_elmt*sizeof(PetscInt), &(ibm->Mobility));
    if(periodic_morpho) PetscMalloc(Nseg*sizeof(PetscReal), &(ibm->SedFlux));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->fluid_den));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->fluid_vis));
    // end added


     PetscMalloc(n_v*sizeof(Node2Cell), &ibm->n2c); 
     PetscMalloc(n_elmt*sizeof(Cell2Cell), &ibm->c2c);
     PetscMalloc(n_elmt*sizeof(Cell2Cells), &ibm->c2cs); 
 
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
    
		//seokkoo
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->count);
		PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->local2global_elmt);
		//seokkoo
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->shear);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->mean_shear);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress1);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress2);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->reynolds_stress3);
		PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->pressure);

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
  }
    PetscPrintf(PETSC_COMM_WORLD, "Finished reading ucd file!\n\n");
 
  //MPI_Barrier(PETSC_COMM_WORLD); 
  return(0);
}

PetscErrorCode ibm_read_ucd_old(IBMNodes *ibm)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal	t, dr;
  PetscInt 	temp;
  double xt;
  char string[128];
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata", "r");
    if (!fd) SETERRQ(1, "Cannot open IBM node file")
    n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);

      fscanf(fd, "%i %i %i %i %i\n", &n_v, &n_elmt, &temp, &temp, &temp);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;

      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

      
      PetscReal cl = 1.;

      PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);
      
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &temp, &x_bp[i], &y_bp[i], &z_bp[i]);
	x_bp[i] = (x_bp[i]    ) / cl;
	y_bp[i] = (y_bp[i]+6. ) / cl ;
	z_bp[i] = (z_bp[i]+15.) / cl;
	
/* 	ibm->x_bp[i] = x_bp[i]; */
/* 	ibm->y_bp[i] = y_bp[i]; */
/* 	ibm->z_bp[i] = z_bp[i]; */

	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

	ibm->u[i].x = 0.;
	ibm->u[i].y = 0.;
	ibm->u[i].z = 0.;
      }
      ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp;

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

	
      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      char str[20];
      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &temp, &temp, str, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;

      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

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

      
    }
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

    /*     for (i=0; i<n_elmt; i++) { */
    /*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
    /*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    for (i=0; i<ibm->n_v; i++) {
      ibm->x_bp[i] = ibm->x_bp0[i];
      ibm->y_bp[i] = ibm->y_bp0[i];
      ibm->z_bp[i] = ibm->z_bp0[i];

      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}


PetscErrorCode read_bed_cell_depth(IBMNodes *ibm, PetscInt tistart)
{
  PetscInt	rank;
  PetscReal	*CellDepth;
  PetscInt	i,ii;
  PetscInt      n_elmt=ibm->n_elmt;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ ib_cell_depth\n");
    char filen[80];  
    sprintf(filen,"ib_cell_depth.dat");
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot open cell_depth file")
    
    if(fd) {
       for (i=0; i<n_elmt; i++) {

	   fscanf(fd, "%le %d\n", &(ibm->elmt_depth[i]), &(ibm->Rigidity[i]));
          
	   }
      fclose(fd);
    }
             
    MPI_Bcast(ibm->elmt_depth, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->Rigidity, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    FILE *f;
    
    sprintf(filen, "surface_elmt_depth_%5.5d.dat",tistart);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,elmt_depth,rigidity\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T=\"TRIANGLES\", N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)\n",ibm->n_v, ibm->n_elmt);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<ibm->n_v; i++) {	
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
    }
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->elmt_depth[i]);
    }
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", ibm->Rigidity[i]);
    }

      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }

    fclose(f);

  }
  else if (rank) {
    MPI_Bcast(ibm->elmt_depth, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->Rigidity, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
                 }

  return(0);
}


#define CROSS(dest, v1, v2) \
	dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
	dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
	dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define SUB(dest, v1, v2) \
	dest[0] = v1[0] - v2[0]; \
	dest[1] = v1[1] - v2[1]; \
	dest[2] = v1[2] - v2[2];


PetscErrorCode calc_ibm_volumeFlux(IBMNodes *ibm,PetscReal delti, PetscReal *VolumeFlux )
{
  PetscInt   i,ii, n_elmt=ibm->n_elmt;
  PetscInt n1e, n2e, n3e;
  PetscReal  Vol=0.,p1[3],p2[3],p3[3],p4[3],p5[3],p6[3],sign,nf[3];
  PetscReal  edge1[3],edge2[3],edge3[3],edge2c3[3],volTH,cent[3],cent_o[3],dir[3];

  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];

    nf[0]=ibm->nf_x[i];nf[1]=ibm->nf_y[i];nf[2]=ibm->nf_z[i];

    p1[0]=ibm->x_bp[n1e];p1[1]=ibm->y_bp[n1e];p1[2]=ibm->z_bp[n1e];
    p2[0]=ibm->x_bp[n2e];p2[1]=ibm->y_bp[n2e];p2[2]=ibm->z_bp[n2e];
    p3[0]=ibm->x_bp[n3e];p3[1]=ibm->y_bp[n3e];p3[2]=ibm->z_bp[n3e];

    p4[0]=ibm->x_bp_o[n1e];p4[1]=ibm->y_bp_o[n1e];p4[2]=ibm->z_bp_o[n1e];
    p5[0]=ibm->x_bp_o[n2e];p5[1]=ibm->y_bp_o[n2e];p5[2]=ibm->z_bp_o[n2e];
    p6[0]=ibm->x_bp_o[n3e];p6[1]=ibm->y_bp_o[n3e];p6[2]=ibm->z_bp_o[n3e];

    for (ii=0; ii<3; ii++) {
      cent[ii]  =(p1[ii]+p2[ii]+p3[ii])/3.;
      cent_o[ii]=(p4[ii]+p5[ii]+p6[ii])/3.;
    }
    
    // calculate volume flux
    SUB(dir,cent,cent_o);
    sign=DOT(dir,nf);
    if (fabs(sign)>1e-15) 
      sign /=fabs(sign);
    else
      sign =0.;

    SUB(edge1,p4,p1);
    SUB(edge2,p4,p2);
    SUB(edge3,p4,p3);
    CROSS(edge2c3,edge2,edge3);
    volTH=DOT(edge1,edge2c3);
    
    Vol +=sign*fabs(volTH/6.)/delti; //user->dt;

    SUB(edge1,p5,p4);
    SUB(edge2,p5,p2);
    SUB(edge3,p5,p3);
    CROSS(edge2c3,edge2,edge3);
    volTH=DOT(edge1,edge2c3);
    
    Vol +=sign*fabs(volTH/6.)/delti; //user->dt;

    SUB(edge1,p6,p5);
    SUB(edge2,p6,p4);
    SUB(edge3,p6,p3);
    CROSS(edge2c3,edge2,edge3);
    volTH=DOT(edge1,edge2c3);
    
    Vol +=sign*fabs(volTH/6.)/delti; //user->dt;
  }
  *VolumeFlux = Vol;
  //user->FluxIntpSum=Vol;
  PetscPrintf(PETSC_COMM_WORLD, "Volume Flux %e\n", Vol);

  return(0);
}

PetscErrorCode bed_change_VTKOut(IBMNodes *ibm, PetscInt ti)
{
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "bed_change_%5.5d.vtk",ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
   
    //    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  5523993 double\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibm->n_v);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->x_bp[i],ibm->y_bp[i],ibm->z_bp[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibm->n_elmt, (n_cells+1)*ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n", ibm->n_elmt);
    PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nf float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->Bvel[i].x,ibm->Bvel[i].y,ibm->Bvel[i].z);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS c double\n");    
    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f %f\n",ibm->BShS[i],ibm->BCShS[i],ibm->SCont[i],ibm->Rigidity[i]);
    }
    fclose(f);
  }

  return(0);
}


PetscErrorCode ibm_surface_VTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti)
{
	PetscPrintf(PETSC_COMM_WORLD, "Write surface file for ibm%i from ibm_io for VTK Out\n", ibi);
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surface%3.3d_%2.2d.vtk", ibi,ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");

    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibm->n_v);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->x_bp[i],ibm->y_bp[i],ibm->z_bp[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibm->n_elmt, (n_cells+1)*ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", ibm->n_v);
    PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->u[i].x,ibm->u[i].y,ibm->u[i].z);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n", ibm->n_elmt);
    PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nf float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nf_x[i],ibm->nf_y[i],ibm->nf_z[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS p double\n");    
    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", ibm->nt_y[i]);
    }
    fclose(f);
  }
  return(0);
}


