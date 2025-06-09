#include "variables.h"
#include <algorithm>
#include <math.h>
#include "canopy.h"
//#include "vtk.h"

using namespace std;


CnpyNodes *Canopy;

extern PetscReal  L_dim;
static  PetscReal cnpyCd = 0.15;


void write_output_file2(UserCtx *user, Vec V, char *file)
{
  PetscViewer viewer;
  char s[256];

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file, FILE_MODE_WRITE, &viewer);
  VecView(V, viewer);
  PetscViewerDestroy(viewer);

  if(my_rank) {
    sprintf(s, "%s.info", file);
    unlink(s);
  }

  //sprintf(s, "%s/%s", path, file);

  //PetscViewerBinaryOpen(PETSC_COMM_WORLD, s, FILE_MODE_READ, &viewer);
  //VecLoadIntoVector(viewer, (user->cnpyF));
  //PetscViewerDestroy(viewer);

  //PetscViewer viewer;
  //char filen[128];

  //sprintf(filen, "CnpyForce%06d_%1d.dat", ti, user->_this);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  //PetscInt N;

  //VecGetSize(user->cnpyF, &N);
  //PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  //VecLoadIntoVector(viewer, (user->cnpyF));
  //PetscViewerDestroy(viewer);

  //PetscBarrier(PETSC_NULL);

  /*return 0;*/

}


PetscErrorCode canopy_read_ucd(CnpyNodes *cnpy, int icnpy)
{
  /* This function reads the canopy ucd files (canopy_****.ucd)
   * and compute the center of each element */

  FILE *fd;
  char filen[80];

  int           n_v , n_elmt ;
  int           *nv1 , *nv2 , *nv3 ;
  PetscReal     *x_bp, *y_bp, *z_bp ;
  PetscReal     *nf_x, *nf_y, *nf_z;
  PetscReal     dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
//  PetscReal     *nt_x, *nt_y, *nt_z;
//  PetscReal     *ns_x, *ns_y, *ns_z;
//  PetscReal     *dA;
  int   ii;

  PetscPrintf(PETSC_COMM_WORLD,"\n \n ");
  PetscPrintf(PETSC_COMM_WORLD,"Reading Canopy ucd file %i: \n",icnpy);

  sprintf(filen,"./canopy_%04d.ucd",icnpy);
  PetscPrintf(PETSC_COMM_WORLD,filen);
  PetscPrintf(PETSC_COMM_WORLD,"\n");
//  PetscPrintf(PETSC_COMM_WORLD,"\n \n ");

  fd = fopen(filen,"r");

  PetscReal cl = 1.;
  PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD,"Characteristic Length is %f3: \n",cl);

  if (fd){
    char string[128];
    fgets(string,128,fd);
    fgets(string,128,fd);
    fgets(string,128,fd);
    fscanf(fd, "%i %i %i %i %i",&n_v,&n_elmt,&ii,&ii,&ii);

    cnpy->n_v = n_v;
    cnpy->n_elmt = n_elmt;


    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(cnpy->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(cnpy->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(cnpy->z_bp));

    /*for (int i=0;i<n_v;i++){
      fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);
      
      cnpy->x_bp[i] = x_bp[i]/cl*L_dim;
      cnpy->y_bp[i] = y_bp[i]/cl*L_dim;
      cnpy->z_bp[i] = z_bp[i]/cl*L_dim;
    }*/


    double maxtmp=-1000.;
    for (int i=0;i<n_v;i++){
      fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);
      cnpy->x_bp[i] = x_bp[i]/cl*L_dim;
      cnpy->y_bp[i] = y_bp[i]/cl*L_dim;
      cnpy->z_bp[i] = z_bp[i]/cl*L_dim;
      maxtmp = PetscMax(maxtmp,cnpy->z_bp[i]);
    }
    double rcvbuff;
    MPI_Allreduce(&maxtmp, &rcvbuff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    cnpy->maxHeight=rcvbuff;



    PetscMalloc(n_elmt*sizeof(int), &nv1);
    PetscMalloc(n_elmt*sizeof(int), &nv2);
    PetscMalloc(n_elmt*sizeof(int), &nv3);

    PetscMalloc(n_elmt*sizeof(int), &(cnpy->nv1));
    PetscMalloc(n_elmt*sizeof(int), &(cnpy->nv2));
    PetscMalloc(n_elmt*sizeof(int), &(cnpy->nv3));


    char   ss[20];
    for (int i=0; i<n_elmt; i++) {
      fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, ss, nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
    } 

    for (int i=0; i<n_elmt; i++) {
      cnpy->nv1[i] = nv1[i]; cnpy->nv2[i] = nv2[i]; cnpy->nv3[i] = nv3[i];
    }
    
/*    double maxtmp=-1000.;
    for (int i=0;i<n_v;i++){
      fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);
      cnpy->x_bp[i] = x_bp[i]/cl*L_dim;
      cnpy->y_bp[i] = y_bp[i]/cl*L_dim;
      cnpy->z_bp[i] = z_bp[i]/cl*L_dim;
      maxtmp = PetscMax(maxtmp,cnpy->z_bp[i]);
    }
    double rcvbuff;
    MPI_Allreduce(&maxtmp, &rcvbuff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    cnpy->maxHeight=rcvbuff;*/

    fclose(fd);


    PetscMalloc(n_elmt*sizeof(PetscReal), &(cnpy->cent_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(cnpy->cent_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(cnpy->cent_z));

    PetscMalloc(n_elmt*sizeof(PetscReal), &(cnpy->nf_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(cnpy->nf_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(cnpy->nf_z));

    for (int i=0; i<n_elmt; i++) {
      int n1e = nv1[i]; 
      int n2e = nv2[i]; 
      int n3e = nv3[i];

      cnpy->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      cnpy->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      cnpy->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;

      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      cnpy->nf_x[i] = dy12 * dz13 - dz12 * dy13;
      cnpy->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      cnpy->nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(cnpy->nf_x[i]*cnpy->nf_x[i] + 
                cnpy->nf_y[i]*cnpy->nf_y[i] + 
                cnpy->nf_z[i]*cnpy->nf_z[i]  );

      cnpy->nf_x[i] /=dr; 
      cnpy->nf_y[i] /=dr; 
      cnpy->nf_z[i] /=dr;

    }

    PetscFree(x_bp);
    PetscFree(y_bp);
    PetscFree(z_bp);
    PetscFree(nv1) ;
    PetscFree(nv2) ;
    PetscFree(nv3) ;

  } else{
    PetscPrintf(PETSC_COMM_WORLD,"File %s not found \n",filen);
  }

  PetscPrintf(PETSC_COMM_WORLD,"Canopy file %i has been read\n",icnpy);
  return(0);
}

PetscErrorCode canopy_read_af(CnpyNodes *cnpy, int icnpy, PetscReal cnpyHeightDim, PetscReal cnpyCd)
{
  /* This function reads the canopy af file frontalareadensity.inp */

  FILE *fd;
  char filen[80];
  int   naf;
  PetscReal *aftmp, *yaftmp;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "CnpyHeight and CnpyCd are %le and %le \n:", cnpyHeightDim, cnpyCd);
  PetscPrintf(PETSC_COMM_WORLD,"Reading Frontal Area density \n:",icnpy);

  sprintf(filen,"./frontalareadensity.inp");

  fd = fopen(filen,"r");
    fscanf(fd,"%i",&naf);
    cnpy->naf=naf;
    MPI_Bcast(&naf, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscMalloc(naf*sizeof(PetscReal), &aftmp);
    PetscMalloc(naf*sizeof(PetscReal), &yaftmp);

    PetscMalloc(naf*sizeof(PetscReal), &(cnpy->af));
    PetscMalloc(naf*sizeof(PetscReal), &(cnpy->yaf));


    for (int i = 0; i<naf; i++){
      fscanf(fd, " %le %le", &aftmp[i], &yaftmp[i]);
      cnpy->af[i] = aftmp[i]*cnpyHeightDim;
      cnpy->yaf[i] = yaftmp[i];
    }
  
    PetscFree(aftmp);
    PetscFree(yaftmp);
  

  fclose(fd);
  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Reading Frontal Area density Out\n:",icnpy);

  return 0;
} 

PetscErrorCode search_canopy_nodes(UserCtx *user, CnpyNodes *cnpy, int icnpy, PetscReal cnpyCd)
{
//  PetscReal ts,te,cput;
//  PetscTime(&ts);
  DA    da = user->da, fda = user->fda;

  DALocalInfo info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt      ncx = 40, ncy = 40, ncz = 40;
  LIST          *cell_trg;
  PetscReal     xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscReal     *x_bp = cnpy->x_bp, *y_bp = cnpy->y_bp, *z_bp = cnpy->z_bp;   
  PetscInt      ln_v, n_v = cnpy->n_v;

  PetscInt      i, j, k;

  PetscReal     dcx, dcy, dcz;
  PetscInt      n1e, n2e, n3e;
  PetscReal     xv_min, yv_min, zv_min, xv_max, yv_max, zv_max;
  PetscInt      iv_min, iv_max, jv_min, jv_max, kv_min, kv_max;
  PetscReal     ***nvert;
  PetscInt      ic, jc, kc;
//  Cmpnts        ***csi, ***eta, ***zet;
  PetscReal     ***aj;
  double area, voleu;
  double dhx_, dhy_, dhz_;

//  PetscReal     *aftmp, *yaftmp;
  PetscReal     cnpyHeightmax, cnpyHeightmin, cnpyHeight;
  PetscReal     ***cdy;
  PetscInt   flg_nvert=0;
  

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Identifying Canopy Grid points \n");

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;

  for(i=0; i<n_v; i++) {
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);

    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);

    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }


 PetscPrintf(PETSC_COMM_WORLD,"Canopy Boundaries(Num, x,y,z): %i, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f: \n",icnpy, xbp_min, xbp_max, ybp_min, ybp_max, zbp_min, zbp_max);


  xbp_min -= 0.05; xbp_max += 0.05;
  ybp_min -= 0.05; ybp_max += 0.05;
  zbp_min -= 0.05; zbp_max += 0.05;


  dcx = (xbp_max - xbp_min) / (ncx - 1.);
  dcy = (ybp_max - ybp_min) / (ncy - 1.);
  dcz = (zbp_max - zbp_min) / (ncz - 1.);


  PetscMalloc(ncz * ncy * ncx * sizeof(LIST), &cell_trg);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Initializing List \n");

  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
        initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);
      }
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Bounding Box \n");

  MPI_Barrier(PETSC_COMM_WORLD);
  for (ln_v=0; ln_v < cnpy->n_elmt; ln_v++) {


    n1e = cnpy->nv1[ln_v]; n2e = cnpy->nv2[ln_v]; n3e = cnpy->nv3[ln_v];


    xv_min = PetscMin(PetscMin(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);
    xv_max = PetscMax(PetscMax(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);

    yv_min = PetscMin(PetscMin(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    yv_max = PetscMax(PetscMax(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);

    zv_min = PetscMin(PetscMin(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    zv_max = PetscMax(PetscMax(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);


    iv_min = floor((xv_min - xbp_min) / dcx); 
    iv_max = floor((xv_max - xbp_min) / dcx) +1;

    jv_min = floor((yv_min - ybp_min) / dcy); 
    jv_max = floor((yv_max - ybp_min) / dcy) +1;

    kv_min = floor((zv_min - zbp_min) / dcz); 
    kv_max = floor((zv_max - zbp_min) / dcz) +1;

    iv_min = (iv_min<0) ? 0:iv_min;
    iv_max = (iv_max>ncx) ? ncx:iv_max;

    jv_min = (jv_min<0) ? 0:jv_min;
    jv_max = (jv_max>ncx) ? ncy:jv_max;

    kv_min = (kv_min<0) ? 0:kv_min;
    kv_max = (kv_max>ncz) ? ncz:kv_max;

    for (k=kv_min; k<kv_max; k++) {
      for (j=jv_min; j<jv_max; j++) {
        for (i=iv_min; i<iv_max; i++) {
          insertnode(&(cell_trg[k *ncx*ncy + j*ncx +i]), ln_v);
        }
      }
    }

  }

//PetscPrintf(PETSC_COMM_WORLD,"Box Boundaries(Num, i,j,k): %i, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f: \n",icnpy, iv_min, iv_max, jv_min, jv_max, kv_min, kv_max);



  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"insertnode completed \n");

  Cmpnts ***coor;
  DAVecGetArray(fda, user->lCent, &coor);
  DAVecGetArray(da, user->cnpyNvert, &nvert);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Populating Nvert \n");

  int flg=0;


       int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  for (k=lzs; k<lze; k++)
  for (j=lys; j<lye; j++)
  for (i=lxs; i<lxe; i++) {
      double val=0;

      if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
      coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
      coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {


        ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
        jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
        kc = floor((coor[k][j][i].z - zbp_min )/ dcz);

        val = point_cell_cnpy(coor[k][j][i], ic, jc, kc, cnpy, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg);

        nvert[k][j][i] =  PetscMax(nvert[k][j][i], val);
        if (nvert[k][j][i] < 0) nvert[k][j][i] = 0;

      }
  }


  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Populating Nvert Completed \n");

  int i_min = 1000;
  int j_min = 1000;
  int k_min = 1000;

  int i_max = -1000;
  int j_max = -1000;
  int k_max = -1000;


  for (k=lzs+1; k<lze-1; k++)
  for (j=lys+1; j<lye-1; j++)
  for (i=lxs+1; i<lxe-1; i++) {

     if(nvert[k][j][i]>0.9){ 
//avoid isolated nodes
        if( (nvert[k+1][j][i] >0.9) ||
            (nvert[k-1][j][i] >0.9) ||
            (nvert[k][j+1][i] >0.9) ||
            (nvert[k][j-1][i] >0.9) ||
            (nvert[k][j][i+1] >0.9) ||
            (nvert[k][j][i-1] >0.9) ){

            nvert[k][j][i] = 1.; 

            i_min = PetscMin(i_min, i);
            i_max = PetscMax(i_max, i);
            j_min = PetscMin(j_min, j);
            j_max = PetscMax(j_max, j);
            k_min = PetscMin(k_min, k);
            k_max = PetscMax(k_max, k);

            flg_nvert=1;
         }else{
             nvert[k][j][i] = 0.;
         }

     }else if(nvert[k][j][i] == 1.)
     {

     }
      else {
        nvert[k][j][i] = 0.;
     }
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Nvert Populated \n");
    if (flg_nvert>0) {


      printf ("Canopy %i found by rank %i \n", icnpy, rank);
    }

//  int rank;
//  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  i_min = PetscMax(i_min-2,lxs);
  j_min = PetscMax(j_min-2,lys);
  k_min = PetscMax(k_min-2,lzs);

  i_max = PetscMin(i_max+2,lxe);
  j_max = PetscMin(j_max+2,lye);
  k_max = PetscMin(k_max+2,lze);

  cnpy->i_min = i_min; 
  cnpy->i_max = i_max; 
  cnpy->j_min = j_min; 
  cnpy->j_max = j_max; 
  cnpy->k_min = k_min; 
  cnpy->k_max = k_max; 

  DAVecGetArray(da, user->cnpyCdy, &cdy);

  double  **bufferMax;
  double  **bufferMin;
  double  **bufferTmp1;

  bufferMax =  (double  **) malloc( sizeof(double  *) * mz );
  bufferMin =  (double  **) malloc( sizeof(double  *) * mz );
  bufferTmp1 = (double  **) malloc( sizeof(double  *) * mz );
  for(k=0; k<mz; k++) {
      bufferMax[k] = (double *) malloc( sizeof(double) * mx );
      bufferMin[k] = (double *) malloc( sizeof(double) * mx );
     bufferTmp1[k] = (double *) malloc( sizeof(double) * mx );
    for(i=0; i<mx; i++){
       bufferMax[k][i]=  -10000.;
       bufferMin[k][i]=  10000.;
      bufferTmp1[k][i]= 0.;
    }
  }


  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Finding max height mx, mz, lzs, lze: %i, %i, %i, %i \n",mx, mz, lzs, lze );
  for (k=lzs; k<lze; k++) {
    for (i=lxs; i<lxe; i++) {
      cnpyHeightmax = -10000.;
      cnpyHeightmin =  10000.;

      for (j=lys; j<lye; j++) {
        if(nvert[k][j][i]>0.9){
           cnpyHeightmax = PetscMax(coor[k][j][i].z,cnpyHeightmax); 
           cnpyHeightmin = PetscMin(coor[k][j][i].z,cnpyHeightmin); 
        }
      }
      bufferMax[k][i] = cnpyHeightmax;
      bufferMin[k][i] = cnpyHeightmin;
    }
  }
  

//  MPI_Barrier(PETSC_COMM_WORLD);
//  PetscPrintf(PETSC_COMM_WORLD,"Broadcasting max height %i \n",mx*mz);

   for(k=0; k<mz; k++) {
     for(i=0; i<mx; i++) {
//      MPI_Allreduce(&bufferMax[j][i], &bufferTmp1[j][i], 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
      MPI_Allreduce(&bufferMax[k][i], &bufferTmp1[k][i], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     }
   }

   for(k=0; k<mz; k++){
     for(i=0; i<mx; i++) {
       bufferMax [k][i] = bufferTmp1[k][i];
       bufferTmp1[k][i]=0.;
     }
   }

   for(k=0; k<mz; k++) {
   for(i=0; i<mx; i++) {
//     MPI_Allreduce(&bufferMin[j][i], &bufferTmp1[j][i], 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
     MPI_Allreduce(&bufferMin[k][i], &bufferTmp1[k][i], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   }
   }

   for(k=0; k<mz; k++){
     for(i=0; i<mx; i++) {
       bufferMin[k][i] = bufferTmp1[k][i];
     }
   }


  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Computing Cdy : %le \n",cnpyCd);
  

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

        cnpyHeight =  cnpy->maxHeight - bufferMin[k][i];

        PetscReal adjHeight = (coor[k][j][i].z-bufferMin[k][i])/cnpyHeight;

        if(nvert[k][j][i]>0.9){ 

//          cdy[k][j][i] = adjHeight;
          cdy[k][j][i]=linearinterpaf(adjHeight, cnpy)*cnpyCd;

        } else {
          cdy[k][j][i] = 0.0;
        }

      }
    }
  }


  PetscPrintf(PETSC_COMM_WORLD,"Cdy computed \n");

  DAVecRestoreArray(fda, user->lCent, &coor);
  DAVecRestoreArray(da, user->cnpyNvert, &nvert);
  DAVecRestoreArray(da, user->cnpyCdy, &cdy);

  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
        destroy(&cell_trg[k*ncx*ncy+j*ncx+i]);
      }
    }
  }


  free(bufferMax);
  free(bufferMin);
  free(bufferTmp1);
  PetscFree(cell_trg);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Exiting2  \n");

//  if (flg_nvert==1) {
//  write_output_file2(user, user->cnpyNvert, "CnpyNvert.dat");
//  MPI_Barrier(PETSC_COMM_WORLD);
//  MPI_Abort(PETSC_COMM_WORLD, 1987);
//  }
//
//  flg_nvert=1;

  return 0;
} 

PetscReal linearinterpaf(PetscReal yh, CnpyNodes *cnpy)
{
  PetscReal val=0.;
  int nii;
  int im,ip;

  nii = cnpy->naf;

  for (int i=0; i<nii; i++){ 
    if(cnpy->yaf[i] > yh){
        ip = i;
        im = i-1;
        val = (cnpy->af[ip]-cnpy->af[im])/(cnpy->yaf[ip]-cnpy->yaf[im])*(yh-cnpy->yaf[im]) + cnpy->af[im];
        break;
    }
  }

  return val;
}

PetscErrorCode CanopyForce(UserCtx *user)
{ 
/* Compute the drag force due to the canopy
 * F_i = Cd * Af(y)*u_i*(u_j*u_j)^.5*nvert
 * Cd*Af(y) is precomputed at search_canopy_nodes at initialization
 * nvert is 1 inside the canopy and 0 outside the canopy
*/

  DA              da = user->da, fda = user->fda;
  DALocalInfo info = user->info;
  PetscInt  xs = info.xs, xe = info.xs + info.xm;
  PetscInt    ys = info.ys, ye = info.ys + info.ym;
  PetscInt  zs = info.zs, ze = info.zs + info.zm;
  PetscInt  mx = info.mx, my = info.my, mz = info.mz;
  PetscInt  lxs, lxe, lys, lye, lzs, lze;
  PetscInt        i, j, k, l, icpy;
  
  Cmpnts  ***ucat;
  Cmpnts  ***lFcnpy;
  Cmpnts  ***Fcnpy;
  PetscReal ***cdy, ***nvert;
  
  PetscReal   um ;

  Cmpnts     ***coor, ***csi, ***eta, ***zet;
  PetscReal ***aj;
  double area, voleu;
  double dhx_, dhy_, dhz_;
  
  DAVecGetArray(fda, user->Ucat, &ucat);
  DAVecGetArray(fda, user->cnpyF, &Fcnpy);
  DAVecGetArray(fda, user->lcnpyF, &lFcnpy);
  DAVecGetArray(da, user->cnpyCdy, &cdy);

  DAVecGetArray(da, user->cnpyNvert, &nvert);

  DAVecGetArray(fda, user->lCsi,  &csi);
  DAVecGetArray(fda, user->lEta,  &eta);
  DAVecGetArray(fda, user->lZet,  &zet);
  DAVecGetArray(da,  user->lAj,  &aj);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Computing Canopy Force in Flow Solver \n");

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
        Fcnpy[k][j][i].x = 0.0;
        Fcnpy[k][j][i].y = 0.0;
        Fcnpy[k][j][i].z = 0.0;

        lFcnpy[k][j][i].x = 0.0;
        lFcnpy[k][j][i].y = 0.0;
        lFcnpy[k][j][i].z = 0.0;
      }
    }
  }

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  int g_imin = 100000;
  int g_jmin = 100000;
  int g_kmin = 100000;

  int g_imax = -100000;
  int g_jmax = -100000;
  int g_kmax = -100000;

  for (icpy=0;icpy<NumberOfCnpy;icpy++){
    int kmin = Canopy[icpy].k_min;
    int kmax = Canopy[icpy].k_max;
    int jmin = Canopy[icpy].j_min;
    int jmax = Canopy[icpy].j_max;
    int imin = Canopy[icpy].i_min;
    int imax = Canopy[icpy].i_max;



    g_kmin=PetscMin(g_kmin,kmin);
    g_jmin=PetscMin(g_jmin,jmin);
    g_imin=PetscMin(g_imin,imin);

    g_kmax=PetscMax(g_kmax,kmax);
    g_jmax=PetscMax(g_jmax,jmax);
    g_imax=PetscMax(g_imax,imax);


    for (k=kmin; k<kmax; k++) {
      for (j=jmin; j<jmax; j++) {
        for (i=imin; i<imax; i++) {

          um = sqrt(ucat[k][j][i].x*ucat[k][j][i].x 
                  + ucat[k][j][i].y*ucat[k][j][i].y 
                  + ucat[k][j][i].z*ucat[k][j][i].z);

          lFcnpy[k][j][i].x -= cdy[k][j][i]*ucat[k][j][i].x*nvert[k][j][i]*um;
          lFcnpy[k][j][i].y -= cdy[k][j][i]*ucat[k][j][i].y*nvert[k][j][i]*um;
          lFcnpy[k][j][i].z -= cdy[k][j][i]*ucat[k][j][i].z*nvert[k][j][i]*um;


        }
      }
    }


  }

  DALocalToLocalBegin(fda, user->lcnpyF, INSERT_VALUES, user->lcnpyF);
  DALocalToLocalEnd  (fda, user->lcnpyF, INSERT_VALUES, user->lcnpyF);



  PetscInt kp;
  PetscInt jp;
  PetscInt ip;

  for (k=g_kmin; k<g_kmax; k++) {

        kp=k+1;
        if(kk_periodic && k==mz-2){
          kp = mz+1;
        }

    for (j=g_jmin; j<g_jmax; j++) {
        jp = j+1;
      for (i=g_imin; i<g_imax; i++) {
        ip = i+1;

        PetscReal Fcart_x = 0.5*(lFcnpy[k][j][i].x+lFcnpy[k][j][ip].x);
        PetscReal Fcart_z = 0.5*(lFcnpy[k][j][i].z+lFcnpy[k][jp][i].z);
        PetscReal Fcart_y = 0.5*(lFcnpy[k][j][i].y+lFcnpy[kp][j][i].y);
                                                                        
   
        Fcnpy[k][j][i].x = Fcart_x*csi[k][j][i].x
                         + Fcart_y*csi[k][j][i].y
                         + Fcart_z*csi[k][j][i].z;

        Fcnpy[k][j][i].y = Fcart_x*eta[k][j][i].x
                         + Fcart_y*eta[k][j][i].y
                         + Fcart_z*eta[k][j][i].z;

        Fcnpy[k][j][i].z = Fcart_x*zet[k][j][i].x
                         + Fcart_y*zet[k][j][i].y
                         + Fcart_z*zet[k][j][i].z;

      }
    }
  }



  DAVecRestoreArray(fda, user->Ucat, &ucat);
  DAVecRestoreArray(fda, user->cnpyF, &Fcnpy);
  DAVecRestoreArray(fda, user->lcnpyF, &lFcnpy);
  DAVecRestoreArray(da, user->cnpyCdy, &cdy);
  DAVecRestoreArray(da, user->cnpyNvert, &nvert);

  DAVecRestoreArray(fda, user->lCsi,  &csi);
  DAVecRestoreArray(fda, user->lEta,  &eta);
  DAVecRestoreArray(fda, user->lZet,  &zet);
  DAVecRestoreArray(da,  user->lAj,  &aj);



//    write_output_file2(user, user->cnpyCdy, "CnpyCdy.dat");
    write_output_file2(user, user->cnpyF, "CnpyForce.dat");
//    write_output_file2(user, user->cnpyNvert, "CnpyNvert.dat");
//
//  MPI_Barrier(PETSC_COMM_WORLD);
//  MPI_Abort(PETSC_COMM_WORLD, 1987);

  return 0;
}

PetscErrorCode randomdirection2(Cmpnts p, PetscInt ip, PetscInt jp,
             PetscReal xbp_min, PetscReal ybp_min,
             PetscReal zbp_max, PetscReal dcx, PetscReal dcy,
             PetscReal dir[3],PetscInt seed)
{
  Cmpnts endpoint;
  PetscReal s;

  PetscReal xpc, ypc;

  xpc = dcx * (ip+0.5) + xbp_min;
  ypc = dcy * (jp+0.5) + ybp_min;

  srand(seed);

  s = rand() / ((double)RAND_MAX + 1) - 0.5;
  endpoint.x = xpc + s * dcx;
  endpoint.y = ypc + s * dcy;
  endpoint.z = zbp_max + 0.2;

  dir[0] = endpoint.x - p.x;
  dir[1] = endpoint.y - p.y;
  dir[2] = endpoint.z - p.z;

  s = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  dir[0] /= s;
  dir[1] /= s;
  dir[2] /= s;
  return 0;
}

PetscInt point_cell_cnpy(Cmpnts p, PetscInt ip, PetscInt jp, PetscInt kp,
			     CnpyNodes *cnpy, PetscInt ncx, PetscInt ncy,
			     PetscInt ncz, PetscReal dcx, PetscReal dcy,
			     PetscReal xbp_min, PetscReal ybp_min,
			     PetscReal zbp_max, LIST *cell_trg,
			     PetscInt flg)
{
  /* This subroutine is the same as that in ibm.c */

  int       *nv1  = cnpy->nv1,  *nv2  = cnpy->nv2,  *nv3  = cnpy->nv3;
  PetscReal	*nf_x = cnpy->nf_x, *nf_y = cnpy->nf_y, *nf_z = cnpy->nf_z;
  PetscReal	*x_bp = cnpy->x_bp, *y_bp = cnpy->y_bp, *z_bp = cnpy->z_bp;

  PetscInt	i, j, k, ln_v, n1e, n2e, n3e, nintp;
  
  PetscInt	nvert_l;
  PetscReal	dt[1000], ndotn, dirdotn;
  Cmpnts        dnn[1000], nn;

  PetscReal epsilon = 1.e-8;
  PetscReal eps_tangent=1.e-10;

  PetscTruth	*Element_Searched;
  j = jp; i = ip;

  PetscTruth NotDecided = PETSC_TRUE, Singularity = PETSC_FALSE;
  PetscReal t, u, v;
  PetscReal orig[3], dir[3], vert0[3], vert1[3], vert2[3];

  node *current;
  PetscInt searchtimes=0;
  PetscMalloc(cnpy->n_elmt*sizeof(PetscTruth), &Element_Searched);

  if (flg) PetscPrintf(PETSC_COMM_SELF, " serch itr\n");

  while (NotDecided) {

    searchtimes++;
    nintp = 0 ;
    randomdirection2(p, ip, jp, xbp_min, ybp_min, zbp_max, dcx, dcy, dir, searchtimes);
    Singularity = PETSC_FALSE;

    if (flg) 
      PetscPrintf(PETSC_COMM_SELF, " serch itr, dir %d %le %le %le\n", searchtimes,dir[0],dir[1],dir[2]);

    for (ln_v=0; ln_v<cnpy->n_elmt; ln_v++) {
      Element_Searched[ln_v] = PETSC_FALSE;
    }

    for (k=kp; k<ncz; k++) {
      current = cell_trg[k*ncx*ncy+j*ncx+i].head;
      while (current) {
	      ln_v = current->Node;
	      if (!Element_Searched[ln_v]) {
	        Element_Searched[ln_v] = PETSC_TRUE;
	        n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];
	        nn.x=nf_x[ln_v]; nn.y=nf_y[ln_v]; nn.z=nf_z[ln_v];

	        orig[0] = p.x; orig[1] = p.y, orig[2] = p.z;

	        vert0[0] = x_bp[n1e]; vert0[1] = y_bp[n1e]; vert0[2] = z_bp[n1e];
	        vert1[0] = x_bp[n2e]; vert1[1] = y_bp[n2e]; vert1[2] = z_bp[n2e];
	        vert2[0] = x_bp[n3e]; vert2[1] = y_bp[n3e]; vert2[2] = z_bp[n3e];
                  
	        dirdotn=dir[0]*nn.x+dir[1]*nn.y+dir[2]*nn.z;

	        nvert_l = intsect_triangle(orig, dir, vert0, vert1, vert2, &t, &u, &v);
	  
	        if (nvert_l > 0 && t>0) {
	          dt[nintp] = t;
	          dnn[nintp].x=nn.x;dnn[nintp].y=nn.y;dnn[nintp].z=nn.z;

	          nintp ++;
	          PetscInt temp;
	          for (temp = 0; temp < nintp-1; temp++) {
	            // Two interception points are the same, this leads to huge
	            // trouble for crossing number test
	            // Rather to program for all cases, we use a new line to
	            // repeat the test
	            ndotn=dnn[temp].x*nn.x+dnn[temp].y*nn.y+dnn[temp].z*nn.z;	      
	            
	            if ((fabs(t-dt[temp]) < epsilon && ndotn>-0.95/*-0.97*/)){ 
		            Singularity = PETSC_TRUE;
	            }

	          }

	          if (Singularity) break;
	        }

	      }

	      if (Singularity) {
	        break;
	      }
	      else {
	        current = current->next;
	      }
	    } 
// Search through the list
      if (Singularity) {
        break;
      }

    } // for k

      if (!Singularity) {
        NotDecided = PETSC_TRUE;
        if (nintp%2) { // The interception point number is odd, inside body
        	PetscFree(Element_Searched);
        	return 1;
        }
        else {
        	PetscFree(Element_Searched);
        	return 0;
        }
      }
    }
  PetscFree(Element_Searched);
  return 0;

}

PetscErrorCode CanopyInitialization(UserCtx *user, int NumberOfCnpy, PetscReal cnpyHeightDim, PetscReal cnpyCd)
{

  PetscPrintf(PETSC_COMM_WORLD,"Initializing Canopy Module \n");

  PetscMalloc(NumberOfCnpy*sizeof(CnpyNodes), &Canopy);
  
  for(int icnpy=0; icnpy<NumberOfCnpy; icnpy++){
    canopy_read_ucd(&Canopy[icnpy], icnpy);
    canopy_read_af(&Canopy[icnpy],  icnpy, cnpyHeightDim, cnpyCd);
    search_canopy_nodes(user, &Canopy[icnpy], icnpy, cnpyCd);
  }
    write_output_file2(user, user->cnpyCdy, "CnpyCdy.dat");
    //write_output_file2(user, user->cnpyF, "CnpyForce.dat");
    write_output_file2(user, user->cnpyNvert, "CnpyNvert.dat");
//
  return 0;
}


