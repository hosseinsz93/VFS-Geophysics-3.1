#include "petsc.h"
#include "variables.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <valarray>
//#include <algorithm>

extern PetscInt  NumberOfBodies, rstart, input_ib_depth, periodic_morpho;
extern PetscInt tistart, ti, tiout, bed_roughness, STRONG_COUPLING, Barchans, sed_tio, SimpleCellCheck, osl_inlet_sediment_flux, y_direction;
extern PetscInt projection_method;
extern PetscInt effective_bed_shear_stress, aval_loop, rans, sand_slide;
extern PetscReal Nseg, U_Bulk, sediment_thickness,d50,deltab, FlowDepth, porosity, sbbb, Cs_, Angle_repose;
extern double inlet_sediment_flux;
#define  PI 3.14159265

double E_coeff (double utau, double ks, double nu)
{
  double kplus=utau*ks/nu, dB;
  double kappa=0.41, B=5.5;
  if(kplus<=2.25) dB = 0.0;
  else if(kplus>2.25 && kplus<90.) dB = (B-8.5+1./kappa*log(kplus))*sin(0.4258*(log(kplus)-0.811));
  else if(kplus>=90.) dB = B-8.5+1./kappa*log(kplus);
  return exp(kappa*(B-dB));
}

double u_hydset_roughness(double nu, double y, double utau, double ks)
{
 double y0plus=11.53,f;
 double kappa=0.41, B=5.5;
 double yplus = utau*y/nu;
 
  if(yplus<=y0plus){f= utau * yplus;}
  else {f= utau/kappa*log(E_coeff(utau,ks,nu)*yplus);}
  return f;
}

double f_hydset(double nu, double u, double y, double utau0, double ks) 
{
double y0plus=11.53, f;
double kappa=0.41, B=5.5;
double yplus=utau0*y/nu;
if (yplus<=y0plus) {f= utau0*yplus-u;}
else {f= utau0*(1./kappa*log(E_coeff(utau0,ks,nu)*yplus))-u;}
return f;
}

double df_hydset (double nu, double u, double y, double utau0, double ks)
{
double eps=1.e-7;
return (f_hydset(nu, u, y, utau0 + eps, ks) - f_hydset(nu, u, y, utau0 - eps, ks))/(2.*eps);
}


double find_utau_hydset(double nu,double u, double y, double utau_guess, double ks)
{
 double utau,utau0 = utau_guess;  
 int ii;
 for(ii=0;ii<30;ii++){
 utau=utau0-f_hydset(nu,u,y,utau0,ks)/df_hydset(nu,u,y,utau0,ks);
 if (fabs(utau-utau0)<1.e-7)break;
 utau0=utau;
  }
return utau;
};


void wall_function_roughness_a (UserCtx *user, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	//*ustar = find_utau_Cabot_roughness(1./user->ren, ut_mag, sc, 0.01, 0, ks);
        *ustar = find_utau_hydset(1./user->ren, ut_mag, sc, 0.01, ks);
	//double ut_mag_modeled = u_Cabot_roughness(1./user->ren, sb, *ustar, 0, ks);
	double ut_mag_modeled = u_hydset_roughness(1./user->ren, sb, *ustar, ks);

	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

//Hossein added from NSF petsc3.1
void wall_function_roughness_a_levelset (UserCtx *user, double ks, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz, double fluid_density, double fluid_viscosity)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	//*ustar = find_utau_Cabot_roughness(1./user->ren, ut_mag, sc, 0.01, 0, ks);
        *ustar = find_utau_hydset(fluid_viscosity, ut_mag, sc, 0.01, ks);
	//double ut_mag_modeled = u_Cabot_roughness(1./user->ren, sb, *ustar, 0, ks);
	double ut_mag_modeled = u_hydset_roughness(fluid_viscosity, sb, *ustar, ks);

	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}


PetscErrorCode distance_sed(PetscInt nv, Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
{
  PetscReal xn1, yn1, zn1;
  PetscReal xc, yc, zc;
  PetscInt rank;
  PetscReal dx1, dy1, dz1, dx2, dy2, dz2, r;

  dx1 = p3.x - p1.x;
  dy1 = p3.y - p1.y;
  dz1 = p3.z - p1.z;

  dx2 = p4.x - p2.x;
  dy2 = p4.y - p2.y;
  dz2 = p4.z - p2.z;

  xn1 = dy1 * dz2 - dz1 * dy2;
  yn1 = - (dx1 * dz2 - dz1 * dx2);
  zn1 = dx1 * dy2 - dy1 * dx2;

  r = sqrt(xn1 * xn1 + yn1 * yn1 + zn1 * zn1);
  xn1 /= r; yn1 /= r; zn1 /= r;

  xc = 0.25 * (p1.x + p2.x + p3.x + p4.x);
  yc = 0.25 * (p1.y + p2.y + p3.y + p4.y);
  zc = 0.25 * (p1.z + p2.z + p3.z + p4.z);

  *d = (p.x - xc) * xn1 + (p.y - yc) * yn1 + (p.z - zc) * zn1;
                        // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
                        // if(nv==771 && rank==5) PetscPrintf(PETSC_COMM_SELF, "d %e\n",*d);
  if (fabs(*d)<1.e-6) *d=0.;
  return (0);
}

PetscTruth ISInsideCell_sed(Cmpnts p, Cmpnts cell[8], PetscReal d[6],PetscInt nv)
{
  // k direction
  distance_sed(nv, cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  if (d[4]<0) return(PETSC_FALSE);
  distance_sed(nv, cell[4], cell[7], cell[6], cell[5], p, &(d[5]));
  if (d[5]<0) return(PETSC_FALSE);

  // j direction
  distance_sed(nv ,cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  if (d[2]<0) return(PETSC_FALSE);

  distance_sed(nv, cell[3], cell[2], cell[6], cell[7], p, &(d[3]));
  if (d[3]<0) return(PETSC_FALSE);

  // i direction
  distance_sed(nv, cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  if (d[0]<0) return(PETSC_FALSE);
  
  distance_sed(nv, cell[1], cell[5], cell[6], cell[2], p, &(d[1]));
  if (d[1]<0) return(PETSC_FALSE);
  return(PETSC_TRUE);
}
  

PetscErrorCode ibm_intp_pj_centroid(UserCtx *user, IBMNodes *ibm, PetscInt ti, PetscInt tistart)
{

DA	        da = user->da, fda = user->fda;
DALocalInfo	info = user->info;

PetscInt        i, j, k, nothing;
PetscInt	ncx = 10, ncy = 10, ncz = 10, bag=100000;
PetscInt	cell_trg_idx[ncz*2][ncy*2][ncx*2];
PetscInt        cell_trgi[ncz*2][ncy*2][ncx*2][bag], cell_trgj[ncz*2][ncy*2][ncx*2][bag], cell_trgk[ncz*2][ncy*2][ncx*2][bag];

PetscInt	n_elmt = ibm->n_elmt, n_v = ibm->n_v;
PetscInt	*nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3= ibm->nv3;
PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;

PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;
PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;
PetscReal	*vx = ibm->vx, *vy = ibm->vy, *vz = ibm->vz;

//PetscReal       *shear,*shear_x,*shear_y,*shear_z;

Vec             Coor;
Cmpnts          ***coor, ***ucat;
PetscReal       ***nvert;
PetscReal       dcx, dcy, dcz;

PetscInt        n1e, n2e, n3e;
PetscReal       xc, yc, zc , dx12,dy12,dz12,dx13,dy13,dz13,dr;
Cmpnts          pj;
PetscInt        ivs, jvs, kvs, ive, jve, kve, nv;

FILE            *f;
char            filen[80];
PetscTruth      flag=PETSC_FALSE;

PetscInt        around_cubes, cp, rank, processors;
PetscReal       flag_number_sum, flag_number;

// averaging methods:
PetscInt        TetQuad = 1,TetQuad_search=0, Harmonic = 0, SAver = 0;
PetscReal       Initial_scc=0.05;

//if(ti==tistart || ti==0){ PetscPrintf(PETSC_COMM_WORLD, "TetQuad, TetQuad_search, Harmonic, SAver: %d %d %d %d \n",TetQuad, TetQuad_search, Harmonic, SAver);
//                          PetscPrintf(PETSC_COMM_WORLD, "Initial scc: %e \n",Initial_scc);}


MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
MPI_Comm_size(PETSC_COMM_WORLD, &processors);

PetscInt	xs = info.xs, xe = info.xs + info.xm;
PetscInt  	ys = info.ys, ye = info.ys + info.ym;
PetscInt	zs = info.zs, ze = info.zs + info.zm;
PetscInt	lxs, lxe, lys, lye, lzs, lze;

 lxs = xs; lxe = xe;
 lys = ys; lye = ye;
 lzs = zs; lze = ze;

 if (xs==0) lxs = xs+1;
 if (ys==0) lys = ys+1;
 if (zs==0) lzs = zs+1;

 if (xe==info.mx) lxe = xe-1;
 if (ye==info.my) lye = ye-1;
 if (ze==info.mz) lze = ze-1;

//Get the coordinates of the block
DAGetGhostedCoordinates(da, &Coor);
DAVecGetArray(fda, Coor, &coor);
DAVecGetArray(fda, user->lUcat, &ucat);
//DAVecGetArray(da, user->P, &prs);
DAVecGetArray(da, user->lNvert, &nvert);

/*
//Allocate the memory for shear component
 PetscMalloc(n_elmt*sizeof(PetscReal),&shear);
 PetscMalloc(n_elmt*sizeof(PetscReal),&shear_x);
 PetscMalloc(n_elmt*sizeof(PetscReal),&shear_y);
 PetscMalloc(n_elmt*sizeof(PetscReal),&shear_z);
 PetscMalloc(n_elmt*sizeof(PetscReal),&vx);
 PetscMalloc(n_elmt*sizeof(PetscReal),&vy);
 PetscMalloc(n_elmt*sizeof(PetscReal),&vz);
 PetscMalloc(n_elmt*sizeof(PetscReal),&pressure);
*/

// Import the k values
/*
PetscOptionsGetInt(PETSC_NULL, "-k_min", &k_Export_Min, &flag);
PetscOptionsGetInt(PETSC_NULL, "-k_max", &k_Export_Max, &flag);
PetscOptionsGetReal(PETSC_NULL, "-dpj", &dpj, &flag);
*/

 //Start finding the thing
  xbp_min = 1.e20; xbp_max = -1.e20;
  ybp_min = 1.e20; ybp_max = -1.e20;
  zbp_min = 1.e20; zbp_max = -1.e20;
  
 // for (k=k_Export_Min; k<k_Export_Max; k++) {
    for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
    for (i=lxs; i<lxe; i++) {
	xbp_min = PetscMin(xbp_min, coor[k][j][i].x);
	xbp_max = PetscMax(xbp_max, coor[k][j][i].x);

	ybp_min = PetscMin(ybp_min, coor[k][j][i].y);
	ybp_max = PetscMax(ybp_max, coor[k][j][i].y);

	zbp_min = PetscMin(zbp_min, coor[k][j][i].z);
	zbp_max = PetscMax(zbp_max, coor[k][j][i].z);
      }
    }
  }

  //Adjust 
  xbp_min -= 0.001; xbp_max += 0.001;
  ybp_min -= 0.001; ybp_max += 0.001;
  zbp_min -= 0.001; zbp_max += 0.001;

  //Devide the block in to nc small boxes
    dcx = (xbp_max - xbp_min) / (ncx - 1.);
    dcy = (ybp_max - ybp_min) / (ncy - 1.);
    dcz = (zbp_max - zbp_min) / (ncz - 1.);

  //Assign cell indices = -1
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {

	cell_trg_idx[k][j][i] = -1;
		  
      }
    }
  }

  //Adjust 
  //Identify the box - control volume  that point is now in
  int iv, jv, kv;
    for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
    for (i=lxs; i<lxe; i++) {
	iv = floor((coor[k][j][i].x - xbp_min) / dcx);
	jv = floor((coor[k][j][i].y - ybp_min) / dcy);
	kv = floor((coor[k][j][i].z - zbp_min) / dcz);
	//We can select only the IB and fluid nodes
	 if (nvert[k][j][i] < 2.)
	  {
	    cell_trg_idx[kv][jv][iv] ++;

	    cell_trgi[kv][jv][iv][cell_trg_idx[kv][jv][iv]] = i;
	    cell_trgj[kv][jv][iv][cell_trg_idx[kv][jv][iv]] = j;
	    cell_trgk[kv][jv][iv][cell_trg_idx[kv][jv][iv]] = k;
	  }
       }
     }
   }


for (nv=0; nv<n_elmt; nv++) 
{
ibm->scc[nv]=Initial_scc;
}

// going through the bed cells
for (nv=0; nv<n_elmt; nv++) 
 {
 flag = PETSC_FALSE;
 if(ibm->Rigidity[nv] == 1)
 {
    vx[nv] = 0.; 
    vy[nv] = 0.; 
    vz[nv] = 0.;
    flag=PETSC_TRUE;
 }
 else
 {   
    vx[nv] = 0.; 
    vy[nv] = 0.; 
    vz[nv] = 0.;

    //3 vertice of the triangle
    n1e = nv1[nv]; 
    n2e = nv2[nv]; 
    n3e = nv3[nv];
   

    dx12 = x_bp[ n2e ] - x_bp[ n1e ];
    dy12 = y_bp[ n2e ] - y_bp[ n1e ];
    dz12 = z_bp[ n2e ] - z_bp[ n1e ];
    dx13 = x_bp[ n3e ] - x_bp[ n1e ];
    dy13 = y_bp[ n3e ] - y_bp[ n1e ];
    dz13 = z_bp[ n3e ] - z_bp[ n1e ];
    nf_x[nv] = dy12 * dz13 - dz12 * dy13;
    nf_y[nv] = -dx12 * dz13 + dz12 * dx13;
    nf_z[nv] = dx12 * dy13 - dy12 * dx13;
    dr = sqrt(nf_x[nv] * nf_x[nv] + nf_y[nv] * nf_y[nv]	+ nf_z[nv] * nf_z[nv]);

    nf_x[ nv ] /= dr;
    nf_y[ nv ] /= dr;
    nf_z[ nv ] /= dr;
 
    host:around_cubes = 3;
    // C is the center point of the triangle
          xc = (x_bp[n1e] + x_bp[n2e] + x_bp[n3e]) / 3.;
          yc = (y_bp[n1e] + y_bp[n2e] + y_bp[n3e]) / 3.;
          zc = (z_bp[n1e] + z_bp[n2e] + z_bp[n3e]) / 3.;

    // recontruct a line from C normal to the cell
    // For the aneurysm case we use - normal surface because the normal vector originally outward
   // We want it inward
         //pj.x = xc + 0.125 * nf_x[nv];
         //pj.y = yc + 0.125 * nf_y[nv];
         //pj.z = zc + 0.125 * nf_z[nv];

         pj.x = xc + ibm->scc[nv] * nf_x[nv];
         pj.y = yc + ibm->scc[nv] * nf_y[nv];
         pj.z = zc + ibm->scc[nv] * nf_z[nv];
    // pj is one point inside the fluid domain

    // locate point P into control volume number v
    iv = floor((pj.x - xbp_min) / dcx);
    jv = floor((pj.y - ybp_min) / dcy);
    kv = floor((pj.z - zbp_min) / dcz);
    
    // Total 27 cubic around the current iv,jv, kv
    ivs = PetscMax(iv-3, 0); 
    jvs = PetscMax(jv-3, 0); 
    kvs = PetscMax(kv-3, 0);
    
    ive = PetscMin(iv+4, ncx); 
    jve = PetscMin(jv+4, ncy); 
    kve = PetscMin(kv+4, ncz);

    flag_number=0.;

    // Now interpolate the velocity from box s- e

    for (k=kvs; k<kve; k++) {
    for (j=jvs; j<jve; j++) {
    for (i= ivs; i<ive; i++) {

	  PetscInt        ivh, jvh, kvh, ic;
	  PetscReal       d[6];
	  PetscReal       r[8];
	  PetscReal       pnv[8];
	  Cmpnts          p[8];

	  //If the point (i,j,k) has one cell
	  if (cell_trg_idx[k][j][i]>=0) 
          {

	    for (ic = 0; ic < cell_trg_idx[k][j][i]+1; ic++) 
             {

	      //Go to that cell with number ic
	      ivh = cell_trgi[k][j][i][ic];
	      jvh = cell_trgj[k][j][i][ic];
	      kvh = cell_trgk[k][j][i][ic];
	    
	      if (ivh>=lxs && ivh<lxe && jvh>=lys && jvh<lye && kvh>=lzs && kvh<lze) {

		//Set 8 points around this point
		// p[i] is one point in the bucket
		p[0] = coor[kvh][jvh][ivh];
		p[1] = coor[kvh][jvh][ivh+1];
		p[2] = coor[kvh][jvh+1][ivh+1];
		p[3] = coor[kvh][jvh+1][ivh];

		p[4] = coor[kvh+1][jvh][ivh];
		p[5] = coor[kvh+1][jvh][ivh+1];
		p[6] = coor[kvh+1][jvh+1][ivh+1];
		p[7] = coor[kvh+1][jvh+1][ivh];

		//if point pj - projection 

		if (ISInsideCell_sed(pj, p, d, nv)) {
		  // do the interpolation and then break

		  PetscReal x, y, z;
		  x = d[0] / (d[0]+d[1]);
		  y = d[2] / (d[2]+d[3]);
		  z = d[4] / (d[4]+d[5]);

                 // check if any of 8 points is solid node?
                 if (nvert[kvh  ][jvh  ][ivh  ] >=3.) {pnv[0]=0.;}else{pnv[0]=1.;}
                 if (nvert[kvh  ][jvh  ][ivh+1] >=3.) {pnv[1]=0.;}else{pnv[1]=1.;}
                 if (nvert[kvh  ][jvh+1][ivh+1] >=3.) {pnv[2]=0.;}else{pnv[2]=1.;}
                 if (nvert[kvh  ][jvh+1][ivh  ] >=3.) {pnv[3]=0.;}else{pnv[3]=1.;}
                   
                 if (nvert[kvh+1][jvh  ][ivh  ] >=3.) {pnv[4]=0.;}else{pnv[4]=1.;}
                 if (nvert[kvh+1][jvh  ][ivh+1] >=3.) {pnv[5]=0.;}else{pnv[5]=1.;}
                 if (nvert[kvh+1][jvh+1][ivh+1] >=3.) {pnv[6]=0.;}else{pnv[6]=1.;}
                 if (nvert[kvh+1][jvh+1][ivh  ] >=3.) {pnv[7]=0.;}else{pnv[7]=1.;}
                  
                 if(TetQuad_search)
                 {
                  if(pnv[0] ==0. || pnv[1] ==0. || pnv[2] ==0. || pnv[3]==0. || pnv[4]==0. || pnv[5]==0. || pnv[6]==0. || pnv[7] ==0.)
                    {
                     ibm->scc[nv] +=0.001;
                     if(ibm->scc[nv]>0.1) goto next;
                     goto host;
                    }
                 }
                      
                  
                 if(Harmonic)
                 {
                  // distance of point from each of 8 corner points
                 r[0]=sqrt(pow((pj.x-coor[kvh  ][jvh  ][ivh  ].x),2.)+pow((pj.y-coor[kvh  ][jvh  ][ivh  ].y),2.)+pow((pj.z-coor[kvh  ][jvh  ][ivh  ].z),2.));
                 r[1]=sqrt(pow((pj.x-coor[kvh  ][jvh  ][ivh+1].x),2.)+pow((pj.y-coor[kvh  ][jvh  ][ivh+1].y),2.)+pow((pj.z-coor[kvh  ][jvh  ][ivh+1].z),2.));
                 r[2]=sqrt(pow((pj.x-coor[kvh  ][jvh+1][ivh+1].x),2.)+pow((pj.y-coor[kvh  ][jvh+1][ivh+1].y),2.)+pow((pj.z-coor[kvh  ][jvh+1][ivh+1].z),2.));
                 r[3]=sqrt(pow((pj.x-coor[kvh  ][jvh+1][ivh  ].x),2.)+pow((pj.y-coor[kvh  ][jvh+1][ivh  ].y),2.)+pow((pj.z-coor[kvh  ][jvh+1][ivh  ].z),2.));
		 r[4]=sqrt(pow((pj.x-coor[kvh+1][jvh  ][ivh  ].x),2.)+pow((pj.y-coor[kvh+1][jvh  ][ivh  ].y),2.)+pow((pj.z-coor[kvh+1][jvh  ][ivh  ].z),2.));
                 r[5]=sqrt(pow((pj.x-coor[kvh+1][jvh  ][ivh+1].x),2.)+pow((pj.y-coor[kvh+1][jvh  ][ivh+1].y),2.)+pow((pj.z-coor[kvh+1][jvh  ][ivh+1].z),2.));
                 r[6]=sqrt(pow((pj.x-coor[kvh+1][jvh+1][ivh+1].x),2.)+pow((pj.y-coor[kvh+1][jvh+1][ivh+1].y),2.)+pow((pj.z-coor[kvh+1][jvh+1][ivh+1].z),2.));
                 r[7]=sqrt(pow((pj.x-coor[kvh+1][jvh+1][ivh  ].x),2.)+pow((pj.y-coor[kvh+1][jvh+1][ivh  ].y),2.)+pow((pj.z-coor[kvh+1][jvh+1][ivh  ].z),2.));
                 }
                 
                  if(TetQuad)
                  {
                  // Tetra - quaradtic interpolation
                  vx[nv] = (ucat[kvh  ][jvh  ][ivh  ].x * (1-x) * (1-y) * (1-z) +
			    ucat[kvh  ][jvh  ][ivh+1].x * x * (1-y) * (1-z) +
			    ucat[kvh  ][jvh+1][ivh  ].x * (1-x) * y * (1-z) +
			    ucat[kvh+1][jvh  ][ivh  ].x * (1-x) * (1-y) * z +
			    ucat[kvh+1][jvh  ][ivh+1].x * x * (1-y) * z +
			    ucat[kvh+1][jvh+1][ivh  ].x * (1-x) * y * z +
			    ucat[kvh  ][jvh+1][ivh+1].x * x * y * (1-z) +
			    ucat[kvh+1][jvh+1][ivh+1].x * x * y * z);

		  vy[nv] = (ucat[kvh  ][jvh  ][ivh  ].y * (1-x) * (1-y) * (1-z) +
			    ucat[kvh  ][jvh  ][ivh+1].y * x * (1-y) * (1-z) +
			    ucat[kvh  ][jvh+1][ivh  ].y * (1-x) * y * (1-z) +
			    ucat[kvh+1][jvh  ][ivh  ].y * (1-x) * (1-y) * z +
			    ucat[kvh+1][jvh  ][ivh+1].y * x * (1-y) * z +
			    ucat[kvh+1][jvh+1][ivh  ].y * (1-x) * y * z +
			    ucat[kvh  ][jvh+1][ivh+1].y * x * y * (1-z) +
			    ucat[kvh+1][jvh+1][ivh+1].y * x * y * z);

		  vz[nv] = (ucat[kvh  ][jvh  ][ivh  ].z * (1-x) * (1-y) * (1-z) +
			    ucat[kvh  ][jvh  ][ivh+1].z * x * (1-y) * (1-z) +
			    ucat[kvh  ][jvh+1][ivh  ].z * (1-x) * y * (1-z) +
			    ucat[kvh+1][jvh  ][ivh  ].z * (1-x) * (1-y) * z +
			    ucat[kvh+1][jvh  ][ivh+1].z * x * (1-y) * z +
			    ucat[kvh+1][jvh+1][ivh  ].z * (1-x) * y * z +
			    ucat[kvh  ][jvh+1][ivh+1].z * x * y * (1-z) +
			    ucat[kvh+1][jvh+1][ivh+1].z * x * y * z);
		  }
                  if(SAver)
                  {
                  // simple averaging among none-solid nodes
                  vx[nv] = (ucat[kvh  ][jvh  ][ivh  ].x*pnv[0]+
                            ucat[kvh  ][jvh  ][ivh+1].x*pnv[1]+
                            ucat[kvh  ][jvh+1][ivh+1].x*pnv[2]+
			    ucat[kvh  ][jvh+1][ivh  ].x*pnv[3]+
                            ucat[kvh+1][jvh  ][ivh  ].x*pnv[4]+
			    ucat[kvh+1][jvh  ][ivh+1].x*pnv[5]+
                            ucat[kvh+1][jvh+1][ivh+1].x*pnv[6]+
			    ucat[kvh+1][jvh+1][ivh  ].x*pnv[7])/
                            (1.e-10+pnv[0]+pnv[1]+pnv[2]+pnv[3]+pnv[4]+pnv[5]+pnv[6]+pnv[7]);

                  vy[nv] = (ucat[kvh  ][jvh  ][ivh  ].y*pnv[0]+
                            ucat[kvh  ][jvh  ][ivh+1].y*pnv[1]+
                            ucat[kvh  ][jvh+1][ivh+1].y*pnv[2]+
			    ucat[kvh  ][jvh+1][ivh  ].y*pnv[3]+
                            ucat[kvh+1][jvh  ][ivh  ].y*pnv[4]+
			    ucat[kvh+1][jvh  ][ivh+1].y*pnv[5]+
                            ucat[kvh+1][jvh+1][ivh+1].y*pnv[6]+
			    ucat[kvh+1][jvh+1][ivh  ].y*pnv[7])/ 
                            (1.e-10+pnv[0]+pnv[1]+pnv[2]+pnv[3]+pnv[4]+pnv[5]+pnv[6]+pnv[7]);

                  vz[nv] = (ucat[kvh  ][jvh  ][ivh  ].z*pnv[0]+
                            ucat[kvh  ][jvh  ][ivh+1].z*pnv[1]+
                            ucat[kvh  ][jvh+1][ivh+1].z*pnv[2]+
			    ucat[kvh  ][jvh+1][ivh  ].z*pnv[3]+
                            ucat[kvh+1][jvh  ][ivh  ].z*pnv[4]+
			    ucat[kvh+1][jvh  ][ivh+1].z*pnv[5]+
                            ucat[kvh+1][jvh+1][ivh+1].z*pnv[6]+
			    ucat[kvh+1][jvh+1][ivh  ].z*pnv[7])/
                            (1.e-10+pnv[0]+pnv[1]+pnv[2]+pnv[3]+pnv[4]+pnv[5]+pnv[6]+pnv[7]);
	          }
                 
                  if(Harmonic)
                  {
                  // harmonic averaging among none-solid nodes
                  vx[nv] = (ucat[kvh  ][jvh  ][ivh  ].x*pnv[0]/r[0]+
                            ucat[kvh  ][jvh  ][ivh+1].x*pnv[1]/r[1]+
                            ucat[kvh  ][jvh+1][ivh+1].x*pnv[2]/r[2]+
			    ucat[kvh  ][jvh+1][ivh  ].x*pnv[3]/r[3]+
                            ucat[kvh+1][jvh  ][ivh  ].x*pnv[4]/r[4]+
			    ucat[kvh+1][jvh  ][ivh+1].x*pnv[5]/r[5]+
                            ucat[kvh+1][jvh+1][ivh+1].x*pnv[6]/r[6]+
			    ucat[kvh+1][jvh+1][ivh  ].x*pnv[7]/r[7])/
                            (1.e-10+pnv[0]/r[0]+pnv[1]/r[1]+pnv[2]/r[2]+pnv[3]/r[3]+pnv[4]/r[4]+pnv[5]/r[5]+pnv[6]/r[6]+pnv[7]/r[7]);

                  vy[nv] = (ucat[kvh  ][jvh  ][ivh  ].y*pnv[0]/r[0]+
                            ucat[kvh  ][jvh  ][ivh+1].y*pnv[1]/r[1]+
                            ucat[kvh  ][jvh+1][ivh+1].y*pnv[2]/r[2]+
                            ucat[kvh  ][jvh+1][ivh  ].y*pnv[3]/r[3]+
                            ucat[kvh+1][jvh  ][ivh  ].y*pnv[4]/r[4]+
			    ucat[kvh+1][jvh  ][ivh+1].y*pnv[5]/r[5]+
                            ucat[kvh+1][jvh+1][ivh+1].y*pnv[6]/r[6]+
			    ucat[kvh+1][jvh+1][ivh  ].y*pnv[7]/r[7])/
                            (1.e-10+pnv[0]/r[0]+pnv[1]/r[1]+pnv[2]/r[2]+pnv[3]/r[3]+pnv[4]/r[4]+pnv[5]/r[5]+pnv[6]/r[6]+pnv[7]/r[7]);

                  vz[nv] = (ucat[kvh  ][jvh  ][ivh  ].z*pnv[0]/r[0]+
                            ucat[kvh  ][jvh  ][ivh+1].z*pnv[1]/r[1]+
                            ucat[kvh  ][jvh+1][ivh+1].z*pnv[2]/r[2]+
			    ucat[kvh  ][jvh+1][ivh  ].z*pnv[3]/r[3]+
                            ucat[kvh+1][jvh  ][ivh  ].z*pnv[4]/r[4]+
			    ucat[kvh+1][jvh  ][ivh+1].z*pnv[5]/r[5]+
                            ucat[kvh+1][jvh+1][ivh+1].z*pnv[6]/r[6]+
			    ucat[kvh+1][jvh+1][ivh  ].z*pnv[7]/r[7])/
                            (1.e-10+pnv[0]/r[0]+pnv[1]/r[1]+pnv[2]/r[2]+pnv[3]/r[3]+pnv[4]/r[4]+pnv[5]/r[5]+pnv[6]/r[6]+pnv[7]/r[7]);
		   }
	
                	  flag = PETSC_TRUE;
			   	  
		  goto next;
		}  //is inside the cubic cell "if" 
	      }   //if the ijk of the biger block in the local cpu's ijk "if"
	    }   //loop over all the fluid nodes in the box "for"
	  }   //if there is any fluid cell in the box "if"
	}   // loop over iv  for the boxes from start to end "for"
      }   // loop over jv for the boxes from start to end "for"
    }   // loop over kv for the boxes from start to end "for"
  }  // if the triangle bed cell is not on the side of IB of underbeneath of IB "if"

next: nothing = 0;
} // loop over the bed triangle cells "for"



if(TetQuad_search){
for (nv=0;nv<n_elmt;nv++){
double Maxscc=0.;
double sccc = ibm->scc[nv];
MPI_Allreduce(&sccc, &Maxscc,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
PetscBarrier( PETSC_NULL );
ibm->scc[nv]=Maxscc;
}}



DAVecRestoreArray(fda, Coor, &coor);
DAVecRestoreArray(fda, user->lUcat, &ucat);
//DAVecRestoreArray(da, user->P, &prs);
DAVecRestoreArray(da, user->lNvert, &nvert);

 return(0);
}



PetscErrorCode final_implementing_Projecting( UserCtx * user, IBMNodes * ibm, PetscInt nelmt )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;
        IBMInfo         *ibminfo;
	PetscInt        ibi; 
       
	Cmpnts	        *Bvel;
	PetscReal       *Shvel;
	PetscInt        n_v = ibm->n_v;
	PetscInt        i, j, k;
	PetscInt        nii, kkk;
	PetscReal       ucx, ucy, ucz, riter;
        
        PetscInt	elmt;

		riter = 0.;
                for(kkk = 0; kkk<4; kkk++){
		for(nii = 0; nii < ibm->n_elmt; nii++)
			{
				if(ibm->Rigidity[nii] == 1)
                                  { 
                                  }
                                     else
                                         {
		                           double testi = fabs(ibm->Bvel_u[nii])+fabs(ibm->Bvel_v[nii])+fabs(ibm->Bvel_w[nii]);
                                           if( nii == nelmt || testi < 1.e-7 )
                                           {
                                           }
                                           else
                                           {
					
                    //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
     	if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ]== ibm->nv3[ nii ] || ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ]== ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])  

                         			{
riter = riter + 1.;
ibm->Bvel_u[ nelmt] = ( ibm->Bvel_u[ nii ] + ibm->Bvel_u[ nelmt ] * ( riter - 1. ) ) / riter;
ibm->Bvel_v[ nelmt] = ( ibm->Bvel_v[ nii ] + ibm->Bvel_v[ nelmt ] * ( riter - 1. ) ) / riter;
ibm->Bvel_w[ nelmt] = ( ibm->Bvel_w[ nii ] + ibm->Bvel_w[ nelmt ] * ( riter - 1. ) ) / riter;
						}
					    }
			                 }
			 }
}
return ( 0 );
}


PetscErrorCode Projecting_new( UserCtx * user, IBMNodes * ibm )
{
DA               da = user->da, fda = user->fda;
PetscReal        *NIB;
PetscInt         processors,rank;
PetscInt         number_of_trianglecells, total_trianglecells;
PetscInt         *rec, *dis, *stride, start;
   
PetscReal        uu, vv, ww;
PetscReal        *current_uvel, *current_vvel, *current_wvel;       
PetscReal        *uvel_buffer, *vvel_buffer, *wvel_buffer;

PetscReal       *Bvel_u, *Bvel_v, *Bvel_w;
  
PetscInt         i,j,elmt;
PetscInt         *element_buffer, *current_element;
MPI_Datatype     stype, ltype;

MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
MPI_Comm_size(PETSC_COMM_WORLD, &processors);

// Creating the buffer to send and recieve information

PetscMalloc(processors*sizeof(PetscInt), &rec);
PetscMalloc(processors*sizeof(PetscInt), &dis);
PetscMalloc(processors*sizeof(PetscInt), &stride);

// Finding the total number of trianle cells on each cpu
// already known as:
number_of_trianglecells = ibm->n_elmt;

// Finding total number of triangle cells in all cpu's to all cpu's knowing the total
MPI_Allreduce(&number_of_trianglecells,&total_trianglecells,1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

// All cpu's sending thier number of triangle cells to the root cpu naming "stride" in root
MPI_Gather(&number_of_trianglecells, 1, MPI_INT, stride, 1, MPI_INT, 0, MPI_COMM_WORLD);

// Now root Bcasts nomber of each cpu's triangle cell as stride[i] to all other processors
MPI_Bcast(stride, processors, MPI_INT, 0, MPI_COMM_WORLD);
 
// Allocate the rec and dis (placment) to each and all processors
start = 0;
for (i=0;i<processors;i++)
{ 
  dis[i] = start;
  rec[i] = stride[i];
  start = start + stride[i];
}

PetscMalloc(number_of_trianglecells*sizeof(PetscInt), &current_element);
//PetscMalloc(number_of_trianglecells*sizeof(PetscReal), &current_utau);
PetscMalloc(number_of_trianglecells*sizeof(PetscReal), &current_uvel);
PetscMalloc(number_of_trianglecells*sizeof(PetscReal), &current_vvel);
PetscMalloc(number_of_trianglecells*sizeof(PetscReal), &current_wvel);

PetscMalloc(total_trianglecells*sizeof(PetscInt), &element_buffer);
//PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &utau_buffer);
PetscMalloc(total_trianglecells*sizeof(PetscReal), &uvel_buffer);
PetscMalloc(total_trianglecells*sizeof(PetscReal), &vvel_buffer);
PetscMalloc(total_trianglecells*sizeof(PetscReal), &wvel_buffer);

// Now putting the triangle cells info into a local current array to be sent to the root process later
// point to the first IB node..

for (i=0; i<number_of_trianglecells; i++)
{
         current_element[i] = i;
       // current_utau[i] = ibminfo->utaub;
        current_uvel[i] = ibm->vx[i];
        current_vvel[i] = ibm->vy[i];
        current_wvel[i] = ibm->vz[i];
 
} // End of local triangle cell on each process 

// first sending the number of elements to the ROOT
// then sending the local info from processes to ROOT
MPI_Type_contiguous(number_of_trianglecells, MPI_INT, &stype);
MPI_Type_commit(&stype);
MPI_Gatherv(current_element, 1, stype, element_buffer, rec, dis, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Type_contiguous(number_of_trianglecells, MPI_DOUBLE, &ltype);
MPI_Type_commit(&ltype);
//MPI_Gatherv(current_utau, 1, ltype, utau_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(current_uvel, 1, ltype, uvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(current_vvel, 1, ltype, vvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(current_wvel, 1, ltype, wvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// Now all info are at ROOT, then ROOt will postprocess the info to assemble the info on the bed surface

if(rank == 0)
{
  PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &NIB);
  // clear former info
  for(i=0; i<ibm->n_elmt;i++)
   {
   ibm->Bvel_u[i] = 0.0;
   ibm->Bvel_v[i] = 0.0;
   ibm->Bvel_w[i] = 0.0;
   NIB[i]=0.;
   }

  // Projecting onto the body surface
  for(i=0; i<total_trianglecells; i++)
   {
   elmt = element_buffer[i];
   if(elmt < ibm->n_elmt)
    { 
   //   uutau = utau_buffer[i];
      uu = uvel_buffer[i];
      vv = vvel_buffer[i];
      ww = wvel_buffer[i];
      
      // transforming IB node's info into the body surface element 
      if(fabs(uu) + fabs(vv) + fabs(ww) > 1.e-6)
        {
         ibm->Bvel_u[elmt] += uu;
         ibm->Bvel_v[elmt] += vv;
         ibm->Bvel_w[elmt] += ww;
    //     ibm->Shvel[elmt] += uutau;
         NIB[elmt] += 1.;
        }
     }
   } // End of all cpu triangle cells 
 
   // Now averaging the variables on each triangle cell of body surface
   for(i=0; i<ibm->n_elmt; i++)
   {
       if(NIB[i]>0)
         {
          ibm->Bvel_u[i] /= NIB[i];
          ibm->Bvel_v[i] /= NIB[i];
          ibm->Bvel_w[i] /= NIB[i];
         // ibm->Shvel[i] /= NIB[i];
         }
      //if(fabs(ibm->Bvel_u[i])>0.) PetscPrintf(PETSC_COMM_WORLD, "Bvelx,Bvely,Bvelz %e %e %e\n",ibm->Bvel_u[i],ibm->Bvel_v[i],ibm->Bvel_w[i]);
   }
   
         
PetscFree(NIB);
// Now ROOT Bcast the info to all other processes
MPI_Bcast(ibm->Bvel_u, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
MPI_Bcast(ibm->Bvel_v, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
MPI_Bcast(ibm->Bvel_w, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
//MPI_Bcast(ibm->Shvel, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);

}  // End of ROOT process work (end of rank == 0)


// other cpu's getting the info from ROOT here
if(rank)
  {
    MPI_Bcast(ibm->Bvel_u, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(ibm->Bvel_v, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(ibm->Bvel_w, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
  //  MPI_Bcast(ibm->Shvel, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
   
    //for(i=0;i<ibm->n_elmt; i++)
    // {
    //     ibm->Bvel[i].x = ibm->Bvel_u[i];
    //     ibm->Bvel[i].y = ibm->Bvel_v[i];
    //     ibm->Bvel[i].z = ibm->Bvel_w[i];
    // }   
  } // End of recieving info from ROOT

PetscBarrier(PETSC_NULL);

PetscFree(current_element);
//PetscFree(current_utau);
PetscFree(current_uvel);
PetscFree(current_vvel);
PetscFree(current_wvel);
PetscFree(element_buffer);
//PetscFree(utau_buffer);
PetscFree(uvel_buffer);
PetscFree(vvel_buffer);
PetscFree(wvel_buffer);

// if an elmt doesnt have value on it yet, it should be interpolated from neighbor elmts
 for(i=0;i<ibm->n_elmt;i++)
 {
 double test = fabs(ibm->Bvel_u[i])+fabs(ibm->Bvel_v[i])+fabs(ibm->Bvel_w[i]);
 if(test < 1.e-6) final_implementing_Projecting(user, ibm, i);
 }

return (0);
}


PetscErrorCode flow_variables_ref_level_new (UserCtx *user, IBMNodes *ibm, PetscInt ti, PetscInt tistart)
{
	DA		 da  = user->da, fda  = user->fda;
	DALocalInfo      info  = user->info;

	PetscInt        n_elmt  = ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt	   i, j, k, ni, nv;
        PetscReal	   ref_level, uutau;
        //PetscReal          sbbb = 0.025; // jhoke-crossvane-rockvane  0.025;

if(ti==tistart || ti==0) PetscPrintf(PETSC_COMM_WORLD, "sbb: %e \n",sbbb);

for (ni=0; ni<n_elmt; ni++) 
 {
   if(ibm->Rigidity[ni] == 1)
   {
    ibm->Shvel[ni]=0.;
    ibm->Bvel[ni].x=0.;
    ibm->Bvel[ni].y=0.;
    ibm->Bvel[ni].z=0.;
   }
   else
   {

   ibm->sbb[ni] = sbbb;      
 // ibm->scc[ni] = sbbb;      
   Cmpnts Ua, Uc, Urf;
	if (ni>=0) {
	  Ua.x = (ibm->u[ibm->nv1[ni]].x + ibm->u[ibm->nv2[ni]].x + ibm->u[ibm->nv3[ni]].x )/3.;
	  Ua.y = (ibm->u[ibm->nv1[ni]].y + ibm->u[ibm->nv2[ni]].y + ibm->u[ibm->nv3[ni]].y )/3.;
	  Ua.z = (ibm->u[ibm->nv1[ni]].z + ibm->u[ibm->nv2[ni]].z + ibm->u[ibm->nv3[ni]].z )/3.;
	           }
     	else       {
	Ua.x = Ua.y = Ua.z = 0.;
                   }
			
	Uc.x = ibm->Bvel_u[ni];
	Uc.y = ibm->Bvel_v[ni];
	Uc.z = ibm->Bvel_w[ni];

    if(!bed_roughness)wall_function_roughness_a (user, 0.001, ibm->scc[ni], ibm->sbb[ni], Ua, Uc, &Urf, &uutau, ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	if(bed_roughness)wall_function_roughness_a (user, k_ss, ibm->scc[ni], ibm->sbb[ni], Ua, Uc, &Urf, &uutau, ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	// if(!bed_roughness)wall_function_roughness_a (user, 0.00000000001, ibm->scc[ni], ibm->sbb[ni], Ua, Uc, &Urf, &uutau, ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	//if(bed_roughness)wall_function_roughness_a (user, k_ss, sc, 0.01, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);

    double ucx = Urf.x;
    double ucy = Urf.y;
    double ucz = Urf.z;
    //double ucx = Uc.x;
    //double ucy = Uc.y;
    //double ucz = Uc.z;
    double nfx = ibm->nf_x[ni];
    double nfy = ibm->nf_y[ni];
    double nfz = ibm->nf_z[ni];

    ibm->Bvel[ni].x = ucx - (ucx * nfx + ucy * nfy + ucz * nfz) * nfx;
    ibm->Bvel[ni].y = ucy - (ucx * nfx + ucy * nfy + ucz * nfz) * nfy;
    ibm->Bvel[ni].z = ucz - (ucx * nfx + ucy * nfy + ucz * nfz) * nfz;
    ibm->Shvel[ni] = uutau;
     
    }
 }
	
 return (0);
}

double Half_Depth_Averaged(UserCtx *user, PetscInt i, PetscInt j_current, PetscInt k, PetscInt y_direction, PetscInt index_ave)
{
DA		da = user->da, fda = user->fda, fda2 = user->fda2;
DALocalInfo	info = user->info;

Cmpnts          ***ucat;
PetscReal       ***nvert;
Cmpnts2  	***K_Omega;

PetscInt	xs = info.xs, xe = info.xs + info.xm;
PetscInt  	ys = info.ys, ye = info.ys + info.ym;
PetscInt	zs = info.zs, ze = info.zs + info.zm;
PetscInt	lxs, lxe, lys, lye, lzs, lze;

 lxs = xs; lxe = xe;
 lys = ys; lye = ye;
 lzs = zs; lze = ze;

 if (xs==0) lxs = xs+1;
 if (ys==0) lys = ys+1;
 if (zs==0) lzs = zs+1;

 if (xe==info.mx) lxe = xe-1;
 if (ye==info.my) lye = ye-1;
 if (ze==info.mz) lze = ze-1;

DAVecGetArray(fda, user->lUcat, &ucat);
DAVecGetArray(da, user->lNvert, &nvert);
DAVecGetArray(fda2, user->lK_Omega, &K_Omega);

	PetscInt	jj; 
        PetscInt        jj_mid_depth = j_current + int((lye - j_current)/2);
        if(j_current >= lye-3) jj_mid_depth = lye;
        PetscInt        j_end = PetscMin(lye,jj_mid_depth);
        double ave_half_depth;
        double sum = 0.;
               ave_half_depth = 0.;
       for (jj=j_current; jj<lye; jj++)
       //for (jj=lys+1; jj<j_end; jj++)
             {
		if ( nvert[k][jj][i]<1.1) {
                         sum += 1.;
			 if(!y_direction && index_ave == 3) ave_half_depth += ucat[k][jj][i].z; 
 			 if(y_direction && index_ave == 3) ave_half_depth += ucat[k][jj][i].y; 

			 if(index_ave == 4) ave_half_depth += K_Omega[k][jj][i].x; 
                         }
             }
         if(sum>0.) {ave_half_depth = ave_half_depth/sum;} else {ave_half_depth = 0.;}                
 
DAVecRestoreArray(fda, user->lUcat, &ucat);
DAVecRestoreArray(da, user->lNvert, &nvert);
DAVecRestoreArray(fda2, user->lK_Omega, &K_Omega);

return ave_half_depth;
}


PetscErrorCode flow_variables_ref_level (UserCtx *user, IBMNodes *ibm, PetscInt ti, PetscInt tistart)
{
	DA		 da  = user->da, fda  = user->fda;
	DALocalInfo      info  = user->info;

	PetscInt	 xs  = info.xs, xe  = info.xs + info.xm;
	PetscInt	 ys  = info.ys, ye  = info.ys + info.ym;
	PetscInt	 zs  = info.zs, ze  = info.zs + info.zm;
	PetscInt	 mx  = info.mx, my  = info.my, mz  = info.mz;
	PetscInt	 lxs, lys, lzs, lxe, lye, lze;
	lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

	if( xs == 0 )lxs = xs + 1;
	if( ys == 0 )lys = ys + 1;
	if( zs == 0 )lzs = zs + 1;

	if( xe == mx )lxe = xe - 1;
	if( ye == my )lye = ye - 1;
	if( ze == mz )lze = ze - 1;

        //IBMNodes *ibm = &ibm_ptr[ibi];
	Cmpnts	   *** ucat, *** lucat;
	PetscReal  *** ustar;
	PetscReal  *** conc, ***lnu_t;
	PetscInt       n_elmt  = ibm->n_elmt, 
		       n_v  = ibm->n_v;
	PetscInt	   i, j, k;
        PetscReal	   ref_level;
	//PetscInt	   ip1,	ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
	//PetscReal	   sk1,	sk2, sk3, cv1, cv2, cv3;
	PetscReal  *** nvert, *** nvert_o;
//	PetscReal	cs1, cs2, cs3;
//	PetscInt	nii;
	IBMInfo     *   ibminfo;
//	PetscInt	NumberOfBodies = 0, ibi;
	//	PetscReal	nfx, 
	//		nfy, 
	//		nfz;
          PetscInt   ibi;              
         //PetscReal   sbbb = 0.025; // all stream structures  0.025;
	 double cc; 
	 double w_half_depth_ave; 
	 double k_half_depth_ave; 

if(ti==tistart || ti==0) PetscPrintf(PETSC_COMM_WORLD, "sbb: %e \n",sbbb);

		DAVecGetArray(da, user->lUstar, &ustar);
		DAVecGetArray(fda, user->Ucat,  &ucat);
		DAVecGetArray(fda, user->lUcat,  &lucat);
		if(LiveBed) DAVecGetArray(da, user->lConc, &conc);
		if(LiveBed) DAVecGetArray(da, user->lNu_t, &lnu_t);
		if(effective_bed_shear_stress) DAVecGetArray(da, user->lNvert, &nvert);

		for(ibi=0; ibi<NumberOfBodies; ibi++)
		{
			IBMListNode *current;
			current = user->ibmlist[ibi].head;
			while (current) {
				IBMInfo *ibminfo = &current->ibm_intp;
				current = current->next;
				
				int ni = ibminfo->cell;
				int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
				int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
				int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
				i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
				
				double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
				double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
				double cs1 = ibminfo->cs1, cs2 = ibminfo->cs2, cs3 = ibminfo->cs3;
                                
	//			  if (ti==tistart){ref_level = sb;} else {ref_level = ibm->Deltab;}

				Cmpnts Ua, Uc;

			        int rigid = ibm->Rigidity[ni];
                                	
        			if (ni>=0) {
					Ua.x = ibm->u[ibm->nv1[ni]].x * cs1 + ibm->u[ibm->nv2[ni]].x * cs2 + ibm->u[ibm->nv3[ni]].x * cs3;
					Ua.y = ibm->u[ibm->nv1[ni]].y * cs1 + ibm->u[ibm->nv2[ni]].y * cs2 + ibm->u[ibm->nv3[ni]].y * cs3;
					Ua.z = ibm->u[ibm->nv1[ni]].z * cs1 + ibm->u[ibm->nv2[ni]].z * cs2 + ibm->u[ibm->nv3[ni]].z * cs3;
			                   }
				else       {
					Ua.x = Ua.y = Ua.z = 0.;
	                                   }
	 			
        			Uc.x = (lucat[kp1][jp1][ip1].x * sk1 + lucat[kp2][jp2][ip2].x * sk2 +lucat[kp3][jp3][ip3].x * sk3);
				Uc.y = (lucat[kp1][jp1][ip1].y * sk1 + lucat[kp2][jp2][ip2].y * sk2 + lucat[kp3][jp3][ip3].y * sk3);
				Uc.z = (lucat[kp1][jp1][ip1].z * sk1 + lucat[kp2][jp2][ip2].z * sk2 + lucat[kp3][jp3][ip3].z * sk3);
				if(LiveBed) {
                                             double c1 = (conc[kp1][jp1][ip1] * sk1 + conc[kp2][jp2][ip2] * sk2 + conc[kp3][jp3][ip3] * sk3);
                                             double c  = PetscMax (c1 , conc[k][j][i]);
                                             double nu = 1./user->ren;
                                                    if(levelset) nu = 1.e-6;
                                             double sigma_phi = 0.75;
                                             double nu_t = lnu_t[kp1][jp1][ip1] * sk1 + lnu_t[kp2][jp2][ip2] * sk2 +lnu_t[kp3][jp3][ip3] * sk3;
                                                    nu_t = PetscMax (nu_t, 0.0);
                                             double nu_tt = nu + sigma_phi * nu_t;
                                             double alfaa = fabs (sc - sbbb);
                                                    alfaa = alfaa * w_s / nu_tt;
                                                    alfaa = exp (alfaa);
                                                    alfaa = PetscMin (1.0, alfaa);
                                                    alfaa = 1.0 - alfaa;
                                             if(!rigid) {cc = c + ibm->SCont[ni] * alfaa;} else {cc = c;}
                                            }
                                if(effective_bed_shear_stress) {w_half_depth_ave = Half_Depth_Averaged (user,i,j,k,y_direction,3);
                                                                if(rans) k_half_depth_ave = Half_Depth_Averaged (user,i,j,k,y_direction,4);
                                                               }
				
    if(!bed_roughness)wall_function_roughness_a (user, 0.001, sc, sbbb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	if(bed_roughness)wall_function_roughness_a (user, k_ss, sc, sbbb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
     
 if( i >= lxs && i < lxe && j >= lys && j < lye && k >= lzs && k < lze )
   {
	 ibminfo->ucomp = ucat[k][j][i].x;
	 ibminfo->vcomp = ucat[k][j][i].y;                                                 
         ibminfo->wcomp = ucat[k][j][i].z;
         ibminfo->utaub = ustar[k][j][i];
         if(LiveBed) ibminfo->CC = cc;
         if(effective_bed_shear_stress) {
                           ibminfo->w_v = w_half_depth_ave;
                           if(!y_direction && (w_half_depth_ave = 0.)) ibminfo->w_v = ucat[k][j][i].z;
						   if(y_direction && (w_half_depth_ave = 0.))  ibminfo->w_v = ucat[k][j][i].y;

                           if(rans) ibminfo->k_v = k_half_depth_ave;
                                        }
   }
			}
		}
		
		if(LiveBed) DAVecRestoreArray(da, user->lConc, &conc);
		if(effective_bed_shear_stress) DAVecRestoreArray(da, user->lNvert, &nvert);
		DAVecRestoreArray(da, user->lUstar, &ustar);
		DAVecRestoreArray(fda, user->Ucat,  &ucat);
		DAVecRestoreArray(fda, user->lUcat,  &lucat);
		DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
		DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	
 return (0);
}


PetscErrorCode flow_variables_restoring (UserCtx *user, IBMNodes *ibm)
{
	DA		 da  = user->da, fda  = user->fda;
	DALocalInfo      info  = user->info;

	PetscInt	 xs  = info.xs, xe  = info.xs + info.xm;
	PetscInt	 ys  = info.ys, ye  = info.ys + info.ym;
	PetscInt	 zs  = info.zs, ze  = info.zs + info.zm;
	PetscInt	 mx  = info.mx, my  = info.my, mz  = info.mz;
	PetscInt	 lxs, lys, lzs, lxe, lye, lze;
	lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

	if( xs == 0 )lxs = xs + 1;
	if( ys == 0 )lys = ys + 1;
	if( zs == 0 )lzs = zs + 1;

	if( xe == mx )lxe = xe - 1;
	if( ye == my )lye = ye - 1;
	if( ze == mz )lze = ze - 1;

        //IBMNodes *ibm = &ibm_ptr[ibi];
	Cmpnts	   *** ucat, *** lucat;
	PetscReal  *** ustar;
     	PetscInt       n_elmt  = ibm->n_elmt, 
		       n_v  = ibm->n_v;
	PetscInt	   i, j, k;
        PetscReal	   ref_level;
	//PetscInt	   ip1,	ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
	//PetscReal	   sk1,	sk2, sk3, cv1, cv2, cv3;
	PetscReal  *** nvert, *** nvert_o;
//	PetscReal	cs1, cs2, cs3;
//	PetscInt	nii;
	IBMInfo     *   ibminfo;
//	PetscInt	NumberOfBodies = 0, ibi;
	//	PetscReal	nfx, 
	//		nfy, 
	//		nfz;
          PetscInt   ibi;              
	  
		DAVecGetArray(da, user->lUstar, &ustar);
		DAVecGetArray(fda, user->Ucat,  &ucat);
		DAVecGetArray(fda, user->lUcat,  &lucat);

		for(ibi=0; ibi<NumberOfBodies; ibi++)
		{
			IBMListNode *current;
			current = user->ibmlist[ibi].head;
			while (current) {
				IBMInfo *ibminfo = &current->ibm_intp;
				current = current->next;
				
				int ni = ibminfo->cell;
				int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
				int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
				int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
				i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
				
				double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
				double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
				double cs1 = ibminfo->cs1, cs2 = ibminfo->cs2, cs3 = ibminfo->cs3;
                                
				Cmpnts Ua, Uc;
				
				if (ni>=0) {
					Ua.x = ibm->u[ibm->nv1[ni]].x * cs1 + ibm->u[ibm->nv2[ni]].x * cs2 + ibm->u[ibm->nv3[ni]].x * cs3;
					Ua.y = ibm->u[ibm->nv1[ni]].y * cs1 + ibm->u[ibm->nv2[ni]].y * cs2 + ibm->u[ibm->nv3[ni]].y * cs3;
					Ua.z = ibm->u[ibm->nv1[ni]].z * cs1 + ibm->u[ibm->nv2[ni]].z * cs2 + ibm->u[ibm->nv3[ni]].z * cs3;
				}
				else {
					Ua.x = Ua.y = Ua.z = 0.;
				}
				
				Uc.x = (lucat[kp1][jp1][ip1].x * sk1 + lucat[kp2][jp2][ip2].x * sk2 +lucat[kp3][jp3][ip3].x * sk3);
				Uc.y = (lucat[kp1][jp1][ip1].y * sk1 + lucat[kp2][jp2][ip2].y * sk2 + lucat[kp3][jp3][ip3].y * sk3);
				Uc.z = (lucat[kp1][jp1][ip1].z * sk1 + lucat[kp2][jp2][ip2].z * sk2 + lucat[kp3][jp3][ip3].z * sk3);
				
    if(!bed_roughness)wall_function_roughness_a (user, 0.001, ibm->scc[ni], ibm->sbb[ni], Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	if(bed_roughness)wall_function_roughness_a (user, k_ss, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);		
			        }
		}
		
		DAVecRestoreArray(da, user->lUstar, &ustar);
		DAVecRestoreArray(fda, user->Ucat,  &ucat);
		DAVecRestoreArray(fda, user->lUcat,  &lucat);
		DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
		DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
	
 return (0);
}



PetscInt IB_node_total_each_cpu(UserCtx *user, IBMNodes *ibm)
{
IBMInfo      *ibminfo;
IBMListNode  *current;
PetscInt      number_of_IBnodes_at_this_cpu;
// point to the first node at each cpu
//IBMListNode * current = user->ibmlist[0].head;
current = user->ibmlist[0].head;
number_of_IBnodes_at_this_cpu = 0;
while (current)
{
IBMInfo *ibminfo = &current->ibm_intp;
  current = current->next;
  number_of_IBnodes_at_this_cpu++;
}
return number_of_IBnodes_at_this_cpu;
} 


PetscErrorCode Projecting( UserCtx * user, IBMNodes * ibm )
{
IBMInfo          *ibminfo;
//IBMListNode      *current;
PetscInt         *current_element;
DA               da = user->da, fda = user->fda;

PetscReal         *NIB;
PetscInt         processors,rank;
PetscInt         number_of_IBnodes=0, total_IBnodes;
PetscInt         *rec, *dis, *stride, start;
   
PetscReal        uu, vv, ww, uutau, cc, w_v_, k_v_;
PetscReal        *current_cc, *current_utau, *current_uvel, *current_vvel, *current_wvel, *current_w_v, *current_k_v;       
PetscReal        *cc_buffer, *utau_buffer, *uvel_buffer, *vvel_buffer, *wvel_buffer, *w_v_buffer, *k_v_buffer;

PetscReal       *Bvel_u, *Bvel_v, *Bvel_w;
Cmpnts	        *Bvel;
PetscReal       *Shvel, *C, *w_ave, *k_ave;
  
PetscInt         i,j,elmt;
PetscInt         *element_buffer;
MPI_Datatype     stype, ltype;

MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
MPI_Comm_size(PETSC_COMM_WORLD, &processors);

// Creating the buffer to send and recieve information

PetscMalloc(processors*sizeof(PetscInt), &rec);
PetscMalloc(processors*sizeof(PetscInt), &dis);
PetscMalloc(processors*sizeof(PetscInt), &stride);

// Finding the total number of IB nodes on each cpu
number_of_IBnodes = IB_node_total_each_cpu(user, ibm);

// Finding total number of IB nodes in all cpu's to all cpu's knowing the total
MPI_Allreduce(&number_of_IBnodes,&total_IBnodes,1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

// All cpu's sending thier number of IB nodes to the root cpu naming "stride" in root
MPI_Gather(&number_of_IBnodes, 1, MPI_INT, stride, 1, MPI_INT, 0, MPI_COMM_WORLD);

// Now root Bcasts nomber of each cpu's IB node as stride[i] to all other processors
MPI_Bcast(stride, processors, MPI_INT, 0, MPI_COMM_WORLD);
 
// Allocate the rec and dis (placment) to each and all processors
start = 0;
for (i=0;i<processors;i++)
{ 
  dis[i] = start;
  rec[i] = stride[i];
  start = start + stride[i];
}

PetscMalloc(number_of_IBnodes*sizeof(PetscInt), &current_element);
PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_utau);
if(LiveBed) PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_cc);
if(effective_bed_shear_stress) PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_w_v);
if(effective_bed_shear_stress && rans) PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_k_v);
PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_uvel);
PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_vvel);
PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_wvel);


PetscMalloc(total_IBnodes*sizeof(PetscInt), &element_buffer);
PetscMalloc(total_IBnodes*sizeof(PetscReal), &utau_buffer);
if(LiveBed) PetscMalloc(total_IBnodes*sizeof(PetscReal), &cc_buffer);
if(effective_bed_shear_stress) PetscMalloc(total_IBnodes*sizeof(PetscReal), &w_v_buffer);
if(effective_bed_shear_stress && rans) PetscMalloc(total_IBnodes*sizeof(PetscReal), &k_v_buffer);
PetscMalloc(total_IBnodes*sizeof(PetscReal), &uvel_buffer);
PetscMalloc(total_IBnodes*sizeof(PetscReal), &vvel_buffer);
PetscMalloc(total_IBnodes*sizeof(PetscReal), &wvel_buffer);

// Now putting the IB infor into a local current array to be sent to the root process later
// point to the first IB node..

IBMListNode * current = user->ibmlist[0].head;
//current = user->ibmlist.head;
i = 0;

while(current)
{
 ibminfo = &current->ibm_intp;
 current=current->next;
 
 // point to the related neighbour bed cell
 elmt = ibminfo->cell;

 if (i<number_of_IBnodes)
      {
        current_element[i] = elmt;
        current_utau[i] = ibminfo->utaub;
        if(LiveBed) current_cc[i] = ibminfo->CC;
        if(effective_bed_shear_stress) current_w_v[i] = ibminfo->w_v;
        if(effective_bed_shear_stress && rans) current_k_v[i] = ibminfo->k_v;
        current_uvel[i] = ibminfo->ucomp;
        current_vvel[i] = ibminfo->vcomp;
        current_wvel[i] = ibminfo->wcomp;
      }
 i++;
} // End of local IB node on each process ( end of IBM list nodes)

// first sending the number of elements to the ROOT
// then sending the local info from processes to ROOT
MPI_Type_contiguous(number_of_IBnodes, MPI_INT, &stype);
MPI_Type_commit(&stype);
MPI_Gatherv(current_element, 1, stype, element_buffer, rec, dis, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Type_contiguous(number_of_IBnodes, MPI_DOUBLE, &ltype);
MPI_Type_commit(&ltype);
MPI_Gatherv(current_utau, 1, ltype, utau_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
if(LiveBed) MPI_Gatherv(current_cc, 1, ltype, cc_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
if(effective_bed_shear_stress) MPI_Gatherv(current_w_v, 1, ltype, w_v_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
if(effective_bed_shear_stress && rans) MPI_Gatherv(current_k_v, 1, ltype, k_v_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(current_uvel, 1, ltype, uvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(current_vvel, 1, ltype, vvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gatherv(current_wvel, 1, ltype, wvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// Now all info are at ROOT, then ROOt will postprocess the info to assemble the info on the bed surface


if(rank == 0)
{
  PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &NIB);

  // clear former info
  for(i=0; i<ibm->n_elmt;i++)
   {

   ibm->Bvel[i].x = 0.0;
   ibm->Bvel[i].y = 0.0;
   ibm->Bvel[i].z = 0.0;
   ibm->Shvel[i] = 0.0;
   if(LiveBed) ibm->C[i] = 0.0;
   if(effective_bed_shear_stress) ibm->w_ave[i] = 0.0;
   if(effective_bed_shear_stress && rans) ibm->k_ave[i] = 0.0;
   ibm->Bvel_u[i] = 0.0;
   ibm->Bvel_v[i] = 0.0;
   ibm->Bvel_w[i] = 0.0;
   NIB[i]=0.;
   }

  // Projecting onto the body surface
  for(i=0; i<total_IBnodes; i++)
   {
   elmt = element_buffer[i];
   if(elmt < ibm->n_elmt)
    { 
      uutau = utau_buffer[i];
      if(LiveBed) cc = cc_buffer[i];
      if(effective_bed_shear_stress) w_v_ = w_v_buffer[i];
      if(effective_bed_shear_stress && rans) k_v_ = k_v_buffer[i];
      uu = uvel_buffer[i];
      vv = vvel_buffer[i];
      ww = wvel_buffer[i];
      
      // transforming IB node's info into the body surface element 
      double tte;
      if(LiveBed) { tte = fabs(uutau) + fabs(uu) + fabs(vv) + fabs(ww) + fabs(cc); } else { tte = fabs(uutau) + fabs(uu) + fabs(vv) + fabs(ww); }
      if(tte > 0.0)
        {
         ibm->Bvel[elmt].x += uu;
         ibm->Bvel[elmt].y += vv;
         ibm->Bvel[elmt].z += ww;
         ibm->Shvel[elmt] += uutau;
         if(LiveBed) ibm->C[elmt] += cc;
         if(effective_bed_shear_stress) ibm->w_ave[elmt] += w_v_;
         if(effective_bed_shear_stress && rans) ibm->k_ave[elmt] += k_v_;
         NIB[elmt] += 1.;
        }
     }
   } // End of all IB nodes 
 
   // Now averaging the variables on each triangle cell of body surface
   for(i=0; i<ibm->n_elmt; i++)
   {
       if(NIB[i]>0)
         {
          ibm->Bvel[i].x /= NIB[i];
          ibm->Bvel[i].y /= NIB[i];
          ibm->Bvel[i].z /= NIB[i];
          ibm->Shvel[i] /= NIB[i];
          if(LiveBed) ibm->C[i] /= NIB[i];
          if(effective_bed_shear_stress) ibm->w_ave[i] /= NIB[i];
          if(effective_bed_shear_stress && rans) ibm->k_ave[i] /= NIB[i];
         }
   }
   
   for(i=0; i<ibm->n_elmt; i++)
     {
         ibm->Bvel_u[i] = ibm->Bvel[i].x;
         ibm->Bvel_v[i] = ibm->Bvel[i].y;
         ibm->Bvel_w[i] = ibm->Bvel[i].z;
     }   
         
PetscFree(NIB);
// Now ROOT Bcast the info to all other processes
MPI_Bcast(ibm->Bvel_u, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
MPI_Bcast(ibm->Bvel_v, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
MPI_Bcast(ibm->Bvel_w, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
MPI_Bcast(ibm->Shvel, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
if(LiveBed) MPI_Bcast(ibm->C, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
if(effective_bed_shear_stress) MPI_Bcast(ibm->w_ave, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
if(effective_bed_shear_stress && rans) MPI_Bcast(ibm->k_ave, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);

}  // End of ROOT process work (end of rank == 0)

// other cpu's getting the info from ROOT here
if(rank)
  {
    MPI_Bcast(ibm->Bvel_u, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(ibm->Bvel_v, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(ibm->Bvel_w, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(ibm->Shvel, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    if(LiveBed)MPI_Bcast(ibm->C, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    if(effective_bed_shear_stress)MPI_Bcast(ibm->w_ave, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
    if(effective_bed_shear_stress && rans)MPI_Bcast(ibm->k_ave, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
   
    for(i=0;i<ibm->n_elmt; i++)
     {
         ibm->Bvel[i].x = ibm->Bvel_u[i];
         ibm->Bvel[i].y = ibm->Bvel_v[i];
         ibm->Bvel[i].z = ibm->Bvel_w[i];
     }   
  } // End of recieving info from ROOT


PetscBarrier(PETSC_NULL);

PetscFree(current_element);
PetscFree(current_utau);
if(LiveBed) PetscFree(current_cc);
PetscFree(current_uvel);
PetscFree(current_vvel);
PetscFree(current_wvel);
PetscFree(element_buffer);
PetscFree(utau_buffer);
if(LiveBed) PetscFree(cc_buffer);
PetscFree(uvel_buffer);
PetscFree(vvel_buffer);
PetscFree(wvel_buffer);
if(effective_bed_shear_stress) PetscFree(current_w_v);
if(effective_bed_shear_stress) PetscFree(w_v_buffer);
if(effective_bed_shear_stress && rans) PetscFree(current_k_v);
if(effective_bed_shear_stress && rans) PetscFree(k_v_buffer);


return (0);
}


PetscErrorCode Smoothing_Shear_vel( UserCtx * user, IBMNodes * ibm ,PetscInt ti, PetscInt tistart)
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;
        IBMInfo         *ibminfo;
	PetscInt        ibi; 
       
	PetscReal       *Shvel;
	PetscInt        n_elmt = ibm->n_elmt, 
		        n_v = ibm->n_v;
	PetscInt        i, j, k;
	PetscInt        nii;
	PetscReal       riter;
       
        PetscInt	elmt,nelmt;

        PetscReal	nfx, 
			nfy, 
			nfz;
//Smoothing the variables at all cells to correct the high gradient regions
PetscInt sm_number = 1;
if(ti==tistart || ti==0) PetscPrintf(PETSC_COMM_WORLD, "Number of Smoothening iterations for Effective_Shear_velocity: %d \n",sm_number);

for(i=0;i<sm_number;i++){
	for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )
	{
		if(ibm->Rigidity[nelmt] == 1)
		{
			ibm->Shvel[ nelmt ] = 0.;
			ibm->Shvel_v[ nelmt ] = 0.;
		}
		else
		{
				riter = 0.;
				int ii;
				for(ii=0;ii<200;ii++)  // xiaolei add  // ali corrected
				{
					nii=ibm->c2cs[nelmt].cs[ii]; // xiaolei add  // ali corrected
					if (nii>=0) {   // xiaolei add 
					if(ibm->Rigidity[nelmt] == 1 || nii==nelmt) 
                                         { 
                                         }
                                        else
                                         {

                    //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
     	if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ]== ibm->nv3[ nii ] || ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ]== ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])  

                         			{
riter = riter + 1.;
//ibm->Shvel[ nelmt ] = ( ibm->Shvel[ nii ] + ibm->Shvel[ nelmt ] * ( riter ) ) / (riter+1.); // It is done once in Finalizing_Projecting
ibm->Shvel_v[ nelmt ] = ( ibm->Shvel_v[ nii ] + ibm->Shvel_v[ nelmt ] * ( riter ) ) / (riter+1.);

						}
					   
			              }
				      }
			      }
	             
                }
         }


}
return ( 0 );
}

PetscErrorCode Finalizing_Projecting( UserCtx * user, IBMNodes * ibm ,PetscInt ti, PetscInt tistart)
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;
        IBMInfo         *ibminfo;
	PetscInt        ibi; 
       
	Cmpnts	        *Bvel;
	PetscReal       *Shvel;
	PetscReal       *C,*w_ave;
	PetscInt        n_elmt = ibm->n_elmt, 
		        n_v = ibm->n_v;
	PetscInt        i, j, k;
	PetscInt        nii;
	PetscReal       ucx, ucy, ucz, riter;
       
	PetscReal       w_v, w_v_;
	PetscReal       k_v, k_v_;

        PetscInt	elmt,nelmt;

        PetscReal	nfx, 
			nfy, 
			nfz;
       /* //PetscReal Cs_ = 1.5; //varies between 0.8 to 1.5
        
        if(effective_bed_shear_stress && rans)
          {
             for(elmt=0; elmt<ibm->n_elmt; elmt++)
                {
		nfx = ibm->nf_x[ elmt ];
		nfy = ibm->nf_y[ elmt ];
		nfz = ibm->nf_z[ elmt ];

		ucx = ibm->Bvel[ elmt ].x;
		ucy = ibm->Bvel[ elmt ].y;
		ucz = ibm->Bvel[ elmt ].z;

                  if(ibm->nf_z[elmt] < 1.e-7 || ibm->elmt_depth[nelmt] > 3.)
                   {
                    ibm->Shvel_v[elmt] = 0.;
                   } else {
                        w_v_ = -(ucx * nfx / nfz + ucy * nfy / nfz + ibm->w_ave[elmt])/(1.e-8 + sqrt(1.0 + nfx * nfx /(1.e-8 + nfz * nfz) + nfy * nfy /(1.e-8 + nfz * nfz)));
                        ibm->w_ave [elmt] = fabs(w_v_/(1.e-8 + sqrt(ucx * ucx + ucy * ucy + w_v_ * w_v_)));
                        ibm->Shvel_v[elmt] = Cs_ * 0.456 * sqrt(fabs(ibm->w_ave[elmt] * ibm->k_ave[elmt]));
	                  }
                 }	   
         } */

//	PetscReal ts, te, cputime;     // xiaolei
//        PetscGetTime(&ts);  // xiaolei
	for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )
	{
		if(ibm->Rigidity[nelmt] == 1)
		{
			ibm->Shvel[ nelmt ] = 0.;
			ibm->Bvel[ nelmt ].x = 0.;
			ibm->Bvel[ nelmt ].y = 0.;
			ibm->Bvel[ nelmt ].z = 0.;
			ibm->Shvel[ nelmt ] = 0.;
		}
		else
		{
		 double test;
		 if(LiveBed){ test = fabs(ibm->Shvel[nelmt])+fabs(ibm->Bvel[nelmt].x)+fabs(ibm->Bvel[nelmt].y)+fabs(ibm->Bvel[nelmt].z) + fabs(ibm->C[nelmt]);}
                 else { test = fabs(ibm->Shvel[nelmt])+fabs(ibm->Bvel[nelmt].x)+fabs(ibm->Bvel[nelmt].y)+fabs(ibm->Bvel[nelmt].z);}
		 if(test < 1.e-7)
			{  
				riter = 0.;
			//	for(nii = 0; nii < ibm->n_elmt; nii++)  // xiaolei deactivate 
				//int _nv[3];
				//_nv[0]=ibm->nv1[nelmt];
				//_nv[1]=ibm->nv2[nelmt];
				//_nv[2]=ibm->nv3[nelmt];
	
				//int ii,jj;
				//for(jj=0;jj<3;jj++)
				int ii;
				for(ii=0;ii<200;ii++)
				{
					nii=ibm->c2cs[nelmt].cs[ii];	// xiaolei add // ali corrected
					if (nii>=0) {    // xiaolei add //ali corrected
					if(ibm->Rigidity[nii] == 1)
                                         { 
                                         }
                                        else
                                         {
		 double testi;
		 if(LiveBed){ testi = fabs(ibm->Shvel[nii])+fabs(ibm->Bvel[nii].x)+fabs(ibm->Bvel[nii].y)+fabs(ibm->Bvel[nii].z) + fabs(ibm->C[nii]);}
                 else { testi = fabs(ibm->Shvel[nii])+fabs(ibm->Bvel[nii].x)+fabs(ibm->Bvel[nii].y)+fabs(ibm->Bvel[nii].z);}
                                           if( nii == nelmt || testi < 1.e-7 )
                                           {
                                           }
                                           else
                                           {
					
				//--------- interpolation between around elment (max 3 around elemts)***************
/*
	if(((ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ])
          &&(ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ] == ibm->nv3[ nii ]))

	|| ((ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ])
          &&(ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ]))

	|| ((ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ] == ibm->nv3[ nii ])
         &&(ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])))

*/
                    //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
     	if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ]== ibm->nv3[ nii ] || ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ]== ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])   

                         			{
riter = riter + 1.;
ibm->Shvel[ nelmt ] = ( ibm->Shvel[ nii ] + ibm->Shvel[ nelmt ] * ( riter - 1. ) ) / riter;
if(LiveBed) ibm->C[ nelmt ] = ( ibm->C[ nii ] + ibm->C[ nelmt ] * ( riter - 1. ) ) / riter;
if(effective_bed_shear_stress) ibm->w_ave[ nelmt ] = ( ibm->w_ave[ nii ] + ibm->w_ave[ nelmt ] * ( riter - 1. ) ) / riter;
if(effective_bed_shear_stress && rans) ibm->k_ave[ nelmt ] = ( ibm->k_ave[ nii ] + ibm->k_ave[ nelmt ] * ( riter - 1. ) ) / riter;
ibm->Bvel[ nelmt].x = ( ibm->Bvel[ nii ].x + ibm->Bvel[ nelmt ].x * ( riter - 1. ) ) / riter;
ibm->Bvel[ nelmt].y = ( ibm->Bvel[ nii ].y + ibm->Bvel[ nelmt ].y * ( riter - 1. ) ) / riter;
ibm->Bvel[ nelmt].z = ( ibm->Bvel[ nii ].z + ibm->Bvel[ nelmt ].z * ( riter - 1. ) ) / riter;
						}
					   }
			              }
				      }
			      }
	              }
                }
         }

//        PetscGetTime(&te);  // xiaolei

//	cputime = te-ts;
//        PetscPrintf(PETSC_COMM_WORLD, "Finalizing_Projecting 1 cputime %d %le\n", ti,cputime);

//        PetscGetTime(&ts);  // xiaolei
//Smoothing the variables at all cells to correct the high gradient regions
PetscInt sm_number =1;
if(ti==tistart || ti==0) PetscPrintf(PETSC_COMM_WORLD, "Number of Smoothening iterations for variables projection: %d \n",sm_number);

for(i=0;i<sm_number;i++){
	for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )
	{
		if(ibm->Rigidity[nelmt] == 1)
		{
			ibm->Shvel[ nelmt ] = 0.;
			if(LiveBed)ibm->C[ nelmt ] = 0.;
			ibm->Bvel[ nelmt ].x = 0.;
			ibm->Bvel[ nelmt ].y = 0.;
			ibm->Bvel[ nelmt ].z = 0.;
			if(effective_bed_shear_stress)ibm->w_ave[ nelmt ] = 0.;
			if(effective_bed_shear_stress && rans)ibm->k_ave[ nelmt ] = 0.;
		}
		else
		{
				riter = 0.;
		//		for(nii = 0; nii < ibm->n_elmt; nii++)  // xiaolei deactivate 
				//int _nv[3];
				//_nv[0]=ibm->nv1[nelmt];
				//_nv[1]=ibm->nv2[nelmt];
				//_nv[2]=ibm->nv3[nelmt];
	
				//int ii,jj;
				//for(jj=0;jj<3;jj++)
				int ii;
				for(ii=0;ii<200;ii++)  // xiaolei add  // ali corrected
				{
					nii=ibm->c2cs[nelmt].cs[ii]; // xiaolei add  // ali corrected
					if (nii>=0) {   // xiaolei add 
					if(ibm->Rigidity[nii] == 1 || nii==nelmt) 
                                         { 
                                         }
                                        else
                                         {
		                                     

                    //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
     	if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ]== ibm->nv3[ nii ] || ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ]== ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])  

                         			{
riter = riter + 1.;
if (!effective_bed_shear_stress) ibm->Shvel[ nelmt ] = ( ibm->Shvel[ nii ] + ibm->Shvel[ nelmt ] * ( riter ) ) / (riter+1.);
if(LiveBed) ibm->C[ nelmt ] = ( ibm->C[ nii ] + ibm->C[ nelmt ] * ( riter ) ) / (riter+1.);
//if(effective_bed_shear_stress) ibm->w_ave[ nelmt ] = ( ibm->w_ave[ nii ] + ibm->w_ave[ nelmt ] * ( riter ) ) / (riter+1.);
//if(effective_bed_shear_stress && rans) ibm->k_ave[ nelmt ] = ( ibm->k_ave[ nii ] + ibm->k_ave[ nelmt ] * ( riter ) ) / (riter+1.);
ibm->Bvel[ nelmt].x = ( ibm->Bvel[ nii ].x + ibm->Bvel[ nelmt ].x * ( riter ) ) / (riter+1.);
ibm->Bvel[ nelmt].y = ( ibm->Bvel[ nii ].y + ibm->Bvel[ nelmt ].y * ( riter ) ) / (riter+1.);
ibm->Bvel[ nelmt].z = ( ibm->Bvel[ nii ].z + ibm->Bvel[ nelmt ].z * ( riter ) ) / (riter+1.);
						}
					   
			              }
				      }
			      }
	             
                }
         }


}

//        PetscGetTime(&te);  // xiaolei

//	cputime = te-ts;
//        PetscPrintf(PETSC_COMM_WORLD, "Finalizing_Projecting 2 cputime %d %le\n", ti,cputime);


//        PetscGetTime(&ts);  // xiaolei

// finalizing the projection of velocity vector and bed shear stress on body  surface 
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
		nfx = ibm->nf_x[ elmt ];
		nfy = ibm->nf_y[ elmt ];
		nfz = ibm->nf_z[ elmt ];

		ucx = ibm->Bvel[ elmt ].x;
		ucy = ibm->Bvel[ elmt ].y;
		ucz = ibm->Bvel[ elmt ].z;
                        
		if(ibm->Rigidity[elmt] == 1)
		{
			ibm->Bvel[ elmt ].x = 0.;
			ibm->Bvel[ elmt ].y = 0.;
			ibm->Bvel[ elmt ].z = 0.;
                        ibm->Shvel[ elmt ] = 0.;
                        if(effective_bed_shear_stress)ibm->Shvel_v[ elmt ] = 0.;
			if(effective_bed_shear_stress)ibm->k_ave[elmt] = 0.;
			if(effective_bed_shear_stress)ibm->w_ave[elmt] = 0.;
       		}
		else   // U_par = U - (U.n)n => (U_par)i = Ui - (Uj.nj).ni  : for every i we have j=1,2,3
		{
			ibm->Bvel[ elmt ].x = ucx - (ucx * nfx + ucy * nfy + ucz * nfz) * nfx;
	        	ibm->Bvel[ elmt ].y = ucy - (ucx * nfx + ucy * nfy + ucz * nfz) * nfy;
			//ibm->Bvel[ elmt ].z = ucz - (ucx * nfx + ucy * nfy + ucz * nfz) * nfz;
			//ibm->Bvel[ elmt ].z = ucz;
			ibm->Bvel[ elmt ].z =  - (ucx * nfx + ucy * nfy + ucz * nfz) * nfz;
			if(y_direction){
				        ibm->Bvel[ elmt ].y = - (ucx * nfx + ucy * nfy + ucz * nfz) * nfy;
			                ibm->Bvel[ elmt ].z =  ucz - (ucx * nfx + ucy * nfy + ucz * nfz) * nfz;
			               }
            if(effective_bed_shear_stress) {
            if(!y_direction)w_v_ = -(ucx * nfx / nfz + ucy * nfy / nfz + ibm->w_ave[elmt])/(1.e-8 + sqrt(1.0 + nfx * nfx /(1.e-8 + nfz * nfz) + nfy * nfy /(1.e-8 + nfz * nfz)));
            if(y_direction)w_v_ = -(ucx * nfx / nfy + ucz * nfz / nfy + ibm->w_ave[elmt])/(1.e-8 + sqrt(1.0 + nfx * nfx /(1.e-8 + nfy * nfy) + nfz * nfz /(1.e-8 + nfy * nfy)));
            if(!y_direction)ibm->w_ave [elmt] = fabs(w_v_/(1.e-8 + sqrt(ucx * ucx + ucy * ucy + w_v_ * w_v_)));
            if(y_direction)ibm->w_ave [elmt] = fabs(w_v_/(1.e-8 + sqrt(ucx * ucx + ucz * ucz + w_v_ * w_v_)));
			ibm->Shvel_v[elmt] = Cs_ * 0.456 * sqrt(fabs(ibm->w_ave[elmt] * ibm->k_ave[elmt]));
                                           }
		}
	}
	
if(effective_bed_shear_stress) Smoothing_Shear_vel( user, ibm, ti, tistart );

//  PetscGetTime(&te);  // xiaolei
//	cputime = te-ts;
//  PetscPrintf(PETSC_COMM_WORLD, "Finalizing_Projecting 3 cputime %d %le\n", ti,cputime);


return ( 0 );
}




PetscErrorCode sediment_variables_projection( UserCtx * user, IBMNodes * ibm )
{
	DA		 da  = user->da, fda  = user->fda;
	DALocalInfo      info  = user->info;
	PetscInt	 xs  = info.xs, xe  = info.xs + info.xm;
	PetscInt	 ys  = info.ys, ye  = info.ys + info.ym;
	PetscInt	 zs  = info.zs, ze  = info.zs + info.zm;
	PetscInt	 mx  = info.mx, my  = info.my, mz  = info.mz;
	PetscInt	 lxs, lys, lzs, lxe, lye, lze;
	lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

	if( xs == 0 )lxs = xs + 1;
	if( ys == 0 )lys = ys + 1;
	if( zs == 0 )lzs = zs + 1;

	if( xe == mx )lxe = xe - 1;
	if( ye == my )lye = ye - 1;
	if( ze == mz )lze = ze - 1;

	PetscInt  ibi  = 0; 
       // IBMNodes *ibm = &ibm_ptr[ibi];
	Cmpnts	   *** ucat;
	Cmpnts	   *   Bvel;
	PetscReal  *** ustar;
	PetscReal  *   Shvel;
	PetscInt	   n_elmt  = ibm->n_elmt, 
		           n_v  = ibm->n_v;
	PetscInt	   i, j, k;
        PetscInt  smoother_iterator;
        PetscInt  smoothing = 1;
	PetscReal	   sb, sc;
	PetscInt	   ip1,	ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
	PetscReal	   sk1,	sk2, sk3, cv1, cv2, cv3;
	PetscReal  *** nvert, *** nvert_o;
	PetscReal	cs1, cs2, cs3;
	PetscInt	nii;
	PetscReal	ucx, ucy, ucz;
	IBMInfo     * ibminfo;
	PetscInt	elmt,nelmt;
	PetscReal	n_of_cell[ n_elmt ], 
			n_of_cell_Sum[ n_elmt ], 
			riter;
	PetscReal	Bvel_x_Sum[ n_elmt ], 
			Bvel_y_Sum[ n_elmt ], 
			Bvel_z_Sum[ n_elmt ], 
			Shvel_Sum[ n_elmt],
                        Shvel_old[n_elmt];
        PetscReal	nfx, 
			nfy, 
			nfz;
        PetscReal       Bvel_x_o[n_elmt], Bvel_y_o[n_elmt], Bvel_z_o[n_elmt],Bed_shear_st_o[n_elmt];
//	IBMListNode  *	current;

//	PetscPrintf( PETSC_COMM_WORLD, "intp 001\n" );


	DAVecGetArray( fda, user->lUcat, &ucat );
	DAVecGetArray( da, user->lUstar, &ustar );

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
		Bvel_x_Sum[ elmt ]	  = 0.;
		Bvel_y_Sum[ elmt ]	  = 0.;
		Bvel_z_Sum[ elmt ]	  = 0.;
		n_of_cell[ elmt ]	  = 0.;
		n_of_cell_Sum[ elmt ] = 0.;

		ibm->Bvel[ elmt ].x	  = 0.;
		ibm->Bvel[ elmt ].y	  = 0.;
		ibm->Bvel[ elmt ].z	  = 0.;
		ibm->Shvel[ elmt ]	  = 0.;  
                if(ti==tistart)
                {
                Bed_shear_st_o[elmt] = 0.0;
                Bvel_x_o[elmt] = 0.0;
                Bvel_y_o[elmt] = 0.0;
                Bvel_z_o[elmt] = 0.0;
                }
	}																				    
					    
IBMListNode *current = user->ibmlist[ ibi ].head;

while( current )
{
IBMInfo  *ibminfo	= &current->ibm_intp;
	current = current->next;
	i		= ibminfo->ni;
	j		= ibminfo->nj;
	k		= ibminfo->nk;    
	sb		= ibminfo->d_s;
	sc		= sb + ibminfo->d_i;
        nelmt           = ibminfo->cell;		
        nfx  	= ibm->nf_x[ nelmt ];
        nfy  	= ibm->nf_y[ nelmt ];
	nfz  	= ibm->nf_z[ nelmt ];

	if(ibm->Rigidity[nelmt] == 1)
	{
              n_of_cell[nelmt] = 1.;		
              ibm->Shvel[ nelmt ] = 0.;
              ibm->Bvel[ nelmt ].x = 0.;
              ibm->Bvel[nelmt].y = 0.;
              ibm->Bvel[nelmt].z = 0.;
        }
	else
	{
	      if( nelmt >= 0  && i >= lxs && i < lxe && j >= lys && j < lye && k >= lzs && k < lze )
	     	{
	          n_of_cell[ nelmt ]   += 1.;
		  ibm->Bvel[ nelmt ].x += ucat[ k ][ j ][ i ].x;
		  ibm->Bvel[ nelmt ].y += ucat[ k ][ j ][ i ].y;
		  ibm->Bvel[ nelmt ].z += ucat[ k ][ j ][ i ].z;
		  ibm->Shvel[ nelmt ]   += ustar[ k ][ j ][ i ];
        	}
	}
}	
						    

for( elmt = 0; elmt < n_elmt; elmt++ )
{
   PetscGlobalSum( &n_of_cell[ elmt ], &n_of_cell_Sum[ elmt ], PETSC_COMM_WORLD );
   PetscGlobalSum( &(ibm->Bvel[ elmt ].x), &Bvel_x_Sum[ elmt ], PETSC_COMM_WORLD );
   PetscGlobalSum( &(ibm->Bvel[ elmt ].y), &Bvel_y_Sum[ elmt ], PETSC_COMM_WORLD );
   PetscGlobalSum( &(ibm->Bvel[ elmt ].z), &Bvel_z_Sum[ elmt ], PETSC_COMM_WORLD );
   PetscGlobalSum( &(ibm->Shvel[ elmt ]), &Shvel_Sum[ elmt ], PETSC_COMM_WORLD );
//   PetscPrintf(PETSC_COMM_WORLD, "celno and Glcelno %d %d\n",n_of_cell[elmt],n_of_cell_Sum[elmt]);
}

  PetscBarrier( PETSC_NULL );  								//	stop till all procc to the jobs  			    


	
for( elmt = 0; elmt < n_elmt; elmt++ )
	{

		if( n_of_cell_Sum[ elmt ] < 1.e-6 )
		{
		}
		else
		{
			ibm->Shvel[ elmt ] = Shvel_Sum[ elmt ] / n_of_cell_Sum[ elmt ];
			ibm->Bvel[ elmt ].x = Bvel_x_Sum[ elmt ] / n_of_cell_Sum[ elmt ];
			ibm->Bvel[ elmt ].y = Bvel_y_Sum[ elmt ] / n_of_cell_Sum[ elmt ];
			ibm->Bvel[ elmt ].z = Bvel_z_Sum[ elmt ] / n_of_cell_Sum[ elmt ];
            
		}

	}


///******************************************************
///*****************************************************
/*

	// print out for check and debuging  6 nov 2009
	//---------------------------------------
PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  FILE            *f;
  char            filen[80];
  if(rank==9){
//  if (!rank) {
  //  PetscTruth dyn=PETSC_FALSE,Hg=PETSC_FALSE,OSI=PETSC_FALSE;
  //  sprintf(filen, "stress%3.3d.dat",ti);
    sprintf(filen, "new_bed.dat");
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_SELF, f, "Variables=x,y,z, delta_z\n");
    PetscFPrintf(PETSC_COMM_SELF, f, "ZONE T='TRIANGLES', N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE VARLOCATION=([4]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);
  //  PetscPrintf(PETSC_COMM_WORLD, "new_bed level\n");
    // x - component
    for (i=0; i<ibm->n_v; i++) {
      
      PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->x_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_SELF, f, "\n");
    // y - component
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->y_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_SELF, f, "\n");
    // z - component
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->z_bp[i]);
    }
    
    PetscFPrintf(PETSC_COMM_SELF, f, "\n");
    
    
    //Write out the Delta_z
    for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%le\n",sqrt(ibm->Bvel[i].x*ibm->Bvel[i].x+ibm->Bvel[i].y*ibm->Bvel[i].y+ibm->Bvel[i].z*ibm->Bvel[i].z));
      }
    
    PetscFPrintf(PETSC_COMM_SELF, f, "\n");

    
    //Write out the link nodes
    for (i=0; i<ibm->n_elmt; i++) 
      {
	
	PetscFPrintf(PETSC_COMM_SELF, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      } fclose(f);
  }  


*/

///***************************************************
//***************************************************
     //   for (smoother_iterator=0;smoother_iterator<1;smoother_iterator++){
	for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )
	{
		//if( ibm->nf_z[ nelmt ] < 1.e-6 )
		if(ibm->Rigidity[nelmt] == 1)
		{
			ibm->Shvel[ nelmt ] = 0.;
			ibm->Bvel[ nelmt ].x = 0.;
			ibm->Bvel[ nelmt ].y = 0.;
			ibm->Bvel[ nelmt ].z = 0.;
		}
		else
		{
		//	if( fabs( Shvel_old[ nelmt ] ) < 1.e-6 ) //|| smoothing == 1)
			if( fabs( ibm->Shvel[ nelmt ] ) < 1.e-6 )//|| smoothing == 1)
		//	if( smoothing == 1)
			{  //for(i=0;i<2;i++){
				riter = 0.;
				for( nii = 0; nii < ibm->n_elmt; nii++ )
				{
					if(ibm->Rigidity[nii] == 1)
                                         { 
                                         }
                                        else
                                         {
                                           if( nii == nelmt || ibm->Shvel[nii]<1.e-6 )
                                           //if( ibm->Shvel[nii]<1.e-6 )
                                           {
                                           }
                                           else
                                           {
					
				//--------- interpolation between around elment (max 3 around elemts)***************
/*
	if(((ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ])
          &&(ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ] == ibm->nv3[ nii ]))

	|| ((ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ])
          &&(ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ]))


	|| ((ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ] == ibm->nv3[ nii ])
         &&(ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])))

*/
                    //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
     	if( ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ]== ibm->nv3[ nii ] || ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ]== ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ]   )  

                         			{
riter = riter + 1.;
ibm->Shvel[ nelmt ]	= ( ibm->Shvel[ nii ] + ibm->Shvel[ nelmt ] * ( riter - 1. ) ) / riter;
ibm->Bvel[ nelmt].x = ( ibm->Bvel[ nii ].x + ibm->Bvel[ nelmt ].x * ( riter - 1. ) ) / riter;
ibm->Bvel[ nelmt].y = ( ibm->Bvel[ nii ].y + ibm->Bvel[ nelmt ].y * ( riter - 1. ) ) / riter;
ibm->Bvel[ nelmt].z = ( ibm->Bvel[ nii ].z + ibm->Bvel[ nelmt ].z * ( riter - 1. ) ) / riter;
						}
					}
				}
			}
		//}
    	  }

    }

}
	DAVecRestoreArray( fda, user->lUcat, &ucat );
	DAVecRestoreArray( da, user->lUstar, &ustar );

// computing the projection of velocity vector on bed surface 
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
		nfx = ibm->nf_x[ elmt ];
		nfy = ibm->nf_y[ elmt ];
		nfz = ibm->nf_z[ elmt ];

		ucx = ibm->Bvel[ elmt ].x;
		ucy = ibm->Bvel[ elmt ].y;
		ucz = ibm->Bvel[ elmt ].z;

	//	if( nfz < 1.e-6 )
		if(ibm->Rigidity[elmt] == 1)
		{
			ibm->Bvel[ elmt ].x = 0.;
			ibm->Bvel[ elmt ].y = 0.;
			ibm->Bvel[ elmt ].z = 0.;
       		}
		else
		{
			//ibm->Bvel[ elmt ].x = ucx - ucx * nfx;// + ucz - ucz * nfz;
	        	//ibm->Bvel[ elmt ].y = ucy - ucy * nfy;// + ucz - ucz * nfz;
			//ibm->Bvel[ elmt ].z = ucz - ucz * nfz;
			ibm->Bvel[ elmt ].x = ucx - (ucx * nfx + ucy * nfy + ucz * nfz) * nfx;
	        	ibm->Bvel[ elmt ].y = ucy - (ucx * nfx + ucy * nfy + ucz * nfz) * nfy;
			ibm->Bvel[ elmt ].z = ucz - (ucx * nfx + ucy * nfy + ucz * nfz) * nfz;
		}
	}

return ( 0 );
}        
      

double integrating (double a , double b, double rp)
{
PetscInt i;
double x,N=100;
double delx=(b-a)/(2.*N);
    double sumeven=0.;
    for(i=1;i<N; i++)
       {
        x=a+delx*2.*i;
        sumeven += pow(x,1.5) * exp(-0.5*(x-rp)*(x-rp));
       }

    double sumodd =0.;
    for(i=1;i<N+1; i++)
       {
        x=a+delx*(2.*i-1.);
        sumodd += pow(x,1.5) * exp(-0.5*(x-rp)*(x-rp));
       }
 
return delx*(pow(a,1.5) * exp(-0.5*(a-rp)*(a-rp)) + pow(b,1.5) * exp(-0.5*(b-rp)*(b-rp))+2.*sumeven+4.*sumodd)/3.;
}




PetscErrorCode BCShS_Smoothening(IBMNodes *ibm, PetscInt ti, PetscInt tistart)
{
PetscInt i,nii, nelmt;
double riter;

PetscInt sm_number = 1;
if(ti==tistart || ti==0) PetscPrintf(PETSC_COMM_WORLD, "Number of Smoothening iterations for critical bed shear stress: %d \n",sm_number);

for(i=0;i<sm_number;i++){
	for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )
	{
		if(ibm->Rigidity [nelmt] == 1)
		{
			ibm->BCShS[ nelmt ] = 100.;
		}
		else
		{
				riter = 0.;
				for(nii = 0; nii < ibm->n_elmt; nii++)
				{
					if(ibm->Rigidity[nii] == 1 || nii==nelmt) 
                                         { 
                                         }
                                        else
                                         {
		                                     

                    //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
     	if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ]== ibm->nv3[ nii ] || ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ]== ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])  

                         			{
riter = riter + 1.;
ibm->BCShS[ nelmt ] = ( ibm->BCShS[ nii ] + ibm->BCShS[ nelmt ] * ( riter ) ) / (riter+1.);
						}
					}
			      }
	             }
         }
}

return(0);
}


PetscErrorCode	bed_concentration( UserCtx * user, IBMNodes * ibm, PetscInt ti, PetscInt tistart )
{
	DA	 da = user->da, 
		 fda = user->fda;
	Cmpnts		 *	Bvel;
	PetscReal	 *	BShS, *BShS_v, *Shvel, *Shvel_v, 
			 *	BCShS;
	PetscReal	 *	SCont;
	PetscReal	 *	w_ave;
	PetscInt		n_elmt	= ibm->n_elmt, 
	         		n_v  = ibm->n_v;
	PetscInt		i, j, k;
	PetscReal		sb,sc;
	PetscInt		ni,nii;
	PetscReal		ucx,ucy,ucz;
	IBMInfo  	 *	ibminfo;
	IBMListNode      *	current;
	PetscInt		elmt;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz;
	PetscReal		xc,yc, zc;
	PetscReal		tmp, tt;
      //PetscReal               Angle_repose	 = 45.; //45.;
	//PetscReal               D50		 = 0.0002;///jhoke-crossvane-rockvane 0.0011 bends 0.0005
	PetscReal               D50		 = d50;///jhoke-crossvane-rockvane 0.0011 bends 0.0005
	PetscReal               D90		 = 0.001;
	//PetscReal               particle_prosity    = 0.45;
	
	//ibm->Deltab = 10.* D50;//0.04*0.3;//0.04*0.152;//0.01*0.152;//0.01;//0.04;// 2.* D50 ;//2.0 * D50 * 7./2.5 ;
	ibm->Deltab = deltab;//0.04*0.3;//0.04*0.152;//0.01*0.152;//0.01;//0.04;// 2.* D50 ;//2.0 * D50 * 7./2.5 ;
	PetscReal       fdepth	= FlowDepth,//0.3,//1.15, // jhoke 0.172, //cross-vane 0.179, // 0.312, //0.4											    
	                particle_sg  = 2.65,//1.600,											    
	        	particle_dens  = 1922.0,//1600., tempor changed from 2650.										    
			fluid_density  = 1000.0,						  			    
			particle_ws,											    
			Shields_param,														    
			CShields_param,  												    
			Dstar,														    
			visk  = 1.e-6,											    
			gi	= 9.806,									  			    
			sgd_tmp,				  										    
			Ctau_0,  												    
			Custar_0,												    
			//U_bulk_vel	= 1.5,//1.0,//0.501,//1.25,//former run for vendetti case: 1.5,//0.7, // jhoke  0.355,//cross-vane 0.33, //0.285*1.5, //times 1.						    
			U_bulk_vel	= U_Bulk,//1.0,//0.501,//1.25,//former run for vendetti case: 1.5,//0.7, // jhoke  0.355,//cross-vane 0.33, //0.285*1.5, //times 1.						    
			max_bed_slope,										    
			angle_of_dev,
                        ttau, ttau_v;
       PetscReal        bita=PI*0.51/6.,DD;
       PetscInt         stochastic = 0, Fredsoe=1, VanRijn=0;
       if(stochastic){Fredsoe=0; VanRijn=1;} 
                  						    
      if(ti==tistart || ti==0){ PetscPrintf(PETSC_COMM_WORLD, "Stochastic: %d \n",stochastic);
      PetscPrintf(PETSC_COMM_WORLD, "Fredsoe: %d \n",Fredsoe);
      PetscPrintf(PETSC_COMM_WORLD, "Barchans: %d \n",Barchans);
      if(Barchans)PetscPrintf(PETSC_COMM_WORLD, "sediment_thickness: %e \n",sediment_thickness);
      PetscPrintf(PETSC_COMM_WORLD, "VanRijn: %d \n",VanRijn);
      PetscPrintf(PETSC_COMM_WORLD, "D50: %e \n",D50);
      PetscPrintf(PETSC_COMM_WORLD, "Delta_b: %e \n",ibm->Deltab);
      PetscPrintf(PETSC_COMM_WORLD, "Flow_Depth (m): %e \n",fdepth);
      PetscPrintf(PETSC_COMM_WORLD, "Sediment Particles Density (kg/m3): %e \n",particle_dens);
      PetscPrintf(PETSC_COMM_WORLD, "Flow_Bulk_Velocity (m/sec): %e \n",U_bulk_vel);}
//PetscReal   z_minimum,z_maximum,depth_elmt;
//---------------------------------------------------------------------------
// computing critical bed shear stress on a flat bed, by shields curve critera 
	PetscReal    angle_repose = Angle_repose * PI / 180.;
	//Angle_repose	  /= 180.;

	sgd_tmp   = ( particle_sg - 1. ) * gi * D50;
	particle_ws
	= sqrt( sgd_tmp )
	   * ( sqrt( 2. / 3. + 36. * visk * visk / ( sgd_tmp * D50 * D50 ) )
			- 6. * visk / ( sqrt( sgd_tmp ) * D50 ) );
	Dstar = D50 * pow(((particle_sg-1.) * gi / ( visk * visk )), 0.3333333333 );


      	PetscPrintf(PETSC_COMM_WORLD, "#### Dstar: %le \n",Dstar);

        int shields = 0, soulsby = 1, old_method = 0, whitehouse = 1; // soulseby and Whitehouse methods: see J Geoph. Res. vol 115, 2010, Chou & Fringer

        if(shields){

	if( Dstar <= 4. )
		CShields_param = 0.24 / Dstar;
	if( Dstar > 4. && Dstar <= 10. )
		CShields_param = 0.14 / pow( Dstar, 0.64 );
	if( Dstar > 10. && Dstar <= 20. )
		CShields_param = 0.04 / pow( Dstar, 0.1 );
	if( Dstar > 20. && Dstar <= 150. )
		CShields_param = 0.013 * pow( Dstar, 0.29 );
	if( Dstar > 150. )
		CShields_param = 0.055;
        } 
        if (soulsby) {
        CShields_param = 0.3/ ( 1. + 1.2 * Dstar) + 0.055 * (1. - exp (-0.02 * Dstar)) ;
        }
        
	if(VanRijn)
        {
         Ctau_0 = CShields_param * sgd_tmp * fluid_density;
	 Custar_0	 = sqrt( Ctau_0 / fluid_density );// change Nov 17 2009---canceled nov 18 2009
        }
      
	if(Fredsoe)
        {
         Ctau_0 = CShields_param; 
	}
	
        // non-dimensionalizing
	//Custar_0	/= U_bulk_vel;
	//Ctau_0		 = Custar_0 * Custar_0;
	particle_ws /= U_bulk_vel;
	D50  		/= fdepth;
	D90  		/= fdepth;
	ibm->Deltab		/= fdepth;
	
// computing critical bed shear stress on inclined bed
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
		//xc	= ibm->cent_x[ elmt ];
		//yc	= ibm->cent_y[ elmt ];
		//zc	= ibm->cent_z[ elmt ];

		nfx = ibm->nf_x[ elmt ];
		nfy = ibm->nf_y[ elmt ];
		nfz = ibm->nf_z[ elmt ];

		ucx = ibm->Bvel[ elmt ].x;
		ucy = ibm->Bvel[ elmt ].y;
		ucz = ibm->Bvel[ elmt ].z;

		if(ibm->Rigidity[elmt])
		{
			ibm->BCShS[ elmt ] = 100.0;
		}
		else
		{
			if (old_method) {
                        max_bed_slope =  atan(sqrt((ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8)) + (ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))));
			//= sqrt( ( nfy / nfz ) * ( nfy / nfz ) + ( nfx / nfz ) * ( nfx / nfz ) );

			if(fabs(max_bed_slope)<1.e-6)
                        {
                        tmp=1.0;
                        }
                        else
                        {
                        //max_bed_slope = atan( max_bed_slope );

			angle_of_dev  = ( nfx / (nfz) ) * ucx + ( nfy / (nfz) ) * ucy;

                          //PetscPrintf(PETSC_COMM_WORLD, "angle of dev %le \n",angle_of_dev);
                        angle_of_dev /= ( sqrt( ( nfy / (nfz) ) * ( nfy / (nfz) ) + ( nfx / (nfz) ) * ( nfx / (nfz) )));

			angle_of_dev /= ( sqrt( ucx * ucx + ucy * ucy ));
                        
                        if(fabs(angle_of_dev)>1.)angle_of_dev=angle_of_dev*PetscMax(.99,angle_of_dev)/fabs(angle_of_dev);
			
                        angle_of_dev = acos( angle_of_dev );

			tmp = sin( angle_of_dev ) * sin( angle_of_dev ) * tan( max_bed_slope )
				   * tan( max_bed_slope );
			tmp /= ( tan( angle_repose ) * tan( angle_repose ) );
			tmp = PetscMax((1. - tmp),0.);
			tmp = cos( max_bed_slope ) * sqrt( tmp );
			tmp = tmp - cos( angle_of_dev ) * sin( max_bed_slope ) / tan( angle_repose );
                        }
                 	
                        //ibm->BCShS[ elmt ] = Ctau_0 *tmp*nfz*0.75/tmp; //times a 0.75 to reduce the critical bed shear stress
                 	//ibm->BCShS[ elmt ] = Ctau_0 *tmp; //times a 0.75 to reduce the critical bed shear stress
                 	ibm->BCShS[ elmt ] = Ctau_0; //times a 0.75 to reduce the critical bed shear stress
                // 	ibm->BCShS[ elmt ] = 0.03; //times a 0.75 to reduce the critical bed shear stress
                // 	DD=Ctau_0/tan(angle_repose);
                         }
 
                         if (whitehouse){
                         double tita_x, tita_y, tita_z, tita, tita1, tita2, tita3;
                         if (!y_direction) {
											tita_x = atan ( - nfx / (nfz+1.e-8) );
											tita_y = atan ( - nfy / (nfz+1.e-8) );
											tita1  = ucx * sin (tita_x) + ucy * sin (tita_y); 
											tita2  = sqrt (ucx * ucx + ucy * ucy);
						                   }
										    else {
													tita_x = atan ( - nfx / (nfy+1.e-8) );
													tita_z = atan ( - nfz / (nfy+1.e-8) );
													tita1  = ucx * sin (tita_x) + ucz * sin (tita_z); 
													tita2  = sqrt (ucx * ucx + ucz * ucz);
											 }
						 tita3  = tita1 / tita2;
                         if (fabs(tita3) >= 1. || fabs(tita3)<= -1.) {tmp =1.;} else {
                         tita   = asin (tita3);
                         tmp           = sin ( tita + angle_repose ) / sin (angle_repose);}
                         
                         ibm->BCShS [ elmt ] = Ctau_0 * tmp;
                                       }
		}
	}

// Smoothening the computed critical bed shear stress over entire domian
//BCShS_Smoothening(ibm, ti, tistart);

        PetscReal SCont_max = -10.;
// computing bed concentration by Van Rijn relation for equilibrium condition // Fredsoe method added as well
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
		//xc	= ibm->cent_x[ elmt ];
		//yc	= ibm->cent_y[ elmt ];
		//zc	= ibm->cent_z[ elmt ];
                 
                int rigid = ibm->Rigidity[elmt];
                 
		nfx = ibm->nf_x[ elmt ];
		nfy = ibm->nf_y[ elmt ];
		nfz = ibm->nf_z[ elmt ];

	//	if( nfz < 1.e-6 )
		if(ibm->Rigidity[elmt] == 1)
		{
		        ibm->SCont[ elmt ] = 0.0;
                        ibm->BShS[elmt] = 0.0;
                        ibm->BCShS[elmt] = 100.0;
                       if (effective_bed_shear_stress && rans) ibm->Shvel_v[elmt]=0.0;
                       if (effective_bed_shear_stress && rans) ibm->BShS_v[elmt]=0.0;
		}
		else
		{

                //Calc efective shear velocity
                       if (effective_bed_shear_stress && rans)
                         {
                         //double pz_px =  - nfx / nfz ;
                         //double pz_py = - nfy / nfz ;
                         //double w_p = (ucx * pz_px + ucy * pz_py - ucz * 2.0)/sqrt(pz_px * pz_px + pz_py * pz_py + 1.0);
                         //double w_p = (ucx * pz_px + ucy * pz_py - ibm->w_ave[elmt])/(1.e-8 + sqrt(pz_px * pz_px + pz_py * pz_py + 1.0));
                         //double R_p = fabs(w_p /(1.e-8 + sqrt(ucx * ucx + ucy * ucy + w_p * w_p)));
                         //PetscReal Csyyy 1.5; //varies between 0.8 to 1.5
                         //ibm->Shvel_v[elmt] = Cs_ * 0.456 * sqrt(fabs(ibm->w_ave[elmt] * ibm->k_ave[elmt]));
                         ibm->Shvel[elmt] = ibm->Shvel[elmt] + ibm->Shvel_v[elmt];
                         }
               //end of effective bed shear stress.  See Jia, Yafei, Altinkar and others, JHR, 2015-2016
                             
                if(VanRijn){
                        //ttau=sqrt(ibm->BShS[elmt])*U_bulk_vel*1.; //inlet vel. increase by 1.5 times
                        //ttau=ttau*ttau*fluid_density;    // dimensionalizing bed shear stress
                        //ibm->BShS[elmt]=ttau;
                        ibm->BShS[elmt] = fluid_density*(ibm->Shvel[elmt]*U_bulk_vel)*(ibm->Shvel[elmt]*U_bulk_vel);//+ DD*(nfx*nfz+nfy*nfz+nfz*nfz-1.);
                        

                                  if(!stochastic)
                                              {

			tt = (ibm->BShS[ elmt ] - ibm->BCShS[ elmt ])/ibm->BCShS[elmt];

			tt = PetscMax( tt, 0. );
                                       
                       //*************************** 1. first strategy ********************************************************************
                         ibm->SCont[ elmt ] =  0.015 *(D50/ibm->Deltab)*pow( tt, 1.5 )/pow( Dstar,0.3);
                                                                     
                       //**************************** 1-2. strategy ******************************************************************

		       //      ibm->SCont[ elmt ] = PetscMax(0.,(0.015 * (D50 /ibm->Deltab ) * pow( tt, 1.5 ) / pow( Dstar, 0.3 )+
                       //                            ibm->netflux_old[elmt]*user->dt/(ibm->dA[elmt]*ibm->Deltab)));

                       //**************************** 2. second strategy ******************************************************************
                       // SCont is the portable part of bed cell concent. Depending on the last time step deposition or scour
                       // the current time step cell concentration can be more or less than the Ca by Van Rijn formula here 
                       // but if tt>0 it gonna be always greater than zero:
                       /* if (tt>0.)
                            {
                       
		             ibm->SCont[ elmt ] = PetscMax(0.,(0.015 * (D50 /ibm->Deltab ) * pow( tt, 1.5 ) / pow( Dstar, 0.3 )+
                                                   ibm->netflux_old[elmt]*user->dt/(ibm->dA[elmt]*ibm->Deltab)));
                                                //   ibm->netflux_old[elmt]*3.0/(ibm->dA[elmt]*ibm->Deltab)));
                            }
                        else
                            {

		             ibm->SCont[ elmt ] = 0.015 * (D50 /ibm->Deltab ) * pow( tt, 1.5 ) / pow( Dstar, 0.3 )+
                                                  PetscMax(0.,((ibm->netflux_old[elmt])*(user->dt)/(ibm->dA[elmt]*ibm->Deltab)));
                                                  //PetscMax(0.,((ibm->netflux_old[elmt])*(3.0)/(ibm->dA[elmt]*ibm->Deltab)));
                            }*/
                       //**************************** 3. Third one  ************************************************************************
                       // SCont is the portable part of bed cell sediemnt cont anb the portable sed-concent
                       // at a cell can not be more than the Ca which is the max. portable sed-concent by Van Rijn formula here 
                       // but because of deposition at last time step it could be less than Ca :
                       /* if (tt>0.)
                            {
                                  ibm->SCont[ elmt ] = 0.015 * (D50 /ibm->Deltab ) * pow( tt, 1.5 ) / pow( Dstar, 0.3 )+
                                                   PetsMin(0., (ibm->netflux_old[elmt]*user->dt/(ibm->dA[elmt]*ibm->Deltab)));
                            }
                        else
                            {
                             ibm->SCont[ elmt ] =  0.015 * (D50 /ibm->Deltab ) * pow( tt, 1.5 ) / pow( Dstar, 0.3 )+
                                                  PetscMax(0.,((ibm->netflux_old[elmt])*(user->dt)/(ibm->dA[elmt]*ibm->Deltab)));
                            }*/
                       //****************************************************************************************************************                                            
                                      }
                                      if(stochastic)
                                      {
                                        double Int1, Int2;
                                        //double sigma = 0.4*fabs(ibm->BShS[elmt]);
                                        double sigma = fabs(ibm->BShS[elmt]);
                                        //double Tau_Cr= 1.5*Ctau_0;
                                        //double Tau_Cr1= 1.5 * ibm->BCShS[elmt];
                                        //double Tau_Cr2= 1.5 * ibm->BCShS[elmt];     
                                        double Tau_Cr1= ibm->BCShS[elmt];
                                        double Tau_Cr2= ibm->BCShS[elmt];     
                                        double r = (ibm->BShS[elmt]-Tau_Cr1)/(sigma+1.e-8);
                                        double p = (-ibm->BShS[elmt]-Tau_Cr2)/(sigma+1.e-8);
                                        double a1 = PetscMax(0.0,r - 4.);
                                        double b1 = r + 5.;            
                                        double a2 = PetscMax(0.0,p - 4.);
                                        double b2 = p + 5.;       
                                        if(r <  -3.0)  Int1 = 0.0001;     
                                        if(r >= -3.0)  Int1 = integrating (a1,b1,r); 
                                        if(p <  -3.0)  Int2 = 0.0001;
                                        if(p >= -3.0)  Int2 = integrating (a2,b2,p);
                                        //Int1 = integrating (a1,b1,r); 
                                        //Int2 = integrating (a2,b2,p);
                                        
                                        ibm->SCont[elmt] =  100*0.03 * (D50/ibm->Deltab)*fabs(Int1+Int2)/pow( Dstar,0.3);
                                      }                                             
                         }
                  
        
                if(Fredsoe) // ali 10 March 2010
                         {
                          ttau=(ibm->Shvel[elmt]*U_bulk_vel)*(ibm->Shvel[elmt]*U_bulk_vel); // the inlet vel. increase by 1.5 times
                          ibm->BShS[elmt]=ttau/(sgd_tmp);
                            
                         if(effective_bed_shear_stress){ttau_v = (ibm->Shvel_v[elmt]*U_bulk_vel)*(ibm->Shvel_v[elmt]*U_bulk_vel); 
                          ibm->BShS_v[elmt] = ttau_v/(sgd_tmp);}
                                            
			  tt = ibm->BShS[ elmt ] - ibm->BCShS[ elmt ];

			  tt = PetscMax( tt, 1.e-7 );
                         
                          PetscReal bita=PI*0.51/6.;
                          
                          tt= tt*tt*tt*tt;
                          bita=bita*bita*bita*bita;
                          
                          ibm->SCont[elmt]=pow((1.+ bita/tt),-0.25)*PI/6.;
                
                          //Hitting the rigid bed; if a cell is on the rigid bed then the bed concentration is zero
                         }

		}
                //if(RigidBed && ibm->cent_z[elmt] <0.225) {ibm->SCont[elmt] = 0.0;}                                                                
		if(ibm->SCont[ elmt ] > porosity) ibm->SCont[ elmt ] = porosity;
		if(ibm->SCont[ elmt ] < 1.e-7) ibm->SCont[ elmt ] = 0.0;

        SCont_max = PetscMax(SCont_max,ibm->SCont[elmt]);
	}

    PetscPrintf(PETSC_COMM_WORLD, "Maximum bed sediment concentration: %e\n",SCont_max);
   //for( elmt = 0; elmt < n_elmt; elmt++ )
   //{
    //if((ibm->cent_z[elmt]< 0.2001 || fabs(ibm->cent_z[elmt]-0.2)<1.e-6)&&(ibm->cent_x[elmt] < 4.0 || ibm->cent_x[elmt] > 6.54 || ibm->cent_y[elmt] < 2.0 || ibm->cent_y[elmt] > 4.54)) ibm->SCont[elmt] = 0.0; 
    //if(ibm->cent_z[elmt]<= 0.01) ibm->SCont[elmt] = 0.0; 
    //if(ibm->cent_z[elmt]<= 0.1001) ibm->SCont[elmt] = 0.0;                                   
    //if(Barchans && ibm->cent_x[elmt] <= sediment_thickness) ibm->SCont[elmt] = 0.0;  
   //}

             
return ( 0 );		// Ali completed 2 nov. 2009		
}

//Hossein added from NSF Petsc3.1 - bed_concentration_levelset



PetscErrorCode calc_cell_grad_y_direction( UserCtx * user, IBMNodes * ibm )
{
DA			 da  = user->da, 
			 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
			 * F_flux_13, 
			 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	 n_elmt  = ibm->n_elmt, 
			 n_v  = ibm->n_v;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
	PetscReal	   xc, yc, zc;
	PetscReal	   tmp;
	PetscReal	  Ave_c, Ave_w, 
					Ave_u, 
					nfx12, 
					nfz12, 
					nfx13, 
					nfz13, 
					nfx23, 
					nfz23, 
					nor_vel_to_face12, 
					nor_vel_to_face13, 
					nor_vel_to_face23, 
					le12z, 
					le12x, 
					le13z, 
					le13x, 
					le23z, 
					le23x, 
					le_dot_n12, 
					le_dot_n13, 
					le_dot_n23, ndotIJ, IJz, IJx;
//**************** initiation
	for( elmt = 0; elmt < n_elmt; elmt++ )
	      {
           ibm->w_grad_z[elmt] = 0.;
	       ibm->u_grad_z[elmt] = 0.;
	       ibm->c_grad_z[elmt] = 0.;
	       ibm->w_grad_x[elmt] = 0.;
	       ibm->u_grad_x[elmt] = 0.;
	       ibm->c_grad_x[elmt] = 0.;
	       ibm->w_max[elmt] = ibm->Bvel[elmt].z;
	       ibm->u_max[elmt] = ibm->Bvel[elmt].x;
	       ibm->c_max[elmt] = ibm->SCont[elmt];
              }

// computing gradient for cell elmt with neighbor nbelmt through all three faces: 1-->2 , 1-->3 , 2-->3
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	   if( ibm->Rigidity[elmt]==1)
		{
		}
		else
		{
			n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
			
            xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
			zc	 = ibm->cent_z[ elmt ];

			dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and  
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/


				nbelmt=ibm->c2c[elmt].c1;  // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{

				    if( ibm->Rigidity[nbelmt] == 1)
					{
					}
					else
					{
						Ave_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );

						Ave_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

				        Ave_c = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );
				
						nfz12 = dx12;
						nfx12 = -dz12;
						//ds	  = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr = sqrt( dz12 * dz12 + dx12 * dx12 );
                                                ds = dr;
						nfz12 /= dr;
						nfx12 /= dr;                 
                                                IJz= ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                IJx= ibm->cent_x[nbelmt]- ibm->cent_x[elmt];

						le12z  = 0.5 * ( ibm->z_bp[ n2e ] + ibm->z_bp[ n1e ] ) - zc;
						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - xc;

                                                ndotIJ = IJz * nfz12 + IJx * nfx12;
                                                if(ndotIJ >0.)
                                                {
						ibm->w_grad_z[elmt] += Ave_w *  nfz12 * ds;
						ibm->u_grad_z[elmt] += Ave_u *  nfz12 * ds;
						ibm->c_grad_z[elmt] += Ave_c *  nfz12 * ds;
						ibm->w_grad_x[elmt] += Ave_w *  nfx12 * ds;
						ibm->u_grad_x[elmt] += Ave_u *  nfx12 * ds;
						ibm->c_grad_x[elmt] += Ave_c *  nfx12 * ds;
						} else{                                                
                        ibm->w_grad_z[elmt] -= Ave_w *  nfz12 * ds;
						ibm->u_grad_z[elmt] -= Ave_u *  nfz12 * ds;
						ibm->c_grad_z[elmt] -= Ave_c *  nfz12 * ds;
						ibm->w_grad_x[elmt] -= Ave_w *  nfx12 * ds;
						ibm->u_grad_x[elmt] -= Ave_u *  nfx12 * ds;
						ibm->c_grad_x[elmt] -= Ave_c *  nfx12 * ds;
                                                }

						if(fabs(ibm->w_max[elmt])<fabs(ibm->w_max[nbelmt]))ibm->w_max[elmt] = ibm->w_max[nbelmt];
						if(fabs(ibm->u_max[elmt])<fabs(ibm->u_max[nbelmt]))ibm->u_max[elmt] = ibm->u_max[nbelmt];
						if(fabs(ibm->c_max[elmt])<fabs(ibm->c_max[nbelmt]))ibm->c_max[elmt] = ibm->c_max[nbelmt];
					}
				}
			

				//					1---->3 edge flux 

				/*   // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/


				nbelmt=ibm->c2c[elmt].c2;  // xiaolei add SEDI
	
				if (nbelmt>=0)  // xiaolei add 
				{
	                              	if( ibm->Rigidity[nbelmt] == 1)
					{
					}
					else
					{

						Ave_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );

						Ave_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

                        Ave_c = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfz13 = dx13;
						nfx13 = -dz13;
						//ds	  = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr = sqrt( dz13 * dz13 + dx13 * dx13 );
                                                ds = dr;
						nfz13 /= dr;
						nfx13 /= dr;
                                                
                                                IJz= ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                IJx= ibm->cent_x[nbelmt]- ibm->cent_x[elmt];

						le13z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n1e ] ) - zc;
						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - xc;
						
                                                ndotIJ = IJz * nfz13 + IJx * nfx13;
                                                if(ndotIJ > 0.){
                        ibm->w_grad_z[elmt] += Ave_w *  nfz13 * ds;
						ibm->u_grad_z[elmt] += Ave_u *  nfz13 * ds;
						ibm->c_grad_z[elmt] += Ave_c *  nfz13 * ds;
						ibm->w_grad_x[elmt] += Ave_w *  nfx13 * ds;
						ibm->u_grad_x[elmt] += Ave_u *  nfx13 * ds;
						ibm->c_grad_x[elmt] += Ave_c *  nfx13 * ds;}
                                                else{
                        ibm->w_grad_z[elmt] -= Ave_w *  nfz13 * ds;
						ibm->u_grad_z[elmt] -= Ave_u *  nfz13 * ds;
						ibm->c_grad_z[elmt] -= Ave_c *  nfz13 * ds;
						ibm->w_grad_x[elmt] -= Ave_w *  nfx13 * ds;
						ibm->u_grad_x[elmt] -= Ave_u *  nfx13 * ds;
						ibm->c_grad_x[elmt] -= Ave_c *  nfx13 * ds;}

						if(fabs(ibm->w_max[elmt])<fabs(ibm->w_max[nbelmt]))ibm->w_max[elmt] = ibm->w_max[nbelmt];
						if(fabs(ibm->u_max[elmt])<fabs(ibm->u_max[nbelmt]))ibm->u_max[elmt] = ibm->u_max[nbelmt];
						if(fabs(ibm->c_max[elmt])<fabs(ibm->c_max[nbelmt]))ibm->c_max[elmt] = ibm->c_max[nbelmt];
					 }
				}
				

				//					2---->3 edge  flux

				/*  // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c3;  // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{
	                               	if( ibm->Rigidity[nbelmt] == 1)
					{
					}
					else
					{
						Ave_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );

						Ave_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

                        Ave_c = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfz23 = dx23;
						nfx23 = -dz23;
						//ds	  = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr = sqrt( dz23 * dz23 + dx23 * dx23 );
                                                ds = dr;
						nfz23 /= dr;
						nfx23 /= dr;
                                                
                                                IJz= ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                IJx= ibm->cent_x[nbelmt]- ibm->cent_x[elmt];

						le23z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n2e ] ) - zc;
						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - xc;

                                                ndotIJ = IJz * nfz23 + IJx * nfx23;
                                                if(ndotIJ > 0.){
                        ibm->w_grad_z[elmt] += Ave_w *  nfz23 * ds;
						ibm->u_grad_z[elmt] += Ave_u *  nfz23 * ds;
						ibm->c_grad_z[elmt] += Ave_c *  nfz23 * ds;
						ibm->w_grad_x[elmt] += Ave_w *  nfx23 * ds;
						ibm->u_grad_x[elmt] += Ave_u *  nfx23 * ds;
						ibm->c_grad_x[elmt] += Ave_c *  nfx23 * ds;}
                                                else{
						ibm->w_grad_z[elmt] -= Ave_w *  nfz23 * ds;
						ibm->u_grad_z[elmt] -= Ave_u *  nfz23 * ds;
						ibm->c_grad_z[elmt] -= Ave_c *  nfz23 * ds;
						ibm->w_grad_x[elmt] -= Ave_w *  nfx23 * ds;
						ibm->u_grad_x[elmt] -= Ave_u *  nfx23 * ds;
						ibm->c_grad_x[elmt] -= Ave_c *  nfx23 * ds;}
						if(fabs(ibm->w_max[elmt])<fabs(ibm->w_max[nbelmt]))ibm->w_max[elmt] = ibm->w_max[nbelmt];
						if(fabs(ibm->u_max[elmt])<fabs(ibm->u_max[nbelmt]))ibm->u_max[elmt] = ibm->u_max[nbelmt];
						if(fabs(ibm->c_max[elmt])<fabs(ibm->c_max[nbelmt]))ibm->c_max[elmt] = ibm->c_max[nbelmt];
				       }
				}
			} //end of neighnor cells loop
		} //end of outer empty cells if

           ibm->w_grad_z[elmt] /= ibm->dA[elmt];
	       ibm->u_grad_z[elmt] /= ibm->dA[elmt];
	       ibm->c_grad_z[elmt] /= ibm->dA[elmt];
	       ibm->w_grad_x[elmt] /= ibm->dA[elmt];
	       ibm->u_grad_x[elmt] /= ibm->dA[elmt];
	       ibm->c_grad_x[elmt] /= ibm->dA[elmt];

 //if(ibm->nf_z[elmt]>0. && ibm->elmt_depth[elmt]<0.1 && ibm->cent_z[elmt]>=.9999)PetscPrintf(PETSC_COMM_WORLD, "elmt /// c_grad_x /// c_grad_y: %d %le %le\n",elmt,ibm->c_grad_x[elmt],ibm->c_grad_y[elmt]);
	} //end of elmt cells loop

	return ( 0 );												//	ali completed on 3 nov. 2009				    
}



PetscErrorCode calc_cell_grad( UserCtx * user, IBMNodes * ibm )
{
DA			 da  = user->da, 
			 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
			 * F_flux_13, 
			 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	 n_elmt  = ibm->n_elmt, 
			 n_v  = ibm->n_v;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
	PetscReal	   xc, yc, zc;
	PetscReal	   tmp;
	PetscReal	  Ave_c, Ave_u, 
					Ave_v, 
					nfx12, 
					nfy12, 
					nfx13, 
					nfy13, 
					nfx23, 
					nfy23, 
					nor_vel_to_face12, 
					nor_vel_to_face13, 
					nor_vel_to_face23, 
					le12x, 
					le12y, 
					le13x, 
					le13y, 
					le23x, 
					le23y, 
					le_dot_n12, 
					le_dot_n13, 
					le_dot_n23, ndotIJ, IJx, IJy;
//**************** initiation
	for( elmt = 0; elmt < n_elmt; elmt++ )
	      {
           ibm->u_grad_x[elmt] = 0.;
	       ibm->v_grad_x[elmt] = 0.;
	       ibm->c_grad_x[elmt] = 0.;
	       ibm->u_grad_y[elmt] = 0.;
	       ibm->v_grad_y[elmt] = 0.;
	       ibm->c_grad_y[elmt] = 0.;
	       ibm->u_max[elmt] = ibm->Bvel[elmt].x;
	       ibm->v_max[elmt] = ibm->Bvel[elmt].y;
	       ibm->c_max[elmt] = ibm->SCont[elmt];
              }

// computing gradient for cell elmt with neighbor nbelmt through all three faces: 1-->2 , 1-->3 , 2-->3
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	   if( ibm->Rigidity[elmt] == 1)
		{
		}
		else
		{
			n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
			
            xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
			zc	 = ibm->cent_z[ elmt ];

			dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and  
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/


				nbelmt=ibm->c2c[elmt].c1;  // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{

				    if(ibm->Rigidity[nbelmt] == 1)
					{
					}
					else
					{
						Ave_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

				                Ave_c = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfx12 = dy12;
						nfy12 = -dx12;
						//ds	  = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr = sqrt( dx12 * dx12 + dy12 * dy12 );
                                                ds = dr;
						nfx12 /= dr;
						nfy12 /= dr;
                                                 
                                                IJx= ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                IJy= ibm->cent_y[nbelmt]- ibm->cent_y[elmt];

						le12x  = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - xc;
						le12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] ) - yc;

                                                ndotIJ = IJx * nfx12 + IJy * nfy12;
                                                if(ndotIJ >0.)
                                                {
						ibm->u_grad_x[elmt] += Ave_u *  nfx12 * ds;
						ibm->v_grad_x[elmt] += Ave_v *  nfx12 * ds;
						ibm->c_grad_x[elmt] += Ave_c *  nfx12 * ds;
						ibm->u_grad_y[elmt] += Ave_u *  nfy12 * ds;
						ibm->v_grad_y[elmt] += Ave_v *  nfy12 * ds;
						ibm->c_grad_y[elmt] += Ave_c *  nfy12 * ds;
						} else{                                                
                                                ibm->u_grad_x[elmt] -= Ave_u *  nfx12 * ds;
						ibm->v_grad_x[elmt] -= Ave_v *  nfx12 * ds;
						ibm->c_grad_x[elmt] -= Ave_c *  nfx12 * ds;
						ibm->u_grad_y[elmt] -= Ave_u *  nfy12 * ds;
						ibm->v_grad_y[elmt] -= Ave_v *  nfy12 * ds;
						ibm->c_grad_y[elmt] -= Ave_c *  nfy12 * ds;
                                                }

						if(fabs(ibm->u_max[elmt])<fabs(ibm->u_max[nbelmt]))ibm->u_max[elmt] = ibm->u_max[nbelmt];
						if(fabs(ibm->v_max[elmt])<fabs(ibm->v_max[nbelmt]))ibm->v_max[elmt] = ibm->v_max[nbelmt];
						if(fabs(ibm->c_max[elmt])<fabs(ibm->c_max[nbelmt]))ibm->c_max[elmt] = ibm->c_max[nbelmt];
					}
				}
			

				//					1---->3 edge flux 

				/*   // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/


				nbelmt=ibm->c2c[elmt].c2;  // xiaolei add SEDI
	
				if (nbelmt>=0)  // xiaolei add 
				{
	                              	if( ibm->Rigidity[nbelmt] == 1)
					{
					}
					else
					{

						Ave_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

                                                Ave_c = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfx13 = dy13;
						nfy13 = -dx13;
						//ds	  = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr = sqrt( dx13 * dx13 + dy13 * dy13 );
                                                ds = dr;
						nfx13 /= dr;
						nfy13 /= dr;
                                                
                                                IJx= ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                IJy= ibm->cent_y[nbelmt]- ibm->cent_y[elmt];

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - xc;
						le13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] ) - yc;
						
                                                ndotIJ = IJx * nfx13 + IJy * nfy13;
                                                if(ndotIJ > 0.){
                                                ibm->u_grad_x[elmt] += Ave_u *  nfx13 * ds;
						ibm->v_grad_x[elmt] += Ave_v *  nfx13 * ds;
						ibm->c_grad_x[elmt] += Ave_c *  nfx13 * ds;
						ibm->u_grad_y[elmt] += Ave_u *  nfy13 * ds;
						ibm->v_grad_y[elmt] += Ave_v *  nfy13 * ds;
						ibm->c_grad_y[elmt] += Ave_c *  nfy13 * ds;}
                                                else{
                                                ibm->u_grad_x[elmt] -= Ave_u *  nfx13 * ds;
						ibm->v_grad_x[elmt] -= Ave_v *  nfx13 * ds;
						ibm->c_grad_x[elmt] -= Ave_c *  nfx13 * ds;
						ibm->u_grad_y[elmt] -= Ave_u *  nfy13 * ds;
						ibm->v_grad_y[elmt] -= Ave_v *  nfy13 * ds;
						ibm->c_grad_y[elmt] -= Ave_c *  nfy13 * ds;}

						if(fabs(ibm->u_max[elmt])<fabs(ibm->u_max[nbelmt]))ibm->u_max[elmt] = ibm->u_max[nbelmt];
						if(fabs(ibm->v_max[elmt])<fabs(ibm->v_max[nbelmt]))ibm->v_max[elmt] = ibm->v_max[nbelmt];
						if(fabs(ibm->c_max[elmt])<fabs(ibm->c_max[nbelmt]))ibm->c_max[elmt] = ibm->c_max[nbelmt];
					 }
				}
				

				//					2---->3 edge  flux

				/*  // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c3;  // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{
	                               	if( ibm->Rigidity[nbelmt] == 1)
					{
					}
					else
					{
						Ave_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

                                                Ave_c = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfx23 = dy23;
						nfy23 = -dx23;
						//ds	  = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr = sqrt( dx23 * dx23 + dy23 * dy23 );
                                                ds = dr;
						nfx23 /= dr;
						nfy23 /= dr;
                                                
                                                IJx= ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                IJy= ibm->cent_y[nbelmt]- ibm->cent_y[elmt];

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - xc;
						le23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] ) - yc;

                                                ndotIJ = IJx * nfx23 + IJy * nfy23;
                                                if(ndotIJ > 0.){
                                                ibm->u_grad_x[elmt] += Ave_u *  nfx23 * ds;
						ibm->v_grad_x[elmt] += Ave_v *  nfx23 * ds;
						ibm->c_grad_x[elmt] += Ave_c *  nfx23 * ds;
						ibm->u_grad_y[elmt] += Ave_u *  nfy23 * ds;
						ibm->v_grad_y[elmt] += Ave_v *  nfy23 * ds;
						ibm->c_grad_y[elmt] += Ave_c *  nfy23 * ds;}
                                                else{
                                                ibm->u_grad_x[elmt] -= Ave_u *  nfx23 * ds;
						ibm->v_grad_x[elmt] -= Ave_v *  nfx23 * ds;
						ibm->c_grad_x[elmt] -= Ave_c *  nfx23 * ds;
						ibm->u_grad_y[elmt] -= Ave_u *  nfy23 * ds;
						ibm->v_grad_y[elmt] -= Ave_v *  nfy23 * ds;
						ibm->c_grad_y[elmt] -= Ave_c *  nfy23 * ds;}
						if(fabs(ibm->u_max[elmt])<fabs(ibm->u_max[nbelmt]))ibm->u_max[elmt] = ibm->u_max[nbelmt];
						if(fabs(ibm->v_max[elmt])<fabs(ibm->v_max[nbelmt]))ibm->v_max[elmt] = ibm->v_max[nbelmt];
						if(fabs(ibm->c_max[elmt])<fabs(ibm->c_max[nbelmt]))ibm->c_max[elmt] = ibm->c_max[nbelmt];
				       }
				}
			} //end of neighnor cells loop
		} //end of outer empty cells if

               ibm->u_grad_x[elmt] /= ibm->dA[elmt];
	       ibm->v_grad_x[elmt] /= ibm->dA[elmt];
	       ibm->c_grad_x[elmt] /= ibm->dA[elmt];
	       ibm->u_grad_y[elmt] /= ibm->dA[elmt];
	       ibm->v_grad_y[elmt] /= ibm->dA[elmt];
	       ibm->c_grad_y[elmt] /= ibm->dA[elmt];

 //if(ibm->nf_z[elmt]>0. && ibm->elmt_depth[elmt]<0.1 && ibm->cent_z[elmt]>=.9999)PetscPrintf(PETSC_COMM_WORLD, "elmt /// c_grad_x /// c_grad_y: %d %le %le\n",elmt,ibm->c_grad_x[elmt],ibm->c_grad_y[elmt]);
	} //end of elmt cells loop

	return ( 0 );												//	ali completed on 3 nov. 2009				    
}

PetscErrorCode AngleSkewness(IBMNodes *ibm){
PetscInt	 n_elmt =ibm->n_elmt, n_v  = ibm->n_v;
PetscReal	 nfx, nfy, nfz, dx12, dy12, dz12, dx13, dy13, dz13, dx23, dy23, dz23, dr, ds;
PetscInt         elmt,nbelmt,n1e,n2e,n3e;
PetscReal        amax, over, sum, a12_13,a13_23,a23_12;

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
			n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];

			dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
                   
                       a12_13 =fabs((180./PI)* acos((dx12*dx13+dy12*dy13)/(sqrt(dx12*dx12+dy12*dy12)*sqrt(dx13*dx13+dy13*dy13))));
                       a13_23 =fabs((180./PI)* acos((dx13*dx23+dy13*dy23)/(sqrt(dx13*dx13+dy13*dy13)*sqrt(dx23*dx23+dy23*dy23))));
                       a23_12 =fabs((180./PI)* acos((-dx23*dx12-dy23*dy12)/(sqrt(dx23*dx23+dy23*dy23)*sqrt(dx12*dx12+dy12*dy12))));
                       
                       sum = a12_13 + a13_23 + a23_12;
                       over = fabs(180.-sum);
                       
                       if(ibm->elmt_depth[elmt]<3. && over>1.) PetscPrintf(PETSC_COMM_WORLD, "Error elmt, over, a1, a2, a3::  %d %e %e %e %e \n",elmt,over, a12_13,a13_23,a23_12);
                       
                       amax = PetscMax (a12_13, a13_23);
                       amax = PetscMax (amax, a23_12);
                       
                       ibm->cent_z_AVE[elmt] = fabs(amax-60.0);

                       if(ibm->elmt_depth[elmt]<3.)PetscPrintf(PETSC_COMM_WORLD, "elmt, over, a1, a2, a3, AngleSkewness::  %d %e %e %e %e %e \n",elmt,over, a12_13,a13_23,a23_12,ibm->cent_z_AVE[elmt]);
         }
return(0);
}  
 
  
PetscErrorCode	sediment_flux_y_direction( UserCtx * user, IBMNodes * ibm, PetscInt ti, PetscInt tistart)
{
	DA			   da  = user->da, 
					fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
	PetscReal	    xc, yc, zc;
	PetscReal	    tmp;
	PetscReal	    caver, 
	                Ave_face12_u,
	                Ave_face12_v, 
					Ave_face12_w, 
					Ave_face13_u,
					Ave_face13_v, 
					Ave_face13_w, 
					Ave_face23_u,
					Ave_face23_v,
					Ave_face23_w, 
					nfx12, 
					nfz12, 
					nfx13, 
					nfz13, 
					nfx23, 
					nfz23, 
					nor_vel_to_face12, 
					nor_vel_to_face13, 
					nor_vel_to_face23, 
					le12x, 
					le12z, 					
					le13x, 
					le13z, 
					le23x, 
					le23z,
					e12x, 
					e12z, 
					e13x, 
					e13z,
					e23x, 
					e23z, 
					le_dot_n12, 
					le_dot_n13, 
					le_dot_n23;
      PetscReal  d_z,d_x,delY, delS;
      PetscReal  w_face, u_face, c_face, landaw,landau, landac,phiw, phiu, phic, phi1w, phi2w, phi1u, phi2u, phi1c, phi2c, RL, fR, fL, fx, beta=1./3.;
      PetscReal  cof=0.;
      PetscInt   upwind_v = 0;
      PetscInt   central_v = 0;
      PetscInt   upwind_c = 0;
      PetscInt   central_c = 0;
      PetscInt   gamma_c = 1;
      PetscInt   gamma_v = 1;
      PetscInt   SOU_v = 0;
      if(ti==tistart || ti==0){ PetscPrintf(PETSC_COMM_WORLD, "COF: %e \n",cof);
       PetscPrintf(PETSC_COMM_WORLD, "upwind_v and upwind_c: %d %d \n",upwind_v,upwind_c);
       PetscPrintf(PETSC_COMM_WORLD, "central_v and central_c: %d %d \n",central_v,central_c);
       PetscPrintf(PETSC_COMM_WORLD, "GAMMA_v and GAMMA_c: %d %d \n",gamma_v,gamma_c);}

// computing variables gradient if the higher order schems for velocity reconstruction or scalar asked
if(gamma_c || gamma_v || SOU_v) calc_cell_grad_y_direction(user, ibm);
               

// computing fluxes through all three faces: 1-->2 , 1-->3 , 2-->3
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	//	if( ibm->nf_z[ elmt ] < 1.e-6 )
		if(ibm->Rigidity[elmt])
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
			n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
			
            xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
			zc	 = ibm->cent_z[ elmt ];

			dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c1;  // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{

				//	if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		                        if(ibm->Rigidity[nbelmt] == 1)
					{
						ibm->F_flux_12[ elmt ] = 0.;
					}
					else
					{
					  	Ave_face12_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_face12_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

						Ave_face12_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );
						
                                                caver = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfz12 = dx12;
						nfx12 = -dz12;
						//ds	  = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr = sqrt( dz12 * dz12 + dx12 * dx12 );
                                                ds = dr;
						nfz12 /= dr;
						nfx12 /= dr;

						nor_vel_to_face12  = nfz12 * Ave_face12_w + nfx12 * Ave_face12_u;

						le12z = 0.5 * ( ibm->z_bp[ n2e ] + ibm->z_bp[ n1e ] ) - zc;
						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - xc;
                                                 
						e12z = 0.5 * ( ibm->z_bp[ n2e ] + ibm->z_bp[ n1e ] );
						e12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] );

						le_dot_n12 = nfz12 * le12z + nfx12 * le12x;
                                                 
                                                delS= sqrt((xc-ibm->cent_x[nbelmt])*(xc-ibm->cent_x[nbelmt])+
                                                                  (yc-ibm->cent_y[nbelmt])*(yc-ibm->cent_y[nbelmt])+
                                                                  (zc-ibm->cent_z[nbelmt])*(zc-ibm->cent_z[nbelmt]));
                                                
						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
														if(cof>1.e-6) delY=(ibm->cent_y[elmt]-ibm->cent_y_AVE[elmt])-(ibm->cent_y[nbelmt]-ibm->cent_y_AVE[nbelmt]);
														else delY = ibm->cent_y[elmt] - ibm->cent_y[nbelmt];
														
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
							if(upwind_v) tmp = ibm->Bvel[ elmt ].z * nfz12 + ibm->Bvel[ elmt ].x * nfx12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             
															 phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_z[elmt]*d_z + ibm->u_grad_x[elmt]*d_x;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);
															 
                                                             phi1w = ibm->Bvel[nbelmt].z - ibm->Bvel[elmt].z;
                                                             phi2w = ibm->w_grad_z[elmt]*d_z + ibm->w_grad_x[elmt]*d_x;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
                                                                                                                         
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[elmt].z;
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;

                                                             fR = sqrt((e12z-ibm->cent_z[nbelmt])*(e12z-ibm->cent_z[nbelmt])+
                                                                      (e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt]));
                                                             RL = sqrt((ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt])+
                                                                       (ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landaw = phiw/beta; 
                                                             landau = phiu/beta;
 
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[elmt].z + (1.-fx)* ibm->Bvel[nbelmt].z;
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[elmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[nbelmt].z;
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;

                                                             tmp =  nfz12 * w_face + nfx12 * u_face;
                                                           }

							if(upwind_c) ibm->F_flux_12[elmt]= -tmp * ds * ibm->SCont[ elmt ] 
                                                                                           -fabs(tmp * ds * ibm->SCont[ elmt ])*cof*delY/delS; 
							if(central_c)ibm->F_flux_12[elmt]= -tmp * ds * caver 
                                                                                           -fabs(tmp * ds * caver )*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_z[elmt]*d_z + ibm->c_grad_x[elmt]*d_x;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e12z-ibm->cent_z[nbelmt])*(e12z-ibm->cent_z[nbelmt])+
                                                                      (e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt]));
                                                             RL = sqrt((ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt])+
                                                                       (ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_12[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delY/delS; 
// PetscPrintf(PETSC_COMM_WORLD, "elmt /// phiu /// phiv /// phic: %d %le %le %le\n",elmt,phiu,phiv,phic);

                                                           }
              						}
							else	//the velocity is in negative direction					    
							{
							if(upwind_v) tmp = ibm->Bvel[ nbelmt ].z * nfz12 + ibm->Bvel[ nbelmt ].x * nfx12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             
                                                             phi1w = ibm->Bvel[elmt].z - ibm->Bvel[nbelmt].z;
                                                             phi2w = ibm->w_grad_z[nbelmt]*d_z + ibm->w_grad_x[nbelmt]*d_x;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);

                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_z[nbelmt]*d_z + ibm->u_grad_x[nbelmt]*d_x;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);
                                                             
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[nbelmt].z;
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;

                                                             fR = sqrt((e12z-ibm->cent_z[elmt])*(e12z-ibm->cent_z[elmt])+
                                                                      (e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt]));
                                                             RL = sqrt((ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt])+
                                                                       (ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landaw = phiw/beta; 
                                                             landau = phiu/beta;
 
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[nbelmt].z + (1.-fx)* ibm->Bvel[elmt].z;
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[nbelmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[elmt].z;

                                                             tmp =  nfx12 * u_face + nfz12 * w_face;
                                                           }

							if(upwind_c) ibm->F_flux_12[elmt]= -tmp * ds * ibm->SCont[ nbelmt ] 
                                                                                          -fabs(tmp * ds * ibm->SCont[ nbelmt ])*cof*delY/delS; 
							if(central_c)ibm->F_flux_12[elmt]= -tmp * ds * caver 
                                                                                            -fabs(tmp * ds * caver)*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_z[nbelmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12z-ibm->cent_z[elmt])*(e12z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_12[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
														if(cof>1.e-6) delY=(ibm->cent_y[elmt]-ibm->cent_y_AVE[elmt])-(ibm->cent_y[nbelmt]-ibm->cent_y_AVE[nbelmt]);
														else delY=ibm->cent_y[elmt] - ibm->cent_y[nbelmt];
														
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
							if(upwind_v) tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].z * nfz12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_z[nbelmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1w = ibm->Bvel[elmt].z - ibm->Bvel[nbelmt].z;
                                                             phi2w = ibm->w_grad_x[nbelmt]*d_x + ibm->w_grad_z[nbelmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[nbelmt].z;

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12z-ibm->cent_z[elmt])*(e12z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[nbelmt].z + (1.-fx)* ibm->Bvel[elmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[nbelmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[elmt].z;

                                                             tmp =  nfx12 * u_face + nfz12 * w_face;
                                                           }

							if(upwind_c) ibm->F_flux_12[elmt] = tmp * ds * ibm->SCont[ nbelmt ]
                                                                                            +fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_12[ elmt ] = tmp * ds * caver 
                                                                                              +fabs(tmp * ds *caver)*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_z[nbelmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12z-ibm->cent_z[elmt])*(e12z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_12[elmt]= tmp * ds * c_face 
                                                                  +fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
							else	//the velocity is in negative direction		
                                                        {
							if(upwind_v)tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].z * nfz12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_z[elmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1w = ibm->Bvel[nbelmt].z - ibm->Bvel[elmt].z;
                                                             phi2w = ibm->w_grad_x[elmt]*d_x + ibm->w_grad_z[elmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
												   
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[elmt].z;

                                                             fR = sqrt((e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt])+
                                                                      (e12z-ibm->cent_z[nbelmt])*(e12z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[elmt].z + (1.-fx)* ibm->Bvel[nbelmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[elmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[nbelmt].z;

                                                             tmp =  nfx12 * u_face + nfz12 * w_face;
                                                           }

							if(upwind_c)ibm->F_flux_12[elmt] = tmp * ds * ibm->SCont[elmt]
                                                                                           +fabs(tmp * ds * ibm->SCont[elmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_12[ elmt ] = tmp * ds * caver 
                                                                                                    + fabs(tmp * ds * caver )*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_z[elmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt])+
                                                                      (e12z-ibm->cent_z[nbelmt])*(e12z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_12[elmt]= tmp * ds * c_face 
                                                                                   +fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }

              						
							}
						}
					}
				}

				//					1---->3 edge flux 

				// xiaolei deactivate 
				/*
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{
				//	if( ibm->nf_z[ nbelmt ] < 1.e-6 )
	                              	if(ibm->Rigidity[nbelmt] == 1)
					{
						ibm->F_flux_13[ elmt ] = 0.;
					}
					else
					{

						Ave_face13_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_face13_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

						Ave_face13_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );
                                                
                        caver = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfz13 = dx13;
						nfx13 = -dz13;
						//ds	  = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr = sqrt( dx13 * dx13 + dz13 * dz13 );
                                                ds = dr;
						nfx13 /= dr;
						nfz13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfz13 * Ave_face13_w;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - xc;
						le13z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n1e ] ) - zc;
						e13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] );
						e13z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n1e ] );

						le_dot_n13 = nfx13 * le13x + nfz13 * le13z;
                                                
                                                delS= sqrt((xc-ibm->cent_x[nbelmt])*(xc-ibm->cent_x[nbelmt])+
                                                                  (yc-ibm->cent_y[nbelmt])*(yc-ibm->cent_y[nbelmt])+
                                                                  (zc-ibm->cent_z[nbelmt])*(zc-ibm->cent_z[nbelmt]));


						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
                                                        if(cof>1.e-6) delY=(ibm->cent_y[elmt]-ibm->cent_y_AVE[elmt])-(ibm->cent_y[nbelmt]-ibm->cent_y_AVE[nbelmt]);
														else delY=ibm->cent_y[elmt] - ibm->cent_y[nbelmt];

							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
							if(upwind_v)tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].z * nfz13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_z[elmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1w = ibm->Bvel[nbelmt].z - ibm->Bvel[elmt].z;
                                                             phi2w = ibm->w_grad_x[elmt]*d_x + ibm->w_grad_z[elmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[elmt].z;

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13z-ibm->cent_z[nbelmt])*(e13z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[elmt].z + (1.-fx)* ibm->Bvel[nbelmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[elmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[nbelmt].z;

                                                             tmp =  nfx13 * u_face + nfz13 * w_face;
                                                           }

					          	if(upwind_c)ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ]
                                                                                             -fabs(tmp * ds * ibm->SCont[elmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_13[ elmt ] = -tmp * ds * caver
                                                                                              -fabs(tmp * ds * caver)*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_z[elmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13z-ibm->cent_z[nbelmt])*(e13z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_13[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
							else	//the velocity is in negative direction		
							{
							if(upwind_v)tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].z * nfz13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_z[nbelmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1w = ibm->Bvel[elmt].z - ibm->Bvel[nbelmt].z;
                                                             phi2w = ibm->w_grad_x[nbelmt]*d_x + ibm->w_grad_z[nbelmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[nbelmt].z;

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13z-ibm->cent_z[elmt])*(e13z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[nbelmt].z + (1.-fx)* ibm->Bvel[elmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[nbelmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[elmt].z;

                                                             tmp =  nfx13 * u_face + nfz13 * w_face;
                                                           }

							if(upwind_c)ibm->F_flux_13[elmt] = -tmp * ds * ibm->SCont[nbelmt]
                                                                                           -fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_13[elmt]= -tmp * ds * caver
                                                                                           -fabs(tmp * ds * caver)*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_z[nbelmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13z-ibm->cent_z[elmt])*(e13z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_13[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
                                                        //delZ=(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt])-(ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt]);
														if(cof>1.e-6) delY=(ibm->cent_y[elmt]-ibm->cent_y_AVE[elmt])-(ibm->cent_y[nbelmt]-ibm->cent_y_AVE[nbelmt]);
														else delY=ibm->cent_y[elmt] - ibm->cent_y[nbelmt];

							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
							if(upwind_v)tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].z * nfz13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_z[nbelmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1w = ibm->Bvel[elmt].z - ibm->Bvel[nbelmt].z;
                                                             phi2w = ibm->w_grad_x[nbelmt]*d_x + ibm->w_grad_z[nbelmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[nbelmt].z;

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13z-ibm->cent_z[elmt])*(e13z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[nbelmt].z + (1.-fx)* ibm->Bvel[elmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[nbelmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[elmt].z;

                                                             tmp =  nfx13 * u_face + nfz13 * w_face;
                                                           }

							if(upwind_c)ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ]
                                                                                           + fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_13[ elmt ] = tmp * ds * caver 
                                                                                           + fabs(tmp * ds * caver )*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_z[nbelmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13z-ibm->cent_z[elmt])*(e13z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_13[elmt]= tmp * ds * c_face 
                                                                  +fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
							else	//the velocity is in negative direction	
                                 			{
							if(upwind_v)tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].z * nfz13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_z[elmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1w = ibm->Bvel[nbelmt].z - ibm->Bvel[elmt].z;
                                                             phi2w = ibm->w_grad_x[elmt]*d_x + ibm->w_grad_z[elmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[elmt].z;

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13z-ibm->cent_z[nbelmt])*(e13z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[elmt].z + (1.-fx)* ibm->Bvel[nbelmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[elmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[nbelmt].z;

                                                             tmp =  nfx13 * u_face + nfz13 * w_face;
                                                           }

							if(upwind_c)ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ]
                                                                                           + fabs(tmp * ds * ibm->SCont[elmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_13[ elmt ] = tmp * ds * caver
                                                                                            + fabs(tmp * ds * caver)*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_z[elmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13z-ibm->cent_z[nbelmt])*(e13z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_13[elmt]= tmp * ds * c_face 
                                                                                   +fabs(tmp * ds * c_face )*cof*delY/delS; 
							}}
						}
					}
				}

				//					2---->3 edge  flux

				/* xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c3; // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{
			//		if( ibm->nf_z[ nbelmt ] < 1.e-6 )
	                               	if(ibm->Rigidity[nbelmt] == 1)
					{
						ibm->F_flux_23[ elmt ] = 0.;
					}
					else
					{
						Ave_face23_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_face23_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

						Ave_face23_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );

                                                caver = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfz23 = dx23;
						nfx23 = -dz23;
						//ds	  = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr = sqrt( dx23 * dx23 + dz23 * dz23 );
                                                ds = dr;
						nfx23 /= dr;
						nfz23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfz23 * Ave_face23_w;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - xc;
						le23z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n2e ] ) - zc;
						e23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] );
						e23z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n2e ] );

						le_dot_n23 = nfx23 * le23x + nfz23 * le23z;

                                                delS= sqrt((xc-ibm->cent_x[nbelmt])*(xc-ibm->cent_x[nbelmt])+
                                                                  (yc-ibm->cent_y[nbelmt])*(yc-ibm->cent_y[nbelmt])+
                                                                  (zc-ibm->cent_z[nbelmt])*(zc-ibm->cent_z[nbelmt]));


						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
                                                        //delZ=(ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
														if(cof>1.e-6) delY=(ibm->cent_y[elmt]-ibm->cent_y_AVE[elmt])-(ibm->cent_y[nbelmt]-ibm->cent_y_AVE[nbelmt]);
														else delY=ibm->cent_y[elmt] - ibm->cent_y[nbelmt];
														
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
 							if(upwind_v) tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].z * nfz23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_z[elmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u +1.e-8);

                                                             phi1w = ibm->Bvel[nbelmt].z - ibm->Bvel[elmt].z;
                                                             phi2w = ibm->w_grad_x[elmt]*d_x + ibm->w_grad_z[elmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[elmt].z;

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23z-ibm->cent_z[nbelmt])*(e23z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[elmt].z + (1.-fx)* ibm->Bvel[nbelmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[elmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[nbelmt].z;

                                                             tmp =  nfx23 * u_face + nfz23 * w_face;
                                                           }

							if(upwind_c) ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ]
                                                                                             -fabs(tmp * ds * ibm->SCont[ elmt ])*cof*delY/delS; 
						        if(central_c)ibm->F_flux_23[ elmt ] = -tmp * ds * caver 
                                                                                             -fabs(tmp * ds *caver)*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_z[elmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23z-ibm->cent_z[nbelmt])*(e23z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_23[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
							else	//the velocity is in negative direction				
							{
							if(upwind_v)tmp = ibm->Bvel[nbelmt].x * nfx23 + ibm->Bvel[ nbelmt ].z * nfz23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_z[nbelmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u+1.e-8);

                                                             phi1w = ibm->Bvel[elmt].z - ibm->Bvel[nbelmt].z;
                                                             phi2w = ibm->w_grad_x[nbelmt]*d_x + ibm->w_grad_z[nbelmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[nbelmt].z;

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23z-ibm->cent_z[elmt])*(e23z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[nbelmt].z + (1.-fx)* ibm->Bvel[elmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[nbelmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[elmt].z;

                                                             tmp =  nfx23 * u_face + nfz23 * w_face;
                                                           }

							if(upwind_c)ibm->F_flux_23[elmt]= -tmp * ds * ibm->SCont[ nbelmt ]
                                                                                          -fabs(tmp * ds * ibm->SCont[ nbelmt ])*cof*delY/delS; 
							if(central_c)ibm->F_flux_23[elmt]= -tmp * ds * caver 
                                                                                           -fabs(tmp * ds * caver )*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_z[nbelmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23z-ibm->cent_z[elmt])*(e23z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_23[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
                                                        //delZ=(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt])-(ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt]);
														if(cof>1.e-6) delY=(ibm->cent_y[elmt]-ibm->cent_y_AVE[elmt])-(ibm->cent_y[nbelmt]-ibm->cent_y_AVE[nbelmt]);
														else delY=ibm->cent_y[elmt] - ibm->cent_y[nbelmt];
														
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
							if(upwind_v) tmp = ibm->Bvel[nbelmt].x * nfx23 + ibm->Bvel[ nbelmt ].z * nfz23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_z[nbelmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u+1.e-8);

                                                             phi1w = ibm->Bvel[elmt].z - ibm->Bvel[nbelmt].z;
                                                             phi2w = ibm->w_grad_x[nbelmt]*d_x + ibm->w_grad_z[nbelmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[nbelmt].z;

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23z-ibm->cent_z[elmt])*(e23z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[nbelmt].z + (1.-fx)* ibm->Bvel[elmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[nbelmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[elmt].z;

                                                             tmp =  nfx23 * u_face + nfz23 * w_face;
                                                           }

							if(upwind_c) ibm->F_flux_23[elmt]= tmp * ds * ibm->SCont[ nbelmt ]
                                                                                           +fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delY/delS; 
						        if(central_c)ibm->F_flux_23[elmt] = tmp * ds * caver 
                                                                                            +fabs(tmp * ds * caver )*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_z = ibm->cent_z[elmt]- ibm->cent_z[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_z[nbelmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23z-ibm->cent_z[elmt])*(e23z-ibm->cent_z[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_23[elmt]= tmp * ds * c_face 
                                                                  +fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							}
							else	//the velocity is in negative direction					    
							{
							if(upwind_v) tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].z * nfz23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_z[elmt]*d_z;
                                                             phiu  = 1.-phi1u/(2.*phi2u+1.e-8);

                                                             phi1w = ibm->Bvel[nbelmt].z - ibm->Bvel[elmt].z;
                                                             phi2w = ibm->w_grad_x[elmt]*d_x + ibm->w_grad_z[elmt]*d_z;
                                                             phiw  = 1.-phi1w/(2.*phi2w+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiw <= 0. || phiw >= 1.) w_face = ibm->Bvel[elmt].z;

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23z-ibm->cent_z[nbelmt])*(e23z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landaw = phiw/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiw < 1. && phiw >= beta) w_face = fx * ibm->Bvel[elmt].z + (1.-fx)* ibm->Bvel[nbelmt].z;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiw < beta && phiw > 0.) w_face = (1.-landaw*(1.-fx)) * ibm->Bvel[elmt].z +
                                                                                                  (1.-fx)*landaw * ibm->Bvel[nbelmt].z;

                                                             tmp =  nfx23 * u_face + nfz23 * w_face;
                                                           }

							if(upwind_c) ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ]
                                                                                              +fabs(tmp * ds * ibm->SCont[elmt])*cof*delY/delS; 
							if(central_c)ibm->F_flux_23[ elmt ] = tmp * ds * caver
                                                                                              +fabs(tmp * ds * caver )*cof*delY/delS; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_z = ibm->cent_z[nbelmt]- ibm->cent_z[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_z[elmt]*d_z;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23z-ibm->cent_z[nbelmt])*(e23z-ibm->cent_z[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_z[elmt]-ibm->cent_z[nbelmt])*(ibm->cent_z[elmt]-ibm->cent_z[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_23[elmt]= tmp * ds * c_face 
                                                                                    +fabs(tmp * ds * c_face )*cof*delY/delS; 
                                                           }
							 }
						}
					}
				}
			}
		}
	}
	PetscPrintf(PETSC_COMM_WORLD, "test1 \n"); //Hossein-debug
	return ( 0 );												//	ali completed on 3 nov. 2009				    
}


PetscErrorCode	outlet_sediment_flux_bend_y_direction( UserCtx * user, IBMNodes * ibm )
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscReal	   sb, 
					sc;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
       	PetscReal	   tmp;
	PetscReal	        	nfx12, 
					nfz12, 
					nfx13, 
					nfz13, 
					nfx23, 
					nfz23; 
	PetscReal                       //y_outlet = 3.41,
                                        //y_plus_x_outlet=29.11,                                            
                                        someX_minus_Z_outlet = 39.17196667;//33.07459,                                            
                                        //someX_plus_Y_inlet = 100.05944;                                            
                                        //someX_plus_Y_inlet = 90.16268;                                            
        PetscReal                       zdif1, zdif2, zdif3,xzdif1,xzdif2,xzdif3;
        PetscReal                       sxzdif1,sxzdif2,sxzdif3;
        PetscReal                       ssxzdif1,ssxzdif2,ssxzdif3;
	PetscReal                       Ave_face12_u,
                                        Ave_face13_u,
                                        Ave_face23_u,
                                        Ave_face12_v, 
                                        Ave_face13_v, 
                                        Ave_face23_v,
										Ave_face12_w, 
                                        Ave_face13_w, 
                                        Ave_face23_w,
                                        nor_vel_to_face12,
                                        nor_vel_to_face13,
                                        nor_vel_to_face23,
                                     	le12x,
                                        le12y,
										le12z,
                                        le13x,
                                        le13y,le13z,
                                        le23x,
                                        le23y,le23z,
                                        le_dot_n12,
                                        le_dot_n13,
                                        le_dot_n23,
                                        xc,yc,zc;
                

// computing fluxes through oulet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
		//	zc       = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
			zc	 = ibm->cent_z[ elmt ];
                        
                        //ydif1=fabs(y_outlet-ibm->y_bp[n1e]);
			//ydif2=fabs(y_outlet-ibm->y_bp[n2e]); 
	                //ydif3=fabs(y_outlet-ibm->y_bp[n3e]);
                      
                        //xydif1=fabs(y_plus_x_outlet-ibm->y_bp[n1e]-ibm->x_bp[n1e]);
                        //xydif2=fabs(y_plus_x_outlet-ibm->y_bp[n2e]-ibm->x_bp[n2e]);
                        //xydif3=fabs(y_plus_x_outlet-ibm->y_bp[n3e]-ibm->x_bp[n3e]);
		     
                        sxzdif1=fabs(someX_minus_Z_outlet-(0.85942*ibm->x_bp[n1e]-ibm->z_bp[n1e]));
                        sxzdif2=fabs(someX_minus_Z_outlet-(0.85942*ibm->x_bp[n2e]-ibm->z_bp[n2e]));
                        sxzdif3=fabs(someX_minus_Z_outlet-(0.85942*ibm->x_bp[n3e]-ibm->z_bp[n3e]));

                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
                        
                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
//    Bend 90 degree                                       
	     // if((ydif1<1.e-6 && ydif2<1.e-6)||(ydif1<1.e-6 && ydif3<1.e-6)||(ydif2<1.e-6 && ydif3<1.e-6)) 


 //   Bend 135 degree
	     // if((xydif1<1.e-6 && xydif2<1.e-6)||(xydif1<1.e-6 && xydif3<1.e-6)||(xydif2<1.e-6 && xydif3<1.e-6))  
 

//   OSL-2013
	       //if((sxydif1<1.e-4 && sxydif2<1.e-4)||(sxydif1<1.e-4 && sxydif3<1.e-4)||(sxydif2<1.e-4 && sxydif3<1.e-4)) 
               // if((ssxydif1<1.e-4 && ssxydif2<1.e-4)||(ssxydif1<1.e-4 && ssxydif3<1.e-4)||(ssxydif2<1.e-4 && sxydif3<1.e-4))  
//   OSL-2016
	      if((sxzdif1<1.e-4 && sxzdif2<1.e-4)||(sxzdif1<1.e-4 && sxzdif3<1.e-4)||(sxzdif2<1.e-4 && sxzdif3<1.e-4)) 
	      {


		if(ibm->Rigidity[elmt])
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		    dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				//
				/* // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	                        if( ibm->Rigidity[nbelmt] == 1)
      					{

					       	Ave_face12_u = ibm->Bvel[ elmt ].x;

						Ave_face12_w = ibm->Bvel[ elmt ].z;

						nfz12 = dx12;
						nfx12 = -dz12;
						//ds    = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr    = sqrt( dx12 * dx12 + dz12 * dz12 );
                                                ds = dr;
						nfx12 /= dr;
						nfz12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfz12 * Ave_face12_w;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le12z = 0.5 * ( ibm->z_bp[ n2e ] + ibm->z_bp[ n1e ] ) - ibm->cent_z[ elmt ];

						le_dot_n12 = nfx12 * le12x + nfz12 * le12z;


						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].z * nfz12;
								ibm->F_flux_12[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
              						}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].z * nfz12;
								ibm->F_flux_12[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].z * nfz12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction		
                                                        {
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].z * nfz12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
                                             
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*
				if( nbelmt != elmt   // xiaolei deactivate 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt])
					{

                                                Ave_face13_u = ibm->Bvel[ elmt ].x;

						Ave_face13_w = ibm->Bvel[ elmt ].z;

						nfz13 = dx13;
						nfx13 = -dz13;
						//ds    = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr    = sqrt( dx13 * dx13 + dz13 * dz13 );
                                                ds = dr;
						nfx13 /= dr;
						nfz13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfz13 * Ave_face13_w;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le13z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n1e ] ) - ibm->cent_z[ elmt ];

						le_dot_n13 = nfx13 * le13x + nfz13 * le13z;


						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].z * nfz13;
								ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].z * nfz13;
								ibm->F_flux_13[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].z * nfz13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction	
                                 			{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].z * nfz13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}

					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                            	if( ibm->Rigidity[nbelmt] == 1)
					{
                                              	Ave_face23_u = ibm->Bvel[ elmt ].x;

						Ave_face23_w = ibm->Bvel[ elmt ].z;

						nfz23 = dx23;
						nfx23 = -dz23;
						//ds    = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr    = sqrt( dx23 * dx23 + dz23 * dz23 );
                                                ds = dr;
						nfx23 /= dr;
						nfz23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfz23 * Ave_face23_w;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - ibm->cent_x[ elmt ];
						le23z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n2e ] ) - ibm->cent_z[ elmt ];

						le_dot_n23 = nfx23 * le23x + nfz23 * le23z;


						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].z * nfz23;
								ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction				
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].z * nfz23;
								ibm->F_flux_23[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].z * nfz23;
								ibm->F_flux_23[ elmt ]
								= tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].z * nfz23;
								ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
						
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling
   PetscPrintf(PETSC_COMM_WORLD, "test2:  \n"); //Hossein-debug
	return ( 0 );		//	ali completed on 19 nov. 2009				    
}


PetscErrorCode	outlet_sediment_flux_contra_y_direction( UserCtx * user, IBMNodes * ibm )
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscReal	   sb, 
					sc;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
       	PetscReal	   tmp;
	PetscReal	        	nfx12, 
					nfz12, 
					nfx13, 
					nfz13, 
					nfx23, 
					nfz23; 
	PetscReal                       x_outlet = 165.0;//739.9;//90.0; //total osl//25.; // Jhoke Crossvane Rockvane 50.;//49.71;
	PetscReal                       z_outlet = 10; //35.0; 
        PetscReal                       xdif1,xdif2,xdif3;
        PetscReal                       zdif1,zdif2,zdif3;
	PetscReal                       Ave_face12_u,
                                        Ave_face13_u,
                                        Ave_face23_u,
                                        Ave_face12_w, 
                                        Ave_face13_w, 
                                        Ave_face23_w,
                                        nor_vel_to_face12,
                                        nor_vel_to_face13,
                                        nor_vel_to_face23,
                                     	le12x,
                                        le12z,
                                        le13x,
                                        le13z,
                                        le23x,
                                        le23z,
                                        le_dot_n12,
                                        le_dot_n13,
                                        le_dot_n23;

// computing fluxes through oulet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                                                                   
                        xdif1=fabs(x_outlet-ibm->x_bp[n1e]);
			xdif2=fabs(x_outlet-ibm->x_bp[n2e]); 
	                xdif3=fabs(x_outlet-ibm->x_bp[n3e]);

                        zdif1=fabs(z_outlet-ibm->z_bp[n1e]);
			zdif2=fabs(z_outlet-ibm->z_bp[n2e]); 
	                zdif3=fabs(z_outlet-ibm->z_bp[n3e]);

	      //if((xdif1<1.e-4 && xdif2<1.e-4)||(xdif1<1.e-4 && xdif3<1.e-4)||(xdif2<1.e-4 && xdif3<1.e-4))  
	      if((zdif1<1.e-4 && zdif2<1.e-4)||(zdif1<1.e-4 && zdif3<1.e-4)||(zdif2<1.e-4 && zdif3<1.e-4))  
	      {


		if( ibm->Rigidity[elmt] == 1)
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		        dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ ) // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/* // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	            if( ibm->Rigidity[nbelmt] == 1)
      					{
					    Ave_face12_u = ibm->Bvel[ elmt ].x;
						Ave_face12_w = ibm->Bvel[ elmt ].z;

						nfz12 = dx12;
						nfx12 = -dz12;
						//ds    = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr    = sqrt( dx12 * dx12 + dz12 * dz12 );
                                                ds = dr;
						nfx12 /= dr;
						nfz12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfz12 * Ave_face12_w;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le12z = 0.5 * ( ibm->z_bp[ n2e ] + ibm->z_bp[ n1e ] ) - ibm->cent_z[ elmt ];

						le_dot_n12 = nfx12 * le12x + nfz12 * le12z;


						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].z * nfz12;
								ibm->F_flux_12[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
              						}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].z * nfz12;
								ibm->F_flux_12[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].z * nfz12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction		
                                                        {
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].z * nfz12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
                                             
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*   // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if( ibm->Rigidity[nbelmt] == 1)
					{

                        Ave_face13_u = ibm->Bvel[ elmt ].x;
						Ave_face13_w = ibm->Bvel[ elmt ].z;

						nfz13 = dx13;
						nfx13 = -dz13;
						//ds    = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr    = sqrt( dx13 * dx13 + dz13 * dz13 );
                                                ds = dr;
						nfx13 /= dr;
						nfz13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfz13 * Ave_face13_w;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le13z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n1e ] ) - ibm->cent_z[ elmt ];

						le_dot_n13 = nfx13 * le13x + nfz13 * le13z;


						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].z * nfz13;
								ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].z * nfz13;
								ibm->F_flux_13[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].z * nfz13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction	
                                 			{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].z * nfz13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}

					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                if( ibm->Rigidity[nbelmt] == 1)
					{
                       	Ave_face23_u = ibm->Bvel[ elmt ].x;
						Ave_face23_w = ibm->Bvel[ elmt ].z;

						nfz23 = dx23;
						nfx23 = -dz23;
						//ds    = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr    = sqrt( dx23 * dx23 + dz23 * dz23 );
                                                ds = dr;
						nfx23 /= dr;
						nfz23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfz23 * Ave_face23_w;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - ibm->cent_x[ elmt ];
						le23z = 0.5 * ( ibm->z_bp[ n3e ] + ibm->z_bp[ n2e ] ) - ibm->cent_z[ elmt ];

						le_dot_n23 = nfx23 * le23x + nfz23 * le23z;


						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].z * nfz23;
								ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction				
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].z * nfz23;
								ibm->F_flux_23[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].z * nfz23;
								ibm->F_flux_23[ elmt ]
								= tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].z * nfz23;
								ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
						
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the x_bp value if at outlet
   }  // for loop, elmt cycling
   PetscPrintf(PETSC_COMM_WORLD, "test3 %d \n"); //Hossein-debug
	return ( 0 );		//	ali completed on 4 Dec. 2009				    
}


PetscErrorCode	sediment_flux( UserCtx * user, IBMNodes * ibm, PetscInt ti, PetscInt tistart)
{
	DA			   da  = user->da, 
					fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
	PetscReal	   xc, yc, zc;
	PetscReal	   tmp;
	PetscReal	  caver, Ave_face12_u, 
					Ave_face12_v, 
					Ave_face12_w, 
					Ave_face13_u, 
					Ave_face13_v, 
					Ave_face13_w, 
					Ave_face23_u, 
					Ave_face23_v, 
					Ave_face23_w, 
					nfx12, 
					nfy12, 
					nfx13, 
					nfy13, 
					nfx23, 
					nfy23, 
					nor_vel_to_face12, 
					nor_vel_to_face13, 
					nor_vel_to_face23, 
					le12x, 
					le12y, 
					le13x, 
					le13y, 
					le23x, 
					le23y, 
					e12x, 
					e12y, 
					e13x, 
					e13y, 
					e23x, 
					e23y, 
					le_dot_n12, 
					le_dot_n13, 
					le_dot_n23;
      PetscReal  d_x,d_y,delZ, delX;
      PetscReal  u_face, v_face, c_face, landau,landav, landac,phiu, phiv, phic, phi1u, phi2u, phi1v, phi2v, phi1c, phi2c, RL, fR, fL, fx, beta=1./3.;
      PetscReal  cof=0.;
      PetscInt   upwind_v = 0;
      PetscInt   central_v = 0;
      PetscInt   upwind_c = 0;
      PetscInt   central_c = 0;
      PetscInt   gamma_c = 1;
      PetscInt   gamma_v = 1;
      PetscInt   SOU_v = 0;
      if(ti==tistart || ti==0){ PetscPrintf(PETSC_COMM_WORLD, "COF: %e \n",cof);
       PetscPrintf(PETSC_COMM_WORLD, "upwind_v and upwind_c: %d %d \n",upwind_v,upwind_c);
       PetscPrintf(PETSC_COMM_WORLD, "central_v and central_c: %d %d \n",central_v,central_c);
       PetscPrintf(PETSC_COMM_WORLD, "GAMMA_v and GAMMA_c: %d %d \n",gamma_v,gamma_c);}

// computing variables gradient if the higher order schems for velocity reconstruction or scalar asked
if(gamma_c || gamma_v || SOU_v) calc_cell_grad(user, ibm);
               
// computing fluxes through all three faces: 1-->2 , 1-->3 , 2-->3
	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	//	if( ibm->nf_z[ elmt ] < 1.e-6 )
		if(ibm->Rigidity[elmt] == 1)
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
			n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
			
            xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
			zc	 = ibm->cent_z[ elmt ];

			dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c1;  // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{

				//	if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		                        if(ibm->Rigidity[nbelmt] == 1)
					{
						ibm->F_flux_12[ elmt ] = 0.;
					}
					else
					{
					  	Ave_face12_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_face12_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

						Ave_face12_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );
						
                                                caver = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfx12 = dy12;
						nfy12 = -dx12;
						//ds	  = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr = sqrt( dx12 * dx12 + dy12 * dy12 );
                                                ds = dr;
						nfx12 /= dr;
						nfy12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfy12 * Ave_face12_v;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - xc;
						le12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] ) - yc;
                                                 
						e12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] );
						e12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] );

						le_dot_n12 = nfx12 * le12x + nfy12 * le12y;
                                                 
                                                delX= sqrt((xc-ibm->cent_x[nbelmt])*(xc-ibm->cent_x[nbelmt])+
                                                                  (yc-ibm->cent_y[nbelmt])*(yc-ibm->cent_y[nbelmt])+
                                                                  (zc-ibm->cent_z[nbelmt])*(zc-ibm->cent_z[nbelmt]));
                                                
						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
                                                        if (cof > 1.e-6) delZ = (ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
                                                        else delZ = ibm->cent_z[elmt] - ibm->cent_z[nbelmt];
														
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
							if(upwind_v) tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_y[elmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[nbelmt].y - ibm->Bvel[elmt].y;
                                                             phi2v = ibm->v_grad_x[elmt]*d_x + ibm->v_grad_y[elmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[elmt].y;

                                                             fR = sqrt((e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt])+
                                                                      (e12y-ibm->cent_y[nbelmt])*(e12y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[elmt].y + (1.-fx)* ibm->Bvel[nbelmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[elmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[nbelmt].y;

                                                             tmp =  nfx12 * u_face + nfy12 * v_face;
                                                           }

							if(upwind_c) ibm->F_flux_12[elmt]= -tmp * ds * ibm->SCont[ elmt ] 
                                                                                           -fabs(tmp * ds * ibm->SCont[ elmt ])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_12[elmt]= -tmp * ds * caver 
                                                                                           -fabs(tmp * ds * caver )*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_y[elmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt])+
                                                                      (e12y-ibm->cent_y[nbelmt])*(e12y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_12[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delZ/delX; 
// PetscPrintf(PETSC_COMM_WORLD, "elmt /// phiu /// phiv /// phic: %d %le %le %le\n",elmt,phiu,phiv,phic);

                                                           }

              						}
							else	//the velocity is in negative direction					    
							{
							if(upwind_v) tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].y * nfy12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_y[nbelmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[elmt].y - ibm->Bvel[nbelmt].y;
                                                             phi2v = ibm->v_grad_x[nbelmt]*d_x + ibm->v_grad_y[nbelmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[nbelmt].y;

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12y-ibm->cent_y[elmt])*(e12y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[nbelmt].y + (1.-fx)* ibm->Bvel[elmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[nbelmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[elmt].y;

                                                             tmp =  nfx12 * u_face + nfy12 * v_face;
                                                           }

							if(upwind_c) ibm->F_flux_12[elmt]= -tmp * ds * ibm->SCont[ nbelmt ] 
                                                                                          -fabs(tmp * ds * ibm->SCont[ nbelmt ])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_12[elmt]= -tmp * ds * caver 
                                                                                            -fabs(tmp * ds * caver)*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_y[nbelmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12y-ibm->cent_y[elmt])*(e12y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_12[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }

							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
														if (cof > 1.e-6) delZ = (ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
                                                        else delZ = ibm->cent_z[elmt] - ibm->cent_z[nbelmt];
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
							if(upwind_v) tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].y * nfy12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_y[nbelmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[elmt].y - ibm->Bvel[nbelmt].y;
                                                             phi2v = ibm->v_grad_x[nbelmt]*d_x + ibm->v_grad_y[nbelmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[nbelmt].y;

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12y-ibm->cent_y[elmt])*(e12y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[nbelmt].y + (1.-fx)* ibm->Bvel[elmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[nbelmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[elmt].y;

                                                             tmp =  nfx12 * u_face + nfy12 * v_face;
                                                           }

							if(upwind_c) ibm->F_flux_12[elmt] = tmp * ds * ibm->SCont[ nbelmt ]
                                                                                            +fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_12[ elmt ] = tmp * ds * caver 
                                                                                              +fabs(tmp * ds *caver)*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_y[nbelmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e12x-ibm->cent_x[elmt])*(e12x-ibm->cent_x[elmt])+
                                                                      (e12y-ibm->cent_y[elmt])*(e12y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_12[elmt]= tmp * ds * c_face 
                                                                  +fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
							else	//the velocity is in negative direction		
                                                        {
							if(upwind_v)tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
                                                        if(central_v)tmp = nor_vel_to_face12;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_y[elmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[nbelmt].y - ibm->Bvel[elmt].y;
                                                             phi2v = ibm->v_grad_x[elmt]*d_x + ibm->v_grad_y[elmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[elmt].y;

                                                             fR = sqrt((e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt])+
                                                                      (e12y-ibm->cent_y[nbelmt])*(e12y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[elmt].y + (1.-fx)* ibm->Bvel[nbelmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[elmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[nbelmt].y;

                                                             tmp =  nfx12 * u_face + nfy12 * v_face;
                                                           }

							if(upwind_c)ibm->F_flux_12[elmt] = tmp * ds * ibm->SCont[elmt]
                                                                                           +fabs(tmp * ds * ibm->SCont[elmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_12[ elmt ] = tmp * ds * caver 
                                                                                                    + fabs(tmp * ds * caver )*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_y[elmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e12x-ibm->cent_x[nbelmt])*(e12x-ibm->cent_x[nbelmt])+
                                                                      (e12y-ibm->cent_y[nbelmt])*(e12y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_12[elmt]= tmp * ds * c_face 
                                                                                   +fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }

              						
							}
						}
					}
				}

				//					1---->3 edge flux 

				// xiaolei deactivate 
				/*
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/


				nbelmt=ibm->c2c[elmt].c2; // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{
				//	if( ibm->nf_z[ nbelmt ] < 1.e-6 )
	                              	if(ibm->Rigidity[nbelmt] == 1)
					{
						ibm->F_flux_13[ elmt ] = 0.;
					}
					else
					{

						Ave_face13_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_face13_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

						Ave_face13_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );
                                                
                                                caver = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfx13 = dy13;
						nfy13 = -dx13;
						//ds	  = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr = sqrt( dx13 * dx13 + dy13 * dy13 );
                                                ds = dr;
						nfx13 /= dr;
						nfy13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfy13 * Ave_face13_v;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - xc;
						le13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] ) - yc;
						e13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] );
						e13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] );

						le_dot_n13 = nfx13 * le13x + nfy13 * le13y;
                                                
                                                delX= sqrt((xc-ibm->cent_x[nbelmt])*(xc-ibm->cent_x[nbelmt])+
                                                                  (yc-ibm->cent_y[nbelmt])*(yc-ibm->cent_y[nbelmt])+
                                                                  (zc-ibm->cent_z[nbelmt])*(zc-ibm->cent_z[nbelmt]));


						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
														if (cof > 1.e-6) delZ = (ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
                                                        else delZ = ibm->cent_z[elmt] - ibm->cent_z[nbelmt];
														
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
							if(upwind_v)tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_y[elmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[nbelmt].y - ibm->Bvel[elmt].y;
                                                             phi2v = ibm->v_grad_x[elmt]*d_x + ibm->v_grad_y[elmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[elmt].y;

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13y-ibm->cent_y[nbelmt])*(e13y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[elmt].y + (1.-fx)* ibm->Bvel[nbelmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[elmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[nbelmt].y;

                                                             tmp =  nfx13 * u_face + nfy13 * v_face;
                                                           }

					          	if(upwind_c)ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ]
                                                                                             -fabs(tmp * ds * ibm->SCont[elmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_13[ elmt ] = -tmp * ds * caver
                                                                                              -fabs(tmp * ds * caver)*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_y[elmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13y-ibm->cent_y[nbelmt])*(e13y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_13[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
							else	//the velocity is in negative direction		
							{
							if(upwind_v)tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].y * nfy13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_y[nbelmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[elmt].y - ibm->Bvel[nbelmt].y;
                                                             phi2v = ibm->v_grad_x[nbelmt]*d_x + ibm->v_grad_y[nbelmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[nbelmt].y;

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13y-ibm->cent_y[elmt])*(e13y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[nbelmt].y + (1.-fx)* ibm->Bvel[elmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[nbelmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[elmt].y;

                                                             tmp =  nfx13 * u_face + nfy13 * v_face;
                                                           }

							if(upwind_c)ibm->F_flux_13[elmt] = -tmp * ds * ibm->SCont[nbelmt]
                                                                                           -fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_13[elmt]= -tmp * ds * caver
                                                                                           -fabs(tmp * ds * caver)*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_y[nbelmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13y-ibm->cent_y[elmt])*(e13y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_13[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
														if (cof > 1.e-6) delZ = (ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
                                                        else delZ = ibm->cent_z[elmt] - ibm->cent_z[nbelmt];

							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
							if(upwind_v)tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].y * nfy13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_y[nbelmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[elmt].y - ibm->Bvel[nbelmt].y;
                                                             phi2v = ibm->v_grad_x[nbelmt]*d_x + ibm->v_grad_y[nbelmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[nbelmt].y;

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13y-ibm->cent_y[elmt])*(e13y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[nbelmt].y + (1.-fx)* ibm->Bvel[elmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[nbelmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[elmt].y;

                                                             tmp =  nfx13 * u_face + nfy13 * v_face;
                                                           }

							if(upwind_c)ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ]
                                                                                           + fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_13[ elmt ] = tmp * ds * caver 
                                                                                           + fabs(tmp * ds * caver )*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_y[nbelmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e13x-ibm->cent_x[elmt])*(e13x-ibm->cent_x[elmt])+
                                                                      (e13y-ibm->cent_y[elmt])*(e13y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_13[elmt]= tmp * ds * c_face 
                                                                  +fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
							else	//the velocity is in negative direction	
                                 			{
							if(upwind_v)tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
                                                        if(central_v)tmp = nor_vel_to_face13;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_y[elmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u + 1.e-8);

                                                             phi1v = ibm->Bvel[nbelmt].y - ibm->Bvel[elmt].y;
                                                             phi2v = ibm->v_grad_x[elmt]*d_x + ibm->v_grad_y[elmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v + 1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[elmt].y;

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13y-ibm->cent_y[nbelmt])*(e13y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[elmt].y + (1.-fx)* ibm->Bvel[nbelmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[elmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[nbelmt].y;

                                                             tmp =  nfx13 * u_face + nfy13 * v_face;
                                                           }

							if(upwind_c)ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ]
                                                                                           + fabs(tmp * ds * ibm->SCont[elmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_13[ elmt ] = tmp * ds * caver
                                                                                            + fabs(tmp * ds * caver)*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_y[elmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c + 1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e13x-ibm->cent_x[nbelmt])*(e13x-ibm->cent_x[nbelmt])+
                                                                      (e13y-ibm->cent_y[nbelmt])*(e13y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_13[elmt]= tmp * ds * c_face 
                                                                                   +fabs(tmp * ds * c_face )*cof*delZ/delX; 
							}}
						}
					}
				}

				//					2---->3 edge  flux

				/* xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c3; // xiaolei add SEDI
				if (nbelmt>=0)  // xiaolei add 
				{
			//		if( ibm->nf_z[ nbelmt ] < 1.e-6 )
	                               	if(ibm->Rigidity[nbelmt] == 1)
					{
						ibm->F_flux_23[ elmt ] = 0.;
					}
					else
					{
						Ave_face23_u = 0.5 * ( ibm->Bvel[ elmt ].x + ibm->Bvel[ nbelmt ].x );

						Ave_face23_v = 0.5 * ( ibm->Bvel[ elmt ].y + ibm->Bvel[ nbelmt ].y );

						Ave_face23_w = 0.5 * ( ibm->Bvel[ elmt ].z + ibm->Bvel[ nbelmt ].z );

                                                caver = 0.5 * ( ibm->SCont[ elmt ] + ibm->SCont[ nbelmt ] );

						nfx23 = dy23;
						nfy23 = -dx23;
						//ds	  = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr = sqrt( dx23 * dx23 + dy23 * dy23 );
                                                ds = dr;
						nfx23 /= dr;
						nfy23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfy23 * Ave_face23_v;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - xc;
						le23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] ) - yc;
						e23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] );
						e23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] );

						le_dot_n23 = nfx23 * le23x + nfy23 * le23y;

                                                delX= sqrt((xc-ibm->cent_x[nbelmt])*(xc-ibm->cent_x[nbelmt])+
                                                                  (yc-ibm->cent_y[nbelmt])*(yc-ibm->cent_y[nbelmt])+
                                                                  (zc-ibm->cent_z[nbelmt])*(zc-ibm->cent_z[nbelmt]));


						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
														if (cof > 1.e-6) delZ = (ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
                                                        else delZ = ibm->cent_z[elmt] - ibm->cent_z[nbelmt];
														
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
 							if(upwind_v) tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_y[elmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u +1.e-8);

                                                             phi1v = ibm->Bvel[nbelmt].y - ibm->Bvel[elmt].y;
                                                             phi2v = ibm->v_grad_x[elmt]*d_x + ibm->v_grad_y[elmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[elmt].y;

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23y-ibm->cent_y[nbelmt])*(e23y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[elmt].y + (1.-fx)* ibm->Bvel[nbelmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[elmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[nbelmt].y;

                                                             tmp =  nfx23 * u_face + nfy23 * v_face;
                                                           }

							if(upwind_c) ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ]
                                                                                             -fabs(tmp * ds * ibm->SCont[ elmt ])*cof*delZ/delX; 
						        if(central_c)ibm->F_flux_23[ elmt ] = -tmp * ds * caver 
                                                                                             -fabs(tmp * ds *caver)*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_y[elmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23y-ibm->cent_y[nbelmt])*(e23y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_23[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
							else	//the velocity is in negative direction				
							{
							if(upwind_v)tmp = ibm->Bvel[nbelmt].x * nfx23 + ibm->Bvel[ nbelmt ].y * nfy23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_y[nbelmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u+1.e-8);

                                                             phi1v = ibm->Bvel[elmt].y - ibm->Bvel[nbelmt].y;
                                                             phi2v = ibm->v_grad_x[nbelmt]*d_x + ibm->v_grad_y[nbelmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[nbelmt].y;

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23y-ibm->cent_y[elmt])*(e23y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[nbelmt].y + (1.-fx)* ibm->Bvel[elmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[nbelmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[elmt].y;

                                                             tmp =  nfx23 * u_face + nfy23 * v_face;
                                                           }

							if(upwind_c)ibm->F_flux_23[elmt]= -tmp * ds * ibm->SCont[ nbelmt ]
                                                                                          -fabs(tmp * ds * ibm->SCont[ nbelmt ])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_23[elmt]= -tmp * ds * caver 
                                                                                           -fabs(tmp * ds * caver )*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_y[nbelmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23y-ibm->cent_y[elmt])*(e23y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_23[elmt]= -tmp * ds * c_face 
                                                                  -fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
														if (cof > 1.e-6) delZ = (ibm->cent_z[elmt]-ibm->cent_z_AVE[elmt])-(ibm->cent_z[nbelmt]-ibm->cent_z_AVE[nbelmt]);
                                                        else delZ = ibm->cent_z[elmt] - ibm->cent_z[nbelmt];
														
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
							if(upwind_v) tmp = ibm->Bvel[nbelmt].x * nfx23 + ibm->Bvel[ nbelmt ].y * nfy23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1u = ibm->Bvel[elmt].x - ibm->Bvel[nbelmt].x;
                                                             phi2u = ibm->u_grad_x[nbelmt]*d_x + ibm->u_grad_y[nbelmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u+1.e-8);

                                                             phi1v = ibm->Bvel[elmt].y - ibm->Bvel[nbelmt].y;
                                                             phi2v = ibm->v_grad_x[nbelmt]*d_x + ibm->v_grad_y[nbelmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[nbelmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[nbelmt].y;

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23y-ibm->cent_y[elmt])*(e23y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[nbelmt].x + (1.-fx)* ibm->Bvel[elmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[nbelmt].y + (1.-fx)* ibm->Bvel[elmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[nbelmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[elmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[nbelmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[elmt].y;

                                                             tmp =  nfx23 * u_face + nfy23 * v_face;
                                                           }

							if(upwind_c) ibm->F_flux_23[elmt]= tmp * ds * ibm->SCont[ nbelmt ]
                                                                                           +fabs(tmp * ds * ibm->SCont[nbelmt])*cof*delZ/delX; 
						        if(central_c)ibm->F_flux_23[elmt] = tmp * ds * caver 
                                                                                            +fabs(tmp * ds * caver )*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[elmt]- ibm->cent_x[nbelmt];
                                                             d_y = ibm->cent_y[elmt]- ibm->cent_y[nbelmt];
                                                             
                                                             phi1c = ibm->SCont[elmt] - ibm->SCont[nbelmt];
                                                             phi2c = ibm->c_grad_x[nbelmt]*d_x + ibm->c_grad_y[nbelmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[nbelmt];

                                                             fR = sqrt((e23x-ibm->cent_x[elmt])*(e23x-ibm->cent_x[elmt])+
                                                                      (e23y-ibm->cent_y[elmt])*(e23y-ibm->cent_y[elmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[nbelmt] + (1.-fx)* ibm->SCont[elmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[nbelmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[elmt];

						             ibm->F_flux_23[elmt]= tmp * ds * c_face 
                                                                  +fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							}
							else	//the velocity is in negative direction					    
							{
							if(upwind_v) tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
                                                        if(central_v)tmp = nor_vel_to_face23;
                                                        if(gamma_v)
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1u = ibm->Bvel[nbelmt].x - ibm->Bvel[elmt].x;
                                                             phi2u = ibm->u_grad_x[elmt]*d_x + ibm->u_grad_y[elmt]*d_y;
                                                             phiu  = 1.-phi1u/(2.*phi2u+1.e-8);

                                                             phi1v = ibm->Bvel[nbelmt].y - ibm->Bvel[elmt].y;
                                                             phi2v = ibm->v_grad_x[elmt]*d_x + ibm->v_grad_y[elmt]*d_y;
                                                             phiv  = 1.-phi1v/(2.*phi2v+1.e-8);
                                                             
                                                             if(phiu <= 0. || phiu >= 1.) u_face = ibm->Bvel[elmt].x;
                                                             if(phiv <= 0. || phiv >= 1.) v_face = ibm->Bvel[elmt].y;

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23y-ibm->cent_y[nbelmt])*(e23y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landau = phiu/beta; 
                                                             landav = phiv/beta;
 
                                                             if(phiu < 1. && phiu >= beta) u_face = fx * ibm->Bvel[elmt].x + (1.-fx)* ibm->Bvel[nbelmt].x;
                                                             if(phiv < 1. && phiv >= beta) v_face = fx * ibm->Bvel[elmt].y + (1.-fx)* ibm->Bvel[nbelmt].y;
                                                             
                                                             if(phiu < beta && phiu > 0.) u_face = (1.-landau*(1.-fx)) * ibm->Bvel[elmt].x +
                                                                                                  (1.-fx)*landau * ibm->Bvel[nbelmt].x;
                                                             if(phiv < beta && phiv > 0.) v_face = (1.-landav*(1.-fx)) * ibm->Bvel[elmt].y +
                                                                                                  (1.-fx)*landav * ibm->Bvel[nbelmt].y;

                                                             tmp =  nfx23 * u_face + nfy23 * v_face;
                                                           }

							if(upwind_c) ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ]
                                                                                              +fabs(tmp * ds * ibm->SCont[elmt])*cof*delZ/delX; 
							if(central_c)ibm->F_flux_23[ elmt ] = tmp * ds * caver
                                                                                              +fabs(tmp * ds * caver )*cof*delZ/delX; 
                                                        if(gamma_c) 
                                                          {
                                                             d_x = ibm->cent_x[nbelmt]- ibm->cent_x[elmt];
                                                             d_y = ibm->cent_y[nbelmt]- ibm->cent_y[elmt];
                                                             
                                                             phi1c = ibm->SCont[nbelmt] - ibm->SCont[elmt];
                                                             phi2c = ibm->c_grad_x[elmt]*d_x + ibm->c_grad_y[elmt]*d_y;
                                                             phic  = 1.-phi1c/(2.*phi2c+1.e-8);
                                                            
                                                             if(phic <= 0. || phic >= 1.) c_face = ibm->SCont[elmt];

                                                             fR = sqrt((e23x-ibm->cent_x[nbelmt])*(e23x-ibm->cent_x[nbelmt])+
                                                                      (e23y-ibm->cent_y[nbelmt])*(e23y-ibm->cent_y[nbelmt]));
                                                             RL = sqrt((ibm->cent_x[elmt]-ibm->cent_x[nbelmt])*(ibm->cent_x[elmt]-ibm->cent_x[nbelmt])+
                                                                       (ibm->cent_y[elmt]-ibm->cent_y[nbelmt])*(ibm->cent_y[elmt]-ibm->cent_y[nbelmt]));
                                                             fx = fR/RL;
 
                                                             landac = phic/beta; 
 
                                                             if(phic < 1. && phic >= beta) c_face = fx * ibm->SCont[elmt] + (1.-fx)* ibm->SCont[nbelmt];
                                                             
                                                             if(phic < beta && phic > 0.) c_face = (1.-landac*(1.-fx)) * ibm->SCont[elmt] +
                                                                                                  (1.-fx)*landac * ibm->SCont[nbelmt];

						             ibm->F_flux_23[elmt]= tmp * ds * c_face 
                                                                                    +fabs(tmp * ds * c_face )*cof*delZ/delX; 
                                                           }
							 }
						}
					}
				}
			}
		}
	}

	return ( 0 );												//	ali completed on 3 nov. 2009				    
}


PetscErrorCode	outlet_sediment_flux_bend( UserCtx * user, IBMNodes * ibm )
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscReal	   sb, 
					sc;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
       	PetscReal	   tmp;
	PetscReal	        	nfx12, 
					nfy12, 
					nfx13, 
					nfy13, 
					nfx23, 
					nfy23; 
	PetscReal                       //y_outlet = 3.41,
                                        //y_plus_x_outlet=29.11,                                            
                                        someX_minus_Y_outlet = 39.17196667;//33.07459,                                            
                                        //someX_plus_Y_inlet = 100.05944;                                            
                                        //someX_plus_Y_inlet = 90.16268;                                            
        PetscReal                       ydif1, ydif2, ydif3,xydif1,xydif2,xydif3;
        PetscReal                       sxydif1,sxydif2,sxydif3;
        PetscReal                       ssxydif1,ssxydif2,ssxydif3;
	PetscReal                       Ave_face12_u,
                                        Ave_face13_u,
                                        Ave_face23_u,
                                        Ave_face12_v, 
                                        Ave_face13_v, 
                                        Ave_face23_v,
										Ave_face12_w, 
                                        Ave_face13_w, 
                                        Ave_face23_w,
                                        nor_vel_to_face12,
                                        nor_vel_to_face13,
                                        nor_vel_to_face23,
                                     	le12x,
                                        le12y,
										le12z,
                                        le13x,
                                        le13y,le13z,
                                        le23x,
                                        le23y,le23z,
                                        le_dot_n12,
                                        le_dot_n13,
                                        le_dot_n23,
                                        xc,yc,zc;
                

// computing fluxes through oulet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
		//	zc       = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
			zc	 = ibm->cent_z[ elmt ];
                        
                        //ydif1=fabs(y_outlet-ibm->y_bp[n1e]);
			//ydif2=fabs(y_outlet-ibm->y_bp[n2e]); 
	                //ydif3=fabs(y_outlet-ibm->y_bp[n3e]);
                      
                        //xydif1=fabs(y_plus_x_outlet-ibm->y_bp[n1e]-ibm->x_bp[n1e]);
                        //xydif2=fabs(y_plus_x_outlet-ibm->y_bp[n2e]-ibm->x_bp[n2e]);
                        //xydif3=fabs(y_plus_x_outlet-ibm->y_bp[n3e]-ibm->x_bp[n3e]);
		     
                        sxydif1=fabs(someX_minus_Y_outlet-(0.85942*ibm->x_bp[n1e]-ibm->y_bp[n1e]));
                        sxydif2=fabs(someX_minus_Y_outlet-(0.85942*ibm->x_bp[n2e]-ibm->y_bp[n2e]));
                        sxydif3=fabs(someX_minus_Y_outlet-(0.85942*ibm->x_bp[n3e]-ibm->y_bp[n3e]));

                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
                        
                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
//    Bend 90 degree                                       
	     // if((ydif1<1.e-6 && ydif2<1.e-6)||(ydif1<1.e-6 && ydif3<1.e-6)||(ydif2<1.e-6 && ydif3<1.e-6)) 


 //   Bend 135 degree
	     // if((xydif1<1.e-6 && xydif2<1.e-6)||(xydif1<1.e-6 && xydif3<1.e-6)||(xydif2<1.e-6 && xydif3<1.e-6))  
 

//   OSL-2013
	       //if((sxydif1<1.e-4 && sxydif2<1.e-4)||(sxydif1<1.e-4 && sxydif3<1.e-4)||(sxydif2<1.e-4 && sxydif3<1.e-4)) 
               // if((ssxydif1<1.e-4 && ssxydif2<1.e-4)||(ssxydif1<1.e-4 && ssxydif3<1.e-4)||(ssxydif2<1.e-4 && sxydif3<1.e-4))  
//   OSL-2016
	      if((sxydif1<1.e-4 && sxydif2<1.e-4)||(sxydif1<1.e-4 && sxydif3<1.e-4)||(sxydif2<1.e-4 && sxydif3<1.e-4)) 
	      {


		if(ibm->Rigidity[elmt])
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		    dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				//
				/* // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	                        if( ibm->Rigidity[nbelmt] == 1)
      					{

					       	Ave_face12_u = ibm->Bvel[ elmt ].x;

						Ave_face12_v = ibm->Bvel[ elmt ].y;

						nfx12 = dy12;
						nfy12 = -dx12;
						//ds    = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr    = sqrt( dx12 * dx12 + dy12 * dy12 );
                                                ds = dr;
						nfx12 /= dr;
						nfy12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfy12 * Ave_face12_v;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n12 = nfx12 * le12x + nfy12 * le12y;


						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
              						}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].y * nfy12;
								ibm->F_flux_12[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction		
                                                        {
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
                                             
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*
				if( nbelmt != elmt   // xiaolei deactivate 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt])
					{

                                                Ave_face13_u = ibm->Bvel[ elmt ].x;

						Ave_face13_v = ibm->Bvel[ elmt ].y;

						nfx13 = dy13;
						nfy13 = -dx13;
						//ds    = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr    = sqrt( dx13 * dx13 + dy13 * dy13 );
                                                ds = dr;
						nfx13 /= dr;
						nfy13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfy13 * Ave_face13_v;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n13 = nfx13 * le13x + nfy13 * le13y;


						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].y * nfy13;
								ibm->F_flux_13[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction	
                                 			{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}

					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                            	if( ibm->Rigidity[nbelmt] == 1)
					{
                                              	Ave_face23_u = ibm->Bvel[ elmt ].x;

						Ave_face23_v = ibm->Bvel[ elmt ].y;

						nfx23 = dy23;
						nfy23 = -dx23;
						//ds    = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr    = sqrt( dx23 * dx23 + dy23 * dy23 );
                                                ds = dr;
						nfx23 /= dr;
						nfy23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfy23 * Ave_face23_v;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - ibm->cent_x[ elmt ];
						le23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] ) - ibm->cent_y[ elmt ];

						le_dot_n23 = nfx23 * le23x + nfy23 * le23y;


						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction				
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
						
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling

	return ( 0 );		//	ali completed on 19 nov. 2009				    
}


PetscErrorCode	outlet_sediment_flux_contra( UserCtx * user, IBMNodes * ibm )
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscReal	   sb, 
					sc;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
       	PetscReal	   tmp;
	PetscReal	        	nfx12, 
					nfy12, 
					nfx13, 
					nfy13, 
					nfx23, 
					nfy23; 
	PetscReal                       x_outlet = 165.0;//739.9;//90.0; //total osl//25.; // Jhoke Crossvane Rockvane 50.;//49.71;
	PetscReal                       y_outlet = 85.1; 
        PetscReal                       xdif1,xdif2,xdif3;
        PetscReal                       ydif1,ydif2,ydif3;
	PetscReal                       Ave_face12_u,
                                        Ave_face13_u,
                                        Ave_face23_u,
                                        Ave_face12_v, 
                                        Ave_face13_v, 
                                        Ave_face23_v,
                                        nor_vel_to_face12,
                                        nor_vel_to_face13,
                                        nor_vel_to_face23,
                                     	le12x,
                                        le12y,
                                        le13x,
                                        le13y,
                                        le23x,
                                        le23y,
                                        le_dot_n12,
                                        le_dot_n13,
                                        le_dot_n23;

// computing fluxes through oulet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                                                                   
                        xdif1=fabs(x_outlet-ibm->x_bp[n1e]);
			xdif2=fabs(x_outlet-ibm->x_bp[n2e]); 
	                xdif3=fabs(x_outlet-ibm->x_bp[n3e]);

                        ydif1=fabs(y_outlet-ibm->y_bp[n1e]);
			ydif2=fabs(y_outlet-ibm->y_bp[n2e]); 
	                ydif3=fabs(y_outlet-ibm->y_bp[n3e]);

	      if((xdif1<1.e-4 && xdif2<1.e-4)||(xdif1<1.e-4 && xdif3<1.e-4)||(xdif2<1.e-4 && xdif3<1.e-4))  
	      //if((ydif1<1.e-4 && ydif2<1.e-4)||(ydif1<1.e-4 && ydif3<1.e-4)||(ydif2<1.e-4 && ydif3<1.e-4))  
	      {


		if( ibm->Rigidity[elmt] == 1)
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		        dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ ) // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/* // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	                        if( ibm->Rigidity[nbelmt] == 1)
      					{

					       	Ave_face12_u = ibm->Bvel[ elmt ].x;

						Ave_face12_v = ibm->Bvel[ elmt ].y;

						nfx12 = dy12;
						nfy12 = -dx12;
						//ds    = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr    = sqrt( dx12 * dx12 + dy12 * dy12 );
                                                ds = dr;
						nfx12 /= dr;
						nfy12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfy12 * Ave_face12_v;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n12 = nfx12 * le12x + nfy12 * le12y;


						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
              						}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].y * nfy12;
								ibm->F_flux_12[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
								tmp = ibm->Bvel[ nbelmt ].x * nfx12 + ibm->Bvel[ nbelmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction		
                                                        {
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
                                             
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*   // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if( ibm->Rigidity[nbelmt] == 1)
					{

                                                Ave_face13_u = ibm->Bvel[ elmt ].x;

						Ave_face13_v = ibm->Bvel[ elmt ].y;

						nfx13 = dy13;
						nfy13 = -dx13;
						//ds    = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr    = sqrt( dx13 * dx13 + dy13 * dy13 );
                                                ds = dr;
						nfx13 /= dr;
						nfy13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfy13 * Ave_face13_v;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n13 = nfx13 * le13x + nfy13 * le13y;


						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].y * nfy13;
								ibm->F_flux_13[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
								tmp = ibm->Bvel[ nbelmt ].x * nfx13 + ibm->Bvel[ nbelmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction	
                                 			{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}

					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                            	if( ibm->Rigidity[nbelmt] == 1)
					{
                                              	Ave_face23_u = ibm->Bvel[ elmt ].x;

						Ave_face23_v = ibm->Bvel[ elmt ].y;

						nfx23 = dy23;
						nfy23 = -dx23;
						//ds    = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr    = sqrt( dx23 * dx23 + dy23 * dy23 );
                                                ds = dr;
						nfx23 /= dr;
						nfy23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfy23 * Ave_face23_v;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - ibm->cent_x[ elmt ];
						le23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] ) - ibm->cent_y[ elmt ];

						le_dot_n23 = nfx23 * le23x + nfy23 * le23y;


						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction				
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= -tmp * ds * ibm->SCont[ nbelmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ nbelmt ].x * nfx23 + ibm->Bvel[ nbelmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= tmp * ds * ibm->SCont[ nbelmt ];
							}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
						
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the x_bp value if at outlet
   }  // for loop, elmt cycling

	return ( 0 );		//	ali completed on 4 Dec. 2009				    
}

//Hossein added from NSF Petsc3.1
PetscErrorCode	inlet_sediment_flux_contra_TexWash( UserCtx * user, IBMNodes * ibm )
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscReal	   sb, 
					sc;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
       	PetscReal	   tmp;
	PetscReal	        	nfx12, 
					nfy12, 
					nfx13, 
					nfy13, 
					nfx23, 
					nfy23; 
        PetscReal                       xdif1,xdif2,xdif3;
	PetscReal                       ydif1,ydif2,ydif3;
	PetscReal                       Ave_face12_u,
                                        Ave_face13_u,
                                        Ave_face23_u,
                                        Ave_face12_v, 
                                        Ave_face13_v, 
                                        Ave_face23_v,
                                        nor_vel_to_face12,
                                        nor_vel_to_face13,
                                        nor_vel_to_face23,
                                     	le12x,
                                        le12y,
                                        le13x,
                                        le13y,
                                        le23x,
                                        le23y,
                                        le_dot_n12,
                                        le_dot_n13,
                                        le_dot_n23;

// computing fluxes through oulet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
    if(!X_Limit_Inlet) x_limit_inlet = 1.e8; 
    if(!Y_Limit_Inlet) y_limit_inlet = 1.e8; 

	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                                                                   
            xdif1=fabs(x_inlet-ibm->x_bp[n1e]);
			xdif2=fabs(x_inlet-ibm->x_bp[n2e]); 
	        xdif3=fabs(x_inlet-ibm->x_bp[n3e]);

            ydif1=fabs(y_inlet-ibm->y_bp[n1e]);
			ydif2=fabs(y_inlet-ibm->y_bp[n2e]); 
	        ydif3=fabs(y_inlet-ibm->y_bp[n3e]);

if((XInlet && ((xdif1<1.e-4 && xdif2<1.e-4)||(xdif1<1.e-4 && xdif3<1.e-4)||(xdif2<1.e-4 && xdif3<1.e-4)) && ibm->cent_y[ elmt ] < y_limit_inlet) || 
   (YInlet && ((ydif1<1.e-4 && ydif2<1.e-4)||(ydif1<1.e-4 && ydif3<1.e-4)||(ydif2<1.e-4 && ydif3<1.e-4)) && ibm->cent_x[ elmt ] < x_limit_inlet))
	      {

		  if( ibm->nf_z[elmt] < 1.e-6 || ibm->elmt_depth[elmt] > 5. || ibm->Rigidity[elmt])
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		        dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ ) // xiaolei deactivate 
			{
				//					1---->2 edge flux
				/* // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	                        if( ibm->nf_z[nbelmt] < 1.e-6 || ibm->elmt_depth[nbelmt] > 5. || ibm->Rigidity[nbelmt])
      					{

					       	Ave_face12_u = ibm->Bvel[ elmt ].x;

						Ave_face12_v = ibm->Bvel[ elmt ].y;

						nfx12 = dy12;
						nfy12 = -dx12;
						//ds    = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr    = sqrt( dx12 * dx12 + dy12 * dy12 );
                                                ds = dr;
						nfx12 /= dr;
						nfy12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfy12 * Ave_face12_v;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n12 = nfx12 * le12x + nfy12 * le12y;


                                                
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
              						}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ]
								= -tmp * ds * ibm->SCont[ elmt ];
							}
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*   // xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if( ibm->nf_z[nbelmt] < 1.e-6 || ibm->elmt_depth[nbelmt] > 5. || ibm->Rigidity[nbelmt])
					{

                                                Ave_face13_u = ibm->Bvel[ elmt ].x;

						Ave_face13_v = ibm->Bvel[ elmt ].y;

						nfx13 = dy13;
						nfy13 = -dx13;
						//ds    = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr    = sqrt( dx13 * dx13 + dy13 * dy13 );
                                                ds = dr;
						nfx13 /= dr;
						nfy13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfy13 * Ave_face13_v;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n13 = nfx13 * le13x + nfy13 * le13y;


							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ]
								= -tmp * ds * ibm->SCont[ elmt ];
							}

					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  xiaolei deactivate 
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                            	if( ibm->nf_z[nbelmt] < 1.e-6 || ibm->elmt_depth[nbelmt] > 5. || ibm->Rigidity[nbelmt])
					{
                                              	Ave_face23_u = ibm->Bvel[ elmt ].x;

						Ave_face23_v = ibm->Bvel[ elmt ].y;

						nfx23 = dy23;
						nfy23 = -dx23;
						//ds    = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr    = sqrt( dx23 * dx23 + dy23 * dy23 );
                                                ds = dr;
						nfx23 /= dr;
						nfy23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfy23 * Ave_face23_v;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - ibm->cent_x[ elmt ];
						le23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] ) - ibm->cent_y[ elmt ];

						le_dot_n23 = nfx23 * le23x + nfy23 * le23y;


							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction				
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= -tmp * ds * ibm->SCont[ elmt ];
							}
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the x_bp value if at outlet
   }  // for loop, elmt cycling

	return ( 0 );					    
}


PetscErrorCode	inlet_sediment_flux_osl( UserCtx * user, IBMNodes * ibm )
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	Cmpnts		 * Bvel;
	PetscReal	 * F_flux_12, 
				 * F_flux_13, 
				 * F_flux_23;
	PetscReal	 * SCont;
	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscReal	   sb, 
					sc;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   nfx, 
					nfy, 
					nfz, 
					dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					dr, 
					ds;
       	PetscReal	   tmp;
	PetscReal	        	nfx12, 
					nfy12, 
					nfx13, 
					nfy13, 
					nfx23, 
					nfy23; 
        PetscReal                       Y_plus_someX_inlet = 79.59488;                                            
        PetscReal                       ydif1, ydif2, ydif3,xydif1,xydif2,xydif3;
        PetscReal                       sxydif1,sxydif2,sxydif3;
        PetscReal                       ssxydif1,ssxydif2,ssxydif3;
	PetscReal                       Ave_face12_u,
                                        Ave_face13_u,
                                        Ave_face23_u,
                                        Ave_face12_v, 
                                        Ave_face13_v, 
                                        Ave_face23_v,
                                        nor_vel_to_face12,
                                        nor_vel_to_face13,
                                        nor_vel_to_face23,
                                     	le12x,
                                        le12y,
                                        le13x,
                                        le13y,
                                        le23x,
                                        le23y,
                                        le_dot_n12,
                                        le_dot_n13,
                                        le_dot_n23,
                                        xc,yc,zc;
                

// computing fluxes through oulet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
		//	zc       = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
			zc	 = ibm->cent_z[ elmt ];
                        
                        //ydif1=fabs(y_outlet-ibm->y_bp[n1e]);
			//ydif2=fabs(y_outlet-ibm->y_bp[n2e]); 
	                //ydif3=fabs(y_outlet-ibm->y_bp[n3e]);
                      
                        //xydif1=fabs(y_plus_x_outlet-ibm->y_bp[n1e]-ibm->x_bp[n1e]);
                        //xydif2=fabs(y_plus_x_outlet-ibm->y_bp[n2e]-ibm->x_bp[n2e]);
                        //xydif3=fabs(y_plus_x_outlet-ibm->y_bp[n3e]-ibm->x_bp[n3e]);
		     
                        sxydif1=fabs(Y_plus_someX_inlet - (0.60098*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        sxydif2=fabs(Y_plus_someX_inlet - (0.60098*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        sxydif3=fabs(Y_plus_someX_inlet - (0.60098*ibm->x_bp[n3e]+ibm->y_bp[n3e]));

                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
                        
                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
//    Bend 90 degree                                       
	     // if((ydif1<1.e-6 && ydif2<1.e-6)||(ydif1<1.e-6 && ydif3<1.e-6)||(ydif2<1.e-6 && ydif3<1.e-6)) 


 //   Bend 135 degree
	     // if((xydif1<1.e-6 && xydif2<1.e-6)||(xydif1<1.e-6 && xydif3<1.e-6)||(xydif2<1.e-6 && xydif3<1.e-6))  
 

//   OSL-2013
	       //if((sxydif1<1.e-4 && sxydif2<1.e-4)||(sxydif1<1.e-4 && sxydif3<1.e-4)||(sxydif2<1.e-4 && sxydif3<1.e-4)) 
               // if((ssxydif1<1.e-4 && ssxydif2<1.e-4)||(ssxydif1<1.e-4 && ssxydif3<1.e-4)||(ssxydif2<1.e-4 && sxydif3<1.e-4))  
//   OSL-2016
	      if((sxydif1<1.e-4 && sxydif2<1.e-4)||(sxydif1<1.e-4 && sxydif3<1.e-4)||(sxydif2<1.e-4 && sxydif3<1.e-4)) 
	      {


		if(ibm->Rigidity[elmt])
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		        dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				//
				/* // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	                        if(ibm->Rigidity[nbelmt])
      					{

					       	Ave_face12_u = ibm->Bvel[ elmt ].x;

						Ave_face12_v = ibm->Bvel[ elmt ].y;

						nfx12 = dy12;
						nfy12 = -dx12;
						//ds    = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );
						dr    = sqrt( dx12 * dx12 + dy12 * dy12 );
                                                ds = dr;
						nfx12 /= dr;
						nfy12 /= dr;

						nor_vel_to_face12  = nfx12 * Ave_face12_u + nfy12 * Ave_face12_v;

						le12x = 0.5 * ( ibm->x_bp[ n2e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le12y = 0.5 * ( ibm->y_bp[ n2e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n12 = nfx12 * le12x + nfy12 * le12y;

        //PetscPrintf(PETSC_COMM_WORLD, "entering the loop  \n");

						if( le_dot_n12 > 0. )	//  the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
              						}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ]
								= -tmp * ds * ibm->SCont[ elmt ];
							}
						}
						else    	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"		
						{
							if( nor_vel_to_face12 > 0. )	//the vel is in positive direction    
                                                      	{
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
                                                        {
								tmp = ibm->Bvel[ elmt ].x * nfx12 + ibm->Bvel[ elmt ].y * nfy12;
								ibm->F_flux_12[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
                                             
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*
				if( nbelmt != elmt   // xiaolei deactivate 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt])
					{

                                                Ave_face13_u = ibm->Bvel[ elmt ].x;

						Ave_face13_v = ibm->Bvel[ elmt ].y;

						nfx13 = dy13;
						nfy13 = -dx13;
						//ds    = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );
						dr    = sqrt( dx13 * dx13 + dy13 * dy13 );
                                                ds = dr;
						nfx13 /= dr;
						nfy13 /= dr;

						nor_vel_to_face13  = nfx13 * Ave_face13_u + nfy13 * Ave_face13_v;

						le13x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n1e ] ) - ibm->cent_x[ elmt ];
						le13y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n1e ] ) - ibm->cent_y[ elmt ];

						le_dot_n13 = nfx13 * le13x + nfy13 * le13y;

        //PetscPrintf(PETSC_COMM_WORLD, "entering the loop  \n");

						if( le_dot_n13 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction		
							{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ]
								= -tmp * ds * ibm->SCont[ elmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell  "nbelmt"
                                                {
							if( nor_vel_to_face13 > 0. )	//the vel is in positive direction    
                                                        {
								tmp = ibm->Bvel[elmt].x * nfx13 + ibm->Bvel[elmt].y * nfy13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[elmt];
							}
							else	//the velocity is in negative direction	
                                 			{
								tmp = ibm->Bvel[ elmt ].x * nfx13 + ibm->Bvel[ elmt ].y * nfy13;
								ibm->F_flux_13[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}

					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                            	if(ibm->Rigidity[nbelmt])
					{
                                              	Ave_face23_u = ibm->Bvel[ elmt ].x;

						Ave_face23_v = ibm->Bvel[ elmt ].y;

						nfx23 = dy23;
						nfy23 = -dx23;
						//ds    = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );
						dr    = sqrt( dx23 * dx23 + dy23 * dy23 );
                                                ds = dr;
						nfx23 /= dr;
						nfy23 /= dr;

						nor_vel_to_face23  = nfx23 * Ave_face23_u + nfy23 * Ave_face23_v;

						le23x = 0.5 * ( ibm->x_bp[ n3e ] + ibm->x_bp[ n2e ] ) - ibm->cent_x[ elmt ];
						le23y = 0.5 * ( ibm->y_bp[ n3e ] + ibm->y_bp[ n2e ] ) - ibm->cent_y[ elmt ];

						le_dot_n23 = nfx23 * le23x + nfy23 * le23y;

        //PetscPrintf(PETSC_COMM_WORLD, "entering the loop  \n");

						if( le_dot_n23 > 0. )	// the current cell "elmt" is on the L side or U/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = -tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction				
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= -tmp * ds * ibm->SCont[ elmt ];
							}
						}
						else	// the current cell "elmt" is on the R side or D/S of neighbor cell "nbelmt"
                                                {
							if( nor_vel_to_face23 > 0. )	//the vel is in positive direction    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ]
								= tmp * ds * ibm->SCont[ elmt ];
							}
							else	//the velocity is in negative direction					    
							{
								tmp = ibm->Bvel[ elmt ].x * nfx23 + ibm->Bvel[ elmt ].y * nfy23;
								ibm->F_flux_23[ elmt ] = tmp * ds * ibm->SCont[ elmt ];
							}
						}
						
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling
	
return ( 0 );		//	ali Feb 2016				    
}


PetscErrorCode	inlet_sediment_flux_bend( UserCtx * user, IBMNodes * ibm, double Q_total)
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					ds;
        PetscReal     channle_width = 7.7118;   
	PetscReal                       //y_outlet = 3.41,
                                        //y_plus_x_outlet=29.11,                                            
                                        //someX_minus_Y_outlet = 33.07459,                                            
                                        //someX_plus_Y_inlet = 100.05944;                                            
                                        someX_plus_Y_inlet = 90.16268;                                            
        PetscReal                       ssxydif1,ssxydif2,ssxydif3;
        PetscReal                       xc,yc,zc;
                

// prescribing inlet fluxes at inlet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
			zc	 = ibm->cent_z[ elmt ];
                        
                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
                        ssxydif1=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        ssxydif2=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        ssxydif3=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
//    Bend 90 degree                                       
	     // if((ydif1<1.e-6 && ydif2<1.e-6)||(ydif1<1.e-6 && ydif3<1.e-6)||(ydif2<1.e-6 && ydif3<1.e-6)) 


 //   Bend 135 degree
	     // if((xydif1<1.e-6 && xydif2<1.e-6)||(xydif1<1.e-6 && xydif3<1.e-6)||(xydif2<1.e-6 && xydif3<1.e-6))  
 

//   OSL inlet
                if((ssxydif1<1.e-4 && ssxydif2<1.e-4)||(ssxydif1<1.e-4 && ssxydif3<1.e-4)||(ssxydif2<1.e-4 && ssxydif3<1.e-4))  
	      {


		if(ibm->Rigidity[elmt])
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		        dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  xiaolei deactivate 
			{
				//					1---->2 edge flux
				/*   // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt] && ibm->nf_z[nbelmt]>0.8)
      					{
						ds    = sqrt( dx12 * dx12 + dy12 * dy12 );
						ibm->F_flux_12[ elmt ] += (ds/channle_width) * Q_total;
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/*  // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt] && ibm->nf_z[nbelmt]>0.8)
					{

						ds    = sqrt( dx13 * dx13 + dy13 * dy13 );
						ibm->F_flux_13[ elmt ] += (ds/channle_width) * Q_total;
					}
					else
               				{
					}
				}

				//					2---->3 edge  flux

				/* // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/

				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt] && ibm->nf_z[nbelmt]>0.8)
					{
						ds    = sqrt( dx23 * dx23 + dy23 * dy23 );
						ibm->F_flux_23[ elmt ] += (ds/channle_width) * Q_total;
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling

	return ( 0 );		//	ali completed on 19 nov. 2009				    
}


PetscErrorCode	inlet_sediment_flux_contra( UserCtx * user, IBMNodes * ibm, double Q_total)
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	IBMInfo  	 * ibminfo;
	IBMListNode  * current;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
	PetscReal	   riter;
	PetscReal	   dx12, 
					dy12, 
					dz12, 
					dx13, 
					dy13, 
					dz13, 
					dx23, 
					dy23, 
					dz23, 
					ds;
        PetscReal     channle_width = 2.75;   
	PetscReal                       //y_outlet = 3.41,
                                        //y_plus_x_outlet=29.11,                                            
                                        //someX_minus_Y_outlet = 33.07459,                                            
                                        //someX_plus_Y_inlet = 100.05944;                                            
                                        someX_plus_Y_inlet = 90.16268;                                            
        PetscReal                       ssxydif1,ssxydif2,ssxydif3;
        PetscReal                       xc,yc,zc;
                

// prescribing inlet fluxes at inlet faces which could be 1->2,3 or 2->3
// the outlet charachteristics is that the y is 3.41 there, later I must find a better characteristics for this, 18 Nov. 2009, ali
      
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        xc	 = ibm->cent_x[ elmt ];
			yc	 = ibm->cent_y[ elmt ];
			zc	 = ibm->cent_z[ elmt ];
                        
                        //ssxydif1=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        //ssxydif2=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        //ssxydif3=fabs(someX_plus_Y_inlet-(0.95754*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
                        ssxydif1=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n1e]+ibm->y_bp[n1e]));
                        ssxydif2=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n2e]+ibm->y_bp[n2e]));
                        ssxydif3=fabs(someX_plus_Y_inlet-(0.91932*ibm->x_bp[n3e]+ibm->y_bp[n3e]));
//    Bend 90 degree                                       
	     // if((ydif1<1.e-6 && ydif2<1.e-6)||(ydif1<1.e-6 && ydif3<1.e-6)||(ydif2<1.e-6 && ydif3<1.e-6)) 


 //   Bend 135 degree
	     // if((xydif1<1.e-6 && xydif2<1.e-6)||(xydif1<1.e-6 && xydif3<1.e-6)||(xydif2<1.e-6 && xydif3<1.e-6))  
 

//   OSL inlet
                if((ssxydif1<1.e-4 && ssxydif2<1.e-4)||(ssxydif1<1.e-4 && ssxydif3<1.e-4)||(ssxydif2<1.e-4 && ssxydif3<1.e-4))  
	      {


		if( ibm->Rigidity[elmt] == 1)
		{
			ibm->F_flux_12[ elmt ] = 0.;
			ibm->F_flux_13[ elmt ] = 0.;
			ibm->F_flux_23[ elmt ] = 0.;
		}
		else
		{
		        dx12 = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
			dy12 = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
			dz12 = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

			dx13 = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
			dy13 = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
			dz13 = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

			dx23 = ibm->x_bp[ n3e ] - ibm->x_bp[ n2e ];
			dy23 = ibm->y_bp[ n3e ] - ibm->y_bp[ n2e ];
			dz23 = ibm->z_bp[ n3e ] - ibm->z_bp[ n2e ];

        //PetscPrintf(PETSC_COMM_WORLD, "entering the loop  \n");

// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
	//		for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{
				//					1---->2 edge flux
				//
				/*  // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	      	                        if(ibm->Rigidity[nbelmt] && ibm->nf_z[nbelmt]>0.8)
      					{
						ds    = sqrt( dx12 * dx12 + dy12 * dy12 );
						ibm->F_flux_12[ elmt ] += (ds/channle_width) * Q_total;
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 

				/* // xiaolei deactivate  
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(ibm->Rigidity[nbelmt] && ibm->nf_z[nbelmt]>0.8)
					{

						ds    = sqrt( dx13 * dx13 + dy13 * dy13 );
						ibm->F_flux_13[ elmt ] += (ds/channle_width) * Q_total;
					}
					else
               				{
					}
				}

				//					2---->3 edge  flux

				/*  /// xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
	                            	if(ibm->Rigidity[nbelmt] && ibm->nf_z[nbelmt]>0.8)
					{
						ds    = sqrt( dx23 * dx23 + dy23 * dy23 );
						ibm->F_flux_23[ elmt ] += (ds/channle_width) * Q_total;
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling

	return ( 0 );		//	ali completed on 19 nov. 2009				    
}



PetscErrorCode	PeriodicSedimentFlux( UserCtx * user, IBMNodes * ibm)
{
	DA	 da  = user->da, 
		 fda  = user->fda;

	PetscInt	   n_elmt  = ibm->n_elmt, 
					n_v  = ibm->n_v;
	PetscInt   elmt,nbelmt,	n1e,n2e,n3e;
        PetscReal     channle_width = 30.0;   
	PetscReal                       Xin = 0.,
                                        Xout = 739.9;                                            
        PetscReal                       sxdif1,sxdif2,sxdif3;

// finding the outlet Fluxes (of each outlet cell) to be fed at the corresponding inlet cell, April 4, 2012, ali
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        sxdif1=fabs(Xout - ibm->x_bp[n1e]);
                        sxdif2=fabs(Xout - ibm->x_bp[n2e]);
                        sxdif3=fabs(Xout - ibm->x_bp[n3e]);
// Outlet Cells Only 
                if((sxdif1<1.e-4 && sxdif2<1.e-4)||(sxdif1<1.e-4 && sxdif3<1.e-4)||(sxdif2<1.e-4 && sxdif3<1.e-4))  
	      {


		if(ibm->Rigidity[elmt] == 1)
		{
		}
		else
		{
// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{

				//					1---->2 edge flux
				/* // xiaolei deactivate 
				
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(( ibm->elmt_depth[nbelmt] > 0.3 || ibm->Rigidity[nbelmt]) && fabs(ibm->cent_x[nbelmt]- Xout)<1.e-4)
      					{
						int Nbr = int (0.5 * (ibm->y_bp[n1e]+ibm->y_bp[n2e]) / (channle_width/Nseg));
						ibm->SedFlux[Nbr] = ibm->F_flux_12[ elmt ];
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 
				/*
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(( ibm->elmt_depth[nbelmt] > 0.3 || ibm->Rigidity[nbelmt]) && fabs(ibm->cent_x[nbelmt]- Xout)<1.e-4)
					{

						int Nbr = int (0.5 * (ibm->y_bp[n1e]+ibm->y_bp[n3e]) / (channle_width/Nseg));
						ibm->SedFlux[Nbr] = ibm->F_flux_13[ elmt ];
					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  // xiaolei deactivate 
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(( ibm->elmt_depth[nbelmt] > 0.3 || ibm->Rigidity[nbelmt]) && fabs(ibm->cent_x[nbelmt]- Xout)<1.e-4)
					{
						int Nbr = int (0.5 * (ibm->y_bp[n2e]+ibm->y_bp[n3e]) / (channle_width/Nseg));
						ibm->SedFlux[Nbr] = ibm->F_flux_23[ elmt ];
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling


// finding the inlet cell and assigning the flux from outlet cells, April 4, 2012, ali
	for( elmt = 0; elmt < n_elmt; elmt++ )
	   {
	         	n1e  = ibm->nv1[ elmt ];
			n2e  = ibm->nv2[ elmt ];
			n3e  = ibm->nv3[ elmt ];
                        
                        sxdif1=fabs(Xin - ibm->x_bp[n1e]);
                        sxdif2=fabs(Xin - ibm->x_bp[n2e]);
                        sxdif3=fabs(Xin - ibm->x_bp[n3e]);
// Outlet Cells Only 
                if((sxdif1<1.e-4 && sxdif2<1.e-4)||(sxdif1<1.e-4 && sxdif3<1.e-4)||(sxdif2<1.e-4 && sxdif3<1.e-4))  
	      {


		if(ibm->Rigidity[elmt] == 1)
		{
		}
		else
		{
// finding neighbor cells and related inflow and outflow sediment flux for current elmt 
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate 
			{

				//					1---->2 edge flux
				/* // xiaolei deactivate
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(( ibm->elmt_depth[nbelmt] > 0.3 || ibm->Rigidity[nbelmt]) && fabs(ibm->cent_x[nbelmt]- Xin)<1.e-4)
      					{
						int Nbr = int (0.5 * (ibm->y_bp[n1e]+ibm->y_bp[n2e]) / (channle_width/Nseg));
						ibm->F_flux_12[ elmt ] = -ibm->SedFlux[Nbr];
					}
					else
					{		
					}
				}

				//					1---->3 edge flux 
				/* xiaolei deactivate
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(( ibm->elmt_depth[nbelmt] > 0.3 || ibm->Rigidity[nbelmt]) && fabs(ibm->cent_x[nbelmt]- Xin)<1.e-5)
					{

						int Nbr = int (0.5 * (ibm->y_bp[n1e]+ibm->y_bp[n3e]) / (channle_width/Nseg));
						ibm->F_flux_13[ elmt ] = -ibm->SedFlux[Nbr];
					}
					else
               				{
					}
				}

				//					2---->3 edge  flux
				/*  // xiaolei deactivate
				if( nbelmt != elmt 
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
				       	if(( ibm->elmt_depth[nbelmt] > 0.3 || ibm->Rigidity[nbelmt]) && fabs(ibm->cent_x[nbelmt]- Xin)<1.e-5)
					{
						int Nbr = int (0.5 * (ibm->y_bp[n2e]+ibm->y_bp[n3e]) / (channle_width/Nseg));
						ibm->F_flux_23[ elmt ] = -ibm->SedFlux[Nbr];
					}
					else
					{
					}
				}
			} // for loop, nbelmt cycling
		}  // if,to check if elmt is on bed  
	 } // if,to check the cell face if at outlet
   }  // for loop, elmt cycling

	return ( 0 );						    
}


PetscErrorCode	check_correct_new_elev_y_direction( UserCtx * user, IBMNodes * ibm, PetscInt ti )
{
	PetscInt		n_elmt	= ibm->n_elmt, 
					 n_v  = ibm->n_v;
	IBMInfo  	 *	ibminfo;
	IBMListNode  *	current;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx12, dy12, dz12, dx13, dy13, dz13, dx23, dy23, dz23, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
        PetscReal  angle_repose= Angle_repose*PI/180.;
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	//	if( ibm->nf_z[ elmt ] < 1.e-6 )
		if( ibm->Rigidity[elmt] == 1)
		{
		}
		else
		{
	        	//xc	= ibm->cent_x[ elmt ];
		        //yc	= ibm->cent_y[ elmt ];
	                //zc      = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
			//zc	= ibm->cent_z[ elmt ];

			n1e = ibm->nv1[ elmt ];
			n2e = ibm->nv2[ elmt ];
			n3e = ibm->nv3[ elmt ];

// finding neighbor cells  
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ ) // xiaolei deactivate

			{
				/* xiaolei deactivate
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
				            )
                                  )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
					//				  abgle between elmt and  nbelmt beyond 1---->2 edge
				//	if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		if( ibm->Rigidity[nbelmt] == 1 || fabs(ibm->Delz[nbelmt]) < 1.e-6)
					{
					}
					else
					{
						dx12 = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
						dy12 = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
						dz12 = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
						dr	  = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );

						angle = asin( dy12 / dr );

                                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
						if( fabs( angle )
							 >= ( 0.99 * angle_repose ) )
						{
                   PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activated for ti: %d, elmt: %d, nbelmt: %d, angle: %le, phi: %le\n",ti, elmt, nbelmt, angle*180./PI, Angle_repose);
							ibm->cent_y[ nbelmt ]
							= ibm->cent_y[ elmt ]
							   + dr * sin( 0.99 * angle_repose ) * sign;
						}
						else
						{
						}

					}
				}

				//				  abgle between elmt and  nbelmt beyond 1---->3 edge  
				/* // xiaolei deactivate
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
			//		if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		if( ibm->Rigidity[nbelmt] == 1 || fabs(ibm->Delz[nbelmt]) < 1.e-6)
					{
					}
					else
					{
						dx13 = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
						dy13 = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
						dz13 = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];

						dr	  = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );

						angle = asin( dy13 / dr );

                                                if(angle<0.){sign=-1.;}else{sign=1.;}

						if( fabs( angle )
							 >= ( 0.99 * angle_repose ) )
						{
                   PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activated for ti: %d, elmt: %d, nbelmt: %d, angle: %le, phi: %le\n",ti, elmt, nbelmt, angle*180./PI, Angle_repose);
							ibm->cent_y[ nbelmt ]
							= ibm->cent_y[ elmt ]
							   + dr * sin( 0.99 * angle_repose ) * sign;
						}
						else
						{
						}
					}
				}

				//				  abgle between elmt and  nbelmt beyond 2---->3 edge 
				//
				/* // xiaolei deactivate
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
			//		if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		if( ibm->Rigidity[nbelmt] == 1 || fabs(ibm->Delz[nbelmt]) < 1.e-6)
					{
					}
					else

					{
						dx23 = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
						dy23 = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
						dz23 = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];

						dr	  = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );

						angle = asin( dy23 / dr );

                                                if(angle<0.){sign=-1.;}else{sign=1.;}

						if( fabs( angle )
							 >= ( 0.99 * angle_repose ) )
						{
                   PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activated for ti: %d, elmt: %d, nbelmt: %d, angle: %le, phi: %le\n",ti, elmt, nbelmt, angle*180./PI, Angle_repose);
							ibm->cent_y[ nbelmt ]
							= ibm->cent_y[ elmt ]
							   + dr * sin( 0.99 * angle_repose ) * sign;
						}
						else
						{
						}

					}
				}
			}
		}
	}

	return ( 0 );												//	ali completed on 3 nov. 2009				    
}


PetscErrorCode	avalanche_first_sweep_y_direction( UserCtx * user, IBMNodes * ibm, PetscInt aval_loop )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;

        if(!aval_loop) return(0);

	PetscInt		n_elmt	= ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx, dy, dz, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
        PetscReal  angle_repose= Angle_repose * PI / 180.;
// check and correct z-direction angle between element's centerpoint
    

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
                 nfx = ibm->nf_x[elmt];
                 nfy = ibm->nf_y[elmt];
                 nfz = ibm->nf_z[elmt];
           
         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt((ibm->nf_z[elmt]/(ibm->nf_y[elmt]+1.e-8))*(ibm->nf_z[elmt]/(ibm->nf_y[elmt]+1.e-8)) 
                                                                                               + (ibm->nf_x[elmt]/(ibm->nf_y[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_y[elmt]+1.e-8))));
         if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if(ibm->Rigidity[elmt] == 0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{
                 ibm->deltz_p_us[elmt] = 0.;
                 ibm->deltz_p_ds[elmt] = 0.;
                 ibm->A_us[elmt] = 0.;
                 ibm->A_ds[elmt] = 0.;

PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_first_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate
			///  xiaolei add
				//int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				int jj;
                                int ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
			// end add

			{

					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt] ==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
	//&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
         //    &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
 //            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
 //            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dz * dz );

				angle = atan( dy / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
        PetscPrintf(PETSC_COMM_WORLD, "neighbour-elmts with angle larger than PHI, angle: %d %le\n",nbelmt, angle*180./PI);
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[elmt] += ibm->dA[nbelmt]*(ibm->cent_y[nbelmt] - ibm->cent_y[elmt] + dr * tan(angle_repose));
                                      ibm->A_us[elmt] +=ibm->dA[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[elmt] += ibm->dA[nbelmt]*(ibm->cent_y[nbelmt] - ibm->cent_y[elmt] - dr * tan(angle_repose));
                                      ibm->A_ds[elmt] +=ibm->dA[nbelmt];}
				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                                      ii ++;
			 }  // neighbour cell  -FOR
                     
                    ibm->deltz_p_us[elmt] /= (ibm->dA[elmt]+ibm->A_us[elmt]); 
                    ibm->deltz_p_ds[elmt] /= (ibm->dA[elmt]+ibm->A_ds[elmt]); 
                 
                         
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate
			///  xiaolei add
	//			int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				 ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
		                 	{
					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	

					//nbelmt=ibm->n2c[_nv[jj]].c[ii];	// xiaolei add
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt]==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
                        
//	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))

			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dz * dz );

				angle = atan( dy / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
                                   if(sign < 0.) {
                                     ibm->deltz_p_us[nbelmt] = ibm->deltz_p_us[elmt] + ibm->cent_y[elmt] - ibm->cent_y[nbelmt] - dr * tan(angle_repose);
                                     ibm->cent_y[nbelmt] += ibm->deltz_p_us[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[nbelmt] = ibm->deltz_p_ds[elmt] + ibm->cent_y[elmt] - ibm->cent_y[nbelmt] + dr * tan(angle_repose);
                                      ibm->cent_y[nbelmt] += ibm->deltz_p_ds[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                               ii ++;
			 }  // neighbour cell  -FOR

                   ibm->cent_y[elmt] += ibm->deltz_p_us[elmt] + ibm->deltz_p_ds[elmt];

		  } // bed elements with big bed-slope  -IF
	} // bed elements -FOR

	return ( 0 ); // July 2, 2010															    
}


PetscErrorCode	avalanche_second_sweep_y_direction( UserCtx * user, IBMNodes * ibm, PetscInt aval_loop )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;

        if(!aval_loop) return(0);

	PetscInt		n_elmt	= ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx, dy, dz, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
        PetscReal  angle_repose= Angle_repose * PI / 180.;
	//PetscReal		deltz_p_us[n_elmt+2], deltz_p_ds[n_elmt+2], A_us[n_elmt+2], A_ds[n_elmt+2];
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = n_elmt-1; elmt >= 0; elmt-- )
	{
                 nfx = ibm->nf_x[elmt];
                 nfy = ibm->nf_y[elmt];
                 nfz = ibm->nf_z[elmt];
         
         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt((ibm->nf_z[elmt]/(ibm->nf_y[elmt]+1.e-8))*(ibm->nf_z[elmt]/(ibm->nf_y[elmt]+1.e-8)) 
                                                                                               + (ibm->nf_x[elmt]/(ibm->nf_y[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_y[elmt]+1.e-8))));
         if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if(ibm->Rigidity[elmt]==0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{
                 ibm->deltz_p_us[elmt] = 0.;
                 ibm->deltz_p_ds[elmt] = 0.;
                 ibm->A_us[elmt] = 0.;
                 ibm->A_ds[elmt] = 0.;

        PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_second_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			// for( nbelmt = n_elmt; nbelmt >= 0; nbelmt-- ) // xiaolei deactivate
		///  xiaolei add
				//int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				int jj;
				int ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
			// end add
			{
					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	

					//nbelmt=ibm->n2c[_nv[jj]].c[ii];	// xiaolei add
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt]==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
//	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dz * dz );

				angle = atan( dy / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
        PetscPrintf(PETSC_COMM_WORLD, "neighbour-elmts with angle larger than PHI, anlgle: %d %le\n",nbelmt,angle*180./PI);
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[elmt] += ibm->dA[nbelmt]*(ibm->cent_y[nbelmt] - ibm->cent_y[elmt] + dr * tan(angle_repose));
                                      ibm->A_us[elmt] +=ibm->dA[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[elmt] += ibm->dA[nbelmt]*(ibm->cent_y[nbelmt] - ibm->cent_y[elmt] - dr * tan(angle_repose));
                                      ibm->A_ds[elmt] +=ibm->dA[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                                      ii ++;
			 }  // neighbour cell  -FOR
                    
 
                    ibm->deltz_p_us[elmt] /= (ibm->dA[elmt]+ibm->A_us[elmt]); 
                    ibm->deltz_p_ds[elmt] /= (ibm->dA[elmt]+ibm->A_ds[elmt]); 
                 
                         
			// for( nbelmt = n_elmt; nbelmt >=0; nbelmt-- )  // xiaolei deactivate 
		///  xiaolei add
//				int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
			// end add
			        {
					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	

					//nbelmt=ibm->n2c[_nv[jj]].c[ii];	// xiaolei add
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt]==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
//	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
//              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dz * dz );

				angle = atan( dy / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[nbelmt] = ibm->deltz_p_us[elmt] + ibm->cent_y[elmt] - ibm->cent_y[nbelmt] - dr * tan(angle_repose);
                                      ibm->cent_y[nbelmt] += ibm->deltz_p_us[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[nbelmt] = ibm->deltz_p_ds[elmt] + ibm->cent_y[elmt] - ibm->cent_y[nbelmt] + dr * tan(angle_repose);
                                      ibm->cent_y[nbelmt] += ibm->deltz_p_ds[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                                      ii ++;
			 }  // neighbour cell  -FOR

                   ibm->cent_y[elmt] += ibm->deltz_p_us[elmt] + ibm->deltz_p_ds[elmt];

		  } // bed elements with big bed-slope  -IF
	} // bed elements -FOR

	return ( 0 ); // July 2, 2010															    
}


PetscErrorCode	check_correct_new_elev( UserCtx * user, IBMNodes * ibm, PetscInt ti )
{
	PetscInt		n_elmt	= ibm->n_elmt, 
					 n_v  = ibm->n_v;
	IBMInfo  	 *	ibminfo;
	IBMListNode  *	current;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx12, dy12, dz12, dx13, dy13, dz13, dx23, dy23, dz23, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
        PetscReal  angle_repose= Angle_repose*PI/180.;
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	//	if( ibm->nf_z[ elmt ] < 1.e-6 )
		if( ibm->Rigidity[elmt] == 1)
		{
		}
		else
		{
	        	//xc	= ibm->cent_x[ elmt ];
		        //yc	= ibm->cent_y[ elmt ];
	                //zc      = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
			//zc	= ibm->cent_z[ elmt ];

			n1e = ibm->nv1[ elmt ];
			n2e = ibm->nv2[ elmt ];
			n3e = ibm->nv3[ elmt ];

// finding neighbor cells  
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ ) // xiaolei deactivate

			{
				/* xiaolei deactivate
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
				            )
                                  )
				*/
				nbelmt=ibm->c2c[elmt].c1; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
					//				  abgle between elmt and  nbelmt beyond 1---->2 edge
				//	if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		if( ibm->Rigidity[nbelmt] == 1 || fabs(ibm->Delz[nbelmt]) < 1.e-6)
					{
					}
					else
					{
						dx12 = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
						dy12 = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
						dz12 = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
						dr	  = sqrt( dx12 * dx12 + dy12 * dy12 + dz12 * dz12 );

						angle = asin( dz12 / dr );

                                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
						if( fabs( angle )
							 >= ( 0.99 * angle_repose ) )
						{
                   PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activated for ti: %d, elmt: %d, nbelmt: %d, angle: %le, phi: %le\n",ti, elmt, nbelmt, angle*180./PI, Angle_repose);
							ibm->cent_z[ nbelmt ]
							= ibm->cent_z[ elmt ]
							   + dr * sin( 0.99 * angle_repose ) * sign;
						}
						else
						{
						}

					}
				}

				//				  abgle between elmt and  nbelmt beyond 1---->3 edge  
				/* // xiaolei deactivate
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n1e
						   || ibm->nv2[ nbelmt ] == n1e
						   || ibm->nv3[ nbelmt ] == n1e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c2; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
			//		if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		if( ibm->Rigidity[nbelmt] == 1 || fabs(ibm->Delz[nbelmt]) < 1.e-6)
					{
					}
					else
					{
						dx13 = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
						dy13 = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
						dz13 = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];

						dr	  = sqrt( dx13 * dx13 + dy13 * dy13 + dz13 * dz13 );

						angle = asin( dz13 / dr );

                                                if(angle<0.){sign=-1.;}else{sign=1.;}

						if( fabs( angle )
							 >= ( 0.99 * angle_repose ) )
						{
                   PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activated for ti: %d, elmt: %d, nbelmt: %d, angle: %le, phi: %le\n",ti, elmt, nbelmt, angle*180./PI, Angle_repose);
							ibm->cent_z[ nbelmt ]
							= ibm->cent_z[ elmt ]
							   + dr * sin( 0.99 * angle_repose ) * sign;
						}
						else
						{
						}
					}
				}

				//				  abgle between elmt and  nbelmt beyond 2---->3 edge 
				//
				/* // xiaolei deactivate
				if( nbelmt != elmt
					 && ( ibm->nv1[ nbelmt ] == n2e
						   || ibm->nv2[ nbelmt ] == n2e
						   || ibm->nv3[ nbelmt ] == n2e
						 )
					 && ( ibm->nv1[ nbelmt ] == n3e
						   || ibm->nv2[ nbelmt ] == n3e
						   || ibm->nv3[ nbelmt ] == n3e
						 )
				   )
				*/
				nbelmt=ibm->c2c[elmt].c3; // xiaolei add
				if(nbelmt>0) // xiaolei add 
				{
			//		if( ibm->nf_z[ nbelmt ] < 1.e-6 )
		if( ibm->Rigidity[nbelmt] == 1 || fabs(ibm->Delz[nbelmt]) < 1.e-6)
					{
					}
					else

					{
						dx23 = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
						dy23 = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
						dz23 = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];

						dr	  = sqrt( dx23 * dx23 + dy23 * dy23 + dz23 * dz23 );

						angle = asin( dz23 / dr );

                                                if(angle<0.){sign=-1.;}else{sign=1.;}

						if( fabs( angle )
							 >= ( 0.99 * angle_repose ) )
						{
                   PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activated for ti: %d, elmt: %d, nbelmt: %d, angle: %le, phi: %le\n",ti, elmt, nbelmt, angle*180./PI, Angle_repose);
							ibm->cent_z[ nbelmt ]
							= ibm->cent_z[ elmt ]
							   + dr * sin( 0.99 * angle_repose ) * sign;
						}
						else
						{
						}

					}
				}
			}
		}
	}

	return ( 0 );												//	ali completed on 3 nov. 2009				    
}


PetscErrorCode	avalanche_first_sweep_old( UserCtx * user, IBMNodes * ibm, PetscInt aval_loop )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;

    if(!aval_loop) return(0);
	PetscInt		n_elmt	= ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx, dy, dz, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
    PetscReal  angle_repose= Angle_repose * PI / 180.;
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
         nfx = ibm->nf_x[elmt];
         nfy = ibm->nf_y[elmt];
         nfz = ibm->nf_z[elmt];

         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt( ( nfy / ( nfz+1.e-8) ) * ( nfy / ( nfz+1.e-8) ) + ( nfx / ( nfz+1.e-8) ) * ( nfx / ( nfz+1.e-8) ) ));
         if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if( ibm->Rigidity[elmt] == 0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{

         ibm->deltz_p_us[elmt] = 0.;
         ibm->deltz_p_ds[elmt] = 0.;
         ibm->A_us[elmt] = 0.;
         ibm->A_ds[elmt] = 0.;
         
       
PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_first_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )
			{

                 	 if(nbelmt != elmt && ibm->Rigidity[nbelmt] == 0
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			    dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8) );

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
        PetscPrintf(PETSC_COMM_WORLD, "neighbour-elmts with angle larger than PHI, angle: %d %le\n",nbelmt, angle*180./PI);
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] + dr * tan(angle_repose));
                                      ibm->A_us[elmt] +=ibm->dA[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] - dr * tan(angle_repose));
                                      ibm->A_ds[elmt] +=ibm->dA[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
			 }  // neighbour cell  -FOR
                    
 
                    ibm->deltz_p_us[elmt] /= (ibm->dA[elmt]+ibm->A_us[elmt]); 
                    ibm->deltz_p_ds[elmt] /= (ibm->dA[elmt]+ibm->A_ds[elmt]); 
                 
                         
			for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )
			{

                 	 if(nbelmt != elmt && ibm->Rigidity[nbelmt]==0
                        //    && (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                        //    || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                        //    || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                        //    || ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                        //   || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
                        //
	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))

			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8) );

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[nbelmt] = ibm->deltz_p_us[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] - dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_us[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[nbelmt] = ibm->deltz_p_ds[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] + dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_ds[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
			 }  // neighbour cell  -FOR

                   ibm->cent_z[elmt] += ibm->deltz_p_us[elmt] + ibm->deltz_p_ds[elmt];

		  } // bed elements with big bed-slope  -IF
	} // bed elements -FOR

	return ( 0 ); // July 2, 2010															    
}


PetscErrorCode	avalanche_second_sweep_old( UserCtx * user, IBMNodes * ibm, PetscInt aval_loop )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;

    if(!aval_loop) return(0);
	PetscInt		n_elmt	= ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx, dy, dz, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
    PetscReal  angle_repose= Angle_repose * PI / 180.;
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = n_elmt-1; elmt >= 0; elmt-- )
	{
         nfx = ibm->nf_x[elmt];
         nfy = ibm->nf_y[elmt];
         nfz = ibm->nf_z[elmt];

         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt( ( nfy / ( nfz+1.e-8) ) * ( nfy / ( nfz+1.e-8) ) + ( nfx / ( nfz+1.e-8) ) * ( nfx / ( nfz+1.e-8) ) ));
         if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if( ibm->Rigidity[elmt]==0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{

         ibm->deltz_p_us[elmt] = 0.;
         ibm->deltz_p_ds[elmt] = 0.;
         ibm->A_us[elmt] = 0.;
         ibm->A_ds[elmt] = 0.;
         
       
        PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_second_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			for( nbelmt = n_elmt-1; nbelmt >= 0; nbelmt-- )
			{

                 	 if(nbelmt != elmt && ibm->Rigidity[nbelmt]==0
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8) );

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
        PetscPrintf(PETSC_COMM_WORLD, "neighbour-elmts with angle larger than PHI, anlgle: %d %le\n",nbelmt,angle*180./PI);
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] + dr * tan(angle_repose));
                                      ibm->A_us[elmt] +=ibm->dA[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] - dr * tan(angle_repose));
                                      ibm->A_ds[elmt] +=ibm->dA[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
			 }  // neighbour cell  -FOR
                    
 
                    ibm->deltz_p_us[elmt] /= (ibm->dA[elmt]+ibm->A_us[elmt]); 
                    ibm->deltz_p_ds[elmt] /= (ibm->dA[elmt]+ibm->A_ds[elmt]); 
                 
                         
			for( nbelmt = n_elmt-1; nbelmt >=0; nbelmt-- )
			{

                 	 if(nbelmt != elmt && ibm->Rigidity[nbelmt]==0
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))

	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8) );

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[nbelmt] = ibm->deltz_p_us[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] - dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_us[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[nbelmt] = ibm->deltz_p_ds[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] + dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_ds[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
			 }  // neighbour cell  -FOR

                   ibm->cent_z[elmt] += ibm->deltz_p_us[elmt] + ibm->deltz_p_ds[elmt];

		  } // bed elements with big bed-slope  -IF
	} // bed elements -FOR

	return ( 0 ); // July 2, 2010															    
}


PetscErrorCode	avalanche_first_sweep( UserCtx * user, IBMNodes * ibm, PetscInt aval_loop )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;

        if(!aval_loop) return(0);

	PetscInt		n_elmt	= ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx, dy, dz, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
        PetscReal  angle_repose= Angle_repose * PI / 180.;
// check and correct z-direction angle between element's centerpoint
    

	for( elmt = 0; elmt < n_elmt; elmt++ )
	{
                 nfx = ibm->nf_x[elmt];
                 nfy = ibm->nf_y[elmt];
                 nfz = ibm->nf_z[elmt];
           
         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt((ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8)) 
                                                                                               + (ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))));
         if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if(ibm->Rigidity[elmt] == 0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{
                 ibm->deltz_p_us[elmt] = 0.;
                 ibm->deltz_p_ds[elmt] = 0.;
                 ibm->A_us[elmt] = 0.;
                 ibm->A_ds[elmt] = 0.;

PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_first_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate
			///  xiaolei add
				//int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				int jj;
                                int ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
			// end add

			{

					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt] ==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
	//&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
         //    &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
 //            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
 //            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
        PetscPrintf(PETSC_COMM_WORLD, "neighbour-elmts with angle larger than PHI, angle: %d %le\n",nbelmt, angle*180./PI);
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] + dr * tan(angle_repose));
                                      ibm->A_us[elmt] +=ibm->dA[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] - dr * tan(angle_repose));
                                      ibm->A_ds[elmt] +=ibm->dA[nbelmt];}
				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                                      ii ++;
			 }  // neighbour cell  -FOR
                    
 
                    ibm->deltz_p_us[elmt] /= (ibm->dA[elmt]+ibm->A_us[elmt]); 
                    ibm->deltz_p_ds[elmt] /= (ibm->dA[elmt]+ibm->A_ds[elmt]); 
                 
                         
			// for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )  // xiaolei deactivate
			///  xiaolei add
	//			int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				 ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
		                 	{
					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	

					//nbelmt=ibm->n2c[_nv[jj]].c[ii];	// xiaolei add
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt]==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
                        
//	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))

			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
                                   if(sign < 0.) {
                                     ibm->deltz_p_us[nbelmt] = ibm->deltz_p_us[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] - dr * tan(angle_repose);
                                     ibm->cent_z[nbelmt] += ibm->deltz_p_us[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[nbelmt] = ibm->deltz_p_ds[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] + dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_ds[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                               ii ++;
			 }  // neighbour cell  -FOR

                   ibm->cent_z[elmt] += ibm->deltz_p_us[elmt] + ibm->deltz_p_ds[elmt];

		  } // bed elements with big bed-slope  -IF
	} // bed elements -FOR

	return ( 0 ); // July 2, 2010															    
}


PetscErrorCode	avalanche_second_sweep( UserCtx * user, IBMNodes * ibm, PetscInt aval_loop )
{
	DA		da  = user->da, fda  = user->fda;
	DALocalInfo     info  = user->info;

        if(!aval_loop) return(0);

	PetscInt		n_elmt	= ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		elmt, nbelmt, n1e, n2e, n3e;
	PetscReal		riter;
	PetscReal		nfx, nfy, nfz, dx, dy, dz, dr, angle;
	PetscReal		xc, yc, zc;
	PetscReal		tmp;
	PetscReal		rtinyn	= 1.e-7, sign;
        PetscReal  angle_repose= Angle_repose * PI / 180.;
	//PetscReal		deltz_p_us[n_elmt+2], deltz_p_ds[n_elmt+2], A_us[n_elmt+2], A_ds[n_elmt+2];
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = n_elmt-1; elmt >= 0; elmt-- )
	{
                 nfx = ibm->nf_x[elmt];
                 nfy = ibm->nf_y[elmt];
                 nfz = ibm->nf_z[elmt];
         
         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt((ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8)) 
                                                                                               + (ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))));
         if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if(ibm->Rigidity[elmt]==0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{
                 ibm->deltz_p_us[elmt] = 0.;
                 ibm->deltz_p_ds[elmt] = 0.;
                 ibm->A_us[elmt] = 0.;
                 ibm->A_ds[elmt] = 0.;

        PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_second_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			// for( nbelmt = n_elmt; nbelmt >= 0; nbelmt-- ) // xiaolei deactivate
		///  xiaolei add
				//int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				int jj;
				int ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
			// end add
			{
					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	

					//nbelmt=ibm->n2c[_nv[jj]].c[ii];	// xiaolei add
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt]==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
//	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//            &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
        PetscPrintf(PETSC_COMM_WORLD, "neighbour-elmts with angle larger than PHI, anlgle: %d %le\n",nbelmt,angle*180./PI);
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] + dr * tan(angle_repose));
                                      ibm->A_us[elmt] +=ibm->dA[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[elmt] += ibm->dA[nbelmt]*(ibm->cent_z[nbelmt] - ibm->cent_z[elmt] - dr * tan(angle_repose));
                                      ibm->A_ds[elmt] +=ibm->dA[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                                      ii ++;
			 }  // neighbour cell  -FOR
                    
 
                    ibm->deltz_p_us[elmt] /= (ibm->dA[elmt]+ibm->A_us[elmt]); 
                    ibm->deltz_p_ds[elmt] /= (ibm->dA[elmt]+ibm->A_ds[elmt]); 
                 
                         
			// for( nbelmt = n_elmt; nbelmt >=0; nbelmt-- )  // xiaolei deactivate 
		///  xiaolei add
//				int _nv[3];
				//_nv[0]=ibm->nv1[elmt];
				//_nv[1]=ibm->nv2[elmt];
				//_nv[2]=ibm->nv3[elmt];
	
				ii = 1;
				for(jj=0;jj<3;jj++)
				//for(ii=0;ii<100;ii++)
			// end add
			        {
					if(ii == 1) nbelmt=ibm->c2c[elmt].c1;	
					if(ii == 2) nbelmt=ibm->c2c[elmt].c2;	
					if(ii == 3) nbelmt=ibm->c2c[elmt].c3;	

					//nbelmt=ibm->n2c[_nv[jj]].c[ii];	// xiaolei add
                 	 if(nbelmt>=0 && nbelmt != elmt && ibm->Rigidity[nbelmt]==0)
                            //&& (ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ]
                            //|| ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ]
                            //|| ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ]== ibm->nv3[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ]== ibm->nv2[ nbelmt ]
                            //|| ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))	
//	&& (((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//             &&(ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv1[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv1[ elmt ] == ibm->nv3[ nbelmt ])
//              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))
//
//	    ||((ibm->nv2[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv2[ elmt ] == ibm->nv3[ nbelmt ])
//              &&(ibm->nv3[ elmt ] == ibm->nv1[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv2[ nbelmt ] || ibm->nv3[ elmt ] == ibm->nv3[ nbelmt ]))))
			       { 
 			        dx = ibm->cent_x[ nbelmt ] - ibm->cent_x[ elmt ];
				dy = ibm->cent_y[ nbelmt ] - ibm->cent_y[ elmt ];
				dz = ibm->cent_z[ nbelmt ] - ibm->cent_z[ elmt ];
				dr   = sqrt( dx * dx + dy * dy );

				angle = atan( dz / (dr + 1.e-8));

                                if(angle<0.){sign=-1.;}else{sign=1.;}
                                                
			//	if( fabs( angle )> angle_repose || sign < 0.)
				  {
                                   if(sign < 0.) {
                                      ibm->deltz_p_us[nbelmt] = ibm->deltz_p_us[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] - dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_us[nbelmt];}

                                   if(sign > 0.){
                                      ibm->deltz_p_ds[nbelmt] = ibm->deltz_p_ds[elmt] + ibm->cent_z[elmt] - ibm->cent_z[nbelmt] + dr * tan(angle_repose);
                                      ibm->cent_z[nbelmt] += ibm->deltz_p_ds[nbelmt];}

				  } // neighbour cells with big angle of connecting line -IF
			       } // neighbour cell -IF
                                      ii ++;
			 }  // neighbour cell  -FOR

                   ibm->cent_z[elmt] += ibm->deltz_p_us[elmt] + ibm->deltz_p_ds[elmt];

		  } // bed elements with big bed-slope  -IF
	} // bed elements -FOR

	return ( 0 ); // July 2, 2010															    
}


PetscErrorCode recomputing_geometry(UserCtx * user, IBMNodes * ibm, PetscInt tistart, PetscInt ti, PetscInt itr_sc, PetscInt avalanche_check_number)
{

	DA                da  = user->da, fda  = user->fda;
	DALocalInfo  	 info  = user->info;
	PetscReal	 max_delz = 1./double(info.my);
	PetscReal	 max_dely = 1./double(info.my);
	PetscInt	         n_elmt  = ibm->n_elmt, n_v  = ibm->n_v;
	PetscInt		 iter, elmt, vert,n1e,n2e,n3e;
	PetscReal		 xc, yc, zc, dx12, dx13, dy12, dy13, dz12, dz13, dr, nfx, nfy, nfz;
        PetscInt itt=0;        
        
        for( elmt = 0; elmt < n_elmt; elmt++ )
	{
                itt++;
		n1e  			  = ibm->nv1[ elmt ];
		n2e  			  = ibm->nv2[ elmt ];
		n3e  			  = ibm->nv3[ elmt ];
		dx12			  = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
		dy12			  = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
		dz12			  = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

		dx13			  = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
		dy13			  = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
		dz13			  = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

		ibm->nf_x[ elmt ] = dy12 * dz13 - dz12 * dy13;
		ibm->nf_y[ elmt ] = -dx12 * dz13 + dz12 * dx13;
		ibm->nf_z[ elmt ] = dx12 * dy13 - dy12 * dx13;

		dr = sqrt(ibm->nf_x[ elmt ] * ibm->nf_x[ elmt ] + ibm->nf_y[ elmt ] * ibm->nf_y[ elmt ]
				+ ibm->nf_z[ elmt ] * ibm->nf_z[ elmt ] );

		ibm->nf_x[ elmt ] /= dr;
		ibm->nf_y[ elmt ] /= dr;
		ibm->nf_z[ elmt ] /= dr;

		nfx = ibm->nf_x[ elmt ];
		nfy = ibm->nf_y[ elmt ];
		nfz = ibm->nf_z[ elmt ];

                if(ibm->Rigidity[elmt]==0 && y_direction) ibm->max_bed_angle[elmt] = atan(sqrt((ibm->nf_z[elmt]/(ibm->nf_y[elmt]+1.e-8))*(ibm->nf_z[elmt]/(ibm->nf_y[elmt]+1.e-8)) 
                                                                                               + (ibm->nf_x[elmt]/(ibm->nf_y[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_y[elmt]+1.e-8))));		
                if(ibm->Rigidity[elmt]==0 && !y_direction) ibm->max_bed_angle[elmt] = atan(sqrt((ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_y[elmt]/(ibm->nf_z[elmt]+1.e-8)) 
                                                                                               + (ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))*(ibm->nf_x[elmt]/(ibm->nf_z[elmt]+1.e-8))));
                if(ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;

		ibm->dA[ elmt ]    = dr / 2.;

                if(SimpleCellCheck)
                {                
                PetscInt steps = sed_tio;
                if((ti==tistart || ti==0) && itt==1 && avalanche_check_number==0) PetscPrintf(PETSC_COMM_WORLD,"Averaging Cz/y onto z/y_bp every: %d \n",steps);
                if(((ti-tistart)/steps)*steps==(ti-tistart) && ibm->Rigidity[elmt]==0)
                   {
        	//ibm->cent_x[ elmt ] = ( ibm->x_bp[ n1e ] + ibm->x_bp[ n2e ] + ibm->x_bp[ n3e ] ) / 3.;
		//ibm->cent_y[ elmt ] = ( ibm->y_bp[ n1e ] + ibm->y_bp[ n2e ] + ibm->y_bp[ n3e ] ) / 3.;
		if(!y_direction){
		                 if(fabs(ibm->Delz[elmt]) > 1.e-6) ibm->cent_z[ elmt ] = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
						} 
						else {
							  if(fabs(ibm->Dely[elmt]) > 1.e-6) ibm->cent_y[ elmt ] = ( ibm->y_bp[ n1e ] + ibm->y_bp[ n2e ] + ibm->y_bp[ n3e ] ) / 3.;
							 }
	
               // PetscPrintf(PETSC_COMM_WORLD,"Averaging Cz onto z_bp done: ti = %d, cell number = %d, avalanche_check_number = %d\n",ti, elmt, avalanche_check_number);
                //    }                                                               		
                //PetscPrintf(PETSC_COMM_WORLD, "\n***** Fluxin0:%f, Fluxin1:%f, Area:%f\n\n", sumFluxIn0, sumFluxIn1, inletArea);
		//if(ibm->cent_z[elmt] <= 0.000001) ibm->cent_z[ elmt ] = 0.0;
		//ibm->cent_z_AVE[ elmt ] = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
                /*
                double Dis1 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n1e])*(ibm->cent_x[elmt]-ibm->x_bp[n1e])
                                 +(ibm->cent_y[elmt]-ibm->y_bp[n1e])*(ibm->cent_y[elmt]-ibm->y_bp[n1e])),1.e-10);
                            
                double Dis2 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n2e])*(ibm->cent_x[elmt]-ibm->x_bp[n2e])
                                 +(ibm->cent_y[elmt]-ibm->y_bp[n2e])*(ibm->cent_y[elmt]-ibm->y_bp[n2e])),1.e-10);

                double Dis3 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n3e])*(ibm->cent_x[elmt]-ibm->x_bp[n3e])
                                 +(ibm->cent_y[elmt]-ibm->y_bp[n3e])*(ibm->cent_y[elmt]-ibm->y_bp[n3e])),1.e-10);
                            

          	ibm->cent_z[elmt] =  (ibm->z_bp[n1e]/Dis1+ibm->z_bp[n2e]/Dis2+ibm->z_bp[n3e]/Dis3)/(1./Dis1+1./Dis2+1./Dis3);*/
          	                              
                   }
                } 
				else {
                      if((ti==tistart || ti==0) && itt==1 && avalanche_check_number==0) PetscPrintf(PETSC_COMM_WORLD,"Averaging Cz/y onto z/y_bp every time step if needed \n");
                      PetscReal        cent_edge_div = 0.;
		             if(!y_direction){
						PetscReal Ave_cent_z  = (ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
		                PetscReal ratio_z_z = ibm->cent_z[elmt]/(Ave_cent_z+1.e-08);
		                if(ratio_z_z >0.)cent_edge_div = fabs(( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.)-fabs(ibm->cent_z[elmt]);
		                if(ratio_z_z <0.)cent_edge_div = fabs((( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.)-ibm->cent_z[elmt]);
                        if(fabs(cent_edge_div) >= max_delz){
                              PetscPrintf(PETSC_COMM_WORLD,"Averaging: Deviation, max_delz, ti %le %le %d \n", cent_edge_div,max_delz,ti);
                              ibm->cent_z[elmt] = Ave_cent_z;
		                      }
						  }
						  else {
									PetscReal Ave_cent_y  = (ibm->y_bp[ n1e ] + ibm->y_bp[ n2e ] + ibm->y_bp[ n3e ] ) / 3.;
		                            PetscReal ratio_y_y = ibm->cent_y[elmt]/(Ave_cent_y+1.e-08);
		                            if(ratio_y_y >0.)cent_edge_div = fabs(( ibm->y_bp[ n1e ] + ibm->y_bp[ n2e ] + ibm->y_bp[ n3e ] ) / 3.)-fabs(ibm->cent_y[elmt]);
		                            if(ratio_y_y <0.)cent_edge_div = fabs((( ibm->y_bp[ n1e ] + ibm->y_bp[ n2e ] + ibm->y_bp[ n3e ] ) / 3.)-ibm->cent_y[elmt]);
                                        if(fabs(cent_edge_div) >= max_dely)
                                            {
                                             PetscPrintf(PETSC_COMM_WORLD,"Averaging: Deviation, max_dely, ti %le %le %d \n", cent_edge_div,max_dely,ti);
                                             ibm->cent_y[elmt] = Ave_cent_y;
                                            }
						     }
                       }
         }
return (0);
}



PetscErrorCode Scour(UserCtx * user, IBMNodes * ibm, PetscInt tistart, PetscInt ti, PetscInt itr_sc)
{
        PetscTruth SMOOTHING;
	DA                da  = user->da, fda  = user->fda;
	DALocalInfo  	 info  = user->info;
		PetscReal	 max_delz = 1./double(info.my);
        PetscReal	 max_dely = 1./double(info.my);
		
//---------------------------
	Cmpnts		 *** ucat;
	Cmpnts		 *   Bvel;  
//	PetscReal	 *** ustar;       
	PetscInt	     n_elmt  = ibm->n_elmt, n_v  = ibm->n_v;
	PetscReal	 *	 F_flux_12, 
			 *	 F_flux_13, 
			 *	 F_flux_23, *elmt_depth,*netflux_old, *Dflux, *Eflux, *qflux; //Hossein
	PetscReal	 *	 SCont;
	//---------------------------
//	IBMListNode      *	 current;
//	PetscReal	 *** nvert, *** nvert_o;
//	IBMInfo  	 *	 ibminfo;
//--------------------------
        PetscInt rank;
	PetscInt		 i,iter, iteration, elmt, vert,n1e,n2e,n3e;
	PetscReal		 nfx, nfy, nfz;
	PetscReal		 ucx, ucy, ucz;
	PetscReal		 xc, yc, zc, dx12, dx13, dy12, dy13, dz12, dz13, dr, riter;
	PetscReal		 netflux, dzdt, z_minimum,z_maximum;
	PetscReal		 dydt, y_minimum,y_maximum;
        PetscReal                maxz, minz;
        PetscReal                maxy, miny;
        PetscInt                 bend = 0;
        PetscInt                 contra = 1;
        PetscReal                angle_repose = Angle_repose * PI / 180.;
        PetscReal                checoo;
        PetscReal                ts,te,cputime;
        PetscReal                Total_Area, infinity_norm, smoothing_residual;
		PetscInt					avalanche_check_number;
		PetscInt					ijk;
        //PetscInt                 aval_loop = 2;
        //PetscReal                dZ[n_elmt], dZ_old[n_elmt], atke[n_elmt], under_relax[n_elmt];
        extern PetscInt projection_method;
        PetscInt                 AngleSkewness_compute = 0;
		
//---------------------------

// compute new normal vec to elmts, elmt area, elmt center point coordinate

  PetscReal ts_1, te_1;

  PetscGetTime(&ts_1);

 if(ti==tistart || ti==0)
 {
        ibm->dtime[tistart-1]=0.;
        ibm->time_bedchange[tistart-1]=0.;
        PetscInt itt=0;
        
        for( elmt = 0; elmt < n_elmt; elmt++ )
	{
                itt++;
		n1e  			  = ibm->nv1[ elmt ];
		n2e  			  = ibm->nv2[ elmt ];
		n3e  			  = ibm->nv3[ elmt ];
		dx12			  = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
		dy12			  = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
		dz12			  = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

		dx13			  = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
		dy13			  = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
		dz13			  = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

		ibm->nf_x[ elmt ] = dy12 * dz13 - dz12 * dy13;
		ibm->nf_y[ elmt ] = -dx12 * dz13 + dz12 * dx13;
		ibm->nf_z[ elmt ] = dx12 * dy13 - dy12 * dx13;

		dr = sqrt(ibm->nf_x[ elmt ] * ibm->nf_x[ elmt ] + ibm->nf_y[ elmt ] * ibm->nf_y[ elmt ]
				+ ibm->nf_z[ elmt ] * ibm->nf_z[ elmt ] );

		ibm->nf_x[ elmt ] /= dr;
		ibm->nf_y[ elmt ] /= dr;
		ibm->nf_z[ elmt ] /= dr;

		ibm->dA[ elmt ]    = dr / 2.;

		ibm->cent_x[ elmt ] = ( ibm->x_bp[ n1e ] + ibm->x_bp[ n2e ] + ibm->x_bp[ n3e ] ) / 3.;
		ibm->cent_y[ elmt ] = ( ibm->y_bp[ n1e ] + ibm->y_bp[ n2e ] + ibm->y_bp[ n3e ] ) / 3.;
		ibm->cent_z[ elmt ] = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;

		if(itr_sc==1)
                {
                if(!y_direction){
                ibm->cent_z_AVE[ elmt ] = ibm->cent_z[elmt];
		ibm->cent_z_old[ elmt ] = ibm->cent_z[elmt];
		ibm->cent_zl[ elmt ] = ibm->cent_z[elmt];} 
                else {
                ibm->cent_y_AVE[ elmt ] = ibm->cent_y[elmt];
		ibm->cent_y_old[ elmt ] = ibm->cent_y[elmt];
		ibm->cent_yl[ elmt ] = ibm->cent_y[elmt];}
                if(itt==1){
                PetscPrintf(PETSC_COMM_WORLD, "y_direction: %d \n",y_direction);
                PetscPrintf(PETSC_COMM_WORLD, "ti start: %d \n",ti);
                PetscPrintf(PETSC_COMM_WORLD, "Angle of Repose: %e\n",angle_repose*180./PI);
                PetscPrintf(PETSC_COMM_WORLD, "Porosity:  %e\n",porosity);
                PetscPrintf(PETSC_COMM_WORLD, "Avalanche correction loop numbers:  %d\n",aval_loop);
                PetscPrintf(PETSC_COMM_WORLD, "Sand-slide model activation status:  %d\n",sand_slide);
                PetscPrintf(PETSC_COMM_WORLD, "projection_method:  %d\n",projection_method);
                }
	        }
                /*
                double Dis1 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n1e])*(ibm->cent_x[elmt]-ibm->x_bp[n1e])
                                 +(ibm->cent_y[elmt]-ibm->y_bp[n1e])*(ibm->cent_y[elmt]-ibm->y_bp[n1e])), 1.e-10);
                            
                double Dis2 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n2e])*(ibm->cent_x[elmt]-ibm->x_bp[n2e])
                                 +(ibm->cent_y[elmt]-ibm->y_bp[n2e])*(ibm->cent_y[elmt]-ibm->y_bp[n2e])),1.e-10);

                double Dis3 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n3e])*(ibm->cent_x[elmt]-ibm->x_bp[n3e])
                                 +(ibm->cent_y[elmt]-ibm->y_bp[n3e])*(ibm->cent_y[elmt]-ibm->y_bp[n3e])),1.e-10);
                            
          	ibm->cent_z[elmt] =  (ibm->z_bp[n1e]/Dis1+ibm->z_bp[n2e]/Dis2+ibm->z_bp[n3e]/Dis3)/(1./Dis1+1./Dis2+1./Dis3); */

                double ds_12 = sqrt(pow((ibm->x_bp[n1e]-ibm->x_bp[n2e]),2.)+pow((ibm->y_bp[n1e]-ibm->y_bp[n2e]),2.)+pow((ibm->z_bp[n1e]-ibm->z_bp[n2e]),2.));
                double ds_13 = sqrt(pow((ibm->x_bp[n1e]-ibm->x_bp[n3e]),2.)+pow((ibm->y_bp[n1e]-ibm->y_bp[n3e]),2.)+pow((ibm->z_bp[n1e]-ibm->z_bp[n3e]),2.));
                double ds_23 = sqrt(pow((ibm->x_bp[n2e]-ibm->x_bp[n3e]),2.)+pow((ibm->y_bp[n2e]-ibm->y_bp[n3e]),2.)+pow((ibm->z_bp[n2e]-ibm->z_bp[n3e]),2.));
                PetscReal ds_max = PetscMax(ds_12,ds_13);
                PetscReal ds_max_ = PetscMax(ds_max,ds_23);

	        z_minimum=PetscMin(ibm->z_bp[ibm->nv1[elmt]],ibm->z_bp[ibm->nv2[elmt]]);
                z_minimum=PetscMin(z_minimum,ibm->z_bp[ibm->nv3[elmt]]); 
                z_maximum=PetscMax(ibm->z_bp[ibm->nv1[elmt]],ibm->z_bp[ibm->nv2[elmt]]);
                z_maximum=PetscMax(z_maximum,ibm->z_bp[ibm->nv3[elmt]]); 
                
	        y_minimum=PetscMin(ibm->y_bp[ibm->nv1[elmt]],ibm->y_bp[ibm->nv2[elmt]]);
                y_minimum=PetscMin(y_minimum,ibm->y_bp[ibm->nv3[elmt]]); 
                y_maximum=PetscMax(ibm->y_bp[ibm->nv1[elmt]],ibm->y_bp[ibm->nv2[elmt]]);
                y_maximum=PetscMax(y_maximum,ibm->y_bp[ibm->nv3[elmt]]);
 
                if(!y_direction){ibm->elmt_depth[elmt]=fabs(z_maximum-z_minimum);}
                else {ibm->elmt_depth[elmt]=fabs(y_maximum-y_minimum);}
                
                if(!input_ib_depth)
                                  {
				ibm->elmt_depth[elmt]=fabs(z_maximum-z_minimum);
                if(ds_max_ >= 0.78 || (!y_direction && ibm->nf_z[elmt] < 0.8) || (y_direction && ibm->nf_y[elmt] < 0.8)) ibm->elmt_depth[elmt] = 5.;
                if(ibm->elmt_depth[elmt] > 0.39) ibm->elmt_depth[elmt] = 5.;
                //if(ibm->cent_x[elmt]<5.0)ibm->elmt_depth[elmt]=5.;
                if(!y_direction){if(ibm->cent_z[elmt]<0.99)ibm->elmt_depth[elmt]=5.;
                                 if(ibm->cent_z[elmt]>1.01)ibm->elmt_depth[elmt]=5.;}
                                 else {if(ibm->cent_y[elmt]<-0.1)ibm->elmt_depth[elmt]=5.;
								 if(ibm->cent_x[elmt]<0.0)ibm->elmt_depth[elmt]=5.;
								 if(ibm->cent_x[elmt]>2.75)ibm->elmt_depth[elmt]=5.;
                                             if(ibm->cent_y[elmt]>0.3)ibm->elmt_depth[elmt]=5.;}
                //if(ibm->elmt_depth[elmt]>1.5)ibm->elmt_depth[elmt]=5.;
                // added just for Indoor_Rock_vane case for z over 1.0 :
                // added just for Cross_vane case for z over 0.0 :
                // added just for Jhook case for z over 0.5 :
                //if(ibm->cent_z[elmt]>=1.3 && ibm->cent_x[elmt]>275.3 && ibm->cent_x[elmt]<294.6)ibm->elmt_depth[elmt]=5.;
                // if(ibm->nf_z[elmt]<1.e-6)ibm->elmt_depth[elmt]=5.;
                // if(ibm->cent_x[elmt]<0.1)ibm->elmt_depth[elmt]=5;
                   
                //PetscReal line_1 = 100.05944;
                //PetscReal line_1 = 90.16268;
                //PetscReal line_2 = 33.07459;
                //PetscReal elmt_coord_1 = 0.95754 * ibm->cent_x[elmt] + ibm->cent_y[elmt]; 
                //PetscReal elmt_coord_1 = 0.91932 * ibm->cent_x[elmt] + ibm->cent_y[elmt]; 
                //PetscReal elmt_coord_2 = 1.13087 * ibm->cent_x[elmt] - ibm->cent_y[elmt]; 
                //if(elmt_coord_1 > line_1 && elmt_coord_2 < line_2)
                //if(elmt_coord_1 > line_1)
                // {
                //  ibm->Rigidity[elmt] = 0;
                // } else { 
                //  ibm->Rigidity[elmt] = 1;
                // }

                ibm->Mobility[elmt] = 1;
                ibm->Rigidity[elmt] = 0;
                if(ibm->elmt_depth[elmt] > 4.5) ibm->Rigidity[elmt] = 1;               
                if(ibm->elmt_depth[elmt] > 4.5) ibm->Mobility[elmt] = 0;               
                                   }
                ibm->netflux_old[elmt]=0.;
				ibm->Dflux[elmt]=0.;  //Hossein
				ibm->Eflux[elmt]=0.; //Hossein
				ibm->qflux[elmt]=0.; //Hossein
              }
  }

  PetscGetTime(&te_1);

  PetscPrintf(PETSC_COMM_WORLD, "#### IB_Depth  cputime %le\n", te_1-ts_1);



  PetscReal ts_2, te_2;

  PetscGetTime(&ts_2);

        
if(STRONG_COUPLING)
{
  for(vert = 0; vert < ibm->n_v; vert++)
     {       
      if(!y_direction){ibm->z_bp_l[vert]=ibm->z_bp[vert];}
      else {ibm->y_bp_l[vert]=ibm->y_bp[vert];}
     }
}          

//Computing Angle of Skewness
if(AngleSkewness_compute) AngleSkewness(ibm);


// Read and/or Write the "ibm->elmt_depth[elmt]"
    if(input_ib_depth && (ti==tistart || ti==0))
    {  
    read_bed_cell_depth(ibm, tistart);
    PetscBarrier( PETSC_NULL );  								//	stop till all procc to the jobs  			   
    } 
     else
    { 
    if (ti==tistart || ti==0) {  
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(!rank){
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "WRITE bed cell depth data\n");
    char string[128];
    char filen[80];
    sprintf(filen, "ib_cell_depth.dat");
    fd = fopen(filen, "w");
   
    if (fd) {
     
    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %d\n", ibm->elmt_depth[i], ibm->Rigidity[i]);
                                  }
            } fclose(fd);
            }
            }
      
     PetscBarrier( PETSC_NULL );}  								//	stop till all procc to the jobs  			   
// End of Read and/or Write "ibm->elmt_depth[elmt]"

/*    if (ti==tistart) {  
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ bed cell concentration data\n");
    char filen[90];  
    sprintf(filen,"bed_elmt_concntrtn_%5.5d" , ti-1);
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file")
    for (i=0; i<ibm->n_elmt; i++) {
	fscanf(fd, "%e\n", ibm->SCont[i]);
	                          }
    MPI_Bcast(ibm->SCont, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
              }
    else  
                                  {
    MPI_Bcast(ibm->SCont, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
                                  }
                                      }
	PetscBarrier( PETSC_NULL );  								//	stop till all procc to the jobs  			   
*/


  PetscGetTime(&te_2);

  PetscPrintf(PETSC_COMM_WORLD, "#### Write/Read IB_Depth  cputime %le\n", te_2-ts_2);


  PetscReal ts_3, te_3;

  PetscGetTime(&ts_3);



if(projection_method == 1)  {
// call for computing flow at reference level
        PetscGetTime(&ts);
        flow_variables_ref_level (user, ibm, ti, tistart);
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_WALL-Function_ref_level  cputime %d %le\n", ti,cputime);
    
// call for computing Bvel, and BShS on the bed surface
        PetscGetTime(&ts);
	//sediment_variables_projection( user, ibm );
	Projecting( user, ibm );
	PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "Projecting  cputime %d %le\n", ti,cputime);

        PetscGetTime(&ts);
	Finalizing_Projecting( user, ibm, ti, tistart );
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "Finalizing_Projecting  cputime %d %le\n", ti,cputime);
        //PetscPrintf(PETSC_COMM_WORLD, "EXNER_COMM  cputime %d %le\n", ti,cputime);
                            }


if(projection_method == 2)  {
// call for computing flow to a point above bed load layer normal from centroid of bed cell 
        PetscGetTime(&ts);
        ibm_intp_pj_centroid (user, ibm, ti, tistart);
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_ibm_intp_pj_centroid  cputime %d %le\n", ti,cputime);
    
// call for comunicating velocities above bed load layer 
        PetscGetTime(&ts);
	Projecting_new( user, ibm );
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_COMM through Projecting_new  cputime %d %le\n", ti,cputime);


// call for wallfunction to construct the velociies and shear stress on the ref level 
        PetscGetTime(&ts);
        flow_variables_ref_level_new (user, ibm, ti, tistart);
	PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_wallfunction  cputime %d %le\n", ti,cputime);
                          }



  PetscGetTime(&te_3);

  PetscPrintf(PETSC_COMM_WORLD, "#### Total Projecting cputime %le\n", te_3-ts_3);



  PetscReal ts_4, te_4;

  PetscGetTime(&ts_4);


iteration=0;
SMOOTHING = PETSC_TRUE;
while(SMOOTHING)
{           
        iteration++;
// call for Critical bed shear stress computation
        PetscGetTime(&ts);
        bed_concentration( user, ibm, ti, tistart );
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_BED-Concentration  cputime %d %le\n", ti,cputime);

// computing sediment fluxes for every elemnts on bed surface
        PetscGetTime(&ts);
	    if(y_direction) sediment_flux_y_direction(user, ibm, ti, tistart);
	    else sediment_flux(user, ibm, ti, tistart);
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_SED-Flux-Inside_domain  cputime %d %le\n", ti,cputime);

// outlet flux  must be adjusted for each case study, ali 18 Nov. 2009
        PetscGetTime(&ts);
        if(bend) {
			      if(y_direction) outlet_sediment_flux_bend_y_direction(user, ibm );
				  else outlet_sediment_flux_bend(user, ibm );
		         }
        if(contra) {
			        if(y_direction) outlet_sediment_flux_contra_y_direction(user, ibm);
					else outlet_sediment_flux_contra(user, ibm);
		           }
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_SED-Flux-Outlet_Boundary  cputime %d %le\n", ti,cputime);

// inlet flux to be self-determined at the inlet Feb. 2016
       if(osl_inlet_sediment_flux)
        {
        PetscGetTime(&ts);
	inlet_sediment_flux_osl( user, ibm );
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_SED-Flux-Inlet_Boundary_OSL  cputime %d %le\n", ti,cputime);
        }
// inlet flux  must be adjusted for each case study, ali Jan 4, 2012
       if(inlet_sediment_flux > 0.)
       { 
       //double inlet_sediment_flux = 4.0 /(1650. * 60. * 0.32 * 0.3 * 0.3);
        PetscGetTime(&ts);
        if(bend) inlet_sediment_flux_bend( user, ibm, inlet_sediment_flux );
        //if(contra) inlet_sediment_flux_contra( user, ibm, inlet_sediment_flux); // xiaolei deactivate  RTD 
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "EXNER_SED-Flux-Inlet_Boundary  cputime %d %le\n", ti,cputime);
       }
     
// Inlet_sediment_flux to be set for Periodic B.C. Apr 4, 2012
       if(periodic_morpho)
       { 
        PetscGetTime(&ts);
         //PeriodicSedimentFlux(user, ibm);
        PetscGetTime(&te);
        cputime=te-ts;
        PetscPrintf(PETSC_COMM_WORLD, "Periodicity B.C. for Bed-Load  cputime %d %le\n", ti,cputime);
       }

// defineing the time steps:
        if(itr_sc==1 && iteration ==1)
        {
	/*double Max_dzdt=0.0;
        for (elmt=0;elmt<n_elmt;elmt++)
            {
	      if( ibm->nf_z[elmt] < 1.e-6 || ibm->elmt_depth[elmt] > 0.1)
		{
		}
	      else
		{
	         netflux = ibm->F_flux_12[ elmt ] + ibm->F_flux_13[ elmt ] + ibm->F_flux_23[ elmt ];
                 dzdt = netflux * ( 1. / ibm->dA[ elmt ] )  / ( 1. - porosity );
                 Max_dzdt=PetscMax(Max_dzdt, fabs(dzdt));           
                }
	     }*/
 
        ibm->dtime[ti] = user->dts;//1.0;//1.0;//1.;//user->dt;//3.0; //0.01/Max_dzdt;
        ibm->time_bedchange[ti] = ibm->time_bedchange[ti-1] + ibm->dtime[ti];
        PetscPrintf(PETSC_COMM_WORLD, "dt and time :  %e %e\n",ibm->dtime[ti], ibm->time_bedchange[ti]);
        }
        
// Check for bed smoothening by adjusting the edge-fluxes
           
/*
        double infinity_norm=0.0;
        for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	      if( ibm->nf_z[elmt] < 1.e-6 || ibm->elmt_depth[elmt] > 0.1)
		{
                        ibm->Delz[elmt]=0.; ibm->netflux_old[elmt]=0.;
		}
	      else
		{
			netflux = ibm->F_flux_12[ elmt ] + ibm->F_flux_13[ elmt ] + ibm->F_flux_23[ elmt ];
		        ibm->netflux_old[elmt] = netflux;

                        if((ibm->BShS[elmt]-ibm->BCShS[elmt])<=0.) netflux=PetscMax(0., netflux);                      

                     	dzdt = netflux * ( 1. / ibm->dA[ elmt ] )  / ( 1. - porosity );
			//ibm->Delz[ elmt ] =20.*(user->dt) * dzdt;
			ibm->Delz[ elmt ] = ibm->dtime[ti] * dzdt;
		                 
                        if(ibm->Delz[elmt] >= 0.01)  ibm->Delz[elmt]= 0.01;
                        if(ibm->Delz[elmt] <= -0.01) ibm->Delz[elmt]= -0.01; 
                        ibm->cent_zl[elmt] = ibm->cent_z[elmt];     	
                        ibm->cent_z[ elmt ] = ibm->cent_z_old[elmt] + ibm->Delz[ elmt ];
                        
                        // Under_relaxation
                            atke[elmt]=0.;
                            dZ[elmt] = ibm->cent_zl[elmt]-ibm->cent_z[elmt];
                            if(fabs(dZ[elmt]-dZ_old[elmt])>1.e-8 && ibm->atke_old[elmt]!=0.3)
                            {
                              atke[elmt]+=dZ[elmt]/(dZ_old[elmt]-dZ[elmt]);
                            }
                            atke[elmt]=ibm->atke_old[elmt]+(ibm->atke_old[elmt]-1.)*atke[elmt];
                            if(atke[elmt]>0.9)atke[elmt]=0.9;
                            if(atke[elmt]<-0.2)atke[elmt]=-0.2;
                            under_relax[elmt]=1.-atke[elmt];
                            //ibm->cent_z[elmt] = under_relax[elmt]*ibm->cent_z[elmt] + (1.-under_relax[elmt])*ibm->cent_zl[elmt];
                            ibm->cent_z[elmt] = 0.7*ibm->cent_z[elmt] + (1.-0.7)*ibm->cent_zl[elmt];
                            ibm->atke_old[elmt] = atke[elmt];
                            
                            infinity_norm= PetscMax(infinity_norm, fabs(ibm->cent_z[elmt]-ibm->cent_zl[elmt]));

                }
          // ibm->cent_z_AVE[elmt]=ibm->cent_z[elmt];
	}*/

        smoothing_residual=0.0;
        for( elmt = 0; elmt < n_elmt; elmt++ )
	{
	      if(ibm->Rigidity[elmt]==1)
		{
                        if(!y_direction){ibm->Delz[elmt]=0.; ibm->netflux_old[elmt]=0.; ibm->cent_z[elmt]=ibm->cent_z_old[elmt];}
                        else {ibm->Dely[elmt]=0.; ibm->netflux_old[elmt]=0.; ibm->cent_y[elmt]=ibm->cent_y_old[elmt];}
		}
	      else
		{
		        netflux = (ibm->F_flux_12[ elmt ] + ibm->F_flux_13[ elmt ] + ibm->F_flux_23[ elmt ])*ibm->Deltab;
                         
                        //if(LiveBed) netflux += w_s * (ibm->C[elmt] - ibm->SCont[elmt]) * ibm->dA[elmt]; 

				ibm->qflux[elmt] = netflux;
				ibm->Dflux[elmt] = w_s * (ibm->C[elmt]); //Hossein
				ibm->Eflux[elmt] = w_s * (ibm->SCont[elmt]); //Hossein
                    
		        ibm->netflux_old[elmt] = netflux;

                       // if((ibm->BShS[elmt]-ibm->BCShS[elmt])<=0.) netflux=PetscMax(0., netflux);                      

                     	if (!y_direction)
						{
							dzdt = netflux * ( 1. / ibm->dA[ elmt ] )  / ( 1. - porosity );
	
			                ibm->Delz[ elmt ] = ibm->dtime[ti] * dzdt;
		                 
                            if(fabs(ibm->Delz[elmt]) > max_delz){
                                           PetscPrintf(PETSC_COMM_WORLD, "ti, Delz[elmt], max_delz %d %le %le\n",ti, ibm->Delz[elmt], max_delz);
                                           ibm->Delz[elmt]= max_delz * (fabs(ibm->Delz[elmt])/ibm->Delz[elmt]);
                                                            }
                            // Hitting the rigid bed -- this is OK for OSL computations
                            //if(ibm->cent_z[elmt]>sediment_thickness)
							{ 
							ibm->cent_zl[elmt] = ibm->cent_z[elmt];     	
                        
							ibm->cent_z[ elmt ] = ibm->cent_z_old[elmt] + ibm->Delz[ elmt ];
                        
							// Under_relaxation
							ibm->cent_z[elmt] = 0.7 * ibm->cent_z[elmt] + (1. - 0.7) * ibm->cent_zl[elmt];
							} 
							// Hiting the rigid bed -- it is OK for OSL computations
							if(ibm->cent_x[elmt] < 5.) ibm->cent_z[elmt]= ibm->cent_z_old[elmt];
							if(ibm->cent_z_old[elmt]>=1.999)  ibm->cent_z[elmt] = ibm->cent_z_old[elmt]; 
							//if(ibm->cent_z[elmt] > 0.6) ibm->cent_z[elmt]= ibm->cent_z_old[elmt];
							//if(ibm->cent_x[elmt] > 20.0) ibm->cent_z[elmt]= ibm->cent_z_old[elmt];
							// if(ibm->cent_x[elmt] < 4.0 ) ibm->cent_z[elmt] = 0.2; 
							//if(ibm->cent_z[elmt] < 0.2001)  ibm->cent_z[elmt] = 0.2; 
							//if(ibm->cent_x[elmt] > 6.54 && ibm->cent_z[elmt] < 0.2001)  ibm->cent_z[elmt] = 0.2; 
							//if(ibm->cent_y[elmt] < 2.0  && ibm->cent_z[elmt] < 0.2001)  ibm->cent_z[elmt] = 0.2; 
							//if(ibm->cent_y[elmt] > 4.54 && ibm->cent_z[elmt] < 0.2001)  ibm->cent_z[elmt] = 0.2; 
							//if(ibm->cent_z[elmt] <= 0.01) ibm->cent_z[elmt] = 0.0; 
							//if(ibm->cent_x[elmt] <= 22.5) ibm->cent_z[elmt] = ibm->cent_z_old[elmt]; 
							//ibm->cent_z[elmt] = ibm->cent_z_old[elmt]; 
						    smoothing_residual= PetscMax(smoothing_residual, fabs(ibm->cent_z[elmt]-ibm->cent_zl[elmt]));
						}
					    else{
							dydt = netflux * ( 1. / ibm->dA[ elmt ] )  / ( 1. - porosity );
	
			                ibm->Dely[ elmt ] = ibm->dtime[ti] * dydt;
		                 
                            if(fabs(ibm->Dely[elmt]) > max_dely){
                                           PetscPrintf(PETSC_COMM_WORLD, "ti, Dely[elmt], max_dely %d %le %le\n",ti, ibm->Dely[elmt], max_dely);
                                           ibm->Dely[elmt]= max_dely * (fabs(ibm->Dely[elmt])/ibm->Dely[elmt]);
                                                            }
							ibm->cent_yl[elmt] = ibm->cent_y[elmt];     	
							ibm->cent_y[ elmt ] = ibm->cent_y_old[elmt] + ibm->Dely[ elmt ];
                        
       						// Under_relaxation
							ibm->cent_y[elmt] = 0.7 * ibm->cent_y[elmt] + (1. - 0.7) * ibm->cent_yl[elmt];
							//if(ibm->cent_z[elmt] < 0.) ibm->cent_y[elmt]= ibm->cent_y_old[elmt]; temporary comment for Xiaoli's project
							//if(ibm->cent_y_old[elmt]>=0.999)  ibm->cent_y[elmt] = ibm->cent_y_old[elmt]; temporary comment for Xiaoli's project
						    smoothing_residual= PetscMax(smoothing_residual, fabs(ibm->cent_y[elmt]-ibm->cent_yl[elmt]));}

                }
          
	}
             
//start add Hossein             
// check and correct the new bed elevation based on the material angle of repose
    if(sand_slide) {
		 aval_loop = 0;
		 avalanche_check_number = 0;
		 if (!y_direction){
			check_correct_new_elev(user, ibm, ti);
		 } else {
			check_correct_new_elev_y_direction(user, ibm, ti);
		 }
                 PetscBarrier( PETSC_NULL );  								//	stop till all procc to the jobs  			   
	               } 
    
// check and correct the new bed elevation based on the Avalanche model
    if(aval_loop){
		   avalanche_check_number = 0;
				   
		for(ijk = 0; ijk < aval_loop; ijk++){
			if (!y_direction){
					avalanche_first_sweep(user, ibm, ti);
					avalanche_second_sweep(user, ibm, ti);
			}
			else {
					avalanche_first_sweep_y_direction(user, ibm, ti);
					avalanche_second_sweep_y_direction(user, ibm, ti);
			}
		}
	           double max_slope = 0.;
                   int iel;
                   for(iel=0; iel<ibm->n_elmt; iel++) {max_slope = PetscMax(max_slope, fabs(ibm->max_bed_angle[iel])); }
                   PetscPrintf(PETSC_COMM_WORLD, "max_bed_angle: %le angle_repose: %le ti: %d \n",max_slope, angle_repose, ti);
        PetscBarrier( PETSC_NULL );  								//	stop till all procc to the jobs  			   
	}
	//end add Hossein


	// transfer the bed change from center of elmt to vertices and assign the new z_bp[] valuse
        for(vert = 0; vert < ibm->n_v; vert++)
	{       
             Total_Area = 0.;
             double zb = 0.;             
			 double yb = 0.;

             riter = 0.;

				// xiaole add
	
				int ii,jj;
				for(ii=0;ii<100;ii++)
        		{
	
					elmt=ibm->n2c[vert].c[ii];	// xiaolei add
			        if (elmt>=0) 
					{
                      int rigid = ibm->Rigidity[elmt];
		              if((!y_direction && fabs(ibm->Delz[elmt])<1.e-6) || (y_direction && fabs(ibm->Dely[elmt]) < 1.e-6) || rigid )
			            {}
			            else
			                {
			                 if (ibm->nv1[elmt]==vert || ibm->nv2[elmt]==vert || ibm->nv3[elmt]==vert)
			                    {
								if (!y_direction) {zb +=  ibm->cent_z[elmt] * ibm->dA[elmt];} else {yb +=  ibm->cent_y[elmt] * ibm->dA[elmt];}
                                  riter++;  
                                  Total_Area += ibm->dA[elmt];
			                    }
			                     else
			                        {}
		                    }
			        }
                }
             
			 if(riter > 0.)
			 {
				if(!y_direction)
				{
				                 ibm->z_bp[vert] = zb / Total_Area; 
                                 // Hitting the rigid bed -- it is OK for OSL computations
                                 // if(riter>0. && ibm->z_bp[vert] < sediment_thickness)ibm->z_bp[vert] = sediment_thickness; 
                                 //if(riter>0. && ibm->z_bp[vert] > 0.6)ibm->z_bp[vert] = ibm->z_bp_o[vert]; 
                                 if(ibm->x_bp[vert] <= 5.)ibm->z_bp[vert] = ibm->z_bp_o[vert]; 
                                 if(ibm->z_bp_o[vert] > 1.999)ibm->z_bp[vert] = ibm->z_bp_o[vert]; 
                                 //if(riter>0. && ibm->z_bp[vert] <= -0.599)  ibm->z_bp[vert] = -0.599; 
				}
				 else 
				        {
					    ibm->y_bp[vert] = yb / Total_Area; 
						if(ibm->z_bp[vert] <= 0.)ibm->y_bp[vert] = ibm->y_bp_o[vert]; 
                                                if(ibm->y_bp_o[vert] > 0.999)ibm->y_bp[vert] = ibm->y_bp_o[vert];
 				        }

             }
        }
        

        // fixing inlet shape, preventing bed change at inlet section
        //PetscReal line_1 = 100.05944;
        //PetscReal line_1 = 90.16268;
        //PetscReal line_2 = 32.65227;
        //PetscReal line_2 = 33.07459;

        //for (vert=0; vert<ibm->n_v; vert++)
       // {
                //PetscReal vert_coord_1 = 0.95754 * ibm->x_bp[vert] + ibm->y_bp[vert]; 
                //PetscReal vert_coord_1 = 0.91932 * ibm->x_bp[vert] + ibm->y_bp[vert]; 
                //PetscReal vert_coord_2 = 1.12558 * ibm->x_bp[vert] - ibm->y_bp[vert]; 
                //PetscReal vert_coord_2 = 1.13087 * ibm->x_bp[vert] - ibm->y_bp[vert]; 
                
                //if(vert_coord_1 < line_1 || vert_coord_2 > line_2 ) ibm->z_bp[vert] = ibm->z_bp_o[vert];
                //if(vert_coord_1 < line_1 ) ibm->z_bp[vert] = ibm->z_bp_o[vert];
         //if( ibm->x_bp[vert]< 1.5) ibm->z_bp[vert] = ibm->z_bp_o[vert];
         //if( ibm->z_bp_o[vert]> 1.49) ibm->z_bp[vert] = ibm->z_bp_o[vert];
         //if( ibm->x_bp[vert]< 25.1 && ibm->y_bp[vert]>20.) ibm->z_bp[vert] = ibm->z_bp_o[vert];
         //if( ibm->z_bp[vert]> 1.35) ibm->z_bp[vert] = ibm->z_bp_old[vert];
         //if( ibm->x_bp[vert]<285. && ibm->y_bp[vert]>47.4 && ibm->z_bp[vert]> 1.3) ibm->z_bp[vert] = ibm->z_bp_o[vert];
        // }

        recomputing_geometry( user, ibm, tistart,ti,itr_sc, avalanche_check_number );
        
        /*double max_slope = 0.;
        int iel;

        for(iel=0; iel<ibm->n_elmt; iel++)
        {
         max_slope = PetscMax(max_slope, fabs(ibm->max_bed_angle[iel]));
        }
        

        if(max_slope > angle_repose && avalanche_check_number<aval_loop) goto correct_cell_slope;*/

        PetscPrintf(PETSC_COMM_WORLD, "smoothing_residual ti itr_sc iteration_number %d %d %d %le\n",ti, itr_sc, iteration, smoothing_residual);
        if(smoothing_residual<1.e-6) SMOOTHING = PETSC_FALSE;
        if(iteration>0) SMOOTHING = PETSC_FALSE;

        // computing new nfx nfy nfz after bed changes
       /* for( elmt = 0; elmt < n_elmt; elmt++ )
	{
		n1e  			  = ibm->nv1[ elmt ];
		n2e  			  = ibm->nv2[ elmt ];
		n3e  			  = ibm->nv3[ elmt ];
		dx12			  = ibm->x_bp[ n2e ] - ibm->x_bp[ n1e ];
		dy12			  = ibm->y_bp[ n2e ] - ibm->y_bp[ n1e ];
		dz12			  = ibm->z_bp[ n2e ] - ibm->z_bp[ n1e ];

		dx13			  = ibm->x_bp[ n3e ] - ibm->x_bp[ n1e ];
		dy13			  = ibm->y_bp[ n3e ] - ibm->y_bp[ n1e ];
		dz13			  = ibm->z_bp[ n3e ] - ibm->z_bp[ n1e ];

		ibm->nf_x[ elmt ] = dy12 * dz13 - dz12 * dy13;
		ibm->nf_y[ elmt ] = -dx12 * dz13 + dz12 * dx13;
		ibm->nf_z[ elmt ] = dx12 * dy13 - dy12 * dx13;

		dr = sqrt(ibm->nf_x[ elmt ] * ibm->nf_x[ elmt ] + ibm->nf_y[ elmt ] * ibm->nf_y[ elmt ]
				+ ibm->nf_z[ elmt ] * ibm->nf_z[ elmt ] );

		ibm->nf_x[ elmt ] /= dr;
		ibm->nf_y[ elmt ] /= dr;
		ibm->nf_z[ elmt ] /= dr;

        		//ibm->dA[elmt ]    = dr / 2.;
                 
                if(((ti-tistart)/1)*1== (ti-tistart))
                {
	//	ibm->cent_z[ elmt ] = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
	       //  ibm->cent_z_AVE[ elmt ] = ( ibm->z_bp[ n1e ] + ibm->z_bp[ n2e ] + ibm->z_bp[ n3e ] ) / 3.;
                
          //      double Dis1 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n1e])*(ibm->cent_x[elmt]-ibm->x_bp[n1e])
          //                       +(ibm->cent_y[elmt]-ibm->y_bp[n1e])*(ibm->cent_y[elmt]-ibm->y_bp[n1e])),1.e-10);
          //                  
          //      double Dis2 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n2e])*(ibm->cent_x[elmt]-ibm->x_bp[n2e])
          //                       +(ibm->cent_y[elmt]-ibm->y_bp[n2e])*(ibm->cent_y[elmt]-ibm->y_bp[n2e])),1.e-10);
          // 
          //      double Dis3 = PetscMax(sqrt((ibm->cent_x[elmt]-ibm->x_bp[n3e])*(ibm->cent_x[elmt]-ibm->x_bp[n3e])
          //                       +(ibm->cent_y[elmt]-ibm->y_bp[n3e])*(ibm->cent_y[elmt]-ibm->y_bp[n3e])),1.e-10);
          //                  
          //
          // 	ibm->cent_z[elmt] =  (ibm->z_bp[n1e]/Dis1+ibm->z_bp[n2e]/Dis2+ibm->z_bp[n3e]/Dis3)/(1./Dis1+1./Dis2+1./Dis3);
          	                              
                }
       } */        
       

   }   // end of while SMOOTHIMG
	
  PetscGetTime(&te_4);

  PetscPrintf(PETSC_COMM_WORLD, "#### Smoothing cputime %le\n", te_4-ts_4);


        // computing bed change speed
	for( vert = 0; vert < ibm->n_v; vert++ )
	{
		ibm->u[ vert ].x = 0.;                              
		if(!y_direction){
							ibm->u[ vert ].z = ( ibm->z_bp[ vert ] - ibm->z_bp_o[ vert ] ) / (user->dt);
							ibm->u[ vert ].y = 0.;
						}
						else {
								ibm->u[ vert ].y = ( ibm->y_bp[ vert ] - ibm->y_bp_o[ vert ] ) / (user->dt);
								ibm->u[ vert ].z = 0.;
							
								}
                //if(fabs(ibm->z_bp[vert] - ibm->z_bp_o[vert])>1.e-6) {maxz=PetscMax(ibm->z_bp[vert],maxz); minz=PetscMin(ibm->z_bp[vert],minz);}
	}
        
 //PetscPrintf(PETSC_COMM_WORLD, "min and max bed elevation at time_step is: %d %e %e\n",ti,minz,maxz);
 //PetscPrintf(PETSC_COMM_WORLD, "Infinity Norm of computations for ti, itr_sc , infinity_norm %d %d %e\n",ti,itr_sc,smoothing_residual);

	// print out for check and debuging  6 nov 2009
	//---------------------------------------

if(((ti-tistart)/tiout)*tiout==(ti-tistart)) {  //print out bed elevations for "ibmdata00" and bed cell concentration for "bed_elmt_concntrtn"
 char   ss[20];
 char string[128];
 PetscInt ii=0;    
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  FILE            *fd;
  char            filen[80];
  if(!rank){
     sprintf(filen, "ibmdata%5.5d.dat",ti-1);
    fd = fopen(filen, "w");
   
      PetscReal cl = 1.;
      PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);

    if (fd) {
          PetscFPrintf(PETSC_COMM_WORLD, fd, "#      bed change trainle mesh\n");
          PetscFPrintf(PETSC_COMM_WORLD, fd, "#      bed change trainle mesh\n");
          PetscFPrintf(PETSC_COMM_WORLD, fd, "#      bed change trainle mesh\n");
        
          PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %i %i %i\n",ibm->n_v,ibm->n_elmt,ii,ii,ii);
           
   for (i=0; i<ibm->n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %le %le %le\n", ii,ibm->x_bp[i]*cl, ibm->y_bp[i]*cl, ibm->z_bp[i]*cl);
	          	      }
     
   for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %s %i %i %i\n", ii,ii, "this is for sediment transport module",ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);

	//fscanf(fd, "%i %i %s %i %i %i\n", ii,ii, ss, nv1+i, nv2+i, nv3+i);
	//nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      }
//      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

  //    i=0;
   //   PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);
       }
      fclose(fd);
    }
/*
  if(!rank){
     
    sprintf(filen,"bed_elmt_concntrtn_%5.5d" , ti-1);
    fd = fopen(filen, "w");

   if (fd) {
   for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n",ibm->SCont[i]);
	          	      }
           }
       
      fclose(fd);
    }  */
}

//  PetscPrintf(PETSC_COMM_WORLD, "ti=%d, tistart= %d\n",ti,tistart);
if(((ti-tistart)/sed_tio)*sed_tio==(ti-tistart)) {  //print out bed elevations
//if((ti==40661) || (ti== 40662)){//|| (ti== 14892) || (ti == 14893)) {  //print out bed elevations
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt n_v================== %d %d\n",ti,tistart);

 PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  FILE            *f;
  char            filen[80];
  if(rank==1){
     sprintf(filen, "bed_change_%5.5d.dat",(ti-1));
    f = fopen(filen, "w");
	if(!y_direction){
					if(effective_bed_shear_stress && rans)PetscFPrintf(PETSC_COMM_SELF, f, "Variables=x,y,z,u,v,w,UU,tau0,tau0_v,tau_cr,SCont,flux1,flux2,flux3,delz,w_ave_z,k_ave,angle,rigidity,qflux,Dflux,Eflux\n");
					if(!effective_bed_shear_stress)PetscFPrintf(PETSC_COMM_SELF, f, "Variables=x,y,z,u,v,w,UU,tau0,tau_cr,SCont,flux1,flux2,flux3,delz,cent_z,angle,rigidity,qflux,Dflux,Eflux\n");
					if(effective_bed_shear_stress && rans){PetscFPrintf(PETSC_COMM_SELF, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-22]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);} 
					else{PetscFPrintf(PETSC_COMM_SELF, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-20]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);}
	                }
					 else {
						if(effective_bed_shear_stress && rans)PetscFPrintf(PETSC_COMM_SELF, f, "Variables=x,y,z,u,v,w,UU,tau0,tau0_v,tau_cr,SCont,flux1,flux2,flux3,dely,w_ave_z,k_ave,angle,rigidity,qflux,Dflux,Eflux\n");
						if(!effective_bed_shear_stress)PetscFPrintf(PETSC_COMM_SELF, f, "Variables=x,y,z,u,v,w,UU,tau0,tau_cr,SCont,flux1,flux2,flux3,dely,cent_y,angle,rigidity,qflux,Dflux,Eflux\n");
						if(effective_bed_shear_stress && rans){PetscFPrintf(PETSC_COMM_SELF, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-22]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);} 
						else{PetscFPrintf(PETSC_COMM_SELF, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-20]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);}  
					      }
  //  PetscPrintf(PETSC_COMM_WORLD, "new_bed level\n");
    // x - component
    
    for (i=0; i<ibm->n_v; i++) {
    PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->x_bp[i]);
    }
    //PetscFPrintf(PETSC_COMM_SELF, f, "\n");
    // y - component                                                                 
    for (i=0; i<ibm->n_v; i++) {
		     if(y_direction){
							if(ti==tistart){ PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->y_bp_o[i]);}
							else {PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->y_bp[i]);}
							}
							else {
								  PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->y_bp[i]);
								 } 
    }
    //PetscFPrintf(PETSC_COMM_SELF, f, "\n");
    // z - component
    for (i=0; i<ibm->n_v; i++) {
     if(!y_direction){
					  if(ti==tistart){ PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->z_bp_o[i]);}
					  else {PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->z_bp[i]);}
	                 }
					  else {
						    PetscFPrintf(PETSC_COMM_SELF, f, "%e \n", ibm->z_bp[i]);
						   } 
    }
    
    //Write out the Delta_z
    for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Bvel[i].x);
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Bvel_u[i]);
      }
    for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Bvel[i].y);
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Bvel_v[i]);
      } 
     for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Bvel[i].z);
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Bvel_w[i]);
      }
     for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",sqrt(ibm->Bvel[i].x*ibm->Bvel[i].x+ibm->Bvel[i].y*ibm->Bvel[i].y+ibm->Bvel[i].z*ibm->Bvel[i].z));
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",sqrt(ibm->Bvel_u[i]*ibm->Bvel_u[i]+ibm->Bvel_v[i]*ibm->Bvel_v[i]+ibm->Bvel_w[i]*ibm->Bvel_w[i]));
      } 
     for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->BShS[i]);
      } 
     if(effective_bed_shear_stress && rans) {
      for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->BShS_v[i]);
      }
                                            } 
     for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->BCShS[i]);
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->c_grad_x[i]);
      } 
      for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->SCont[i]);
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->c_grad_y[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) 
      {
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->u_grad_x[i]);
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->F_flux_12[i]);
      }  
     for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->F_flux_13[i]);
      }  
     for (i=0; i<ibm->n_elmt; i++) 
      {
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->cent_z_AVE[i]);
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->F_flux_23[i]);
      }

     for (i=0; i<ibm->n_elmt; i++) 
      {
	if(!y_direction) {PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Delz[i]);} else {PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->Dely[i]);}
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->BCShS[i]);
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->cent_z_AVE[i]);
      }  
     for (i=0; i<ibm->n_elmt; i++) 
      {
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->C[i]);
	if(effective_bed_shear_stress)PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->w_ave[i]);
	if(!effective_bed_shear_stress){if (!y_direction){PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->cent_z[i]);} else {PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->cent_y[i]);}}
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->max_bed_angle[i]*180./3.1415);
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->elmt_depth[i]);
      }
     if(effective_bed_shear_stress && rans){
      for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->k_ave[i]);
      }
      }
     for (i=0; i<ibm->n_elmt; i++) 
      {
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->C[i]);
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->cent_z[i]);
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->max_bed_angle[i]*180./3.1415);
	if(!y_direction)PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",(180./PI)*atan(sqrt((ibm->nf_y[i]/(ibm->nf_z[i]+1.e-8))*(ibm->nf_y[i]/(ibm->nf_z[i]+1.e-8)) + (ibm->nf_x[i]/(ibm->nf_z[i]+1.e-8))*(ibm->nf_x[i]/(ibm->nf_z[i]+1.e-8)))));else PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",(180./PI)*atan(sqrt((ibm->nf_z[i]/(ibm->nf_y[i]+1.e-8))*(ibm->nf_z[i]/(ibm->nf_y[i]+1.e-8)) + (ibm->nf_x[i]/(ibm->nf_y[i]+1.e-8))*(ibm->nf_x[i]/(ibm->nf_y[i]+1.e-8))))); 
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->elmt_depth[i]);
      }
     for (i=0; i<ibm->n_elmt; i++) 
      {
	PetscFPrintf(PETSC_COMM_SELF, f, "%d\n",ibm->Rigidity[i]);
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->cent_z[i]);
//	PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->max_bed_angle[i]*180./3.1415);
	//PetscFPrintf(PETSC_COMM_SELF, f, "%e\n",ibm->elmt_depth[i]);
      }

	  for (i=0; i<ibm->n_elmt; i++) //Hossein
		{
			PetscFPrintf(PETSC_COMM_SELF, f, " %e \n", ibm->qflux[elmt]);
		}
	  
	  for (i=0; i<ibm->n_elmt; i++) //Hossein
		{
			PetscFPrintf(PETSC_COMM_SELF, f, " %e \n", ibm->Dflux[elmt]);
		}

		for (i=0; i<ibm->n_elmt; i++) //Hossein
		{
			PetscFPrintf(PETSC_COMM_SELF, f, " %e \n", ibm->Eflux[elmt]);
		}

 
    //Write out the link nodes
    for (i=0; i<ibm->n_elmt; i++) 
      {
	
	PetscFPrintf(PETSC_COMM_SELF, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }   
    
   
    fclose(f);
  } //end of Comm-wolrd procces
  
} //end if iteration

return ( 0 );
}

// xiaolei add SEDI
PetscErrorCode Connectivity_ib( UserCtx * user, IBMNodes * ibm ) 
{
	int i,j;
	int n1e,n2e,n3e;
	int _n1e,_n2e,_n3e;

	for(i=0;i<ibm->n_elmt;i++) {
		n1e=ibm->nv1[i];	
		n2e=ibm->nv2[i];	
		n3e=ibm->nv3[i];	
		
		ibm->c2c[i].c1=-100;
		ibm->c2c[i].c2=-100;
		ibm->c2c[i].c3=-100;


		for(j=0;j<ibm->n_elmt;j++) {
			_n1e=ibm->nv1[j];
			_n2e=ibm->nv2[j];
			_n3e=ibm->nv3[j];

			if (j!=i && (_n1e==n1e || _n2e==n1e || _n3e==n1e) &&  (_n1e==n2e || _n2e==n2e || _n3e==n2e) ) ibm->c2c[i].c1=j; 
			if (j!=i && (_n1e==n1e || _n2e==n1e || _n3e==n1e) &&  (_n1e==n3e || _n2e==n3e || _n3e==n3e) ) ibm->c2c[i].c2=j; 
			if (j!=i && (_n1e==n2e || _n2e==n2e || _n3e==n2e) &&  (_n1e==n3e || _n2e==n3e || _n3e==n3e) ) ibm->c2c[i].c3=j; 

		}

	}


	for(i=0;i<ibm->n_v;i++) {

		int iii=0;
		int ii;
		
		for(ii=0;ii<100;ii++) {
			ibm->n2c[i].c[ii]=-100;
		}

		for(j=0;j<ibm->n_elmt;j++) {
			_n1e=ibm->nv1[j];
			_n2e=ibm->nv2[j];
			_n3e=ibm->nv3[j];

			if (_n1e==i || _n2e==i || _n3e==i) {
				ibm->n2c[i].c[iii]=j;
				iii++; 
				if (iii>99) PetscPrintf(PETSC_COMM_WORLD, "The numbder of cells sharing the node %d is larger than 100\n",i);
			}

		}

	}


	for(i=0;i<ibm->n_elmt;i++) {

		int iii=0;
		int ii;
		n1e=ibm->nv1[i];	
		n2e=ibm->nv2[i];	
		n3e=ibm->nv3[i];	
		
		for(ii=0;ii<200;ii++) {
			ibm->c2cs[i].cs[ii]=-200;
		}

		for(j=0;j<ibm->n_elmt;j++) {
			_n1e=ibm->nv1[j];
			_n2e=ibm->nv2[j];
			_n3e=ibm->nv3[j];

			if (i!=j && (_n1e==n1e || _n2e==n1e || _n3e==n1e || _n1e==n2e || _n2e==n2e || _n3e==n2e || _n1e==n3e || _n2e==n3e || _n3e==n3e)) {
				ibm->c2cs[i].cs[iii]=j;
				iii++; 
				if (iii>199) PetscPrintf(PETSC_COMM_WORLD, "The numbder of cells around the cell %d is larger than 200\n",i);
			}

		}

	}
}

