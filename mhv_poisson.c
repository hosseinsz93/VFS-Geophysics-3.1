#include "variables.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscmg.h"

extern PetscInt block_number, ti;
extern PetscInt immersed;

PetscErrorCode MyPoissonKSPMonitor(KSP ksp, PetscInt n, PetscReal rnom, void *dummy)
{
  UserCtx *user = (UserCtx*)dummy;
  Vec x;
  PetscInt tttt, mi, i, j, k;
  PetscReal norm;

  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	its;
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

  KSPGetIterationNumber(ksp, &its);
  its++;

  if ((its/20)*20 == its && immersed) {
    Vec TempUcont;
    VecDuplicate(user->Ucont, &TempUcont);
    VecCopy(user->Ucont, TempUcont);
    VecDuplicate(user->Ucont, &user->Ucont_o);
    KSPBuildSolution(ksp, PETSC_NULL, &x);
/*     Projection2(user, TempUcont, x); */
/*     Contra2Cart2(user, TempUcont); */
/*     ibm_interpolation2(user->ibm_intp, user, user->ibm, TempUcont); */
    Cmpnts ***ucont, ***uucont;
    PetscReal ***rhs, ***aj, ***nvert;
    Vec Rhs;
    KSPGetRhs(ksp, &Rhs);
    DAVecGetArray(da, Rhs, &rhs);
    DAVecGetArray(fda, TempUcont, &ucont);
    DAVecGetArray(da, user->lAj, &aj);
    DAVecGetArray(da, user->lNvert, &nvert);
    DAVecGetArray(fda, user->Ucont, &uucont);

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	      nvert[k][j+1][i] + nvert[k][j-1][i] +
	      nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) {
	    rhs[k][j][i] = -(ucont[k][j][i].x - ucont[k][j][i-1].x +
			    ucont[k][j][i].y - ucont[k][j-1][i].y +
			    ucont[k][j][i].z - ucont[k-1][j][i].z) / user->dt * aj[k][j][i] / user->st;

	  }
	  if (nvert[k][j][i] > 0.1) {
	    uucont[k][j][i] = ucont[k][j][i];
	  }
	  else {
	    if (nvert[k][j][i+1] > 0.1) {
	      uucont[k][j][i].x = ucont[k][j][i].x;
	    }
	    if (nvert[k][j+1][i] > 0.1) {
	      uucont[k][j][i].y = ucont[k][j][i].y;
	    }
	    if (nvert[k+1][j][i] > 0.1) {
	      uucont[k][j][i].z = ucont[k][j][i].z;
	    }
	  }
	}
      }
    }

    
    DAVecRestoreArray(fda, user->Ucont, &uucont);
    DAVecRestoreArray(da, Rhs, &rhs);
    DAVecRestoreArray(fda, TempUcont, &ucont);
    DAVecRestoreArray(da, user->lAj, &aj);
    DAVecRestoreArray(da, user->lNvert, &nvert);
    VecDestroy(user->Ucont_o);
    VecDestroy(TempUcont);
  }
  
/*   KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &x); */
/*   VecMax(x, &tttt, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "KSP Max %d %le\n", tttt, norm); */

/*   for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (mi=xs; mi<xe; mi++) { */
/* 	if (lidx(mi,j,k,user)*3 ==tttt) { */
/* 	  PetscPrintf(PETSC_COMM_SELF, "KspMMa %d %d %d %d %le\n", lidx(mi,j,k,user)*3,mi,j, k, norm); */
/* 	} */
/* 	if (lidx(mi,j,k,user)*3+1 ==tttt) { */
/* 	  PetscPrintf(PETSC_COMM_SELF, "KspMMa1 %d %d %d %d %le\n", lidx(mi,j,k,user)*3,mi,j, k, norm); */
/* 	} */
/* 	if (lidx(mi,j,k,user)*3+2 ==tttt) { */
/* 	  PetscPrintf(PETSC_COMM_SELF, "KspMMa2 %d %d %d %d %le\n", lidx(mi,j,k,user)*3,mi,j, k, norm); */
/* /\* 	  PetscPrintf(PETSC_COMM_SELF, "MMa2 %d %d %d %d %le %le %le %le %le %le %le\n", tttt, mi,j, k, ucat[k][j][mi].x, ucat[k][j][mi].y, ucat[k][j][mi].z, nvert[k][j][mi], rct[k][j][mi].x, rct[k][j][mi].y, rct[k][j][mi].z); *\/ */
/* 	} */

/*       } */
/*     } */
/*   } */

/*   KSPBuildSolution(ksp, PETSC_NULL, &x); */
/*   VecMax(x, &tttt, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "KSP Max1 %d %le\n", tttt, norm); */
  return 0;
}

PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user)
{
  if (*(user->isc)) {
    *ih = i;
  }
  else {
    *ih = 2 * i;
  }

  if (*(user->jsc)) {
    *jh = j;
  }
  else {
    *jh = 2 * j;
  }

  if (*(user->ksc)) {
    *kh = k;
  }
  else {
    *kh = 2 * k;
  }

  return 0;
}

#define GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user) \
  if (*(user->isc)) { \
    ic = i; \
    ia = 0; \
  } \
  else { \
    ic = (i+1) / 2; \
    ia = (i - 2 * (ic)) == 0 ? 1 : -1; \
    if (i==1 || i==mx-2) ia = 0; \
  }\
  if (*(user->jsc)) { \
    jc = j; \
    ja = 0; \
  } \
  else { \
    jc = (j+1) / 2; \
    ja = (j - 2 * (jc)) == 0 ? 1 : -1; \
    if (j==1 || j==my-2) ja = 0; \
  } \
  if (*(user->ksc)) { \
    kc = k; \
    ka = 0; \
  } \
  else { \
    kc = (k+1) / 2; \
    ka = (k - 2 * (kc)) == 0 ? 1 : -1; \
    if (k==1 || k==mz-2) ka = 0; \
  } \
  if (ka==-1 && nvert_c[kc-1][jc][ic] > 0.1) ka = 0; \
  else if (ka==1 && nvert_c[kc+1][jc][ic] > 0.1) ka = 0; \
  if (ja==-1 && nvert_c[kc][jc-1][ic] > 0.1) ja = 0; \
  else if (ja==1 && nvert_c[kc][jc+1][ic] > 0.1) ja = 0; \
  if (ia==-1 && nvert_c[kc][jc][ic-1] > 0.1) ia = 0; \
  else if (ia==1 && nvert_c[kc][jc][ic+1] > 0.1) ia = 0;

PetscInt lidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)
{
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;


  if (!(user->aotopetsc)) {
/*     PetscPrintf(PETSC_COMM_SELF, "ttt"); */
/*     PetscMalloc( (info.xm * info.ym * info.zm)*sizeof(int), &(user->idx_from)); */
/*     PetscInt kk, jj, ii, count; */
/*     count = 0; */
/*     for (kk=zs; kk<ze; kk++) { */
/*       for (jj=ys; jj<ye; jj++) { */
/* 	for (ii=xs; ii<xe; ii++) { */
/* 	  user->idx_from[count] = kk*mx*my + jj*mx + ii; count++; */
/* 	} */
/*       } */
/*     } */
/*     AO ao; */
/*     DAGetAO(user->da, &ao); */
/*     AOApplicationToPetsc(ao, count, user->idx_from); */
    user->aotopetsc = PETSC_TRUE;

    DAGetGlobalIndices(user->da, PETSC_NULL, &user->idx_from);
  }
  return (user->idx_from[(k-gzs) * (info.gxm*info.gym) + (j-gys)*(info.gxm) + (i-gxs)]);
  //  return (i + j * mx + k * mx * my);
}

PetscErrorCode mymatmultadd(Mat mat, Vec v1, Vec v2, Vec v3)
{

  Vec vt;
  VecDuplicate(v3, &vt);
  MatMult(mat, v1, vt);
  VecWAXPY(v3, 1., v2, vt);
  VecDestroy(vt);
  return(0);
}


PetscErrorCode mymult(Mat J, Vec X, Vec F)
{
  UserCtx *user;
  MatShellGetContext(J, &user);

  DA	da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;

  Vec F1;
  PetscReal ***f, ***f1, ***x;
  Cmpnts       ***icsi, ***ieta, ***izet, ***jcsi, ***jeta, ***jzet;
  Cmpnts       ***kcsi, ***keta, ***kzet;
  PetscReal    ***iaj, ***jaj, ***kaj, ***aj;

  PetscReal	dpdc, dpde, dpdz;

  Vec lX;

  DAGetLocalVector(da, &lX);
  DAGlobalToLocalBegin(da, X, INSERT_VALUES, lX);
  VecDuplicate(lX, &F1);
  DAVecGetArray(da, F,  &f);

  VecSet(F1, 0.);
  DAVecGetArray(da, F1, &f1);


  DAVecGetArray(fda, user->lICsi, &icsi);
  DAVecGetArray(fda, user->lIEta, &ieta);
  DAVecGetArray(fda, user->lIZet, &izet);

  DAVecGetArray(fda, user->lJCsi, &jcsi);
  DAVecGetArray(fda, user->lJEta, &jeta);
  DAVecGetArray(fda, user->lJZet, &jzet);

  DAVecGetArray(fda, user->lKCsi, &kcsi);
  DAVecGetArray(fda, user->lKEta, &keta);
  DAVecGetArray(fda, user->lKZet, &kzet);

  DAVecGetArray(da, user->lAj, &aj);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);

  lxs = xs-1; lxe = xe;
  lys = ys-1; lye = ye;
  lzs = zs-1; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  DAGlobalToLocalEnd(da, X, INSERT_VALUES, lX);
  DAVecGetArray(da, lX,  &x);

  PetscReal ***nvert;
  DAVecGetArray(da, user->lNvert, &nvert);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (i<mx-2) {
	  dpdc = x[k][j][i+1] - x[k][j][i];

	  if (j==1 || nvert[k][j-1][i] > 0.1
	      || nvert[k][j-1][i+1] > 0.1) {
	    dpde = (x[k][j+1][i] + x[k][j+1][i+1] -
		    x[k][j  ][i] - x[k][j  ][i+1]) * 0.5;
	  }
	  else if (j==my-2 || nvert[k][j+1][i] > 0.1
	      || nvert[k][j+1][i+1] > 0.1) {
	    dpde = (x[k][j  ][i] + x[k][j  ][i+1] -
		    x[k][j-1][i] - x[k][j-1][i+1]) * 0.5;
	  }
	  else {
	    dpde = (x[k][j+1][i] + x[k][j+1][i+1] -
		    x[k][j-1][i] - x[k][j-1][i+1]) * 0.25;
	  }

	  if (k==1 || nvert[k-1][j][i] > 0.1
	      || nvert[k-1][j][i+1] > 0.1) {
	    dpdz = (x[k+1][j][i] + x[k+1][j][i+1] -
		    x[k  ][j][i] - x[k  ][j][i+1]) * 0.5;
	  }
	  else if (k==mz-2 || nvert[k+1][j][i] > 0.1
	      || nvert[k+1][j][i+1] > 0.1) {
	    dpdz = (x[k  ][j][i] + x[k  ][j][i+1] -
		    x[k-1][j][i] - x[k-1][j][i+1]) * 0.5;
	  }
	  else {
	    dpdz = (x[k+1][j][i] + x[k+1][j][i+1] -
		    x[k-1][j][i] - x[k-1][j][i+1]) * 0.25;
	  }

	  f1[k][j][i] = 
	    dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		    icsi[k][j][i].y * icsi[k][j][i].y +
		    icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	    dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		    ieta[k][j][i].y * icsi[k][j][i].y +
		    ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] + 
	    dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		    izet[k][j][i].y * icsi[k][j][i].y +
		    izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	  if (nvert[k][j][i] + nvert[k][j][i+1] > 0.1) {
	    f1[k][j][i] = 0.;
	  }
	  
	}
	else {
	  f1[k][j][i] = 0.;
	}
      }
    }
  }

/*   DALocalToLocalBegin(da, F1, INSERT_VALUES, F1); */
/*   DALocalToLocalEnd(da, F1, INSERT_VALUES, F1); */

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (i==0) {
	  f[k][j][i] = x[k][j][i+1] - x[k][j][i];
	}
	else if (i==mx-1) {
	  f[k][j][i] = x[k][j][i-1] - x[k][j][i];
	}
	else if (j==0) {
	  f[k][j][i] = x[k][j+1][i] - x[k][j][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = x[k][j-1][i] - x[k][j][i];
	}
	else if (k==0) {
	  f[k][j][i] = x[k+1][j][i] - x[k][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = x[k-1][j][i] - x[k][j][i];
	}
	else {
	  f[k][j][i] = f1[k][j][i] - f1[k][j][i-1];
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (j<my-2) {
 	  if (i==1 || nvert[k][j][i-1] > 0.1
	      || nvert[k][j+1][i-1] > 0.1) {
	    dpdc = (x[k][j][i+1] + x[k][j+1][i+1] -
		    x[k][j][i  ] - x[k][j+1][i  ]) * 0.5;
	  }
	  else if (i==mx-2 || nvert[k][j][i+1] > 0.1
	      || nvert[k][j+1][i+1] > 0.1) {
	    dpdc = (x[k][j][i  ] + x[k][j+1][i  ] -
		    x[k][j][i-1] - x[k][j+1][i-1]) * 0.5;
	  }
	  else {
	    dpdc = (x[k][j][i+1] + x[k][j+1][i+1] -
		    x[k][j][i-1] - x[k][j+1][i-1]) * 0.25;
	  }

	  dpde = x[k][j+1][i] - x[k][j][i];

	  if (k==1 || nvert[k-1][j][i] > 0.1
	      || nvert[k-1][j+1][i] > 0.1) {
	    dpdz = (x[k+1][j][i] + x[k+1][j+1][i] -
		    x[k  ][j][i] - x[k  ][j+1][i]) * 0.5;
	  }
	  else if (k==mz-2 || nvert[k+1][j][i] > 0.1
	      || nvert[k+1][j+1][i] > 0.1) {
	    dpdz = (x[k  ][j][i] + x[k  ][j+1][i] -
		    x[k-1][j][i] - x[k-1][j+1][i]) * 0.5;
	  }
	  else {
	    dpdz = (x[k+1][j][i] + x[k+1][j+1][i] -
		    x[k-1][j][i] - x[k-1][j+1][i]) * 0.25;
	  }
	  f1[k][j][i] =
	    dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		    jcsi[k][j][i].y * jeta[k][j][i].y +
		    jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	    dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		    jeta[k][j][i].y * jeta[k][j][i].y +
		    jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	    dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		    jzet[k][j][i].y * jeta[k][j][i].y +
		    jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	  if (nvert[k][j+1][i] + nvert[k][j][i] > 0.1) {
	    f1[k][j][i] = 0.;
	  }
	}
	else {
	  f1[k][j][i] = 0.;
	}
      }
    }
  }

/*   DALocalToLocalBegin(da, F1, INSERT_VALUES, F1); */
/*   DALocalToLocalEnd(da, F1, INSERT_VALUES, F1); */
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (i==0) {
	  f[k][j][i] = x[k][j][i+1] - x[k][j][i];
	}
	else if (i==mx-1) {
	  f[k][j][i] = x[k][j][i-1] - x[k][j][i];
	}
	else if (j==0) {
	  f[k][j][i] = x[k][j+1][i] - x[k][j][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = x[k][j-1][i] - x[k][j][i];
	}
	else if (k==0) {
	  f[k][j][i] = x[k+1][j][i] - x[k][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = x[k-1][j][i] - x[k][j][i];
	}
	else {
	  f[k][j][i] += f1[k][j][i] - f1[k][j-1][i];
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (k < mz-2) {
	  if (i==1 || nvert[k][j][i-1] > 0.1
	      || nvert[k+1][j][i-1] > 0.1) {
	    dpdc = (x[k][j][i+1] + x[k+1][j][i+1] -
		    x[k][j][i  ] - x[k+1][j][i  ]) * 0.5;
	  }
	  else if (i==mx-2 || nvert[k][j][i+1] > 0.1
	      || nvert[k+1][j][i+1] > 0.1) {
	    dpdc = (x[k][j][i  ] + x[k+1][j][i  ] -
		    x[k][j][i-1] - x[k+1][j][i-1]) * 0.5;
	  }
	  else {
	    dpdc = (x[k][j][i+1] + x[k+1][j][i+1] -
		    x[k][j][i-1] - x[k+1][j][i-1]) * 0.25;
	  }

	  if (j==1 || nvert[k][j-1][i] > 0.1
	      || nvert[k+1][j-1][i] > 0.1) {
	    dpde = (x[k][j+1][i] + x[k+1][j+1][i] -
		    x[k][j  ][i] - x[k+1][j  ][i]) * 0.5;
	  }
	  else if (j==my-2 || nvert[k][j+1][i] > 0.1
	      || nvert[k+1][j+1][i] > 0.1) {
	    dpde = (x[k][j  ][i] + x[k+1][j  ][i] -
		    x[k][j-1][i] - x[k+1][j-1][i]) * 0.5;
	  }
	  else {
	    dpde = (x[k][j+1][i] + x[k+1][j+1][i] -
		    x[k][j-1][i] - x[k+1][j-1][i]) * 0.25;
	  }

	  dpdz = x[k+1][j][i] - x[k][j][i];
	  f1[k][j][i] =
	    dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		    kcsi[k][j][i].y * kzet[k][j][i].y +
		    kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	    dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		    keta[k][j][i].y * kzet[k][j][i].y +
		    keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	    dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		    kzet[k][j][i].y * kzet[k][j][i].y +
		    kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	  if (nvert[k][j][i] + nvert[k+1][j][i] > 0.1) {
	    f1[k][j][i] = 0.;
	  }
	}
	else {
	  f1[k][j][i] = 0.;
	}
      }
    }
  }
/*   DALocalToLocalBegin(da, F1, INSERT_VALUES, F1); */
/*   DALocalToLocalEnd(da, F1, INSERT_VALUES, F1); */

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (i==0) {
	  f[k][j][i] = x[k][j][i];//0.;// x[k][j][i+1] - x[k][j][i];
	}
	else if (i==mx-1) {
	  f[k][j][i] = x[k][j][i];//0.; //x[k][j][i-1] - x[k][j][i];
	}
	else if (j==0) {
	  f[k][j][i] = x[k][j][i];//0.;//x[k][j+1][i] - x[k][j][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = x[k][j][i];//0.;//x[k][j-1][i] - x[k][j][i];
	}
	else if (k==0) {
	  f[k][j][i] = x[k][j][i];//0.;//x[k+1][j][i] - x[k][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = x[k][j][i];//0.;//x[k-1][j][i] - x[k][j][i];
	}
	else {
	  f[k][j][i] += f1[k][j][i] - f1[k-1][j][i];
	  f[k][j][i] *= (-aj[k][j][i]);
	}
      }
    }
  }

/*   PetscReal norm; */
/*   VecMax(F, &i, &norm); */

/*   PetscPrintf(PETSC_COMM_WORLD, "VecMax %d %le\n", i, norm); */
/*   MatMult(user->A, X, F); */
/*   VecMax(F, &i, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "VecMax2 %d %le\n", i, norm); */

  DAVecRestoreArray(da, user->lNvert, &nvert);

  DAVecRestoreArray(da, F,  &f);
  DAVecRestoreArray(da, F1, &f1);
  DAVecRestoreArray(da, lX,  &x);

  DAVecRestoreArray(fda, user->lICsi, &icsi);
  DAVecRestoreArray(fda, user->lIEta, &ieta);
  DAVecRestoreArray(fda, user->lIZet, &izet);

  DAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DAVecRestoreArray(fda, user->lJEta, &jeta);
  DAVecRestoreArray(fda, user->lJZet, &jzet);

  DAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DAVecRestoreArray(fda, user->lKEta, &keta);
  DAVecRestoreArray(fda, user->lKZet, &kzet);

  DAVecRestoreArray(da, user->lAj, &aj);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);

  DARestoreLocalVector(da, &lX);
  VecDestroy(F1);
  return 0;
}
                                                                                                                                            
PetscErrorCode MyInterpolation(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, &user);


  
  DA	da = user->da, fda = user->fda;

  DA	da_f = user->da_f, da_c = user->da_c;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert, ***nvert_c;
  PetscInt i, j, k, ic, jc, kc, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


  //  ibm_interpolation_advanced_mg(user->user_c, user->user_c->ibm, X);

  DAVecGetArray(da,   F, &f);


  Vec lX;
  DACreateLocalVector(da_c, &lX);
  //  DAGetLocalVector(da_c, &lX);
  DAGlobalToLocalBegin(da_c, X, INSERT_VALUES, lX);
  DAGlobalToLocalEnd(da_c, X, INSERT_VALUES, lX);  
  DAVecGetArray(da_c, lX, &x);

  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(da_c, *(user->lNvert_c), &nvert_c);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);
/* 	  kc = (k + 1) / 2; */
/* 	  jc = (j + 1) / 2; */
/* 	  ic = (i + 1) / 2; */
/* 	  ka = (k - 2 * kc)==0 ? 1 : -1; */
/* 	  ja = (j - 2 * jc)==0 ? 1 : -1; */
/* 	  ia = (i - 2 * ic)==0 ? 1 : -1; */
/* 	  if (ka==-1 &&(k==1 || nvert_c[kc-1][jc][ic]>0.1)) ka = 0; */
/* 	  else if (ka==1 && (k==mz-2 || nvert_c[kc+1][jc][ic]>0.1)) ka=0; */

/* 	  if (ja==-1 &&(j==1 || nvert_c[kc][jc-1][ic]>0.1)) ja = 0; */
/* 	  else if (ja==1 &&(j==my-2 || nvert_c[kc][jc+1][ic]>0.1)) ja=0; */

/* 	  if (ia==-1 &&(i==1 || nvert_c[kc][jc][ic-1]>0.1)) ia = 0; */
/* 	  else if (ia==1 && (i==mx-2 || nvert_c[kc][jc][ic+1]>0.1)) ia=0; */

	  f[k][j][i] = (x[kc   ][jc   ][ic   ] * 9 +
			x[kc   ][jc+ja][ic   ] * 3 +
			x[kc   ][jc   ][ic+ia] * 3 +
			x[kc   ][jc+ja][ic+ia]) * 3./64. +
	    (x[kc+ka][jc   ][ic   ] * 9 +
	     x[kc+ka][jc+ja][ic   ] * 3 +
	     x[kc+ka][jc   ][ic+ia] * 3 +
	     x[kc+ka][jc+ja][ic+ia]) /64.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {

	if (i==0) {
	  f[k][j][i] = 0.;//-f[k][j][i+1];
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;//-f[k][j][i-1];
	}
	else if (j==0) {
	  f[k][j][i] = 0.;//-f[k][j+1][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;//-f[k][j-1][i];
	}
	else if (k==0) {
	  f[k][j][i] = 0.;//-f[k+1][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;//-f[k-1][j][i];
	}
	if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;
/* 	  f[k][j][i] = 0.125 * */
/* 	    (x[kc  ][jh  ][ih  ] + */
/* 	     x[kc  ][jh  ][ih-1] + */
/* 	     x[kh  ][jh-1][ih  ] + */
/* 	     x[kh-1][jh  ][ih  ] + */
/* 	     x[kh  ][jh-1][ih-1] + */
/* 	     x[kh-1][jh-1][ih  ] + */
/* 	     x[kh-1][jh  ][ih-1] + */
/* 	     x[kh-1][jh-1][ih-1]); */
/* 	} */
      }
    }
  }

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(da_c, *(user->lNvert_c), &nvert_c);

  DAVecRestoreArray(da_c, lX, &x);
  //  DARestoreLocalVector(da_c, &lX);
  VecDestroy(lX);
  DAVecRestoreArray(da,   F,  &f);

/*   PetscReal sum; */
/*   PetscInt N; */
/*   VecSum(F, &sum); */
/*   VecGetSize(F, &N); */
/*   sum = sum / (-1.*N); */
/*   VecShift(F, sum); */

/*   PetscReal max, min; */
/*   PetscInt  maxi, mini; */
/*   VecMax(X, &maxi, &max); */
/*   VecMin(X, &mini, &min); */
/*   PetscPrintf(PETSC_COMM_WORLD, "MM %d %le %d %le\n", maxi, max, mini, min); */
/*   PetscReal norm; */
/*   VecNorm(F, NORM_2, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Restriction Norm %le\n", norm); */

  //  MatNullSpaceRemove(da.nullsp, F, PETSC_NULL);
}

PetscErrorCode MyRestriction(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, &user);


  
  DA	da = user->da, fda = user->fda;

  DA	da_f = user->da_f, da_c = user->da_c;

  DALocalInfo	info;
  DAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert;
  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  DAVecGetArray(da,   F, &f);

  Vec lX;
  //  DAGetLocalVector(da_f, &lX);
  DACreateLocalVector(da_f, &lX);
  DAGlobalToLocalBegin(da_f, X, INSERT_VALUES, lX);
  DAGlobalToLocalEnd(da_f, X, INSERT_VALUES, lX);  
  DAVecGetArray(da_f, lX, &x);

  DAVecGetArray(da, user->lNvert, &nvert);

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (k==0) {
	  f[k][j][i] = 0.;
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;
	}
	else if (j==0) {
	  f[k][j][i] = 0.;
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;
	}
	else if (i==0) {
	  f[k][j][i] = 0.;
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;
	}
	else {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user);
	  f[k][j][i] = 0.125 *
	    (x[kh   ][jh   ][ih   ] +
	     x[kh   ][jh   ][ih-ia] +
	     x[kh   ][jh-ja][ih   ] +
	     x[kh-ka][jh   ][ih   ] +
	     x[kh   ][jh-ja][ih-ia] +
	     x[kh-ka][jh-ja][ih   ] +
	     x[kh-ka][jh   ][ih-ia] +
	     x[kh-ka][jh-ja][ih-ia]);
	  if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;
	}
      }
    }
  }


  DAVecRestoreArray(da_f, lX, &x);
  VecDestroy(lX);
  //  DARestoreLocalVector(da_f, &lX);
  DAVecRestoreArray(da,   F,  &f);
  DAVecRestoreArray(da, user->lNvert, &nvert);
/*   PetscReal norm; */
/*   VecNorm(F, NORM_2, &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Interpolation Norm %le\n", norm); */
  return 0;
}

PetscErrorCode MyNvertRestriction(UserCtx *user_h, UserCtx *user_c)
{
  DA		da = user_c->da, fda = user_c->fda;

  DALocalInfo	info = user_c->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i,j,k;
  PetscInt ih, jh, kh, ia, ja, ka;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert, ***nvert_h;

  DAVecGetArray(user_h->da, user_h->lNvert, &nvert_h);
  DAVecGetArray(user_c->da, user_c->Nvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if (*(user_c->isc)) ia = 0;
  else ia = 1;

  if (*(user_c->jsc)) ja = 0;
  else ja = 1;

  if (*(user_c->ksc)) ka = 0;
  else ka = 1;
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
	if (nvert_h[kh   ][jh   ][ih   ] +
	    nvert_h[kh   ][jh   ][ih-ia] +
	    nvert_h[kh   ][jh-ja][ih   ] +
	    nvert_h[kh-ka][jh   ][ih   ] +
	    nvert_h[kh   ][jh-ja][ih-ia] +
	    nvert_h[kh-ka][jh   ][ih-ia] +
	    nvert_h[kh-ka][jh-ja][ih   ] +
	    nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
	  nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
	}
      }
    }
  }
  DAVecRestoreArray(user_h->da, user_h->lNvert, &nvert_h);
  DAVecRestoreArray(user_c->da, user_c->Nvert, &nvert);

  DAGlobalToLocalBegin(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DAGlobalToLocalEnd(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DAVecGetArray(user_c->da, user_c->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j][i-1] > 1.1 &&
	      nvert[k][j+1][i] + nvert[k][j-1][i] > 1.1 &&
	      nvert[k+1][j][i] + nvert[k-1][j][i] > 1.1) {
	    nvert[k][j][i] = 1.;
	  }
	}
      }
    }
  }
/*   for (k=lzs; k<lze; k++) { */
/*     for (j=lys; j<lye; j++) { */
/*       for (i=lxs; i<lxe; i++) { */
/* 	if (i==1 && nvert[k][j][i+1] >0.1) nvert[k][j][i] = 1.; */
/* 	if (i==mx-2 && nvert[k][j][i-1] >0.1) nvert[k][j][i] = 1.; */
/* 	if (j==1 && nvert[k][j+1][i] >0.1) nvert[k][j][i] = 1.; */
/* 	if (j==my-2 && nvert[k][j-1][i] >0.1) nvert[k][j][i] = 1.; */
/* 	if (abs(nvert[k][j][i]) < 0.1) { */
/* 	  if (nvert[k][j][i+1]+nvert[k][j][i-1]>1.1 || */
/* 	      nvert[k][j+1][i]+nvert[k][j-1][i]>1.1 || */
/* 	      nvert[k+1][j][i]+nvert[k-1][j][i]>1.1) */
/* 	    nvert[k][j][i] = 1.; */
/* 	} */
/*       } */
/*     } */
/*   } */
  DAVecRestoreArray(user_c->da, user_c->lNvert, &nvert);
  DALocalToGlobal(user_c->da, user_c->lNvert, INSERT_VALUES, user_c->Nvert);
  return 0;
}

PetscErrorCode FullyBlocked(UserCtx *user)
{
  DA da = user->da, fda = user->fda;
  Vec nNvert;
  DALocalInfo info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;

  PetscInt *KSKE = user->KSKE;
  PetscReal ***nvert;
  PetscTruth *Blocked;

  DACreateNaturalVector(da, &nNvert);
  DAGlobalToNaturalBegin(da, user->Nvert, INSERT_VALUES, nNvert);
  DAGlobalToNaturalEnd(da, user->Nvert, INSERT_VALUES, nNvert);

  VecScatter ctx;
  Vec Zvert;
  VecScatterCreateToZero(nNvert, &ctx, &Zvert);

  VecScatterBegin(nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD, ctx);
  VecScatterEnd(nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD, ctx);

  VecScatterDestroy(ctx);
  VecDestroy(nNvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {

    VecGetArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    PetscMalloc(mx*my*sizeof(PetscTruth), &Blocked);
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	Blocked[j*mx+i] = PETSC_FALSE;
	for (k=0; k<mz; k++) {
	  if (nvert[k][j][i] > 0.1) {
	    if (!Blocked[j*mx+i]) {
	      KSKE[2*(j*mx+i)] = k;
	      Blocked[j*mx+i] = PETSC_TRUE;
	    }
	    else {
	      KSKE[2*(j*mx+i)] = PetscMin(KSKE[2*(j*mx+i)], k);
	    }
	  }
	}
      }
    }


    user->multinullspace = PETSC_TRUE;
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	if (!Blocked[j*mx+i]) {
	  user->multinullspace = PETSC_FALSE;
	  break;
	}
      }
    }
    PetscFree(Blocked);
    VecRestoreArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);

    }
  }
  else {
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }

  DACreateNaturalVector(da, &nNvert);
  VecDestroy(nNvert);

  VecDestroy(Zvert);
  return 0;
}



PetscErrorCode PoissonNullSpaceFunction(Vec X, void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;

  DA da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	***x, ***nvert;
  MatNullSpace 	nullsp;
  PetscInt	i, j, k;

/*   /\* First remove a constant from the Vec field X *\/ */
/*   MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp); */
/*   MatNullSpaceRemove(nullsp, X, PETSC_NULL); */
/*   MatNullSpaceDestroy(nullsp); */

  /* Then apply boundary conditions */
  DAVecGetArray(da, X, &x);
  DAVecGetArray(da, user->lNvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal lsum, sum;
  PetscInt  lnum, num;

  if (user->multinullspace) PetscPrintf(PETSC_COMM_WORLD, "MultiNullSpace!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  if (!user->multinullspace) {
    lsum = 0;
    lnum = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum = sum / (-1.0 * num);

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    x[k][j][i] +=sum;
	  }
	}
      }
    }
  }
  else {
    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }

    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
    PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD);
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }
			   
  }
  if (zs == 0) {
    k = 0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k = mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ys == 0) {
    j = 0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ye == my) {
    j = my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xs == 0) {
    i = 0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xe == mx) {
    i = mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert[k][j][i] > 0.1)
	  x[k][j][i] = 0.;
      }
    }
  }
  DAVecRestoreArray(da, X, &x);
  DAVecRestoreArray(da, user->lNvert, &nvert);

  return 0;
}

PetscErrorCode PoissonRHS(UserCtx *user, Vec B)
{
  DALocalInfo info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze= info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
  PetscReal	***nvert, ***aj, ***rb, dt = user->dt;
  Cmpnts	***ucont;

  DAVecGetArray(user->da, B, &rb);
  DAVecGetArray(user->fda, user->lUcont, &ucont);
  DAVecGetArray(user->da, user->lNvert, &nvert);
  DAVecGetArray(user->da, user->lAj, &aj);

  PetscReal coef = 1.;
  if (ti<2) coef = 1.5;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (i==0 || i==mx-1 || j==0 || j==my-1 ||
	    k==0 || k==mz-1) {
	  rb[k][j][i] = 0.;
	}
	else if ((int)(nvert[k][j][i]+0.5) !=0) {
	  rb[k][j][i] = 0;
	}
	else {
	  rb[k][j][i] = (-(ucont[k][j][i].x - ucont[k][j][i-1].x +
			   ucont[k][j][i].y - ucont[k][j-1][i].y +
			   ucont[k][j][i].z - ucont[k-1][j][i].z)
			/*  + user->FluxSumIntp */)
	    / dt * aj[k][j][i] / user->st * 1.5 / coef;
	    //- user->FluxIntpSum * 1.5/ dt;
	}

      }
    }
  }

  DAVecRestoreArray(user->fda, user->lUcont, &ucont);
  DAVecRestoreArray(user->da, user->lNvert, &nvert);
  DAVecRestoreArray(user->da, user->lAj, &aj);
  DAVecRestoreArray(user->da, B, &rb);

  return 0;
}

#define CP 0
#define EP 1
#define WP 2
#define NP 3
#define SP 4
#define TP 5
#define BP 6
#define NE 7
#define SE 8
#define NW 9
#define SW 10
#define TN 11
#define BN 12
#define TS 13
#define BS 14
#define TE 15
#define BE 16
#define TW 17
#define BW 18

PetscErrorCode PoissonLHSNew(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;

  PetscInt lxs, lxe, lys, lye, lzs, lze;
  PetscInt gxs, gxe, gys, gye, gzs, gze;

  Vec		G11, G12, G13, G21, G22, G23, G31, G32, G33;
  PetscReal	***g11, ***g12, ***g13, ***g21, ***g22, ***g23;
  PetscReal	***g31, ***g32, ***g33;

  PetscReal	***nvert, ***nvert_o;
  PetscScalar	vol[19];
  PetscInt	idx[19], row;

  PetscInt	i, j, k, N;
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (!user->assignedA) {
    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(user->P, &M);
    MatCreateMPIAIJ(PETSC_COMM_WORLD, M, M, N, N, 19, PETSC_NULL, 19, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  MatZeroEntries(user->A);

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

  DAVecGetArray(da, user->lAj, &aj);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);

  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(da, user->lNvert_o, &nvert_o);

  VecDuplicate(user->lAj, &G11);
  VecDuplicate(user->lAj, &G12);
  VecDuplicate(user->lAj, &G13);
  VecDuplicate(user->lAj, &G21);
  VecDuplicate(user->lAj, &G22);
  VecDuplicate(user->lAj, &G23);
  VecDuplicate(user->lAj, &G31);
  VecDuplicate(user->lAj, &G32);
  VecDuplicate(user->lAj, &G33);

  DAVecGetArray(da, G11, &g11);
  DAVecGetArray(da, G12, &g12);
  DAVecGetArray(da, G13, &g13);
  DAVecGetArray(da, G21, &g21);
  DAVecGetArray(da, G22, &g22);
  DAVecGetArray(da, G23, &g23);
  DAVecGetArray(da, G31, &g31);
  DAVecGetArray(da, G32, &g32);
  DAVecGetArray(da, G33, &g33);

  for (k=gzs; k<gze; k++) {
    for (j=gys; j<gye; j++) {
      for (i=gxs; i<gxe; i++) {
	g11[k][j][i] = (icsi[k][j][i].x * icsi[k][j][i].x +
			icsi[k][j][i].y * icsi[k][j][i].y +
			icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g12[k][j][i] = (ieta[k][j][i].x * icsi[k][j][i].x +
			ieta[k][j][i].y * icsi[k][j][i].y +
			ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g13[k][j][i] = (izet[k][j][i].x * icsi[k][j][i].x +
			izet[k][j][i].y * icsi[k][j][i].y +
			izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
					  			 
	g21[k][j][i] = (jcsi[k][j][i].x * jeta[k][j][i].x +
			jcsi[k][j][i].y * jeta[k][j][i].y +
			jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g22[k][j][i] = (jeta[k][j][i].x * jeta[k][j][i].x +
			jeta[k][j][i].y * jeta[k][j][i].y +
			jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g23[k][j][i] = (jzet[k][j][i].x * jeta[k][j][i].x +
			jzet[k][j][i].y * jeta[k][j][i].y +
			jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
					  			 
	g31[k][j][i] = (kcsi[k][j][i].x * kzet[k][j][i].x +
			kcsi[k][j][i].y * kzet[k][j][i].y +
			kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g32[k][j][i] = (keta[k][j][i].x * kzet[k][j][i].x +
			keta[k][j][i].y * kzet[k][j][i].y +
			keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g33[k][j][i] = (kzet[k][j][i].x * kzet[k][j][i].x +
			kzet[k][j][i].y * kzet[k][j][i].y +
			kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];	
      }
    }
  }

  PetscInt m;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = lidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  vol[CP] = 1.;	idx[CP] = lidx(i, j, k, user);
	  MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	}
	else {
	  if (nvert[k][j][i] > 0.1) { // i, j, k is not a fluid point
	    vol[CP] = 1.; idx[CP] = lidx(i, j, k, user);
	    MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	  }
	  else { // i, j, k is a fluid point
	    for (m=0; m<19; m++) {
	      vol[m] = 0.;
	    }
	    /* Contribution from i+1 - i */
	    if (nvert[k][j][i+1] < 2.1 && i != mx-2) { // i+1, j, k is a fluid point
	      /* dpdc{i} = (p_{i+1} - p_{i}) * g11_{i} */
	      vol[CP] -= g11[k][j][i]; //i, j, k
	      vol[EP] += g11[k][j][i]; // i+1, j, k

	      /* dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1})
	       * 0.25 * g12[k][j][i] */
	      if (j == my-2 || nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && j!=1) {
		  vol[CP] += g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] += g12[k][j][i] * 0.5; // i+1, j, k
		  vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
		  vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
		}
	      }
	      else if (j == 1 || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5;  //i, j+1, k
		  vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
		  vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
		}
	      }
	      else {
		vol[NP] += g12[k][j][i] * 0.25; // i, j+1, k
		vol[NE] += g12[k][j][i] * 0.25; // i+1, j+1, k
		vol[SP] -= g12[k][j][i] * 0.25; // i, j-1, k
		vol[SE] -= g12[k][j][i] * 0.25; // i+1, j-1, k
	      }
	      
	      /* dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if (k == mz-2 || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && k!=1) {
		  vol[CP] += g13[k][j][i] * 0.5; // i, j, k
		  vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
		  vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
		  vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
		}
	      }
	      else if (k == 1 || nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
		  vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
		  vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
		  vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
		}
	      }
	      else {
		vol[TP] += g13[k][j][i] * 0.25; //i, j, k+1
		vol[TE] += g13[k][j][i] * 0.25; //i+1, j, k+1
		vol[BP] -= g13[k][j][i] * 0.25; //i, j, k-1
		vol[BE] -= g13[k][j][i] * 0.25; //i+1, j, k-1
	      }
	    }  // end of i+1 - i

	    /* Contribution from i - i-1 */
	    if (nvert[k][j][i-1] < 2.1 && i != 1) { // i-1, j, k is a fluid point
	      /* -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i} */
	      vol[CP] -= g11[k][j][i-1];  //i, j, k
	      vol[WP] += g11[k][j][i-1];  //i-1, j, k

	      /* -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1})
	       * 0.25 * g12[k][j][i-1] */
	      if (j == my-2 || nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1 && j!=1) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
		  vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
		  vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
		}
	      }
	      else if (j == 1 || nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; // i, j+1, k
		  vol[NW] -= g12[k][j][i-1] * 0.5; // i-1, j+1, k
		  vol[CP] += g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g12[k][j][i-1] * 0.5; // i-1, j, k
		}
	      }
	      else {
		vol[NP] -= g12[k][j][i-1] * 0.25; // i, j+1, k
		vol[NW] -= g12[k][j][i-1] * 0.25; //i-1, j+1, k
		vol[SP] += g12[k][j][i-1] * 0.25; // i, j-1, k
		vol[SW] += g12[k][j][i-1] * 0.25; // i-1, j-1, k
	      }

	      /* -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if (k == mz-2 || nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1 && k!=1) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
		  vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
		  vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
		}
	      }
	      else if (k == 1 || nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; // i, j, k+1
		  vol[TW] -= g13[k][j][i-1] * 0.5; //i-1, j, k+1
		  vol[CP] += g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g13[k][j][i-1] * 0.5; //i-1, j, k
		}
	      }
	      else {
		vol[TP] -= g13[k][j][i-1] * 0.25;  // i, j, k+1
		vol[TW] -= g13[k][j][i-1] * 0.25; // i-1, j, k+1
		vol[BP] += g13[k][j][i-1] * 0.25;  // i, j, k-1
		vol[BW] += g13[k][j][i-1] * 0.25; // i-1, j, k-1
	      }
	    } // end of i - i-1

	    /* Contribution from j+1 - j */
	    if (nvert[k][j+1][i] < 2.1 && j != my-2) {
	      /* dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) *
		 0.25 */
	      if (i == mx-2 || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && i!=1) {
		  vol[CP] += g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
		  vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
		  vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
		}
	      }
	      else if (i == 1 || nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; // i+1, j, k
		  vol[NE] += g21[k][j][i] * 0.5; // i+1, j+1, k
		  vol[CP] -= g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] -= g21[k][j][i] * 0.5; // i, j+1, k
		}
	      }
	      else {
		vol[EP] += g21[k][j][i] * 0.25; //i+1, j, k
		vol[NE] += g21[k][j][i] * 0.25; //i+1, j+1, k
		vol[WP] -= g21[k][j][i] * 0.25; //i-1, j, k
		vol[NW] -= g21[k][j][i] * 0.25; //i-1, j+1, k
	      }
	     
	      /* dpde{j} = (p{j+1} - p{j}) * g22[k][j][i] */
	      vol[CP] -= g22[k][j][i];
	      vol[NP] += g22[k][j][i];

	      /* dpdz{j} = (p{j, k+1}+p{j+1,k+1} - p{j,k-1}-p{j+1,k-1}) *0.25*/
	      if (k == mz-2 || nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && k!=1) {
		  vol[CP] += g23[k][j][i] * 0.5; //i,j,k
		  vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
		  vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
		  vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
		}
	      }
	      else if (k == 1 || nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; //i, j, k+1
		  vol[TN] += g23[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g23[k][j][i] * 0.5;//i, j, k
		  vol[NP] -= g23[k][j][i] * 0.5;//i, j+1, k
		}
	      }
	      else {
		vol[TP] += g23[k][j][i] * 0.25; // i, j, k+1
		vol[TN] += g23[k][j][i] * 0.25; // i, j+1, k+1
		vol[BP] -= g23[k][j][i] * 0.25; // i, j, k-1
		vol[BN] -= g23[k][j][i] * 0.25; // i, j+1, k-1
	      }
	    } // End of j+1 - j

	    /* Contribution j - j-1 */
	    if (nvert[k][j-1][i] < 2.1 && j != 1) {
	      /* -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) *
		 0.25 * g21[k][j-1][i] */
	      if (i == mx-2 || nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1 && i!=1) {
		  vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
		  vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
		  vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
		  vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
		}
	      }
	      else if (i == 1 || nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5;//i+1, j, k
		  vol[SE] -= g21[k][j-1][i] * 0.5;//i+1, j-1, k
		  vol[CP] += g21[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g21[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[EP] -= g21[k][j-1][i] * 0.25;// i+1, j, k
		vol[SE] -= g21[k][j-1][i] * 0.25;// i+1, j-1, k
		vol[WP] += g21[k][j-1][i] * 0.25;// i-1, j, k
		vol[SW] += g21[k][j-1][i] * 0.25;// i-1, j-1, k
	      }
	      
	      /* -dpde{j-1} = -(p{j} - p{j-1}) * g22[k][j-1][i] */
	      vol[CP] -= g22[k][j-1][i];
	      vol[SP] += g22[k][j-1][i];

	      /* -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) *
		 0.25 * g23[k][j-1][i] */
	      if (k == mz-2 || nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 && k!=1) {
		  vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
		  vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
		  vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if (k == 1 || nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5;//i, j, k+1
		  vol[TS] -= g23[k][j-1][i] * 0.5;//i, j-1, k+1
		  vol[CP] += g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g23[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[TP] -= g23[k][j-1][i] * 0.25;//i, j, k+1
		vol[TS] -= g23[k][j-1][i] * 0.25;//i, j-1, k+1
		vol[BP] += g23[k][j-1][i] * 0.25;//i, j, k-1
		vol[BS] += g23[k][j-1][i] * 0.25;//i, j-1, k-1
	      }
	    } // End of j - j-1

	    /* contribution from k+1 - k */
	    if (nvert[k+1][j][i] < 2.1 && k != mz-2) {
	      /* dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) *
		 0.25 * g31[k][j][i] */
	      if (i == mx-2 || nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && i!=1) {
		  vol[CP] += g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
		  vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
		  vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
		}
	      }
	      else if (i == 1 || nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5;//i+1, j, k
		  vol[TE] += g31[k][j][i] * 0.5;//i+1, j, k+1
		  vol[CP] -= g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g31[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[EP] += g31[k][j][i] * 0.25;//i+1, j, k
		vol[TE] += g31[k][j][i] * 0.25;//i+1, j, k+1
		vol[WP] -= g31[k][j][i] * 0.25;//i-1, j, k
		vol[TW] -= g31[k][j][i] * 0.25;//i-1, j, k+1
	      }

	      /* dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 
		 0.25 * g32[k][j][i] */
	      if (j == my-2 || nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && j!=1) {
		  vol[CP] += g32[k][j][i] * 0.5;//i, j,k
		  vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
		  vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
		  vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
		}
	      }
	      else if (j == 1 || nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5;//i, j+1, k
		  vol[TN] += g32[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g32[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g32[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[NP] += g32[k][j][i] * 0.25;//i, j+1, k
		vol[TN] += g32[k][j][i] * 0.25;//i, j+1, k+1
		vol[SP] -= g32[k][j][i] * 0.25;//i, j-1, k
		vol[TS] -= g32[k][j][i] * 0.25;//i, j-1, k+1
	      }

	      /* dpdz{k} = p{k+1} - p{k} */
	      vol[CP] -= g33[k][j][i]; //i, j, k
	      vol[TP] += g33[k][j][i]; //i, j, k+1
	    } // End of k+1 - k

	    /* Contribution from k - k-1 */
	    if (nvert[k-1][j][i] < 2.1 && k != 1) {
	      /* -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) *
		 0.25 * g31[k-1][j][i] */
	      if (i == mx-2 || nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1 && i!=1) {
		  vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
		  vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
		  vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
		}
	      }
	      else if (i == 1 || nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5;//i+1, j, k
		  vol[BE] -= g31[k-1][j][i] * 0.5;//i+1, j, k-1
		  vol[CP] += g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g31[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[EP] -= g31[k-1][j][i] * 0.25;//i+1, j, k
		vol[BE] -= g31[k-1][j][i] * 0.25;//i+1, j, k-1
		vol[WP] += g31[k-1][j][i] * 0.25;//i-1, j, k
		vol[BW] += g31[k-1][j][i] * 0.25;//i-1, j, k-1
	      }
	      
	      /* -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) * 
		 0.25 * g32[k-1][j][i] */
	      if (j == my-2 || nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1 && j!=1) {
		  vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
		  vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
		  vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
		  vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if (j == 1 || nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5;//i, j+1, k
		  vol[BN] -= g32[k-1][j][i] * 0.5;//i, j+1, k-1
		  vol[CP] += g32[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g32[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[NP] -= g32[k-1][j][i] * 0.25;//i, j+1, k
		vol[BN] -= g32[k-1][j][i] * 0.25;//i, j+1, k-1
		vol[SP] += g32[k-1][j][i] * 0.25;//i, j-1, k
		vol[BS] += g32[k-1][j][i] * 0.25;//i, j-1, k-1
	      }
	      
	      /* -dpdz{k-1} = -(p{k} - p{k-1}) * g33[k-1][j][i] */
	      vol[CP] -= g33[k-1][j][i]; // i, j, k
	      vol[BP] += g33[k-1][j][i]; //i, j, k-1
	    } // End of k - k-1
	    for (m=0; m<19; m++) {
	      vol[m] *= -aj[k][j][i];
	    }

	    idx[CP] = lidx(i  , j  , k  , user);
	    idx[EP] = lidx(i+1, j  , k  , user);
	    idx[WP] = lidx(i-1, j  , k  , user);
	    idx[NP] = lidx(i  , j+1, k  , user);
	    idx[SP] = lidx(i  , j-1, k  , user);
	    idx[TP] = lidx(i  , j  , k+1, user);
	    idx[BP] = lidx(i  , j  , k-1, user);
	    idx[NE] = lidx(i+1, j+1, k  , user);
	    idx[SE] = lidx(i+1, j-1, k  , user);
	    idx[NW] = lidx(i-1, j+1, k  , user);
	    idx[SW] = lidx(i-1, j-1, k  , user);
	    idx[TN] = lidx(i  , j+1, k+1, user);
	    idx[BN] = lidx(i  , j+1, k-1, user);
	    idx[TS] = lidx(i  , j-1, k+1, user);
	    idx[BS] = lidx(i  , j-1, k-1, user);
	    idx[TE] = lidx(i+1, j  , k+1, user);
	    idx[BE] = lidx(i+1, j  , k-1, user);
	    idx[TW] = lidx(i-1, j  , k+1, user);
	    idx[BW] = lidx(i-1, j  , k-1, user);
	    MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);

/* 	    if (user->thislevel == 2) { */
/* 	      if (row==179161) { */
/* 		PetscInt tpa, tpb, tpc; */
/* 		for (tpa=k-1; tpa<k+2; tpa++) { */
/* 		  for (tpb=j-1; tpb<j+2; tpb++) { */
/* 		    for (tpc=i-1; tpc<i+2; tpc++) { */
/* 		      PetscPrintf(PETSC_COMM_SELF, "TT %i %i %i %e\n", tpa, tpb, tpc, nvert[tpa][tpb][tpc]); */
/* 		    } */
/* 		  } */
/* 		} */
/* 	      } */
/* 	    } */
	  } // End of fluid point
	} // End of interial points
      }
    }
  }
  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  DAVecRestoreArray(da, G11, &g11);
  DAVecRestoreArray(da, G12, &g12);
  DAVecRestoreArray(da, G13, &g13);
  DAVecRestoreArray(da, G21, &g21);
  DAVecRestoreArray(da, G22, &g22);
  DAVecRestoreArray(da, G23, &g23);
  DAVecRestoreArray(da, G31, &g31);
  DAVecRestoreArray(da, G32, &g32);
  DAVecRestoreArray(da, G33, &g33);
  
  VecDestroy(G11);
  VecDestroy(G12);
  VecDestroy(G13);
  VecDestroy(G21);
  VecDestroy(G22);
  VecDestroy(G23);
  VecDestroy(G31);
  VecDestroy(G32);
  VecDestroy(G33);


  //  VecCopy(user->Phi, user->P);
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

  DAVecRestoreArray(da, user->lAj, &aj);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);

  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(da, user->lNvert_o, &nvert_o);
  return 0;
}


/* PetscErrorCode PoissonSolver_MG(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo) */
/* { */
/*   PetscInt l; */
/*   PetscInt bi; */

/*   MGCtx *mgctx = usermg->mgctx; */

/*   KSP	*mgksp, subksp, csksp; */
/*   PC	*mgpc, subpc; */
/*   UserCtx	*user; */

/*   PetscInt	m_c, m_f, M_c, M_f; */

/*   PetscMalloc(block_number*sizeof(KSP), &mgksp); */
/*   PetscMalloc(block_number*sizeof(PC), &mgpc); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     /\* Create ksp for Multigrid *\/ */
/*     KSPCreate(PETSC_COMM_WORLD, &mgksp[bi]); */
    
/*     /\* Use multigrid as preconditioner *\/ */
/*     KSPGetPC(mgksp[bi], &mgpc[bi]); */
/*     PCSetType(mgpc[bi], PCMG); */
/*     /\* Setup MG levels from usergm->mglevels *\/ */
/*     PCMGSetLevels(mgpc[bi], usermg->mglevels, PETSC_NULL); */
/*     /\* V cycle *\/ */
/*     PCMGSetCycles(mgpc[bi], 1); */
    
/*     PCMGSetType(mgpc[bi], PC_MG_MULTIPLICATIVE); */
/* /\*     PCMGSetType(mgpc[bi], PC_MG_FULL); *\/ */

/*     /\* Create Restriction and Interpolate schemes  */
/*        This is needed for all levels other than the coarsest one *\/ */
/*     for (l=usermg->mglevels-1; l>0; l--) { */
/*       user = mgctx[l].user; */
/*       m_c = (usermg->mgctx[l-1].user[bi].info.xm * */
/* 	     usermg->mgctx[l-1].user[bi].info.ym * */
/* 	     usermg->mgctx[l-1].user[bi].info.zm); */

/*       m_f = (usermg->mgctx[l].user[bi].info.xm * */
/* 	     usermg->mgctx[l].user[bi].info.ym * */
/* 	     usermg->mgctx[l].user[bi].info.zm); */

/*       M_c = (usermg->mgctx[l-1].user[bi].info.mx * */
/* 	     usermg->mgctx[l-1].user[bi].info.my * */
/* 	     usermg->mgctx[l-1].user[bi].info.mz); */

/*       M_f = (usermg->mgctx[l].user[bi].info.mx * */
/* 	     usermg->mgctx[l].user[bi].info.my * */
/* 	     usermg->mgctx[l].user[bi].info.mz); */

/*       MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, (void*)&mgctx[l-1].user[bi], &user[bi].MR); */
/*       MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, (void*)&user[bi], &user[bi].MP); */

/*       PCMGSetRestriction(mgpc[bi], l, user[bi].MR); */
/*       PCMGSetInterpolate(mgpc[bi], l, user[bi].MP); */

/*       /\* Use subroutine MyRestriction and MyInterpolation for  */
/* 	 Mat * Vec operation *\/ */
/*       MatShellSetOperation(user[bi].MR, MATOP_MULT, (void(*)(void))MyRestriction); */
/*       MatShellSetOperation(user[bi].MP, MATOP_MULT, (void(*)(void))MyInterpolation); */

/*       MatShellSetOperation(user[bi].MR, MATOP_MULT_ADD,(void(*)(void))mymatmultadd); */
/*       MatShellSetOperation(user[bi].MP, MATOP_MULT_ADD,(void(*)(void))mymatmultadd); */

/*     } */
/*   } */

/*   if (immersed) { */
/*     for (l=usermg->mglevels-1; l>0; l--) { */
/*       for (bi=0; bi<block_number; bi++) { */
/* 	mgctx[l].user[bi].multinullspace = PETSC_FALSE; */
/* 	MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]); */
/*       } */
/*     } */
/*     /\* At the corsest level, check whether the grid is separated into sevearal  */
/*        blockes by the immersed body *\/ */
/*     l = 0; */
/*     user = mgctx[l].user; */
/*     for (bi=0; bi<block_number; bi++) { */
/*       PetscMalloc(user[bi].info.mx*user[bi].info.my*2*sizeof(PetscInt), &user[bi].KSKE); */
/*       FullyBlocked(&user[bi]); */

/*     } */
/*   } */

/*   for (l=usermg->mglevels-1; l>=0; l--) { */
/*     user = mgctx[l].user; */
/*     for (bi=0; bi<block_number; bi++) { */
/*       PoissonLHSNew(&(user[bi]), ibm, user[bi].ibm_intp); */


/* /\*       if (l<usermg->mglevels-1) { *\/ */
/* /\* 	VecDuplicate(user[bi].P, &user[bi].Rhsp); *\/ */
/* /\* 	PCMGSetRhs(mgpc[bi], l, user[bi].Rhsp); *\/ */
/* /\*       } *\/ */
/* /\*       if (l) { *\/ */
/* /\* 	VecDuplicate(user[bi].P, &user[bi].R); *\/ */
/* /\* 	PCMGSetR(mgpc[bi], l, user[bi].R); *\/ */
/* /\*       } *\/ */

      
/*       if (l) { /\* Not the coarset grid level *\/ */
/* 	PCMGGetSmoother(mgpc[bi], l, &subksp); */
/* 	/\* Set the left hand side for every KSP at each grid level *\/ */
/* 	KSPSetOperators(subksp, user[bi].A, user[bi].A, */
/* 			DIFFERENT_NONZERO_PATTERN); */
      
/* 	KSPSetFromOptions(subksp); */
/* 	KSPSetUp(subksp); */
/* 	KSPGetPC(subksp, &subpc); */
/* 	PCSetType(subpc,PCBJACOBI); */
/* 	PCSetUp(subpc); */

/* 	KSP *subsubksp; */
/* 	PC subsubpc; */
/* 	PetscInt abi, nlocal; */
/* 	PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp); */
/* 	for (abi = 0; abi<nlocal; abi++) { */
/* 	  KSPGetPC(subsubksp[abi], &subsubpc); */
/* 	  PCFactorSetShiftNonzero(subsubpc, 1.e-9); */
/* 	} */

/* 	PCFactorSetShiftNonzero(subpc, 1.e-9); */
/*       } */
/*       else {  /\* Coarsest grid *\/ */

/* 	/\* The default solver for the coarset grid is */
/* 	   KSPPreonly and PCLU. */
/* 	   One can choose other solvers, such as a KSP solver by changing the */
/* 	   following lines. *\/ */
/* 	PCMGGetCoarseSolve(mgpc[bi], &subksp); */
/* /\* 	KSPSetType(subksp, KSPBCGS); *\/ */
/* 	KSPGetPC(subksp, &subpc); */
/* 	PCSetType(subpc, PCBJACOBI); */
/* 	PetscInt cs_max_it=25; */
/* 	PetscOptionsGetInt(PETSC_NULL, "-mg-coarse_ksp_max_it", &cs_max_it, PETSC_NULL); */
/* 	KSPSetTolerances(subksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, cs_max_it); */
	

/* 	KSPSetOperators(subksp, user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN); */
/* 	KSPSetUp(subksp); */
/* 	PCFactorSetShiftNonzero(subpc, 1.e-9); */

/* /\* 	KSPGetPC(subksp, &subpc); *\/ */
/* /\* 	PCFactorSetShiftNonzero(subpc, 1.e-9); *\/ */
/* /\* 	PCSetUp(subpc); *\/ */
/*       } */

/*       /\* The Poisson equation has Neumann boundary conditions, thus */
/* 	 need to use NullSpace to solve the equation. */
/* 	 We use the subroutine PoissonNullSpaceFunction to  */
/* 	 (1) reduce a constant value from the solution */
/* 	 (2) Set the solution values at boundaries and blanking nodes to */
/* 	 zero *\/ */
/*       MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp); */
/*       MatNullSpaceSetFunction(user[bi].nullsp, PoissonNullSpaceFunction, &user[bi]); */
    
/*       KSPSetNullSpace(subksp, user[bi].nullsp); */

/*       PCMGSetResidual(mgpc[bi], l, PCMGDefaultResidual, user[bi].A); */

/*       KSPSetUp(subksp); */
/*     } */
/*   } */

/*   l = usermg->mglevels-1; */
/*   user = mgctx[l].user; */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecDuplicate(user[bi].P, &user[bi].B); */
/*     InterpolationFlux(&user[bi]); */
/*     PoissonRHS(&(user[bi]), user[bi].B); */
/*   } */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     KSPSetOperators(mgksp[bi], user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN); */
/*     KSPSetFromOptions(mgksp[bi]); */
/*     KSPSetUp(mgksp[bi]); */

/*     PetscReal norm; */
/*     VecNorm(user[bi].B, NORM_2, &norm); */
/*     PetscPrintf(PETSC_COMM_WORLD, "KSP RHS Norm %e\n", norm); */
/*     //    KSPSetMonitor(mgksp[bi], MyKSPMonitor, &user[bi], PETSC_NULL); */
/*     KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi); */
/*     Vec Temp; */
/*     VecDuplicate(user[bi].B, &Temp); */
/*     MatMult(user[bi].A, user[bi].Phi, Temp); */
/*     VecAXPY(Temp, -1, user[bi].B); */
/*     VecNorm(Temp, NORM_2, &norm); */
/*     PetscInt itn; */
/*     KSPGetIterationNumber(mgksp[bi], &itn); */
/*     PetscPrintf(PETSC_COMM_WORLD, "PoissonSolver Iteration Number %i Error %e\n", itn, norm); */
/*     VecDestroy(Temp); */

/* /\*     PetscTruth flg = PETSC_FALSE; *\/ */
/*  /\*     PetscOptionsGetTruth(PETSC_NULL, "-after_mg", &flg, PETSC_NULL); *\/ */
/* /\*     if (flg) { *\/ */
/* /\*       PCSetType(mgpc[bi], PCBJACOBI); *\/ */
/* /\*       KSPSetInitialGuessNonzero(mgksp[bi], PETSC_TRUE); *\/ */
/* /\*       KSPSetTolerances(mgksp[bi], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10); *\/ */
/* /\*       KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi); *\/ */
/* /\*     } *\/ */
/*   } */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DAGlobalToLocalBegin(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi); */
/*     DAGlobalToLocalEnd(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi); */
/*   } */

/*   /\* Release the allocated spaces *\/ */
/*   for (l=usermg->mglevels-1; l>=0; l--) { */
/*     user = mgctx[l].user; */
/*     for (bi=0; bi<block_number; bi++) { */
/* /\*       if (l<usermg->mglevels-1){ *\/ */
/* /\* 	VecDestroy(user[bi].Rhsp); *\/ */
/* /\*       } *\/ */
/* /\*       if (l) { *\/ */
/* /\* 	VecDestroy(user[bi].R); *\/ */
/* /\*       } *\/ */
/*       MatNullSpaceDestroy(user[bi].nullsp); */

/*       MatDestroy(user[bi].A); */
/*       user[bi].assignedA = PETSC_FALSE; */
/*       if (l) { /\* Grid level > 0 (coarsest) *\/ */
/* 	MatDestroy(user[bi].MR); */
/* 	MatDestroy(user[bi].MP); */
/*       } */
/*       else { */
/* 	PetscFree(user[bi].KSKE); */
/*       } */
/*     } */
/*   } */
/*   for (bi=0; bi<block_number; bi++) { */
/*     //    PCDestroy(mgpc[bi]); */
/*     KSPDestroy(mgksp[bi]); */
/*     VecDestroy(mgctx[usermg->mglevels-1].user[bi].B); */
/*   } */
/* /\*   PetscFree(mgksp); *\/ */
/* /\*   PetscFree(mgpc); *\/ */

/*   return 0; */
/* } */

PetscErrorCode PoissonSolver_MG(UserMG *usermg, IBMNodes *ibm, IBMInfo *ibminfo)
{
  PetscInt l;
  PetscInt bi;

  MGCtx *mgctx = usermg->mgctx;

  KSP	*mgksp, subksp, csksp;
  PC	*mgpc, subpc;
  UserCtx	*user;

  PetscInt	m_c, m_f, M_c, M_f;

  PetscMalloc(block_number*sizeof(KSP), &mgksp);
  PetscMalloc(block_number*sizeof(PC), &mgpc);

  l = usermg->mglevels-1;
  user = mgctx[l].user;
  for (bi=0; bi<block_number; bi++) {
    VecDuplicate(user[bi].P, &user[bi].B);
    PoissonRHS(&(user[bi]), user[bi].B);
  }

  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      PoissonLHSNew(&(user[bi]), ibm, user[bi].ibm_intp);
    }
  }

/*   l = usermg->mglevels-1; */
/*   user = mgctx[l].user; */
/*   for (bi=0; bi<block_number; bi++) { */
/*     KSPCreate(PETSC_COMM_WORLD, &mgksp[bi]); */
/*     KSPSetOperators(mgksp[bi], user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN); */
/*     KSPSetFromOptions(mgksp[bi]); */
/*     KSPSetUp(mgksp[bi]); */
/*     KSPSetTolerances(mgksp[bi], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10); */

/*     KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi); */

/*     KSPDestroy(mgksp[bi]); */
/*   } */
  for (bi=0; bi<block_number; bi++) {
    /* Create ksp for Multigrid */
    KSPCreate(PETSC_COMM_WORLD, &mgksp[bi]);
    KSPAppendOptionsPrefix(mgksp[bi], "ps_");
/*     KSPSetInitialGuessNonzero(mgksp[bi], PETSC_TRUE); */
    
    /* Use multigrid as preconditioner */
    KSPGetPC(mgksp[bi], &mgpc[bi]);
    PCSetType(mgpc[bi], PCMG);
    /* Setup MG levels from usergm->mglevels */
    PCMGSetLevels(mgpc[bi], usermg->mglevels, PETSC_NULL);
    /* V cycle */
    PCMGSetCycles(mgpc[bi], 1);
    
    PCMGSetType(mgpc[bi], PC_MG_MULTIPLICATIVE);
/*     PCMGSetType(mgpc[bi], PC_MG_FULL); */

    /* Create Restriction and Interpolate schemes 
       This is needed for all levels other than the coarsest one */
    for (l=usermg->mglevels-1; l>0; l--) {
      user = mgctx[l].user;
      m_c = (usermg->mgctx[l-1].user[bi].info.xm *
	     usermg->mgctx[l-1].user[bi].info.ym *
	     usermg->mgctx[l-1].user[bi].info.zm);

      m_f = (usermg->mgctx[l].user[bi].info.xm *
	     usermg->mgctx[l].user[bi].info.ym *
	     usermg->mgctx[l].user[bi].info.zm);

      M_c = (usermg->mgctx[l-1].user[bi].info.mx *
	     usermg->mgctx[l-1].user[bi].info.my *
	     usermg->mgctx[l-1].user[bi].info.mz);

      M_f = (usermg->mgctx[l].user[bi].info.mx *
	     usermg->mgctx[l].user[bi].info.my *
	     usermg->mgctx[l].user[bi].info.mz);

      MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, (void*)&mgctx[l-1].user[bi], &user[bi].MR);
      MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, (void*)&user[bi], &user[bi].MP);

      PCMGSetRestriction(mgpc[bi], l, user[bi].MR);
      PCMGSetInterpolate(mgpc[bi], l, user[bi].MP);

      /* Use subroutine MyRestriction and MyInterpolation for 
	 Mat * Vec operation */
      MatShellSetOperation(user[bi].MR, MATOP_MULT, (void(*)(void))MyRestriction);
      MatShellSetOperation(user[bi].MP, MATOP_MULT, (void(*)(void))MyInterpolation);

      MatShellSetOperation(user[bi].MR, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);
      MatShellSetOperation(user[bi].MP, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);

    }
  }


  if (immersed) {
    for (l=usermg->mglevels-1; l>0; l--) {
      for (bi=0; bi<block_number; bi++) {
	mgctx[l].user[bi].multinullspace = PETSC_FALSE;
	MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]);
      }
    }
    /* At the corsest level, check whether the grid is separated into sevearal 
       blockes by the immersed body */
    l = 0;
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      PetscMalloc(user[bi].info.mx*user[bi].info.my*2*sizeof(PetscInt), &user[bi].KSKE);
      FullyBlocked(&user[bi]);

    }
  }

  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      
      if (l) { /* Not the coarset grid level */
	PCMGGetSmoother(mgpc[bi], l, &subksp);
	/* Set the left hand side for every KSP at each grid level */
	KSPSetOperators(subksp, user[bi].A, user[bi].A,
			DIFFERENT_NONZERO_PATTERN);
      
/* 	KSPSetFromOptions(subksp); */
/* 	KSPSetUp(subksp); */
/* 	KSPGetPC(subksp, &subpc); */
/* 	PCFactorSetShiftNonzero(subpc, 1.e-9); */
      }
      else {  /* Coarsest grid */

	/* The default solver for the coarset grid is
	   KSPPreonly and PCLU.
	   One can choose other solvers, such as a KSP solver by changing the
	   following lines. */
	PCMGGetCoarseSolve(mgpc[bi], &subksp);
/* 	KSPSetType(subksp, KSPBCGS); */
	KSPGetPC(subksp, &subpc);
	PCSetType(subpc, PCBJACOBI);
	KSPSetTolerances(subksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 40);
	

	KSPSetOperators(subksp, user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);
	KSPSetUp(subksp);
/* 	KSPGetPC(subksp, &subpc); */
/* 	PCFactorSetShiftNonzero(subpc, 1.e-9); */
/* 	PCSetUp(subpc); */
      }

      /* The Poisson equation has Neumann boundary conditions, thus
	 need to use NullSpace to solve the equation.
	 We use the subroutine PoissonNullSpaceFunction to 
	 (1) reduce a constant value from the solution
	 (2) Set the solution values at boundaries and blanking nodes to
	 zero */
      MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp);
      MatNullSpaceSetFunction(user[bi].nullsp, PoissonNullSpaceFunction, &user[bi]);
    
/*       KSPSetNullSpace(subksp, user[bi].nullsp); */

      PCMGSetResidual(mgpc[bi], l, PCMGDefaultResidual, user[bi].A);

      KSPSetUp(subksp);

      if (l<usermg->mglevels-1) {
	MatGetVecs(user[bi].A, &user[bi].R, PETSC_NULL);
	//	VecDuplicate(user[bi].P, &user[bi].R);
	PCMGSetRhs(mgpc[bi], l, user[bi].R);
      }

    }
  }

  l = usermg->mglevels-1;
  user = mgctx[l].user;
  
  for (bi=0; bi<block_number; bi++) {
    KSPSetOperators(mgksp[bi], user[bi].A, user[bi].A, DIFFERENT_NONZERO_PATTERN);
    KSPSetFromOptions(mgksp[bi]);
    KSPSetUp(mgksp[bi]);

    PetscReal norm;
    VecNorm(user[bi].B, NORM_2, &norm);
    PetscPrintf(PETSC_COMM_WORLD, "KSP RHS Norm %e\n", norm);
    //    KSPSetMonitor(mgksp[bi], MyKSPMonitor, &user[bi], PETSC_NULL);
    KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi);

/*     PetscTruth flg = PETSC_FALSE; */
/*     PetscOptionsGetTruth(PETSC_NULL, "-after_mg", &flg, PETSC_NULL); */
/*     if (flg) { */
/*       PCSetType(mgpc[bi], PCBJACOBI); */
/*       KSPSetInitialGuessNonzero(mgksp[bi], PETSC_TRUE); */
/*       KSPSetTolerances(mgksp[bi], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10); */
/*       KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi); */
/*     } */
  }

  for (bi=0; bi<block_number; bi++) {
    DAGlobalToLocalBegin(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
    DAGlobalToLocalEnd(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
  }

  for (bi=0; bi<block_number; bi++) {
    //    PCDestroy(mgpc[bi]);
    for (l=0; l<usermg->mglevels-1; l++) {
      VecDestroy(mgctx[l].user[bi].R);
    }

    KSPDestroy(mgksp[bi]);


    VecDestroy(mgctx[usermg->mglevels-1].user[bi].B);
  }

  /* Release the allocated spaces */
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      MatNullSpaceDestroy(user[bi].nullsp);

      MatDestroy(user[bi].A);
      user[bi].assignedA = PETSC_FALSE;
      if (l) { /* Grid level > 0 (coarsest) */
	MatDestroy(user[bi].MR);
	MatDestroy(user[bi].MP);
      }
      else {
	PetscFree(user[bi].KSKE);
      }
    }
  }

  PetscFree(mgksp);
  PetscFree(mgpc);

  return 0;
}


PetscErrorCode Projection(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;
  PetscReal	***p;
  PetscReal	dt = user->dt;

  Cmpnts	***ucont;
  PetscReal	***nvert;

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

  DAVecGetArray(fda, user->lICsi, &icsi);
  DAVecGetArray(fda, user->lIEta, &ieta);
  DAVecGetArray(fda, user->lIZet, &izet);

  DAVecGetArray(fda, user->lJCsi, &jcsi);
  DAVecGetArray(fda, user->lJEta, &jeta);
  DAVecGetArray(fda, user->lJZet, &jzet);

  DAVecGetArray(fda, user->lKCsi, &kcsi);
  DAVecGetArray(fda, user->lKEta, &keta);
  DAVecGetArray(fda, user->lKZet, &kzet);

  DAVecGetArray(da, user->lAj, &aj);
  DAVecGetArray(da, user->lIAj, &iaj);
  DAVecGetArray(da, user->lJAj, &jaj);
  DAVecGetArray(da, user->lKAj, &kaj);

  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(da, user->lPhi, &p);
  
  DAVecGetArray(fda, user->Ucont, &ucont);


/*   for (k=zs; k<ze-1; k++) { */
/*     PetscPrintf(PETSC_COMM_WORLD, "b4ucont%le %le %le\n", ucont[k][1][1].x, ucont[k][1][1].y, ucont[k][1][1].z); */
/*   } */

  PetscReal dpdc, dpde, dpdz;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (i<mx-2) {
	  dpdc = p[k][j][i+1] - p[k][j][i];

	  dpde = 0.;
	  dpdz = 0.;
	  if (j==my-2 || nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	    if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && j!=1) {
	      dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		      p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	    }
	  }
	  else if (j == 1 || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	    if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	      dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		      p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	    }
	  }
	  else {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	  }

	  if (k == mz-2 || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	    if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && k!=1) {
	      dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		      p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	    }
	  }
	  else if (k == 1 || nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	    if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	      dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		      p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	    }
	  }
	  else {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	  }

	  if (!(nvert[k][j][i] + nvert[k][j][i+1])) {
	    ucont[k][j][i].x -= 
	      (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		       icsi[k][j][i].y * icsi[k][j][i].y +
		       icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	       dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		       ieta[k][j][i].y * icsi[k][j][i].y +
		       ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] + 
	       dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		       izet[k][j][i].y * icsi[k][j][i].y +
		       izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i]) * user->dt * user->st * 2/3.;
	  }
	}

	if (j<my-2) {
	  dpdc = 0.;
	  dpdz = 0.;
	  if (i == mx-2 || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	    if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && i!=1) {
	      dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		      p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	    }
	  }
	  else if (i == 1 || nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	    if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	      dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		      p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	    }
	  }
	  else {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	  }

	  dpde = p[k][j+1][i] - p[k][j][i];

	  if (k == mz-2 || nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	    if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && k!=1) {
	      dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		      p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	    }
	  }
	  else if (k == 1 || nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	    if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	      dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		      p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	    }
	  }
	  else {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	  }

	  if (!(nvert[k][j][i] + nvert[k][j+1][i])) {
	    ucont[k][j][i].y -=
	      (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		       jcsi[k][j][i].y * jeta[k][j][i].y +
		       jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	       dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		       jeta[k][j][i].y * jeta[k][j][i].y +
		       jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	       dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		       jzet[k][j][i].y * jeta[k][j][i].y +
		       jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i]) * user->dt * user->st * 2/3.;
	  }
	}
	
	if (k < mz-2) {
	  dpdc = 0.;
	  dpde = 0.;
	  if (i == mx-2 || nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	    if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && i!=1) {
	      dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		      p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	    }
	  }
	  else if (i == 1 || nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	    if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	      dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		      p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	    }
	  }
	  else {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	  }
	  
	  if (j == my-2 || nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	    if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && j!=1) {
	      dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		      p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	    }
	  }
	  else if (j == 1 || nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	    if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	      dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		      p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	    }
	  }
	  else {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	  }

	  dpdz = p[k+1][j][i] - p[k][j][i];
	  if (!(nvert[k][j][i] + nvert[k+1][j][i])) {
	    
	    ucont[k][j][i].z -=
	      (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		       kcsi[k][j][i].y * kzet[k][j][i].y +
		       kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	       dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		       keta[k][j][i].y * kzet[k][j][i].y +
		       keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	       dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		       kzet[k][j][i].y * kzet[k][j][i].y +
		       kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i])*user->dt*user->st * 2/3.;
	  }
	}
      }
    }
  }


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

  DAVecRestoreArray(da, user->lAj, &aj);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);
  DAVecRestoreArray(da, user->lPhi, &p);
  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAVecRestoreArray(da, user->lNvert, &nvert);

/*   VecAssemblyBegin(user->Ucont); */
/*   VecAssemblyEnd(user->Ucont); */
  
  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);


  PetscBarrier(PETSC_NULL);
  InflowFlux(user);
  OutflowFlux(user);
  PetscPrintf(PETSC_COMM_WORLD, "Fluxin %le\n", user->FluxOutSum);

  FormBCS(user);
  
  return(0);
}

PetscErrorCode InterpolationFlux(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  Vec		Ucont = user->Ucont;
  Cmpnts	***ucont;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
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

  PetscReal ***nvert, ***aj;
  PetscReal lFluxIntp, FluxSumIntp;
  PetscReal lInsideBody, InsideBodySum;
  PetscReal lFluidVolume, FluidVolumeSum;
  DAVecGetArray(da, user->lNvert, &nvert);
  DAVecGetArray(fda, user->Ucont, &ucont);
  DAVecGetArray(da, user->lAj, &aj);
  lFluxIntp = 0.;
  lInsideBody = 0.;
  lFluidVolume = 0.;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  lFluidVolume += 1./aj[k][j][i];
	  if (nvert[k][j][i+1] > 0.1) {
	    lFluxIntp -= ucont[k][j][i].x;
	  }
	  if (nvert[k][j+1][i] > 0.1) {
	    lFluxIntp -= ucont[k][j][i].y;
	  }
	  if (nvert[k+1][j][i] > 0.1) {
	    lFluxIntp -= ucont[k][j][i].z;
	  }
	}
	else {
	  lInsideBody += 1.;
	  if (nvert[k][j][i] < 1.1) {
	    if (nvert[k][j][i+1] < 0.1) {
	      lFluxIntp += ucont[k][j][i].x;
	    }
	    if (nvert[k][j+1][i] < 0.1) {
	      lFluxIntp += ucont[k][j][i].y;
	    }
	    if (nvert[k+1][j][i] < 0.1) {
	      lFluxIntp += ucont[k][j][i].z;
	    }
	  }
	}
      }
    }
  }

  PetscGlobalSum(&lFluxIntp, &user->FluxIntpSum, PETSC_COMM_WORLD);
  //  PetscGlobalSum(&lInsideBody, &InsideBodySum, PETSC_COMM_WORLD);
  PetscGlobalSum(&lFluidVolume, &FluidVolumeSum, PETSC_COMM_WORLD);

  user->FluxIntpSum /= FluidVolumeSum;

  PetscPrintf(PETSC_COMM_WORLD, "FluxIntp %e %e\n", user->FluxIntpSum, FluidVolumeSum);
  DAVecRestoreArray(da, user->lAj, &aj);  
  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAVecRestoreArray(da, user->lNvert, &nvert);
  return 0;
}
