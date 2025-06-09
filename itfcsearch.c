static char help[] = "Interface Searching!";

#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct {
  DA	da, fda;
  PetscInt IM, JM, KM;
  Vec	Coor, Cent, BCS;
  /* bctype is used to specify boundary conditions
     if bctype == 0 then the whole side is interface points */
  PetscInt bctype[6];
  PetscInt itfcptsnumber;
  PetscInt *itfcI, *itfcJ, *itfcK;
  PetscInt *itfchostI, *itfchostJ, *itfchostK, *itfchostB;
  PetscReal *itfchostx, *itfchosty, *itfchostz;

  
} UserCtx;

typedef struct {
  PetscInt *is, *ie, *js, *je, *ks, *ke;
} SearchRange;

struct list_node {
  PetscInt i, j, k, bi;
  struct list_node *next;
};

typedef struct list_node llnode;

llnode *list_add(llnode **p, int i, int j, int k, int bi) {
  llnode *n = malloc(sizeof(llnode));
  n->next = *p;
  *p = n;
  n->i = i;
  n->j = j;
  n->k = k;
  n->bi = bi;
  return n;  
}

PetscErrorCode distance(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
{
  PetscReal xn1, yn1, zn1;
  PetscReal xc, yc, zc;
  
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
  if (PetscAbsReal(*d)<1.e-6) *d=0.;
  return (0);
}

PetscTruth ISInsideCell(Cmpnts p, Cmpnts cell[8], PetscReal d[6])
{
  // k direction
  distance(cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  if (d[4]<0) return(PETSC_FALSE);
  distance(cell[4], cell[7], cell[6], cell[5], p, &(d[5]));
  if (d[5]<0) return(PETSC_FALSE);

  // j direction
  distance(cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  if (d[2]<0) return(PETSC_FALSE);

  distance(cell[3], cell[2], cell[6], cell[7], p, &(d[3]));
  if (d[3]<0) return(PETSC_FALSE);

  // i direction
  distance(cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  if (d[0]<0) return(PETSC_FALSE);
  
  distance(cell[1], cell[5], cell[6], cell[2], p, &(d[1]));
  if (d[1]<0) return(PETSC_FALSE);
  return(PETSC_TRUE);
}

int main(int argc, char**argv)
{
  PetscInt block_number;
  UserCtx *user;

  PetscInt bi, i, j, k;
  PetscReal	***bcs;

  PetscReal  Max_X,Max_Y,Max_Z;
  PetscReal  Min_X,Min_Y,Min_Z;
  PetscReal  x,y,z;

  PetscInt dIM=1, dJM=30, dKM=30;
  PetscInt dI, dJ, dK;
  PetscReal xmin=1.e10, xmax=-1.e10, ymin=1.e10, ymax=-1.e10, 
    zmin=1.e10, zmax=-1.e10, tmp;
  PetscReal ddx, ddy, ddz;

  llnode *ptsincell[dKM][dJM][dIM], *head[dKM][dJM][dIM], *current;

  Cmpnts ***coor, ***cent;

  PetscInitialize(&argc, &argv, (char *)0, help);

  for (k=0; k<dKM; k++) {
    for (j=0; j<dJM; j++) {
      for (i=0; i<dIM; i++) {
	head[k][j][i] = NULL;
      }
    }
  }

  // Read in Grid data
  FILE *fd;
  fd = fopen("grid.dat", "r");
  fscanf(fd, "%i\n", &block_number);
  PetscPrintf(PETSC_COMM_WORLD, "%i\n", block_number);
  PetscMalloc(block_number*sizeof(UserCtx), &user);
  
  //for (bi=block_number-1; bi>-1; bi--) {
  for (bi=0; bi<block_number; bi++) {
    fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
  }
  for (bi=0; bi<block_number; bi++) {
	      DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1, 1, PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &(user[bi].da));
    DASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DAGetCoordinateDA(user[bi].da, &(user[bi].fda));
    DACreateGlobalVector(user[bi].fda, &(user[bi].Coor));
    DAVecGetArray(user[bi].fda, user[bi].Coor, &coor);
	  
    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i\n", user[bi].IM, user[bi].JM, user[bi].KM);
    for(k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
        for(i=0; i<user[bi].IM; i++) {
          fscanf(fd, "%le", &(coor[k][j][i].x));
        }
      }
    }
  
    for(k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
        for(i=0; i<user[bi].IM; i++) {
          fscanf(fd, "%le", &(coor[k][j][i].y));
        }
      }
    }
  
    for(k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
        for(i=0; i<user[bi].IM; i++) {
          fscanf(fd, "%le", &(coor[k][j][i].z));
        }
      }
    }
    
    DAVecRestoreArray(user[bi].fda, user[bi].Coor, &coor);
    
    VecStrideMin(user[bi].Coor, 0, PETSC_NULL, &tmp);
    xmin = PetscMin(xmin, tmp);
    VecStrideMax(user[bi].Coor, 0, PETSC_NULL, &tmp);
    xmax = PetscMax(xmax, tmp);

    VecStrideMin(user[bi].Coor, 1, PETSC_NULL, &tmp);
    ymin = PetscMin(ymin, tmp);
    VecStrideMax(user[bi].Coor, 1, PETSC_NULL, &tmp);
    ymax = PetscMax(ymax, tmp);

    VecStrideMin(user[bi].Coor, 2, PETSC_NULL, &tmp);
    zmin = PetscMin(zmin, tmp);
    VecStrideMax(user[bi].Coor, 2, PETSC_NULL, &tmp);
    zmax = PetscMax(zmax, tmp);
  }

  xmin -= 0.01; xmax += 0.01;
  ymin -= 0.01; ymax += 0.01;
  zmin -= 0.01; zmax += 0.01;

  fclose(fd);

  fd = fopen("bcs.dat", "r");
  if(!fd) printf("cannot open bcs.dat !\n"),exit(0);
  for (bi=0; bi<block_number; bi++) {
    fscanf(fd, "%i %i %i %i %i %i\n", &(user[bi].bctype[0]),
	   &(user[bi].bctype[1]), &(user[bi].bctype[2]), &(user[bi].bctype[3]),
	   &(user[bi].bctype[4]), &(user[bi].bctype[5]));
    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i %i %i\n",bi,user[bi].bctype[1],user[bi].bctype[2],user[bi].bctype[3],user[bi].bctype[4],user[bi].bctype[5] );
  }
  fclose(fd);

  ddx = (xmax - xmin) / (dIM-1.);
  ddy = (ymax - ymin) / (dJM-1.);
  ddz = (zmax - zmin) / (dKM-1.);

  /* Calculate the locations of grid center and
     locate the center points into a dIM * dJM * dKM cell system
  */
  for (bi=0; bi<block_number; bi++) {
    VecDuplicate(user[bi].Coor, &(user[bi].Cent));
    DAVecGetArray(user[bi].fda, user[bi].Coor, &coor);
    DAVecGetArray(user[bi].fda, user[bi].Cent, &cent);

    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
	for (i=0; i<user[bi].IM; i++) {
	  dI = (int)((coor[k][j][i].x - xmin) / ddx);
	  dJ = (int)((coor[k][j][i].y - ymin) / ddy);
	  dK = (int)((coor[k][j][i].z - zmin) / ddz);
	  head[dK][dJ][dI] = list_add(&(head[dK][dJ][dI]), i, j, k, bi);
	}
      }
    }
    for (k=0; k<user[bi].KM-1; k++) {
      for (j=0; j<user[bi].JM-1; j++) {
	for (i=0; i<user[bi].IM-1; i++) {
	  cent[k][j][i].x = 0.125 * (coor[k  ][j  ][i  ].x +
				     coor[k  ][j  ][i+1].x +
				     coor[k  ][j+1][i  ].x +
				     coor[k  ][j+1][i+1].x +
				     coor[k+1][j  ][i  ].x +
				     coor[k+1][j  ][i+1].x +
				     coor[k+1][j+1][i  ].x +
				     coor[k+1][j+1][i+1].x);

	  cent[k][j][i].y = 0.125 * (coor[k  ][j  ][i  ].y +
				     coor[k  ][j  ][i+1].y +
				     coor[k  ][j+1][i  ].y +
				     coor[k  ][j+1][i+1].y +
				     coor[k+1][j  ][i  ].y +
				     coor[k+1][j  ][i+1].y +
				     coor[k+1][j+1][i  ].y +
				     coor[k+1][j+1][i+1].y);

	  cent[k][j][i].z = 0.125 * (coor[k  ][j  ][i  ].z +
				     coor[k  ][j  ][i+1].z +
				     coor[k  ][j+1][i  ].z +
				     coor[k  ][j+1][i+1].z +
				     coor[k+1][j  ][i  ].z +
				     coor[k+1][j  ][i+1].z +
				     coor[k+1][j+1][i  ].z +
				     coor[k+1][j+1][i+1].z);

	}
      }
    }

    if (bi==1) {
      i=user[bi].IM-2;
      j=user[bi].JM-2;
      k=user[bi].KM-2;
      Max_X=10;//cent[k  ][j  ][2  ].x;
      Max_Y=cent[k  ][j-10][i  ].y;
      Max_Z=cent[k-10][j  ][i  ].z;
     
      Min_X=-10;//cent[k  ][j  ][  1].x;
      Min_Y=cent[k  ][ 10][i  ].y;
      Min_Z=cent[ 10][j  ][i  ].z;

      PetscPrintf(PETSC_COMM_WORLD, "Max x y z %d %d %d\n",i,j,k);
      PetscPrintf(PETSC_COMM_WORLD, "Max x y z %le %le %le\n",Max_X,Max_Y,Max_Z);
      PetscPrintf(PETSC_COMM_WORLD, "Min x y z %le %le %le\n",Min_X,Min_Y,Min_Z);
  
    }	
    DAVecRestoreArray(user[bi].fda, user[bi].Cent, &cent);
  }


  /* Decide wheter a point is interface point */
  for (bi = 0; bi<block_number; bi++) {
    DACreateGlobalVector(user[bi].da, &(user[bi].BCS));
    VecSet(user[bi].BCS, 1.);
    DAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    DAVecGetArray(user[bi].fda, user[bi].Cent, &cent);

    if (user[bi].bctype[0] == 0) { // imin boundary
      for (k=0; k<user[bi].KM; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  for (i=1; i<2; i++) {
	    bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    if (user[bi].bctype[1] == 0) { // imax boundary
      for (k=0; k<user[bi].KM; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  for (i=user[bi].IM-2; i<user[bi].IM; i++) {
	    bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    if (user[bi].bctype[2] == 0) { // jmin boundary
      for (k=0; k<user[bi].KM; k++) {
	for (j=0; j<2; j++) {
	  for (i=1; i<user[bi].IM; i++) {
	    bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0) { // jmax boundary
      for (k=0; k<user[bi].KM; k++) {
	for (j=user[bi].JM-2; j<user[bi].JM; j++) {
	  for (i=1; i<user[bi].IM; i++) {
	    bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    if (user[bi].bctype[4] == 0) { // kmin boundary
      for (k=0; k<2; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  for (i=1; i<user[bi].IM; i++) {
	    bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    if (user[bi].bctype[5] == 0) { // kmax boundary
      for (k=user[bi].KM-2; k<user[bi].KM; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  for (i=1; i<user[bi].IM; i++) {
	    bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    
    // inside another domain
    if (bi==0) {
      for (k=1; k<user[bi].KM-1; k++) {
	for (j=1; j<user[bi].JM-1; j++) {
	  for (i=1; i<user[bi].IM-1; i++) {
	    x=cent[k][j][i].x;
	    y=cent[k][j][i].y;
	    z=cent[k][j][i].z;
	    
	    if (x<=Max_X && x>=Min_X &&
		y<=Max_Y && y>=Min_Y &&
		z<=Max_Z && z>=Min_Z)
	      bcs[k][j][i] = 0.;
	  }
	}
      }
    }

    DAVecRestoreArray(user[bi].fda, user[bi].Cent, &cent);
    DAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
  }

  Cmpnts ***host[10];
  for (bi=0; bi<block_number; bi++) {
    DAVecGetArray(user[bi].fda, user[bi].Coor, &(host[bi]));
  }

  /* For each cell in the dividing cell system, find out the range
     of center grid nodes from each block */
  PetscReal xc, yc, zc, xh[8], yh[8], zh[8];
  PetscInt bh, sb, si, sj, sk;
  Cmpnts cell[8];
/*   SearchRange sr[dKM][dJM][dIM]; */
  PetscInt ids[dKM][dJM][dIM][10], ide[dKM][dJM][dIM][10],
    jds[dKM][dJM][dIM][10], jde[dKM][dJM][dIM][10],
    kds[dKM][dJM][dIM][10], kde[dKM][dJM][dIM][10];
  for (k=0; k<dKM; k++) {
    for (j=0; j<dJM; j++) {
      for (i=0; i<dIM; i++) {
/*         PetscMalloc(block_number*sizeof(PetscInt), &(sr[k][j][i].is)); */
/*         PetscMalloc(block_number*sizeof(PetscInt), &(sr[k][j][i].ie)); */
/*         PetscMalloc(block_number*sizeof(PetscInt), &(sr[k][j][i].js)); */
/*         PetscMalloc(block_number*sizeof(PetscInt), &(sr[k][j][i].je)); */
/*         PetscMalloc(block_number*sizeof(PetscInt), &(sr[k][j][i].ks)); */
/*         PetscMalloc(block_number*sizeof(PetscInt), &(sr[k][j][i].ke)); */

        PetscInt tt;
	for (tt=0; tt<block_number; tt++) {
	  ids[k][j][i][tt] = 1000;
	  ide[k][j][i][tt] = 0;

	  jds[k][j][i][tt] = 1000;
	  jde[k][j][i][tt] = 0;

	  kds[k][j][i][tt] = 1000;
	  kde[k][j][i][tt] = 0;

/* 	  sr[k][j][i].ie[tt] = 0; */

/* 	  sr[k][j][i].js[tt] = 1000; */
/* 	  sr[k][j][i].je[tt] = 0; */

/* 	  sr[k][j][i].ks[tt] = 1000; */
/* 	  sr[k][j][i].ke[tt] = 0; */
	}
	current = head[k][j][i];
	while (current) {
	  bh = current->bi;
//	  if (bh!=1) PetscPrintf(PETSC_COMM_WORLD, "%i %i %i DDDD\n", i,j ,k);
	  tt = PetscMax(current->i-1,1);
	  
	  tt = PetscMin(tt, ids[k][j][i][bh]);
	  ids[k][j][i][current->bi] = tt;

	  tt = PetscMin(current->i+1,user[bh].IM-1);
	  
	  tt = PetscMax(tt, ide[k][j][i][bh]);
	  ide[k][j][i][bh] = tt;


	  tt = PetscMax(current->j-1,0);
	  
	  tt = PetscMin(tt, jds[k][j][i][bh]);
	  jds[k][j][i][bh] = tt;

	  tt = PetscMin(current->j+1,user[bh].JM-1);
	  
	  tt = PetscMax(tt, jde[k][j][i][bh]);
	  jde[k][j][i][bh] = tt;

	  tt = PetscMax(current->k-1,0);
	  
	  tt = PetscMin(tt, kds[k][j][i][bh]);
	  kds[k][j][i][bh] = tt;

	  tt = PetscMin(current->k+1,user[bh].KM-1);
	  
	  tt = PetscMax(tt, kde[k][j][i][bh]);
	  kde[k][j][i][bh] = tt;


/* 	  sr[k][j][i].ie[bh] = PetscMax(PetscMin(current->i, user[bh].IM-2), */
/* 					sr[k][j][i].ie[bh]); */
/* 	  sr[k][j][i].js[bh] = PetscMin(PetscMax(current->j-1,0), */
/* 					sr[k][j][i].js[bh]); */
/* 	  sr[k][j][i].je[bh] = PetscMax(PetscMin(current->j, user[bh].JM-2), */
/* 					sr[k][j][i].je[bh]); */
/* 	  sr[k][j][i].ks[bh] = PetscMin(PetscMax(current->k-1,0), */
/* 					sr[k][j][i].ks[bh]); */
/* 	  sr[k][j][i].ke[bh] = PetscMax(PetscMin(current->k, user[bh].KM-2), */
/* 					sr[k][j][i].ke[bh]); */
	  current = current->next;
	}
	//       PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i %i %i\n", i,j,k,bh,  sr[k][j][i].ks[0], sr[k][j][i].ke[0]);

      }
    }
  }

  for (bi=0; bi<block_number; bi++) {
    DAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    user[bi].itfcptsnumber=0;
    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
	for (i=0; i<user[bi].IM; i++) {
	  if(bcs[k][j][i]<1.e-6) user[bi].itfcptsnumber++;
	}
      }
    }
    DAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
  }

  PetscReal d[6];
  Cmpnts pc;
  PetscTruth found;
  
  for (bi=0; bi<block_number; bi++) {
    DAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    user[bi].itfcptsnumber=0;
    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
        for (i=0; i<user[bi].IM; i++) {
          if (bcs[k][j][i]<1.e-6) {
            user[bi].itfcptsnumber++;
          }
        }
      }
    }
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfcI);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfcJ);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfcK);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostI);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostJ);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostK);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostB);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscReal), &user[bi].itfchostx);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscReal), &user[bi].itfchosty);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscReal), &user[bi].itfchostz);
    DAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
  }
  
  PetscInt number, is, ie, js, je, ks, ke;
  for (bi=0; bi<block_number; bi++) {
    DAVecGetArray(user[bi].fda, user[bi].Coor, &cent);
    DAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    number = -1;
    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i\n", user[bi].IM, user[bi].JM, user[bi].KM, bi);
    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
	for (i=0; i<user[bi].IM; i++) {
	  if(bcs[k][j][i]<1.e-6) {
	    /* 	    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i\n", i, j, k, bi); */
	    number++;
	    pc.x = cent[k][j][i].x;
	    pc.y = cent[k][j][i].y;
	    pc.z = cent[k][j][i].z;
	    dI = (int)((pc.x - xmin) / ddx);
	    dJ = (int)((pc.y - ymin) / ddy);
	    dK = (int)((pc.z - zmin) / ddz);
	    found = PETSC_FALSE;
	    for (sb=0; sb<block_number; sb++) {
	      if (sb!=bi) {
/* 		ks = sr[dK][dJ][dI].ks[sb]; */
/* 		ke = sr[dK][dJ][dI].ke[sb]; */
/* 		PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i %i %i\n", dK, dJ, dI, sb,  ks, ke); */
         
		for (sk=kds[dK][dJ][dI][sb]; sk<kde[dK][dJ][dI][sb]; sk++) {
		  for (sj=jds[dK][dJ][dI][sb]; sj<jde[dK][dJ][dI][sb]; sj++) {
		    for (si=ids[dK][dJ][dI][sb]; si<ide[dK][dJ][dI][sb]; si++) {
		      //PetscPrintf(PETSC_COMM_WORLD, "sk %i %i %i\n", sk, sj, si);
		      cell[0] = host[sb][sk  ][sj  ][si  ];
		      cell[1] = host[sb][sk  ][sj  ][si+1];
		      cell[2] = host[sb][sk  ][sj+1][si+1];
		      cell[3] = host[sb][sk  ][sj+1][si  ];

		      cell[4] = host[sb][sk+1][sj  ][si  ];
		      cell[5] = host[sb][sk+1][sj  ][si+1];
		      cell[6] = host[sb][sk+1][sj+1][si+1];
		      cell[7] = host[sb][sk+1][sj+1][si  ];

		      if(ISInsideCell(pc, cell, d)) {
			found = PETSC_TRUE;
			user[bi].itfcI[number] = i;
			user[bi].itfcJ[number] = j;
			user[bi].itfcK[number] = k;
			user[bi].itfchostI[number] = si;
			user[bi].itfchostJ[number] = sj;
			user[bi].itfchostK[number] = sk;
			user[bi].itfchostB[number] = sb;
			user[bi].itfchostx[number] = d[0] / (d[0] + d[1]);
			user[bi].itfchosty[number] = d[2] / (d[2] + d[3]);
			user[bi].itfchostz[number] = d[4] / (d[4] + d[5]);

			goto nextp;
		      }
		    }
		  }
		}
	      }
	    }
	    nextp: number=number;
	    if (!found) {
	      PetscPrintf(PETSC_COMM_WORLD, "Can't find the host cell at %i,%i,%i,%i!\n", i, j, k, bi);
	      //PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i\n", user[bi].IM, user[bi].JM, user[bi].KM, bi);
	    }
	    //PetscPrintf(PETSC_COMM_WORLD, "%i\n", number);
	  }
        }
      }
    }
    DAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
    DAVecRestoreArray(user[bi].fda, user[bi].Cent, &cent);
  }
  
  fd = fopen("interface.dat", "w");
  for (bi=0; bi<block_number; bi++) {
    PetscFPrintf(PETSC_COMM_WORLD, fd, "%i\n", user[bi].itfcptsnumber);
    PetscPrintf(PETSC_COMM_WORLD, "bi=%i  ifcp#=%i\n", bi,user[bi].itfcptsnumber);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %i\n", user[bi].itfcI[i], user[bi].itfcJ[i], user[bi].itfcK[i]);
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %i %i\n", user[bi].itfchostI[i], user[bi].itfchostJ[i], user[bi].itfchostK[i], user[bi].itfchostB[i]);
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le %le\n", user[bi].itfchostx[i], user[bi].itfchosty[i], user[bi].itfchostz[i]);
    }
  }
  fclose(fd);
  PetscFinalize();
  return(0);
}

