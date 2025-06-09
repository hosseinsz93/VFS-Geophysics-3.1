//canopy.h
//#ifndef CANOPY_H // include guard
//#define CANOPY_H
#include "variables.h"
#include "petsc.h"
#include "petscvec.h"
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>

extern PetscInt cnpy;
extern PetscInt NumberOfCnpy;
extern PetscReal cnpyDispZ;
extern PetscReal cnpyHeightMlt;

// Define the structure of variables for each canopy model in domain
typedef struct CnpyNodes{ 
	//Canopy
	int		n_v, n_elmt;            	// Number of points and elements
	int		*nv1, *nv2, *nv3;       	// Element to point index
  PetscReal	*x_bp, *y_bp, *z_bp;    	// Coordinates of surface nodes
	PetscReal	*nf_x, *nf_y, *nf_z;    	// Normal direction
//	PetscReal	*nt_x, *nt_y, *nt_z;    	// Tangent dir
//	PetscReal	*ns_x, *ns_y, *ns_z;    	// Azimuthal dir
//	PetscReal	*dA;                    	// Area of an element
	PetscReal	*cent_x,*cent_y,*cent_z;	// Center of Element
	int   i_min, i_max;			// Max and Min indeces	
	int		j_min, j_max;			// Max and Min indeces 
	int		k_min, k_max;			// Max and Min indeces

  int naf;
  PetscReal *af;
  PetscReal *yaf;
  PetscReal maxHeight;
//  PetscReal       *cnpyVert;			// Nvert homologous variable for canopy
//	PetscReal	*CDy;				// Drag coefficient as function of height
//	Cmpnts          *Fcnpy;				// Body Force to Cell

} CnpyNodes;

PetscErrorCode CanopyInitialization(UserCtx *user, int NumberOfCnpy, PetscReal cnpyHeightDim, PetscReal cnpyCd);
PetscErrorCode CanopyForce(UserCtx *user);

PetscErrorCode canopy_read_ucd(CnpyNodes *cnpy, int icnpy);
PetscErrorCode canopy_read_af(CnpyNodes *cnpy, int icnpy, PetscReal cnpyHeightDim);
PetscErrorCode search_canopy_nodes(UserCtx *user, CnpyNodes *cnpy, int icnpy);
PetscReal linearinterpaf(PetscReal yh, CnpyNodes *cnpy);
PetscInt point_cell_cnpy(Cmpnts p, PetscInt ip, PetscInt jp, PetscInt kp,
           CnpyNodes *cnpy, PetscInt ncx, PetscInt ncy,
           PetscInt ncz, PetscReal dcx, PetscReal dcy,
           PetscReal xbp_min, PetscReal ybp_min,
           PetscReal zbp_max, LIST *cell_trg,
           PetscInt flg);

PetscErrorCode randomdirection2(Cmpnts p, PetscInt ip, PetscInt jp,
             PetscReal xbp_min, PetscReal ybp_min,
             PetscReal zbp_max, PetscReal dcx, PetscReal dcy,
             PetscReal dir[3],PetscInt seed);

//#endif
