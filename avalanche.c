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

         ibm->deltz_p_us[elmt] = 0.;
         ibm->deltz_p_ds[elmt] = 0.;
         ibm->A_us[elmt] = 0.;
         ibm->A_ds[elmt] = 0.;
         
         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt( ( nfy / ( nfz+1.e-8) ) * ( nfy / ( nfz+1.e-8) ) + ( nfx / ( nfz+1.e-8) ) * ( nfx / ( nfz+1.e-8) ) ));
         if(nfz  = 0. || ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if( ibm->Rigidity[elmt] == 0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{

       
PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_first_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			for( nbelmt = 0; nbelmt < n_elmt; nbelmt++ )
			{

                 	 if(nbelmt != elmt && ibm->nf_z[nbelmt] > 1.e-7 && ibm->elmt_depth[nbelmt] < 0.3 && ibm->Rigidity[nbelmt] ==0
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
				dr   = sqrt( dx * dx + dy * dy + dz * dz);

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
				dr   = sqrt( dx * dx + dy * dy + dz * dz );

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
    
// check and correct z-direction angle between element's centerpoint

	for( elmt = n_elmt; elmt >= 0; elmt-- )
	{
         nfx = ibm->nf_x[elmt];
         nfy = ibm->nf_y[elmt];
         nfz = ibm->nf_z[elmt];

         ibm->deltz_p_us[elmt] = 0.;
         ibm->deltz_p_ds[elmt] = 0.;
         ibm->A_us[elmt] = 0.;
         ibm->A_ds[elmt] = 0.;
         
         if(ibm->Rigidity[elmt]==0) ibm->max_bed_angle[elmt] = atan(sqrt( ( nfy / ( nfz+1.e-8) ) * ( nfy / ( nfz+1.e-8) ) + ( nfx / ( nfz+1.e-8) ) * ( nfx / ( nfz+1.e-8) ) ));
         if(nfz  = 0. || ibm->Rigidity[elmt]==1) ibm->max_bed_angle[elmt] = 0.;   
            

		if( ibm->Rigidity[elmt]==0 && fabs(ibm->max_bed_angle[elmt]) > angle_repose)
		{

       
        PetscPrintf(PETSC_COMM_WORLD, "avalanch model activated_second_sweep: elmt & Maxangle &  repose: %d %le %le\n",elmt, ibm->max_bed_angle[elmt]*180./PI, angle_repose*180./PI);
			for( nbelmt = n_elmt; nbelmt >= 0; nbelmt-- )
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
				dr   = sqrt( dx * dx + dy * dy + dz * dz );

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
                 
                         
			for( nbelmt = n_elmt; nbelmt >=0; nbelmt-- )
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
				dr   = sqrt( dx * dx + dy * dy + dz * dz );

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
