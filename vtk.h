
PetscErrorCode DAView_3DVTK_StructuredGrid_appended_old(DM da,DM fda,  Vec FIELD,const char file_prefix[])
{
  char           vtk_filename[PETSC_MAX_PATH_LEN];
  PetscMPIInt    rank;
  MPI_Comm       comm;
  FILE           *vtk_fp = NULL;
  PetscInt       si,sj,sk,nx,ny,nz,i;
  PetscInt       f,n_fields,N;
  DM             cda;
  Vec            coords;
  Vec            l_FIELD;
  PetscScalar    *_L_FIELD;
  PetscInt       memory_offset;
  PetscScalar    *buffer;
  PetscLogDouble t0,t1;
  Vec Coor;
  Cmpnts ***coor;

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Writting VTK \n");
  PetscPrintf(PETSC_COMM_WORLD,"!!!!!!!!!!!!!!!!!!!!!!!! \n");

  PetscTime(&t0);

    DMDALocalInfo info;

    DMDAGetLocalInfo(da, &info);
    PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscInt xs = info.xs, xe = xs + info.xm;
    PetscInt ys = info.ys, ye = ys + info.ym;
    PetscInt zs = info.zs, ze = zs + info.zm;


    PetscInt lxs = xs, lxe = xe;
    PetscInt lys = ys, lye = ye;
    PetscInt lzs = zs, lze = ze;

  /* create file name */
  PetscObjectGetComm((PetscObject)da,&comm);
  MPI_Comm_rank(comm,&rank);
  PetscSNPrintf(vtk_filename,sizeof(vtk_filename),"subdomain-%s-p%1.4d.vts",file_prefix,rank);

  /* open file and write header */
  vtk_fp = fopen(vtk_filename,"w");
  if (!vtk_fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SYS,"Cannot open file = %s \n",vtk_filename);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<?xml version=\"1.0\"?>\n");

  /* coords */
//  N  = (ze-1-zs)*(ye-1-ys)*(xe-1-xs);
  N  = (xe-2-xs)*(ye-2-ys)*(ze-2-zs);


  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

  PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n", xs, xe-2, ys, ye-2, zs, ze-2);
  PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "    <Piece Extent=\"%D %D %D %D %D %D\">\n", xs, xe-2, ys, ye-2, zs, ze-2);

  memory_offset = 0;

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <CellData></CellData>\n");

  
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <Points>\n");
  /* copy coordinates */

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",memory_offset);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      </Points>\n");

  memory_offset = memory_offset + sizeof(PetscInt) + sizeof(PetscScalar)*N*3;

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <PointData> \n ");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"        <DataArray type=\"Float64\" Name=\"Field\" format=\"appended\" offset=\"%d\"/>\n", memory_offset);
//memory_offset = memory_offset + sizeof(PetscInt) + sizeof(PetscScalar)*N*3;
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      </PointData>\n");

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </Piece>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  </StructuredGrid>\n");


  PetscMalloc(sizeof(PetscScalar)*N,&buffer);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <AppendedData encoding=\"raw\">\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"_");

  /* write coordinates */

//  DMGetCoordinateDM(da,&cda);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Getting Coordinates \n");


  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Writting Coordinates \n");

  int M  = (ze-1-zs)*(ye-1-ys)*(xe-1-xs);
  {
    int  length = sizeof(PetscScalar)*M*3;
    fwrite(&length,sizeof(int),1,vtk_fp);
    double value[3];
  MPI_Barrier(PETSC_COMM_WORLD);
      printf("xs=%d, xe=%d\n", xs, xe);
  MPI_Barrier(PETSC_COMM_WORLD);
      printf("ys=%d, ye=%d\n", ys, ye);
  MPI_Barrier(PETSC_COMM_WORLD);
      printf("zs=%d, ze=%d\n", zs, ze);
  MPI_Barrier(PETSC_COMM_WORLD);

    for (int k=zs; k<ze-1; k++) {
      for (int j=ys; j<ye-1; j++) {
        for (int i=xs; i<xe-1; i++) {

          value[0] = coor[k][j][i].x;
          value[1] = coor[k][j][i].y;
          value[2] = coor[k][j][i].z;

          fwrite(value,sizeof(PetscScalar),3,vtk_fp);

        }
      }
    }

  }
  DMDAVecRestoreArray(fda, Coor, &coor);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"DONE Writting Coordinates \n");


//  /* write fields */
//
  DMGetLocalVector(da,&l_FIELD);
  DMGlobalToLocalBegin(da, FIELD,INSERT_VALUES,l_FIELD);
  DMGlobalToLocalEnd(da,FIELD,INSERT_VALUES,l_FIELD);
  VecGetArray(l_FIELD,&_L_FIELD);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"Writting Data \n");

  {
    int length = sizeof(PetscScalar)*N;
    fwrite(&length,sizeof(int),1,vtk_fp);
    /* load */
    for (i=0; i<N; i++) buffer[i] = _L_FIELD[i];
    /* write */
    fwrite(buffer,sizeof(PetscScalar),N,vtk_fp);
  }

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"\n  </AppendedData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"</VTKFile>\n");

  PetscFree(buffer);

  if (vtk_fp) {
    fclose(vtk_fp);
    vtk_fp = NULL;
  }

  PetscTime(&t1);

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"DONE With VTK \n");
  return(0);
}


PetscErrorCode DAView_3DVTK_StructuredGrid_appended(DM da,Vec FIELD,const char file_prefix[])
{
  char           vtk_filename[PETSC_MAX_PATH_LEN];
  PetscMPIInt    rank;
  MPI_Comm       comm;
  FILE           *vtk_fp = NULL;
  PetscInt       si,sj,sk,nx,ny,nz,i;
  PetscInt       f,n_fields,N;
  DM             cda;
  Vec            coords;
  Vec            l_FIELD;
  PetscScalar    *_L_FIELD;
  PetscInt       memory_offset;
  PetscScalar    *buffer;
  PetscLogDouble t0,t1;

  PetscTime(&t0);

  /* create file name */
  PetscObjectGetComm((PetscObject)da,&comm);
  MPI_Comm_rank(comm,&rank);
  PetscSNPrintf(vtk_filename,sizeof(vtk_filename),"subdomain-%s-p%1.4d.vts",file_prefix,rank);

  /* open file and write header */
  vtk_fp = fopen(vtk_filename,"w");
  if (!vtk_fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SYS,"Cannot open file = %s \n",vtk_filename);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<?xml version=\"1.0\"?>\n");

  /* coords */
  DMDAGetGhostCorners(da,&si,&sj,&sk,&nx,&ny,&nz);
  N    = nx * ny * nz;

#if defined(PETSC_WORDS_BIGENDIAN)
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",si,si+nx-1,sj,sj+ny-1,sk,sk+nz-1);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",si,si+nx-1,sj,sj+ny-1,sk,sk+nz-1);

  memory_offset = 0;

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <CellData></CellData>\n");

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <Points>\n");

  /* copy coordinates */
  DMGetCoordinateDM(da,&cda);
  DMGetCoordinatesLocal(da,&coords);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",memory_offset);
  memory_offset = memory_offset + sizeof(PetscInt) + sizeof(PetscScalar)*N*3;

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      </Points>\n");

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <PointData Scalars=\" ");
  DMDAGetInfo(da,0,0,0,0,0,0,0,&n_fields,0,0,0,0,0);
  for (f=0; f<n_fields; f++) {
    const char *field_name;
    DMDAGetFieldName(da,f,&field_name);
    PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"%s ",field_name);
  }
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"\">\n");

  for (f=0; f<n_fields; f++) {
    const char *field_name;

    DMDAGetFieldName(da,f,&field_name);
    PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"        <DataArray type=\"Float64\" Name=Data format=\"appended\" offset=\"%d\"/>\n", memory_offset);
    memory_offset = memory_offset + sizeof(PetscInt) + sizeof(PetscScalar)*N;
  }

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      </PointData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </Piece>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  </StructuredGrid>\n");

  PetscMalloc(sizeof(PetscScalar)*N,&buffer);
  DMGetLocalVector(da,&l_FIELD);
  DMGlobalToLocalBegin(da, FIELD,INSERT_VALUES,l_FIELD);
  DMGlobalToLocalEnd(da,FIELD,INSERT_VALUES,l_FIELD);
  VecGetArray(l_FIELD,&_L_FIELD);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <AppendedData encoding=\"raw\">\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"_");

  /* write coordinates */
  {
    int         length = sizeof(PetscScalar)*N*3;
    PetscScalar *allcoords;

    fwrite(&length,sizeof(int),1,vtk_fp);
    VecGetArray(coords,&allcoords);
    fwrite(allcoords,sizeof(PetscScalar),3*N,vtk_fp);
    VecRestoreArray(coords,&allcoords);
  }
  /* write fields */
  for (f=0; f<n_fields; f++) {
    int length = sizeof(PetscScalar)*N;
    fwrite(&length,sizeof(int),1,vtk_fp);
    /* load */
    for (i=0; i<N; i++) buffer[i] = _L_FIELD[n_fields*i + f];

    /* write */
    fwrite(buffer,sizeof(PetscScalar),N,vtk_fp);
  }
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"\n  </AppendedData>\n");

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"</VTKFile>\n");

  PetscFree(buffer);
  VecRestoreArray(l_FIELD,&_L_FIELD);
  DMRestoreLocalVector(da,&l_FIELD);

  if (vtk_fp) {
    fclose(vtk_fp);
    vtk_fp = NULL;
  }

  PetscTime(&t1);
  return(0);
}
