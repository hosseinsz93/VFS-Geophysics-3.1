#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string>
#include <unistd.h>
#include <time.h>
#include <vector>

//#include "petscvec.h"
//#include "petscda.h"
//#include "petscksp.h"
//#include "petscsnes.h"
// FOr Generating inflow from one grid to the the other grid 

#include "TECIO.h"

using namespace std;

double PI=3.14159265359;

int Nt_g1, Nt_g2;
double dt_g1, dt_g2;
int Nx_g1, Ny_g1, Nx_g2, Ny_g2;
int bc[4];


int save_inflow_period, read_inflow_period, inflow_recycle_perioid;
char path_g1[256]; 
char path_g2[256]; 

typedef struct {
	//PetscScalar x, y, z;	
        double x, y, z;
} Cmpnts;

int levelset=0;
int sediment=1;


double Scale_inflow=1.53; //3.75; // 1.0; 1.47 sedi

int main(int argc, char **argv) 
{
  
	sprintf(path_g1, "./inflow_0.000888");
	sprintf(path_g2, "./inflow_0.001");
	int t_g1,t_g2,j,i;
	
	Nt_g1=12000; //40000;
	dt_g1=0.000888;
	Nt_g2=10000;
	dt_g2=0.001;
	Nx_g1=662;
	Ny_g1=393;
	Nx_g2=232;
	Ny_g2=287;
	save_inflow_period=1000;
	inflow_recycle_perioid=12000;
	bc[0]=0;
	bc[1]=0;
	bc[2]=0;
	bc[3]=1;

	read_inflow_period=save_inflow_period;
	
	std::vector< std::vector<double> > x_g1 (Ny_g1);
	std::vector< std::vector<double> > y_g1 (Ny_g1);
	std::vector< std::vector<double> > z_g1 (Ny_g1);
	std::vector< std::vector<double> > level_g1 (Ny_g1);
	std::vector< std::vector<double> > nvert_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) x_g1[j].resize(Nx_g1);
	for( j=0; j<Ny_g1; j++) y_g1[j].resize(Nx_g1);
	for( j=0; j<Ny_g1; j++) z_g1[j].resize(Nx_g1);
	for( j=0; j<Ny_g1; j++) level_g1[j].resize(Nx_g1);
	for( j=0; j<Ny_g1; j++) nvert_g1[j].resize(Nx_g1);
	
	std::vector< std::vector<double> > x_g2 (Ny_g2);
	std::vector< std::vector<double> > y_g2 (Ny_g2);
	std::vector< std::vector<double> > z_g2 (Ny_g2);
	std::vector< std::vector<double> > level_g2 (Ny_g2);
	std::vector< std::vector<double> > nvert_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) x_g2[j].resize(Nx_g2);
	for( j=0; j<Ny_g2; j++) y_g2[j].resize(Nx_g2);
	for( j=0; j<Ny_g2; j++) z_g2[j].resize(Nx_g2);
	for( j=0; j<Ny_g2; j++) level_g2[j].resize(Nx_g2);
	for( j=0; j<Ny_g2; j++) nvert_g2[j].resize(Nx_g2);
	

	std::vector< std::vector<int> > ii_g2 (Ny_g2);
	std::vector< std::vector<int> > jj_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) ii_g2[j].resize(Nx_g2);
	for( j=0; j<Ny_g2; j++) jj_g2[j].resize(Nx_g2);
	

	std::vector< std::vector<Cmpnts> > ucat1_g1 (Ny_g1);
	std::vector< std::vector<Cmpnts> > ucat2_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) ucat1_g1[j].resize(Nx_g1);
	for( j=0; j<Ny_g1; j++) ucat2_g1[j].resize(Nx_g1);
	

	std::vector< std::vector<Cmpnts> > U_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) U_g1[j].resize(Nx_g1);
	
	std::vector< std::vector<double> > uu_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) uu_g1[j].resize(Nx_g1);
						
	std::vector< std::vector<double> > vv_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) vv_g1[j].resize(Nx_g1);

	std::vector< std::vector<double> > ww_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) ww_g1[j].resize(Nx_g1);

	std::vector< std::vector<double> > uv_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) uv_g1[j].resize(Nx_g1);

	std::vector< std::vector<double> > vw_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) vw_g1[j].resize(Nx_g1);

	std::vector< std::vector<double> > wu_g1 (Ny_g1);
	for( j=0; j<Ny_g1; j++) wu_g1[j].resize(Nx_g1);


	std::vector< std::vector<Cmpnts> > ucat_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) ucat_g2[j].resize(Nx_g2);
	

	// Read Grid g1
    	FILE *fd;
    	char filen[80];  
 	char string[128]; 
        sprintf(filen, "InflowPlane_g1.dat");
        fd = fopen(filen, "r");
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	if (levelset) fgets(string, 128, fd);	

	printf("read coor, levelset and immersed body information for grid g1\n");

	double tmp;
       	for (j=1; j<Ny_g1; j++)  
       	for (i=1; i<Nx_g1; i++) { 
        	if (levelset) fscanf(fd, "%le %le %le %le %le", &x_g1[j][i], &y_g1[j][i], &z_g1[j][i], &nvert_g1[j][i], &level_g1[j][i]);
        	if (!levelset) fscanf(fd, "%le %le %le %le", &x_g1[j][i], &y_g1[j][i], &z_g1[j][i], &nvert_g1[j][i]);
		
		y_g1[j][i]+=0.4;
		x_g1[j][i]+=1.375;

		/*
		if(sediment) {
			double tmpx=x_g1[j][i];
			double tmpy=y_g1[j][i];
			double tmpz=z_g1[j][i];

			x_g1[j][i]=tmpy;
			y_g1[j][i]=tmpz;
			z_g1[j][i]=tmpx;
			
		}
		*/

		z_g1[j][i]=0.0;
        }
        fclose(fd);
	
	// Read Grid g2
        sprintf(filen, "InflowPlane_g2.dat");
        fd = fopen(filen, "r");
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	fgets(string, 128, fd);	
	if (levelset) fgets(string, 128, fd);	

	printf("read coor, levelset and immersed body information for grid g2\n");

       	for (j=1; j<Ny_g2; j++)  
       	for (i=1; i<Nx_g2; i++) { 
        	if (levelset) fscanf(fd, "%le %le %le %le %le", &x_g2[j][i], &y_g2[j][i], &z_g2[j][i], &nvert_g2[j][i], &level_g2[j][i]);
        	if (!levelset) fscanf(fd, "%le %le %le %le", &x_g2[j][i], &y_g2[j][i], &z_g2[j][i], &nvert_g2[j][i]);

		if(sediment) {
			double tmpx=x_g2[j][i];
			double tmpy=y_g2[j][i];
			double tmpz=z_g2[j][i];

			x_g2[j][i]=tmpy;
			y_g2[j][i]=tmpz;
			z_g2[j][i]=tmpx;
			
		}

		z_g2[j][i]=0.0;
        }
        fclose(fd);

	printf("find the closest grid node of g1 to g2, Nx_g1=%d, Ny_g1=%d\n", Nx_g1, Ny_g1);
	printf("find the closest grid node of g1 to g2, Nx_g2=%d, Ny_g2=%d\n", Nx_g2, Ny_g2);
	int i1, j1;
	// find the closest grid node of g1 to g2
	for (j=1; j<Ny_g2; j++)  
       	for (i=1; i<Nx_g2; i++) { 
		double xg2=x_g2[j][i];
		double yg2=y_g2[j][i];
		double zg2=z_g2[j][i];

		double dist_min=10e6;
		double dist;

		ii_g2[j][i]=0; 
		jj_g2[j][i]=0; 
	

		for (j1=1; j1<Ny_g1; j1++)  
       		for (i1=1; i1<Nx_g1; i1++) { 
			double xg1=x_g1[j1][i1];
			double yg1=y_g1[j1][i1];
			double zg1=z_g1[j1][i1];
			dist=sqrt(pow((xg1-xg2),2)+pow((yg1-yg2),2)+pow((zg1-zg2),2));
			if (dist<dist_min) {
				dist_min=dist;
				jj_g2[j][i]=j1; 
				ii_g2[j][i]=i1; 
			}
		}

	}		

	      		printf("output tecplot file\n");
    			FILE *f;
	    		char filen_tec[80];
                	sprintf(filen_tec, "%s/grid_g1.plt", path_g2);
//	    		sprintf(filen_tec, "InFlowField%06d.dat",t_LES);
	    		f = fopen(filen_tec, "w");
	    		fprintf(f, "Variables=x,y,z\n");

			fprintf(f, "ZONE T='QUADRILATERAL', N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([1-6]=NODAL)\n", (Nx_g1-1)*(Ny_g1-1), (Nx_g1-2)*(Ny_g1-2));

	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
      				fprintf(f, "%le\n", x_g1[j][i]);
			}

	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
	                        fprintf(f, "%le\n", y_g1[j][i]);
	                }

        	        for (j=1;j<Ny_g1;j++)
                	for (i=1;i<Nx_g1;i++) {
	                       	fprintf(f, "%le\n", z_g1[j][i]);
	                }


                	for (j=1;j<Ny_g1-1;j++)
	                for (i=1;i<Nx_g1-1;i++) {
			     	fprintf(f, "%d %d %d %d\n", (j-1)*(Nx_g1-1)+i, (j-1)*(Nx_g1-1)+i+1, (j)*(Nx_g1-1)+i+1, (j)*(Nx_g1-1)+i);
			}
    
		    	fclose(f);
	




	      		printf("output tecplot file\n");
                	sprintf(filen_tec, "%s/iijj_dt=%g.plt", path_g2, dt_g2);
//	    		sprintf(filen_tec, "InFlowField%06d.dat",t_LES);
	    		f = fopen(filen_tec, "w");
	    		fprintf(f, "Variables=x,y,z,ii,jj\n");

			fprintf(f, "ZONE T='QUADRILATERAL', N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([1-6]=NODAL)\n", (Nx_g2-1)*(Ny_g2-1), (Nx_g2-2)*(Ny_g2-2));

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
      				fprintf(f, "%le\n", x_g2[j][i]);
			}

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
	                        fprintf(f, "%le\n", y_g2[j][i]);
	                }

        	        for (j=1;j<Ny_g2;j++)
                	for (i=1;i<Nx_g2;i++) {
	                       	fprintf(f, "%le\n", z_g2[j][i]);
	                }


	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
        	                fprintf(f, "%d\n", ii_g2[j][i]);
	                }

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
        	                fprintf(f, "%d\n", jj_g2[j][i]);
	                }


                	for (j=1;j<Ny_g2-1;j++)
	                for (i=1;i<Nx_g2-1;i++) {
			     	fprintf(f, "%d %d %d %d\n", (j-1)*(Nx_g2-1)+i, (j-1)*(Nx_g2-1)+i+1, (j)*(Nx_g2-1)+i+1, (j)*(Nx_g2-1)+i);
			}
    
		    	fclose(f);
	



	// read the inflow data on grid g1
	
	FILE *f_g1;
	FILE *f_g2;
    	char fname_g1[80];  
    	char fname_g2[80];  
	int ti;
	int ti2;
	double Time_p, Time_n;
	double Time_;

	int ti_;

	int ti_name;

       	for (j=1; j<Ny_g1; j++)  
       	for (i=1; i<Nx_g1; i++) { 

		U_g1[j][i].x=0.0; 		
		U_g1[j][i].y=0.0; 		
		U_g1[j][i].z=0.0; 		

		uu_g1[j][i]=0.0; 		
		vv_g1[j][i]=0.0; 		
		ww_g1[j][i]=0.0; 		
		uv_g1[j][i]=0.0; 		
		vw_g1[j][i]=0.0; 		
		wu_g1[j][i]=0.0; 		
	}


	for (ti=0;ti<Nt_g1;ti++) {


		printf("Average inflow at ti=%d\n", ti);
		ti2=ti;
		
		if(ti2>inflow_recycle_perioid) {
			ti2 -= (ti2/inflow_recycle_perioid) * inflow_recycle_perioid;
		}
		int ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

		sprintf(fname_g1, "%s/inflow_%06d_dt=%g.dat", path_g1, ti_name, dt_g1);

		if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
		
			if(ti!=0) fclose(f_g1);
			f_g1=fopen(fname_g1, "rb");
			if(!f_g1) printf("\n******************* Cannot open %s ! *******************\n", fname_g1),exit(0);
		}
		for(j=0; j<Ny_g1; j++) fread(&ucat2_g1[j][0], sizeof(Cmpnts), Nx_g1, f_g1);

       		for (j=1; j<Ny_g1; j++)  
	       	for (i=1; i<Nx_g1; i++) { 
 			U_g1[j][i].x+=ucat2_g1[j][i].x;
 			U_g1[j][i].y+=ucat2_g1[j][i].y;
 			U_g1[j][i].z+=ucat2_g1[j][i].z;
			uu_g1[j][i]+=ucat2_g1[j][i].x*ucat2_g1[j][i].x;
			vv_g1[j][i]+=ucat2_g1[j][i].y*ucat2_g1[j][i].y;
			ww_g1[j][i]+=ucat2_g1[j][i].z*ucat2_g1[j][i].z;
			uv_g1[j][i]+=ucat2_g1[j][i].x*ucat2_g1[j][i].y;
			vw_g1[j][i]+=ucat2_g1[j][i].y*ucat2_g1[j][i].z;
			wu_g1[j][i]+=ucat2_g1[j][i].z*ucat2_g1[j][i].x;
		}

	}


       	for (j=1; j<Ny_g1; j++)  
	for (i=1; i<Nx_g1; i++) { 
 		U_g1[j][i].x/=(double)Nt_g1;
 		U_g1[j][i].y/=(double)Nt_g1;
 		U_g1[j][i].z/=(double)Nt_g1;
		uu_g1[j][i]/=(double)Nt_g1;
		uu_g1[j][i]-=U_g1[j][i].x*U_g1[j][i].x;

		vv_g1[j][i]/=(double)Nt_g1;
		vv_g1[j][i]-=U_g1[j][i].y*U_g1[j][i].y;

		ww_g1[j][i]/=(double)Nt_g1;
		ww_g1[j][i]-=U_g1[j][i].z*U_g1[j][i].z;

		uv_g1[j][i]/=(double)Nt_g1;
		uv_g1[j][i]-=U_g1[j][i].x*U_g1[j][i].y;

		vw_g1[j][i]/=(double)Nt_g1;
		vw_g1[j][i]-=U_g1[j][i].y*U_g1[j][i].z;

		wu_g1[j][i]/=(double)Nt_g1;
		wu_g1[j][i]-=U_g1[j][i].z*U_g1[j][i].x;
	}

	FILE *f_avg;
 	char filen_avg[80];

	printf("output tecplot file\n");
        sprintf(filen_avg, "%s/inflow_avg.plt",path_g1);

	f_avg = fopen(filen_avg, "w");
	fprintf(f_avg, "Variables=x,y,z,U,V,W,uu,vv,ww,uv,vw,wu\n");

	fprintf(f_avg, "ZONE T='QUADRILATERAL', N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([1-12]=NODAL)\n", (Nx_g1-1)*(Ny_g1-1), (Nx_g1-2)*(Ny_g1-2));

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", x_g1[j][i]);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", y_g1[j][i]);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", z_g1[j][i]);
	}


	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", U_g1[j][i].x);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", U_g1[j][i].y);
	}


	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", U_g1[j][i].z);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", uu_g1[j][i]);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", vv_g1[j][i]);
	}


	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", ww_g1[j][i]);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", uv_g1[j][i]);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", vw_g1[j][i]);
	}

	for (j=1;j<Ny_g1;j++)
	for (i=1;i<Nx_g1;i++) {
      		fprintf(f_avg, "%le\n", wu_g1[j][i]);
	}


        for (j=1;j<Ny_g1-1;j++)
	for (i=1;i<Nx_g1-1;i++) {
	     	fprintf(f_avg, "%d %d %d %d\n", (j-1)*(Nx_g1-1)+i, (j-1)*(Nx_g1-1)+i+1, (j)*(Nx_g1-1)+i+1, (j)*(Nx_g1-1)+i);
	}
    
    	fclose(f_avg);
	

	//exit(0);

	
	for (ti=0;ti<Nt_g1-1;ti++) {

		if (ti==0) {
			printf("read inflow at ti=%d\n", ti);
			ti2=ti;
		
			if(ti2>inflow_recycle_perioid) {
				ti2 -= (ti2/inflow_recycle_perioid) * inflow_recycle_perioid;
			}
			ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

			sprintf(fname_g1, "%s/inflow_%06d_dt=%g.dat", path_g1, ti_name, dt_g1);

			if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
			
				if(ti!=0) fclose(f_g1);
				f_g1=fopen(fname_g1, "rb");
				if(!f_g1) printf("\n******************* Cannot open %s ! *******************\n", fname_g1),exit(0);
			}
			for(j=0; j<Ny_g1; j++) fread(&ucat1_g1[j][0], sizeof(Cmpnts), Nx_g1, f_g1);

			printf("read inflow at ti=%d\n", ti+1);

			ti2=ti+1;
		
			if(ti2>inflow_recycle_perioid) {
				ti2 -= (ti2/inflow_recycle_perioid) * inflow_recycle_perioid;
			}
			int ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

			sprintf(fname_g1, "%s/inflow_%06d_dt=%g.dat", path_g1, ti_name, dt_g1);

			if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
			
				if(ti!=0) fclose(f_g1);
				f_g1=fopen(fname_g1, "rb");
				if(!f_g1) printf("\n******************* Cannot open %s ! *******************\n", fname_g1),exit(0);
			}
			for(j=0; j<Ny_g1; j++) fread(&ucat2_g1[j][0], sizeof(Cmpnts), Nx_g1, f_g1);
		}
		else {

			printf("read inflow at ti=%d\n", ti+1);
			ti2=ti+1;
		
			if(ti2>inflow_recycle_perioid) {
				ti2 -= (ti2/inflow_recycle_perioid) * inflow_recycle_perioid;
			}
			int ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

			sprintf(fname_g1, "%s/inflow_%06d_dt=%g.dat", path_g1, ti_name, dt_g1);

			if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
			
				if(ti!=0) fclose(f_g1);
				f_g1=fopen(fname_g1, "rb");
				if(!f_g1) printf("\n******************* Cannot open %s ! *******************\n", fname_g1),exit(0);
			}
			for(j=0; j<Ny_g1; j++) fread(&ucat2_g1[j][0], sizeof(Cmpnts), Nx_g1, f_g1);

		}

		Time_p=(double)ti*dt_g1;
		Time_n=(double)(ti+1)*dt_g1;
		
		for(ti_=0; ti_<Nt_g2; ti_++) {

			Time_=(double)ti_*dt_g2;
			if ((Time_-Time_p)>=0.0 && (Time_-Time_n)<0.0) {

				printf("###### interpolate to g2 at ti=%d\n", ti_);
				double fac1, fac2;

				fac1=(Time_n-Time_)/(Time_n-Time_p+1.e-12);
				fac2=(Time_-Time_p)/(Time_n-Time_p+1.e-12);
			
				for(j=1;j<Ny_g2;j++)
				for(i=1;i<Nx_g2;i++) {

					double xx_g2=x_g2[j][i];	
					double yy_g2=y_g2[j][i];	
					double zz_g2=z_g2[j][i];	

					int ii=ii_g2[j][i];
					int jj=jj_g2[j][i];

					double sum_u=0.0;
					double sum_v=0.0;
					double sum_w=0.0;
					double sum_coef=0.0;

					for(j1=jj-2;j1<=jj+2;j1++)
					for(i1=ii-2;i1<=ii+2;i1++) {
						if (i1>0 && i1<Nx_g1 && j1>0 && j1<Ny_g1) {
							double xx_g1=x_g1[j1][i1];	
							double yy_g1=y_g1[j1][i1];	
							double zz_g1=z_g1[j1][i1];	

							double dist=1.0/(sqrt(pow((xx_g2-xx_g1),2)+pow((yy_g2-yy_g1),2)+pow((zz_g2-zz_g1),2))+1.e-12);

							sum_u+=fac1*dist*ucat1_g1[j1][i1].x+fac2*dist*ucat2_g1[j1][i1].x;
							sum_v+=fac1*dist*ucat1_g1[j1][i1].y+fac2*dist*ucat2_g1[j1][i1].y;
							sum_w+=fac1*dist*ucat1_g1[j1][i1].z+fac2*dist*ucat2_g1[j1][i1].z;

							//sum_u+=fac1*dist+fac2*dist;
							//sum_v+=fac1*dist+fac2*dist;
							//sum_w+=fac1*dist+fac2*dist;

							sum_coef+=fac1*dist+fac2*dist;
						}
					}
					/*
					ucat_g2[j][i].x=Scale_inflow*sum_u/(sum_coef+1.e-12);
					ucat_g2[j][i].y=Scale_inflow*sum_v/(sum_coef+1.e-12);
					ucat_g2[j][i].z=Scale_inflow*sum_w/(sum_coef+1.e-12);
					*/
					
					// sediment
					//
					ucat_g2[j][i].y=Scale_inflow*sum_u/(sum_coef+1.e-12);
					ucat_g2[j][i].z=Scale_inflow*sum_v/(sum_coef+1.e-12);
					ucat_g2[j][i].x=Scale_inflow*sum_w/(sum_coef+1.e-12);

				}
		
				for(j=1;j<Ny_g2;j++)
				for(i=1;i<Nx_g2;i++) {
					if (nvert_g2[j][i]>1.0) {
	                                        ucat_g2[j][i].x=0.0;
        	                                ucat_g2[j][i].y=0.0;
                	                        ucat_g2[j][i].z=0.0;
					}
				}

				/*	
				for(j=0;j<Ny_g2;j++)
				for(i=0;i<Nx_g2;i++) {
					// noslip
					if (bc[0]==0 && i==0) {
						ucat_g2[j][i].x=0.0;
						ucat_g2[j][i].y=0.0;
						ucat_g2[j][i].z=0.0;
					}
					if (bc[1]==0 && i==Nx_g2-1) {
						ucat_g2[j][i].x=0.0;
						ucat_g2[j][i].y=0.0;
						ucat_g2[j][i].z=0.0;
					}
					if (bc[2]==0 && j==0) {
						ucat_g2[j][i].x=0.0;
						ucat_g2[j][i].y=0.0;
						ucat_g2[j][i].z=0.0;
					}
					if (bc[3]==0 && j==Ny_g2-1) {
						ucat_g2[j][i].x=0.0;
						ucat_g2[j][i].y=0.0;
						ucat_g2[j][i].z=0.0;
					}

					// neumann
					if (bc[0]==1 && i==0) {
						ucat_g2[j][i].x=0.0;
						ucat_g2[j][i].y=ucat_g2[j][1].y;
						ucat_g2[j][i].z=ucat_g2[j][1].z;
					}

					if (bc[1]==1 && i==Nx_g2-1) {
						ucat_g2[j][i].x=0.0;
						ucat_g2[j][i].y=ucat_g2[j][Nx_g2-2].y;
						ucat_g2[j][i].z=ucat_g2[j][Nx_g2-2].z;
					}
					if (bc[2]==1 && j==0) {
						ucat_g2[j][i].x=ucat_g2[1][i].x;
						ucat_g2[j][i].y=0.0;
						ucat_g2[j][i].z=ucat_g2[1][i].z;
					}
					if (bc[3]==1 && j==Ny_g2-1) {
						ucat_g2[j][i].x=ucat_g2[Ny_g2-2][i].x;
						ucat_g2[j][i].y=0.0;
						ucat_g2[j][i].z=ucat_g2[Ny_g2-2][i].z;
					}

					// Periodic
					if (bc[0]==2 && i==0) {
						ucat_g2[j][i]=ucat_g2[j][Nx_g2-2];
					}

					if (bc[1]==2 && i==Nx_g2-1) {
						ucat_g2[j][i]=ucat_g2[j][1];
					}
					if (bc[2]==2 && j==0) {
						ucat_g2[j][i]=ucat_g2[Ny_g2-2][i];
					}
					if (bc[3]==2 && j==Ny_g2-1) {
						ucat_g2[j][i]=ucat_g2[1][i];
					}


				}

				*/

				printf("write inflow of g2 at ti=%d\n", ti_);
	        	        ti_name = ( (ti_-1) / save_inflow_period ) * save_inflow_period + 1;

        	        	sprintf(fname_g2, "%s/inflow_%06d_dt=%g.dat", path_g2, ti_name, dt_g2);
	                	if(ti_==0 || ti_%save_inflow_period==1) unlink(fname_g2);    // delete existing file

        	        	f_g2=fopen(fname_g2, "ab");
                		if(!f_g2) printf("\n******************* Cannot open %s ! *******************\n", fname_g2),exit(0);

	                	for(j=0; j<Ny_g2; j++) fwrite(&ucat_g2[j][0], sizeof(Cmpnts), Nx_g2, f_g2);
	                	fclose(f_g2);


			}

		}

		for(j=0;j<Ny_g1;j++)
		for(i=0;i<Nx_g1;i++) {
			ucat1_g1[j][i]=ucat2_g1[j][i];
		}
		


	}



	for(ti_=0; ti_<Nt_g2; ti_++) {

		ti2=ti_;
		
		if(ti2>Nt_g2) {
			ti2 -= (ti2/Nt_g2) * Nt_g2;
		}
		ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

		sprintf(fname_g2, "%s/inflow_%06d_dt=%g.dat", path_g2, ti_name, dt_g2);

		if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
			
			if(ti_!=0) fclose(f_g2);
			f_g2=fopen(fname_g2, "rb");
			if(!f_g2) printf("\n******************* Cannot open %s ! *******************\n", fname_g2),exit(0);
		}
		for(j=0; j<Ny_g2; j++) fread(&ucat_g2[j][0], sizeof(Cmpnts), Nx_g2, f_g2);


		if(ti_%50==0) {
	      		printf("output tecplot file\n");
    			FILE *f;
	    		char filen_tec[80];
                	sprintf(filen_tec, "%s/inflowcheck_%06d_dt=%g.plt", path_g2, ti_, dt_g2);
//	    		sprintf(filen_tec, "InFlowField%06d.dat",t_LES);
	    		f = fopen(filen_tec, "w");
	    		fprintf(f, "Variables=x,y,z,U,V,W\n");

			fprintf(f, "ZONE T='QUADRILATERAL', N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([1-6]=NODAL)\n", (Nx_g2-1)*(Ny_g2-1), (Nx_g2-2)*(Ny_g2-2));

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
      				fprintf(f, "%le\n", x_g2[j][i]);
			}

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
	                        fprintf(f, "%le\n", y_g2[j][i]);
	                }

        	        for (j=1;j<Ny_g2;j++)
                	for (i=1;i<Nx_g2;i++) {
	                       	fprintf(f, "%le\n", z_g2[j][i]);
	                }


	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
	                        fprintf(f, "%le\n", ucat_g2[j][i].x);
	                }

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
        	                fprintf(f, "%le\n", ucat_g2[j][i].y);
	                }

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
        	                fprintf(f, "%le\n", ucat_g2[j][i].z);
	                }

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
        	        //        fprintf(f, "%d\n", ii_g2[j][i]);
	                }

	                for (j=1;j<Ny_g2;j++)
	                for (i=1;i<Nx_g2;i++) {
        	        //        fprintf(f, "%d\n", jj_g2[j][i]);
	                }



                	for (j=1;j<Ny_g2-1;j++)
	                for (i=1;i<Nx_g2-1;i++) {
			     	fprintf(f, "%d %d %d %d\n", (j-1)*(Nx_g2-1)+i, (j-1)*(Nx_g2-1)+i+1, (j)*(Nx_g2-1)+i+1, (j)*(Nx_g2-1)+i);
			}
    
		    	fclose(f);
	
		}


	}




	for(ti_=0; ti_<Nt_g1; ti_++) {

		ti2=ti_;
		
		if(ti2>inflow_recycle_perioid) {
			ti2 -= (ti2/inflow_recycle_perioid) * inflow_recycle_perioid;
		}
		ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

		sprintf(fname_g1, "%s/inflow_%06d_dt=%g.dat", path_g1, ti_name, dt_g1);

		if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
			
			if(ti_!=0) fclose(f_g1);
			f_g1=fopen(fname_g1, "rb");
			if(!f_g1) printf("\n******************* Cannot open %s ! *******************\n", fname_g1),exit(0);
		}
		for(j=0; j<Ny_g1; j++) fread(&ucat1_g1[j][0], sizeof(Cmpnts), Nx_g1, f_g1);


		if(ti_%save_inflow_period==0) {
	      		printf("output tecplot file\n");
    			FILE *f;
	    		char filen_tec[80];
                	sprintf(filen_tec, "%s/inflowcheck_%06d_dt=%g.plt", path_g1, ti_, dt_g1);
//	    		sprintf(filen_tec, "InFlowField%06d.dat",t_LES);
	    		f = fopen(filen_tec, "w");
	    		fprintf(f, "Variables=x,y,z,U,V,W\n");

			fprintf(f, "ZONE T='QUADRILATERAL', N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([1-6]=NODAL)\n", (Nx_g1-1)*(Ny_g1-1), (Nx_g1-2)*(Ny_g1-2));

	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
      				fprintf(f, "%le\n", x_g1[j][i]);
			}

	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
	                        fprintf(f, "%le\n", y_g1[j][i]);
	                }

        	        for (j=1;j<Ny_g1;j++)
                	for (i=1;i<Nx_g1;i++) {
	                       	fprintf(f, "%le\n", z_g1[j][i]);
	                }


	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
	                        fprintf(f, "%le\n", ucat1_g1[j][i].x);
	                }

	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
        	                fprintf(f, "%le\n", ucat1_g1[j][i].y);
	                }

	                for (j=1;j<Ny_g1;j++)
	                for (i=1;i<Nx_g1;i++) {
        	                fprintf(f, "%le\n", ucat1_g1[j][i].z);
	                }


                	for (j=1;j<Ny_g1-1;j++)
	                for (i=1;i<Nx_g1-1;i++) {
			     	fprintf(f, "%d %d %d %d\n", (j-1)*(Nx_g1-1)+i, (j-1)*(Nx_g1-1)+i+1, (j)*(Nx_g1-1)+i+1, (j)*(Nx_g1-1)+i);
			}
    
		    	fclose(f);
	
		}


	}

	// ----------------

	// read the inflow data on grid g2
	
	std::vector< std::vector<Cmpnts> > U_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) U_g2[j].resize(Nx_g2);
	
	std::vector< std::vector<double> > uu_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) uu_g2[j].resize(Nx_g2);
						
	std::vector< std::vector<double> > vv_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) vv_g2[j].resize(Nx_g2);

	std::vector< std::vector<double> > ww_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) ww_g2[j].resize(Nx_g2);

	std::vector< std::vector<double> > uv_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) uv_g2[j].resize(Nx_g2);

	std::vector< std::vector<double> > vw_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) vw_g2[j].resize(Nx_g2);

	std::vector< std::vector<double> > wu_g2 (Ny_g2);
	for( j=0; j<Ny_g2; j++) wu_g2[j].resize(Nx_g2);



       	for (j=1; j<Ny_g2; j++)  
       	for (i=1; i<Nx_g2; i++) { 

		U_g2[j][i].x=0.0; 		
		U_g2[j][i].y=0.0; 		
		U_g2[j][i].z=0.0; 		

		uu_g2[j][i]=0.0; 		
		vv_g2[j][i]=0.0; 		
		ww_g2[j][i]=0.0; 		
		uv_g2[j][i]=0.0; 		
		vw_g2[j][i]=0.0; 		
		wu_g2[j][i]=0.0; 		
	}


	for (ti=0;ti<Nt_g2;ti++) {


		printf("Average inflow at ti=%d\n", ti);
		ti2=ti;
		
		if(ti2>inflow_recycle_perioid) {
			ti2 -= (ti2/inflow_recycle_perioid) * inflow_recycle_perioid;
		}
		int ti_name = ( ti2 / read_inflow_period ) * read_inflow_period + 1;

		sprintf(fname_g2, "%s/inflow_%06d_dt=%g.dat", path_g2, ti_name, dt_g2);

		if(ti2==0 || (ti2>90 && ti2%read_inflow_period==1) ) {
		
			if(ti!=0) fclose(f_g2);
			f_g2=fopen(fname_g2, "rb");
			if(!f_g2) printf("\n******************* Cannot open %s ! *******************\n", fname_g2),exit(0);
		}
		for(j=0; j<Ny_g2; j++) fread(&ucat_g2[j][0], sizeof(Cmpnts), Nx_g2, f_g2);

       		for (j=1; j<Ny_g2; j++)  
	       	for (i=1; i<Nx_g2; i++) { 
 			U_g2[j][i].x+=ucat_g2[j][i].x;
 			U_g2[j][i].y+=ucat_g2[j][i].y;
 			U_g2[j][i].z+=ucat_g2[j][i].z;
			uu_g2[j][i]+=ucat_g2[j][i].x*ucat_g2[j][i].x;
			vv_g2[j][i]+=ucat_g2[j][i].y*ucat_g2[j][i].y;
			ww_g2[j][i]+=ucat_g2[j][i].z*ucat_g2[j][i].z;
			uv_g2[j][i]+=ucat_g2[j][i].x*ucat_g2[j][i].y;
			vw_g2[j][i]+=ucat_g2[j][i].y*ucat_g2[j][i].z;
			wu_g2[j][i]+=ucat_g2[j][i].z*ucat_g2[j][i].x;
		}

	}


       	for (j=1; j<Ny_g2; j++)  
	for (i=1; i<Nx_g2; i++) { 
 		U_g2[j][i].x/=(double)Nt_g2;
 		U_g2[j][i].y/=(double)Nt_g2;
 		U_g2[j][i].z/=(double)Nt_g2;
		uu_g2[j][i]/=(double)Nt_g2;
		uu_g2[j][i]-=U_g2[j][i].x*U_g2[j][i].x;

		vv_g2[j][i]/=(double)Nt_g2;
		vv_g2[j][i]-=U_g2[j][i].y*U_g2[j][i].y;

		ww_g2[j][i]/=(double)Nt_g2;
		ww_g2[j][i]-=U_g2[j][i].z*U_g2[j][i].z;

		uv_g2[j][i]/=(double)Nt_g2;
		uv_g2[j][i]-=U_g2[j][i].x*U_g2[j][i].y;

		vw_g2[j][i]/=(double)Nt_g2;
		vw_g2[j][i]-=U_g2[j][i].y*U_g2[j][i].z;

		wu_g2[j][i]/=(double)Nt_g2;
		wu_g2[j][i]-=U_g2[j][i].z*U_g2[j][i].x;
	}

	//FILE *f_avg;
 	//char filen_avg[80];

	printf("output tecplot file\n");
        sprintf(filen_avg, "%s/inflow_avg_g2.plt",path_g2);

	f_avg = fopen(filen_avg, "w");
	fprintf(f_avg, "Variables=x,y,z,U,V,W,uu,vv,ww,uv,vw,wu\n");

	fprintf(f_avg, "ZONE T='QUADRILATERAL', N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([1-12]=NODAL)\n", (Nx_g2-1)*(Ny_g2-1), (Nx_g2-2)*(Ny_g2-2));

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", x_g2[j][i]);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", y_g2[j][i]);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", z_g2[j][i]);
	}


	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", U_g2[j][i].x);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", U_g2[j][i].y);
	}


	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", U_g2[j][i].z);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", uu_g2[j][i]);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", vv_g2[j][i]);
	}


	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", ww_g2[j][i]);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", uv_g2[j][i]);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", vw_g2[j][i]);
	}

	for (j=1;j<Ny_g2;j++)
	for (i=1;i<Nx_g2;i++) {
      		fprintf(f_avg, "%le\n", wu_g2[j][i]);
	}


        for (j=1;j<Ny_g2-1;j++)
	for (i=1;i<Nx_g2-1;i++) {
	     	fprintf(f_avg, "%d %d %d %d\n", (j-1)*(Nx_g2-1)+i, (j-1)*(Nx_g2-1)+i+1, (j)*(Nx_g2-1)+i+1, (j)*(Nx_g2-1)+i);
	}
    
    	fclose(f_avg);
	


}





