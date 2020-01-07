/* ml9_twin.c: Twin domain for J. Chem. Phys.  */

/* Chain Parameters */
#define N       50		
#define b0	31.84e-10	/*22.5e-10; 3.4e-10; 20.4e-10*/   
#define R_hyd	15.92e-10	/*11.25e-10; 10.0e-10; 10.2e-10*/ 
#define DELTA	0.008*b0  
#define PSI	sqrt(b0/500.0e-10)	/* bending parameter (in radiants) */ 
#define ZI	sqrt(b0*kB*temp/2.6e-28)	/* torsion parameter (C_t=2.6e-28 J m) */
#define EPSILON 100.0		/*1.9e-6*/ 
#define SIGMA	b0/pow(2.0,1.0/6.0)	/*3.03e-10; 18.17e-10*/	
#define CUTOFF	b0	/*pow(2.0,1.0/6.0)*SIGMA*/

/* Solvent Parameters */
#define kB	1.3806503e-23	
#define temp	310.0   	
#define ETA	0.01	/*0.007, 0.01*/	
#define D_rot   kB*temp/(Pi*ETA*pow(R_hyd,2.0)*b0)
#define D_trans kB*temp/(6.0*Pi*ETA*R_hyd) 

/* Numerical Parameters */
#define dt      5.0e-12		/*0.1e-12; 5.0e-12*/	
#define Nt      20000000		/*5000000*/   
#define SEED	-3    

#define Pi	3.141592653589793238462643	

#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h" 

/* Subroutines for Stretching Force Components */
double F_s_x_sub(int i, int k, double **r_ij, double **u_norm_x)
{
        double f_s_x;

if(i<=25){
        if(i==1){
                f_s_x=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_x[i][k];  
        }
	if(i==25){
                f_s_x=-( kB*temp/pow(DELTA,2.0) )*(r_ij[25][25-1]-b0)*u_norm_x[25-1][k];
        }
	if(i!=1 && i!=25){
		f_s_x=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_x[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_x[i-1][k] );
        }
}
if(i>25){
        if(i==26){
                f_s_x=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_x[i][k];  
        }
	if(i==N){
                f_s_x=-( kB*temp/pow(DELTA,2.0) )*(r_ij[N][N-1]-b0)*u_norm_x[N-1][k];
        }
	if(i!=26 && i!=N){
		f_s_x=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_x[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_x[i-1][k] );
        }
}
        return f_s_x;
}
double F_s_y_sub(int i, int k, double **r_ij, double **u_norm_y)
{
        double f_s_y;

if(i<=25){
        if(i==1){
                f_s_y=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_y[i][k];  
        }
	if(i==25){
                f_s_y=-( kB*temp/pow(DELTA,2.0) )*(r_ij[25][25-1]-b0)*u_norm_y[25-1][k];
        }
        if(i!=1 && i!=25){
                f_s_y=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_y[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_y[i-1][k] );
        }
}
if(i>25){
        if(i==26){
                f_s_y=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_y[i][k];  
        }
	if(i==N){
                f_s_y=-( kB*temp/pow(DELTA,2.0) )*(r_ij[N][N-1]-b0)*u_norm_y[N-1][k];
        }
        if(i!=26 && i!=N){
                f_s_y=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_y[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_y[i-1][k] );
        }
}
        return f_s_y;
}
double F_s_z_sub(int i, int k, double **r_ij, double **u_norm_z)
{
        double f_s_z;

if(i<=25){
        if(i==1){
		f_s_z=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_z[i][k];
        }
        if(i==25){
                f_s_z=-( kB*temp/pow(DELTA,2.0) )*(r_ij[25][25-1]-b0)*u_norm_z[25-1][k];
        }
        if(i!=1 && i!=25){
                f_s_z=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_z[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_z[i-1][k] );
        }
}
if(i>25){
        if(i==26){
                f_s_z=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_z[i][k];  
        }
	if(i==N){
                f_s_z=-( kB*temp/pow(DELTA,2.0) )*(r_ij[N][N-1]-b0)*u_norm_z[N-1][k];
        }
        if(i!=26 && i!=N){
                f_s_z=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_z[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_z[i-1][k] );
        }
}
        return f_s_z;
}

/* Subroutines for Bending Force Components */
double F_b_x_sub(int i, int k, double *beta, double **r_ij, double **u_norm_x)
{
        double f_b_x;
	double *A_x,*B_x; 

	A_x=dvector(1,N),B_x=dvector(1,N);

if(i<=25){
        if(i==1){
		A_x[1]=( beta[1]/(r_ij[1][2]*sin(beta[1])) )*
			( u_norm_x[2][k]-u_norm_x[1][k]*cos(beta[1]) );

                f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[1] ); 
        }
      	if(i==2){
		A_x[2]=( beta[2]/(r_ij[3][2]*sin(beta[2])) )*
                        ( u_norm_x[3][k]-u_norm_x[2][k]*cos(beta[2]) );
                B_x[2]=( beta[1]/(r_ij[3][2]*sin(beta[1])) )*
                        ( u_norm_x[1][k]-u_norm_x[2][k]*cos(beta[1]) );
                A_x[1]=( beta[1]/(r_ij[1][2]*sin(beta[1])) )*
                        ( u_norm_x[2][k]-u_norm_x[1][k]*cos(beta[1]) );

		f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[2]-A_x[1]+B_x[2] );
	}
	if(i!=1 && i!=2 && i!=25){
          	A_x[i]=( beta[i]/(r_ij[i][i+1]*sin(beta[i])) )*
                        ( u_norm_x[i+1][k]-u_norm_x[i][k]*cos(beta[i]) );
                B_x[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_x[i-1][k]-u_norm_x[i][k]*cos(beta[i-1]) );
                A_x[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_x[i][k]-u_norm_x[i-1][k]*cos(beta[i-1]) );
                B_x[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_x[i-2][k]-u_norm_x[i-1][k]*cos(beta[i-2]) );    

		f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[i]-A_x[i-1]+B_x[i]-B_x[i-1] ); 
        }
        if(i==25){
                A_x[25-1]=( beta[25-1]/(r_ij[25-1][25]*sin(beta[25-1])) )*
                        ( u_norm_x[25][k]-u_norm_x[25-1][k]*cos(beta[25-1]) );
                B_x[25-1]=( beta[25-2]/(r_ij[25-1][25]*sin(beta[25-2])) )*
                        ( u_norm_x[25-2][k]-u_norm_x[25-1][k]*cos(beta[25-2]) );

                f_b_x=( kB*temp/pow(PSI,2.0) )*( A_x[25-1]+B_x[25-1] );
        }
}
if(i>25){
        if(i==26){
		A_x[26]=( beta[26]/(r_ij[26][27]*sin(beta[26])) )*
			( u_norm_x[27][k]-u_norm_x[26][k]*cos(beta[26]) );

                f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[26] ); 
        }
      	if(i==27){
		A_x[27]=( beta[27]/(r_ij[28][27]*sin(beta[27])) )*
                        ( u_norm_x[28][k]-u_norm_x[27][k]*cos(beta[27]) );
                B_x[27]=( beta[26]/(r_ij[28][27]*sin(beta[26])) )*
                        ( u_norm_x[26][k]-u_norm_x[27][k]*cos(beta[26]) );
                A_x[26]=( beta[26]/(r_ij[26][27]*sin(beta[26])) )*
                        ( u_norm_x[27][k]-u_norm_x[26][k]*cos(beta[26]) );

		f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[27]-A_x[26]+B_x[27] );
	}
	if(i!=26 && i!=27 && i!=N){
          	A_x[i]=( beta[i]/(r_ij[i][i+1]*sin(beta[i])) )*
                        ( u_norm_x[i+1][k]-u_norm_x[i][k]*cos(beta[i]) );
                B_x[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_x[i-1][k]-u_norm_x[i][k]*cos(beta[i-1]) );
                A_x[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_x[i][k]-u_norm_x[i-1][k]*cos(beta[i-1]) );
                B_x[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_x[i-2][k]-u_norm_x[i-1][k]*cos(beta[i-2]) );    

		f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[i]-A_x[i-1]+B_x[i]-B_x[i-1] ); 
        }
        if(i==N){
                A_x[N-1]=( beta[N-1]/(r_ij[N-1][N]*sin(beta[N-1])) )*
                        ( u_norm_x[N][k]-u_norm_x[N-1][k]*cos(beta[N-1]) );
                B_x[N-1]=( beta[N-2]/(r_ij[N-1][N]*sin(beta[N-2])) )*
                        ( u_norm_x[N-2][k]-u_norm_x[N-1][k]*cos(beta[N-2]) );

                f_b_x=( kB*temp/pow(PSI,2.0) )*( A_x[N-1]+B_x[N-1] );
        }
}
       	free_dvector(A_x,1,N);
        free_dvector(B_x,1,N); 

	return f_b_x;
}
double F_b_y_sub(int i, int k, double *beta, double **r_ij, double **u_norm_y)
{
	double f_b_y;
        double *A_y,*B_y;

        A_y=dvector(1,N),B_y=dvector(1,N);

if(i<=25){
        if(i==1){
                A_y[1]=( beta[1]/(r_ij[1][2]*sin(beta[1])) )*
                        ( u_norm_y[2][k]-u_norm_y[1][k]*cos(beta[1]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[1] );                     
        }
        if(i==2){
                A_y[2]=( beta[2]/(r_ij[3][2]*sin(beta[2])) )*
                        ( u_norm_y[3][k]-u_norm_y[2][k]*cos(beta[2]) );
                B_y[2]=( beta[1]/(r_ij[3][2]*sin(beta[1])) )*
                        ( u_norm_y[1][k]-u_norm_y[2][k]*cos(beta[1]) );
                A_y[1]=( beta[1]/(r_ij[1][2]*sin(beta[1])) )*
                        ( u_norm_y[2][k]-u_norm_y[1][k]*cos(beta[1]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[2]-A_y[1]+B_y[2] );
        }
        if(i!=1 && i!=2 && i!=25){
                A_y[i]=( beta[i]/(r_ij[i][i+1]*sin(beta[i])) )*
                        ( u_norm_y[i+1][k]-u_norm_y[i][k]*cos(beta[i]) );
                B_y[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_y[i-1][k]-u_norm_y[i][k]*cos(beta[i-1]) );
                A_y[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_y[i][k]-u_norm_y[i-1][k]*cos(beta[i-1]) );
                B_y[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_y[i-2][k]-u_norm_y[i-1][k]*cos(beta[i-2]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[i]-A_y[i-1]+B_y[i]-B_y[i-1] );
        }
        if(i==25){
                A_y[25-1]=( beta[25-1]/(r_ij[25-1][25]*sin(beta[25-1])) )*
                        ( u_norm_y[25][k]-u_norm_y[25-1][k]*cos(beta[25-1]) );
                B_y[25-1]=( beta[25-2]/(r_ij[25-1][25]*sin(beta[25-2])) )*
                        ( u_norm_y[25-2][k]-u_norm_y[25-1][k]*cos(beta[25-2]) );

                f_b_y=( kB*temp/pow(PSI,2.0) )*( A_y[25-1]+B_y[25-1] );
        }
}
if(i>25){
        if(i==26){
                A_y[26]=( beta[26]/(r_ij[26][27]*sin(beta[26])) )*
                        ( u_norm_y[27][k]-u_norm_y[26][k]*cos(beta[26]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[26] );                     
        }
        if(i==27){
                A_y[27]=( beta[27]/(r_ij[28][27]*sin(beta[27])) )*
                        ( u_norm_y[28][k]-u_norm_y[27][k]*cos(beta[27]) );
                B_y[27]=( beta[26]/(r_ij[28][27]*sin(beta[26])) )*
                        ( u_norm_y[26][k]-u_norm_y[27][k]*cos(beta[26]) );
                A_y[26]=( beta[26]/(r_ij[26][27]*sin(beta[26])) )*
                        ( u_norm_y[27][k]-u_norm_y[26][k]*cos(beta[26]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[27]-A_y[26]+B_y[27] );
        }
        if(i!=26 && i!=27 && i!=N){
                A_y[i]=( beta[i]/(r_ij[i][i+1]*sin(beta[i])) )*
                        ( u_norm_y[i+1][k]-u_norm_y[i][k]*cos(beta[i]) );
                B_y[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_y[i-1][k]-u_norm_y[i][k]*cos(beta[i-1]) );
                A_y[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_y[i][k]-u_norm_y[i-1][k]*cos(beta[i-1]) );
                B_y[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_y[i-2][k]-u_norm_y[i-1][k]*cos(beta[i-2]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[i]-A_y[i-1]+B_y[i]-B_y[i-1] );
        }
        if(i==N){
                A_y[N-1]=( beta[N-1]/(r_ij[N-1][N]*sin(beta[N-1])) )*
                        ( u_norm_y[N][k]-u_norm_y[N-1][k]*cos(beta[N-1]) );
                B_y[N-1]=( beta[N-2]/(r_ij[N-1][N]*sin(beta[N-2])) )*
                        ( u_norm_y[N-2][k]-u_norm_y[N-1][k]*cos(beta[N-2]) );

                f_b_y=( kB*temp/pow(PSI,2.0) )*( A_y[N-1]+B_y[N-1] );
        }
}
        free_dvector(A_y,1,N);
        free_dvector(B_y,1,N);

        return f_b_y;
}
double F_b_z_sub(int i, int k, double *beta, double **r_ij, double **u_norm_z)
{
	double f_b_z;
        double *A_z,*B_z;

        A_z=dvector(1,N),B_z=dvector(1,N);

if(i<=25){
        if(i==1){
                A_z[1]=( beta[1]/(r_ij[1][2]*sin(beta[1])) )*
                        ( u_norm_z[2][k]-u_norm_z[1][k]*cos(beta[1]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[1] );                     
        }
        if(i==2){
                A_z[2]=( beta[2]/(r_ij[3][2]*sin(beta[2])) )*
                        ( u_norm_z[3][k]-u_norm_z[2][k]*cos(beta[2]) );
                B_z[2]=( beta[1]/(r_ij[3][2]*sin(beta[1])) )*
                        ( u_norm_z[1][k]-u_norm_z[2][k]*cos(beta[1]) );
                A_z[1]=( beta[1]/(r_ij[1][2]*sin(beta[1])) )*
                        ( u_norm_z[2][k]-u_norm_z[1][k]*cos(beta[1]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[2]-A_z[1]+B_z[2] );
        }
        if(i!=1 && i!=2 && i!=25){
                A_z[i]=( beta[i]/(r_ij[i][i+1]*sin(beta[i])) )*
                        ( u_norm_z[i+1][k]-u_norm_z[i][k]*cos(beta[i]) );
                B_z[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_z[i-1][k]-u_norm_z[i][k]*cos(beta[i-1]) );
                A_z[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_z[i][k]-u_norm_z[i-1][k]*cos(beta[i-1]) );
                B_z[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_z[i-2][k]-u_norm_z[i-1][k]*cos(beta[i-2]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[i]-A_z[i-1]+B_z[i]-B_z[i-1] );
        }
        if(i==25){
                A_z[25-1]=( beta[25-1]/(r_ij[25-1][25]*sin(beta[25-1])) )*
                        ( u_norm_z[25][k]-u_norm_z[25-1][k]*cos(beta[25-1]) );
                B_z[25-1]=( beta[25-2]/(r_ij[25-1][25]*sin(beta[25-2])) )*
                        ( u_norm_z[25-2][k]-u_norm_z[25-1][k]*cos(beta[25-2]) );

                f_b_z=( kB*temp/pow(PSI,2.0) )*( A_z[25-1]+B_z[25-1] );
        }
}
if(i>25){
        if(i==26){
                A_z[26]=( beta[26]/(r_ij[26][27]*sin(beta[26])) )*
                        ( u_norm_z[27][k]-u_norm_z[26][k]*cos(beta[26]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[26] );                     
        }
        if(i==27){
                A_z[27]=( beta[27]/(r_ij[28][27]*sin(beta[27])) )*
                        ( u_norm_z[28][k]-u_norm_z[27][k]*cos(beta[27]) );
                B_z[27]=( beta[26]/(r_ij[28][27]*sin(beta[26])) )*
                        ( u_norm_z[26][k]-u_norm_z[27][k]*cos(beta[26]) );
                A_z[26]=( beta[26]/(r_ij[26][27]*sin(beta[26])) )*
                        ( u_norm_z[27][k]-u_norm_z[26][k]*cos(beta[26]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[27]-A_z[26]+B_z[27] );
        }
        if(i!=26 && i!=27 && i!=N){
                A_z[i]=( beta[i]/(r_ij[i][i+1]*sin(beta[i])) )*
                        ( u_norm_z[i+1][k]-u_norm_z[i][k]*cos(beta[i]) );
                B_z[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_z[i-1][k]-u_norm_z[i][k]*cos(beta[i-1]) );
                A_z[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_z[i][k]-u_norm_z[i-1][k]*cos(beta[i-1]) );
                B_z[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_z[i-2][k]-u_norm_z[i-1][k]*cos(beta[i-2]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[i]-A_z[i-1]+B_z[i]-B_z[i-1] );
        }
        if(i==N){
                A_z[N-1]=( beta[N-1]/(r_ij[N-1][N]*sin(beta[N-1])) )*
                        ( u_norm_z[N][k]-u_norm_z[N-1][k]*cos(beta[N-1]) );
                B_z[N-1]=( beta[N-2]/(r_ij[N-1][N]*sin(beta[N-2])) )*
                        ( u_norm_z[N-2][k]-u_norm_z[N-1][k]*cos(beta[N-2]) );

                f_b_z=( kB*temp/pow(PSI,2.0) )*( A_z[N-1]+B_z[N-1] );
        }
}
        free_dvector(A_z,1,N);
        free_dvector(B_z,1,N);

        return f_b_z;
}

/* Subroutines for Twisting Force Components */
double F_t_x_sub(int i, int k, double **r_ij, double **u_norm_x, double **u_norm_y, 
	double **u_norm_z, double *f_norm_x, double *f_norm_y, double *f_norm_z,
	double *v_norm_x, double *v_norm_y, double *v_norm_z, double **a_p_g)
{
        double f_t_x;
        double *A_plus_x,*B_plus_x;

        A_plus_x=dvector(1,N),B_plus_x=dvector(1,N);

if(i<=25){
	if(i==1){
		B_plus_x[1]=( 1.0/( r_ij[1][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*
			( -f_norm_x[1]*(u_norm_x[2][k]*v_norm_x[1]+u_norm_y[2][k]*v_norm_y[1]+
			u_norm_z[2][k]*v_norm_z[1])+v_norm_x[1]*(u_norm_x[2][k]*f_norm_x[1]+
			u_norm_y[2][k]*f_norm_y[1]+u_norm_z[2][k]*f_norm_z[1]) );

                f_t_x=( kB*temp/pow(ZI,2.0) )*a_p_g[1][k]*B_plus_x[1];
        }
        if(i==2){
     		A_plus_x[1]=( 1.0/( r_ij[3][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*
                        ( f_norm_x[2]*(v_norm_x[2]*u_norm_x[1][k]+v_norm_y[2]*u_norm_y[1][k]
                        +v_norm_z[2]*u_norm_z[1][k])-v_norm_x[2]*(f_norm_x[2]*u_norm_x[1][k]
                        +f_norm_y[2]*u_norm_y[1][k]+f_norm_z[2]*u_norm_z[1][k]) );
                B_plus_x[2]=( 1.0/( r_ij[3][2]*(1.0+u_norm_x[3][k]*u_norm_x[2][k]+
                        u_norm_y[3][k]*u_norm_y[2][k]+u_norm_z[3][k]*u_norm_z[2][k]) ) )*
                        ( -f_norm_x[2]*(u_norm_x[3][k]*v_norm_x[2]+u_norm_y[3][k]*v_norm_y[2]+
                        u_norm_z[3][k]*v_norm_z[2])+v_norm_x[2]*(u_norm_x[3][k]*f_norm_x[2]+
                        u_norm_y[3][k]*f_norm_y[2]+u_norm_z[3][k]*f_norm_z[2]) );
                B_plus_x[1]=( 1.0/( r_ij[1][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+      
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*       
                        ( -f_norm_x[1]*(u_norm_x[2][k]*v_norm_x[1]+u_norm_y[2][k]*v_norm_y[1]+  
                        u_norm_z[2][k]*v_norm_z[1])+v_norm_x[1]*(u_norm_x[2][k]*f_norm_x[1]+    
                        u_norm_y[2][k]*f_norm_y[1]+u_norm_z[2][k]*f_norm_z[1]) );       

                f_t_x=( kB*temp/pow(ZI,2.0) )*( a_p_g[2][k]*B_plus_x[2]+
                        a_p_g[1][k]*A_plus_x[1]-a_p_g[1][k]*B_plus_x[1] );
	}
        if(i!=1 && i!=2 && i!=25){
       		A_plus_x[i-1]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( f_norm_x[i]*(v_norm_x[i]*u_norm_x[i-1][k]+v_norm_y[i]*u_norm_y[i-1][k]
                        +v_norm_z[i]*u_norm_z[i-1][k])-v_norm_x[i]*(f_norm_x[i]*u_norm_x[i-1][k]
                        +f_norm_y[i]*u_norm_y[i-1][k]+f_norm_z[i]*u_norm_z[i-1][k]) );
                A_plus_x[i-2]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i-1][k]*u_norm_x[i-2][k]+
                        u_norm_y[i-1][k]*u_norm_y[i-2][k]+u_norm_z[i-1][k]*u_norm_z[i-2][k]) ) )*
                        ( f_norm_x[i-1]*(v_norm_x[i-1]*u_norm_x[i-2][k]+v_norm_y[i-1]*u_norm_y[i-2][k]
                        +v_norm_z[i-1]*u_norm_z[i-2][k])-v_norm_x[i-1]*(f_norm_x[i-1]*u_norm_x[i-2][k]
                        +f_norm_y[i-1]*u_norm_y[i-2][k]+f_norm_z[i-1]*u_norm_z[i-2][k]) );
                B_plus_x[i]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i+1][k]*u_norm_x[i][k]+
                        u_norm_y[i+1][k]*u_norm_y[i][k]+u_norm_z[i+1][k]*u_norm_z[i][k]) ) )*
                        ( -f_norm_x[i]*(u_norm_x[i+1][k]*v_norm_x[i]+u_norm_y[i+1][k]*v_norm_y[i]+
                        u_norm_z[i+1][k]*v_norm_z[i])+v_norm_x[i]*(u_norm_x[i+1][k]*f_norm_x[i]+
                        u_norm_y[i+1][k]*f_norm_y[i]+u_norm_z[i+1][k]*f_norm_z[i]) );
                B_plus_x[i-1]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( -f_norm_x[i-1]*(u_norm_x[i][k]*v_norm_x[i-1]+u_norm_y[i][k]*v_norm_y[i-1]+
                        u_norm_z[i][k]*v_norm_z[i-1])+v_norm_x[i-1]*(u_norm_x[i][k]*f_norm_x[i-1]+
                        u_norm_y[i][k]*f_norm_y[i-1]+u_norm_z[i][k]*f_norm_z[i-1]) );

                f_t_x=( kB*temp/pow(ZI,2.0) )*( a_p_g[i][k]*B_plus_x[i]+
                        a_p_g[i-1][k]*A_plus_x[i-1]-a_p_g[i-1][k]*B_plus_x[i-1]-
                        a_p_g[i-2][k]*A_plus_x[i-2] ); 
	}
        if(i==25){
                A_plus_x[25-2]=( 1.0/( r_ij[25][25-1]*(1.0+u_norm_x[25-1][k]*u_norm_x[25-2][k]+
                        u_norm_y[25-1][k]*u_norm_y[25-2][k]+u_norm_z[25-1][k]*u_norm_z[25-2][k]) ) )*
                        ( f_norm_x[25-1]*(v_norm_x[25-1]*u_norm_x[25-2][k]+v_norm_y[25-1]*u_norm_y[25-2][k]
                        +v_norm_z[25-1]*u_norm_z[25-2][k])-v_norm_x[25-1]*(f_norm_x[25-1]*u_norm_x[25-2][k]
                        +f_norm_y[25-1]*u_norm_y[25-2][k]+f_norm_z[25-1]*u_norm_z[25-2][k]) );
                B_plus_x[25-1]=( 1.0/( r_ij[25][25-1]*(1.0+u_norm_x[25][k]*u_norm_x[25-1][k]+
                        u_norm_y[25][k]*u_norm_y[25-1][k]+u_norm_z[25][k]*u_norm_z[25-1][k]) ) )*
                        ( -f_norm_x[25-1]*(u_norm_x[25][k]*v_norm_x[25-1]+u_norm_y[25][k]*v_norm_y[25-1]+
                        u_norm_z[25][k]*v_norm_z[25-1])+v_norm_x[25-1]*(u_norm_x[25][k]*f_norm_x[25-1]+
                        u_norm_y[25][k]*f_norm_y[25-1]+u_norm_z[25][k]*f_norm_z[25-1]) );

                f_t_x=( kB*temp/pow(ZI,2.0) )*( 
                        -a_p_g[25-1][k]*B_plus_x[25-1]-
                        a_p_g[25-2][k]*A_plus_x[25-2] );        
	}
}
if(i>25){
	if(i==26){
		B_plus_x[26]=( 1.0/( r_ij[26][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*
			( -f_norm_x[26]*(u_norm_x[27][k]*v_norm_x[26]+u_norm_y[27][k]*v_norm_y[26]+
			u_norm_z[27][k]*v_norm_z[26])+v_norm_x[26]*(u_norm_x[27][k]*f_norm_x[26]+
			u_norm_y[27][k]*f_norm_y[26]+u_norm_z[27][k]*f_norm_z[26]) );

                f_t_x=( kB*temp/pow(ZI,2.0) )*a_p_g[26][k]*B_plus_x[26];
        }
        if(i==27){
     		A_plus_x[26]=( 1.0/( r_ij[28][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*
                        ( f_norm_x[27]*(v_norm_x[27]*u_norm_x[26][k]+v_norm_y[27]*u_norm_y[26][k]
                        +v_norm_z[27]*u_norm_z[26][k])-v_norm_x[27]*(f_norm_x[27]*u_norm_x[26][k]
                        +f_norm_y[27]*u_norm_y[26][k]+f_norm_z[27]*u_norm_z[26][k]) );
                B_plus_x[27]=( 1.0/( r_ij[28][27]*(1.0+u_norm_x[28][k]*u_norm_x[27][k]+
                        u_norm_y[28][k]*u_norm_y[27][k]+u_norm_z[28][k]*u_norm_z[27][k]) ) )*
                        ( -f_norm_x[27]*(u_norm_x[28][k]*v_norm_x[27]+u_norm_y[28][k]*v_norm_y[27]+
                        u_norm_z[28][k]*v_norm_z[27])+v_norm_x[27]*(u_norm_x[28][k]*f_norm_x[27]+
                        u_norm_y[28][k]*f_norm_y[27]+u_norm_z[28][k]*f_norm_z[27]) );
                B_plus_x[26]=( 1.0/( r_ij[26][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+      
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*       
                        ( -f_norm_x[26]*(u_norm_x[27][k]*v_norm_x[26]+u_norm_y[27][k]*v_norm_y[26]+  
                        u_norm_z[27][k]*v_norm_z[26])+v_norm_x[26]*(u_norm_x[27][k]*f_norm_x[26]+    
                        u_norm_y[27][k]*f_norm_y[26]+u_norm_z[27][k]*f_norm_z[26]) );       

                f_t_x=( kB*temp/pow(ZI,2.0) )*( a_p_g[27][k]*B_plus_x[27]+
                        a_p_g[26][k]*A_plus_x[26]-a_p_g[26][k]*B_plus_x[26] );
	}
        if(i!=26 && i!=27 && i!=N){
       		A_plus_x[i-1]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( f_norm_x[i]*(v_norm_x[i]*u_norm_x[i-1][k]+v_norm_y[i]*u_norm_y[i-1][k]
                        +v_norm_z[i]*u_norm_z[i-1][k])-v_norm_x[i]*(f_norm_x[i]*u_norm_x[i-1][k]
                        +f_norm_y[i]*u_norm_y[i-1][k]+f_norm_z[i]*u_norm_z[i-1][k]) );
                A_plus_x[i-2]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i-1][k]*u_norm_x[i-2][k]+
                        u_norm_y[i-1][k]*u_norm_y[i-2][k]+u_norm_z[i-1][k]*u_norm_z[i-2][k]) ) )*
                        ( f_norm_x[i-1]*(v_norm_x[i-1]*u_norm_x[i-2][k]+v_norm_y[i-1]*u_norm_y[i-2][k]
                        +v_norm_z[i-1]*u_norm_z[i-2][k])-v_norm_x[i-1]*(f_norm_x[i-1]*u_norm_x[i-2][k]
                        +f_norm_y[i-1]*u_norm_y[i-2][k]+f_norm_z[i-1]*u_norm_z[i-2][k]) );
                B_plus_x[i]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i+1][k]*u_norm_x[i][k]+
                        u_norm_y[i+1][k]*u_norm_y[i][k]+u_norm_z[i+1][k]*u_norm_z[i][k]) ) )*
                        ( -f_norm_x[i]*(u_norm_x[i+1][k]*v_norm_x[i]+u_norm_y[i+1][k]*v_norm_y[i]+
                        u_norm_z[i+1][k]*v_norm_z[i])+v_norm_x[i]*(u_norm_x[i+1][k]*f_norm_x[i]+
                        u_norm_y[i+1][k]*f_norm_y[i]+u_norm_z[i+1][k]*f_norm_z[i]) );
                B_plus_x[i-1]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( -f_norm_x[i-1]*(u_norm_x[i][k]*v_norm_x[i-1]+u_norm_y[i][k]*v_norm_y[i-1]+
                        u_norm_z[i][k]*v_norm_z[i-1])+v_norm_x[i-1]*(u_norm_x[i][k]*f_norm_x[i-1]+
                        u_norm_y[i][k]*f_norm_y[i-1]+u_norm_z[i][k]*f_norm_z[i-1]) );

                f_t_x=( kB*temp/pow(ZI,2.0) )*( a_p_g[i][k]*B_plus_x[i]+
                        a_p_g[i-1][k]*A_plus_x[i-1]-a_p_g[i-1][k]*B_plus_x[i-1]-
                        a_p_g[i-2][k]*A_plus_x[i-2] ); 
	}
        if(i==N){
                A_plus_x[N-2]=( 1.0/( r_ij[N][N-1]*(1.0+u_norm_x[N-1][k]*u_norm_x[N-2][k]+
                        u_norm_y[N-1][k]*u_norm_y[N-2][k]+u_norm_z[N-1][k]*u_norm_z[N-2][k]) ) )*
                        ( f_norm_x[N-1]*(v_norm_x[N-1]*u_norm_x[N-2][k]+v_norm_y[N-1]*u_norm_y[N-2][k]
                        +v_norm_z[N-1]*u_norm_z[N-2][k])-v_norm_x[N-1]*(f_norm_x[N-1]*u_norm_x[N-2][k]
                        +f_norm_y[N-1]*u_norm_y[N-2][k]+f_norm_z[N-1]*u_norm_z[N-2][k]) );
                B_plus_x[N-1]=( 1.0/( r_ij[N][N-1]*(1.0+u_norm_x[N][k]*u_norm_x[N-1][k]+
                        u_norm_y[N][k]*u_norm_y[N-1][k]+u_norm_z[N][k]*u_norm_z[N-1][k]) ) )*
                        ( -f_norm_x[N-1]*(u_norm_x[N][k]*v_norm_x[N-1]+u_norm_y[N][k]*v_norm_y[N-1]+
                        u_norm_z[N][k]*v_norm_z[N-1])+v_norm_x[N-1]*(u_norm_x[N][k]*f_norm_x[N-1]+
                        u_norm_y[N][k]*f_norm_y[N-1]+u_norm_z[N][k]*f_norm_z[N-1]) );

                f_t_x=( kB*temp/pow(ZI,2.0) )*( 
                        -a_p_g[N-1][k]*B_plus_x[N-1]-
                        a_p_g[N-2][k]*A_plus_x[N-2] );        
	}
}
	free_dvector(A_plus_x,1,N);
	free_dvector(B_plus_x,1,N);

        return f_t_x;
}
double F_t_y_sub(int i, int k, double **r_ij, double **u_norm_x, double **u_norm_y,
        double **u_norm_z, double *f_norm_x, double *f_norm_y, double *f_norm_z,
        double *v_norm_x, double *v_norm_y, double *v_norm_z, double **a_p_g)
{
	double f_t_y;
        double *A_plus_y,*B_plus_y;

        A_plus_y=dvector(1,N),B_plus_y=dvector(1,N);

if(i<=25){
	if(i==1){
		B_plus_y[1]=( 1.0/( r_ij[1][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*
			( -f_norm_y[1]*(u_norm_x[2][k]*v_norm_x[1]+u_norm_y[2][k]*v_norm_y[1]+
			u_norm_z[2][k]*v_norm_z[1])+v_norm_y[1]*(u_norm_x[2][k]*f_norm_x[1]+
			u_norm_y[2][k]*f_norm_y[1]+u_norm_z[2][k]*f_norm_z[1]) );

                f_t_y=( kB*temp/pow(ZI,2.0) )*a_p_g[1][k]*B_plus_y[1];
        }
        if(i==2){
     		A_plus_y[1]=( 1.0/( r_ij[3][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*
                        ( f_norm_y[2]*(v_norm_x[2]*u_norm_x[1][k]+v_norm_y[2]*u_norm_y[1][k]
                        +v_norm_z[2]*u_norm_z[1][k])-v_norm_y[2]*(f_norm_x[2]*u_norm_x[1][k]
                        +f_norm_y[2]*u_norm_y[1][k]+f_norm_z[2]*u_norm_z[1][k]) );
                B_plus_y[2]=( 1.0/( r_ij[3][2]*(1.0+u_norm_x[3][k]*u_norm_x[2][k]+
                        u_norm_y[3][k]*u_norm_y[2][k]+u_norm_z[3][k]*u_norm_z[2][k]) ) )*
                        ( -f_norm_y[2]*(u_norm_x[3][k]*v_norm_x[2]+u_norm_y[3][k]*v_norm_y[2]+
                        u_norm_z[3][k]*v_norm_z[2])+v_norm_y[2]*(u_norm_x[3][k]*f_norm_x[2]+
                        u_norm_y[3][k]*f_norm_y[2]+u_norm_z[3][k]*f_norm_z[2]) );
                B_plus_y[1]=( 1.0/( r_ij[1][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+      
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*       
                        ( -f_norm_y[1]*(u_norm_x[2][k]*v_norm_x[1]+u_norm_y[2][k]*v_norm_y[1]+  
                        u_norm_z[2][k]*v_norm_z[1])+v_norm_y[1]*(u_norm_x[2][k]*f_norm_x[1]+    
                        u_norm_y[2][k]*f_norm_y[1]+u_norm_z[2][k]*f_norm_z[1]) );       

                f_t_y=( kB*temp/pow(ZI,2.0) )*( a_p_g[2][k]*B_plus_y[2]+
                        a_p_g[1][k]*A_plus_y[1]-a_p_g[1][k]*B_plus_y[1] );
	}
        if(i!=1 && i!=2 && i!=25){
       		A_plus_y[i-1]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( f_norm_y[i]*(v_norm_x[i]*u_norm_x[i-1][k]+v_norm_y[i]*u_norm_y[i-1][k]
                        +v_norm_z[i]*u_norm_z[i-1][k])-v_norm_y[i]*(f_norm_x[i]*u_norm_x[i-1][k]
                        +f_norm_y[i]*u_norm_y[i-1][k]+f_norm_z[i]*u_norm_z[i-1][k]) );
                A_plus_y[i-2]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i-1][k]*u_norm_x[i-2][k]+
                        u_norm_y[i-1][k]*u_norm_y[i-2][k]+u_norm_z[i-1][k]*u_norm_z[i-2][k]) ) )*
                        ( f_norm_y[i-1]*(v_norm_x[i-1]*u_norm_x[i-2][k]+v_norm_y[i-1]*u_norm_y[i-2][k]
                        +v_norm_z[i-1]*u_norm_z[i-2][k])-v_norm_y[i-1]*(f_norm_x[i-1]*u_norm_x[i-2][k]
                        +f_norm_y[i-1]*u_norm_y[i-2][k]+f_norm_z[i-1]*u_norm_z[i-2][k]) );
                B_plus_y[i]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i+1][k]*u_norm_x[i][k]+
                        u_norm_y[i+1][k]*u_norm_y[i][k]+u_norm_z[i+1][k]*u_norm_z[i][k]) ) )*
                        ( -f_norm_y[i]*(u_norm_x[i+1][k]*v_norm_x[i]+u_norm_y[i+1][k]*v_norm_y[i]+
                        u_norm_z[i+1][k]*v_norm_z[i])+v_norm_y[i]*(u_norm_x[i+1][k]*f_norm_x[i]+
                        u_norm_y[i+1][k]*f_norm_y[i]+u_norm_z[i+1][k]*f_norm_z[i]) );
                B_plus_y[i-1]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( -f_norm_y[i-1]*(u_norm_x[i][k]*v_norm_x[i-1]+u_norm_y[i][k]*v_norm_y[i-1]+
                        u_norm_z[i][k]*v_norm_z[i-1])+v_norm_y[i-1]*(u_norm_x[i][k]*f_norm_x[i-1]+
                        u_norm_y[i][k]*f_norm_y[i-1]+u_norm_z[i][k]*f_norm_z[i-1]) );

                f_t_y=( kB*temp/pow(ZI,2.0) )*( a_p_g[i][k]*B_plus_y[i]+
                        a_p_g[i-1][k]*A_plus_y[i-1]-a_p_g[i-1][k]*B_plus_y[i-1]-
                        a_p_g[i-2][k]*A_plus_y[i-2] ); 
	}
        if(i==25){
                A_plus_y[25-2]=( 1.0/( r_ij[25][25-1]*(1.0+u_norm_x[25-1][k]*u_norm_x[25-2][k]+
                        u_norm_y[25-1][k]*u_norm_y[25-2][k]+u_norm_z[25-1][k]*u_norm_z[25-2][k]) ) )*
                        ( f_norm_y[25-1]*(v_norm_x[25-1]*u_norm_x[25-2][k]+v_norm_y[25-1]*u_norm_y[25-2][k]
                        +v_norm_z[25-1]*u_norm_z[25-2][k])-v_norm_y[25-1]*(f_norm_x[25-1]*u_norm_x[25-2][k]
                        +f_norm_y[25-1]*u_norm_y[25-2][k]+f_norm_z[25-1]*u_norm_z[25-2][k]) );
                B_plus_y[25-1]=( 1.0/( r_ij[25][25-1]*(1.0+u_norm_x[25][k]*u_norm_x[25-1][k]+
                        u_norm_y[25][k]*u_norm_y[25-1][k]+u_norm_z[25][k]*u_norm_z[25-1][k]) ) )*
                        ( -f_norm_y[25-1]*(u_norm_x[25][k]*v_norm_x[25-1]+u_norm_y[25][k]*v_norm_y[25-1]+
                        u_norm_z[25][k]*v_norm_z[25-1])+v_norm_y[25-1]*(u_norm_x[25][k]*f_norm_x[25-1]+
                        u_norm_y[25][k]*f_norm_y[25-1]+u_norm_z[25][k]*f_norm_z[25-1]) );

                f_t_y=( kB*temp/pow(ZI,2.0) )*( 
                        -a_p_g[25-1][k]*B_plus_y[25-1]-
                        a_p_g[25-2][k]*A_plus_y[25-2] );        
	}
}
if(i>25){
	if(i==26){
		B_plus_y[26]=( 1.0/( r_ij[26][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*
			( -f_norm_y[26]*(u_norm_x[27][k]*v_norm_x[26]+u_norm_y[27][k]*v_norm_y[26]+
			u_norm_z[27][k]*v_norm_z[26])+v_norm_y[26]*(u_norm_x[27][k]*f_norm_x[26]+
			u_norm_y[27][k]*f_norm_y[26]+u_norm_z[27][k]*f_norm_z[26]) );

                f_t_y=( kB*temp/pow(ZI,2.0) )*a_p_g[26][k]*B_plus_y[26];
        }
        if(i==27){
     		A_plus_y[26]=( 1.0/( r_ij[28][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*
                        ( f_norm_y[27]*(v_norm_x[27]*u_norm_x[26][k]+v_norm_y[27]*u_norm_y[26][k]
                        +v_norm_z[27]*u_norm_z[26][k])-v_norm_y[27]*(f_norm_x[27]*u_norm_x[26][k]
                        +f_norm_y[27]*u_norm_y[26][k]+f_norm_z[27]*u_norm_z[26][k]) );
                B_plus_y[27]=( 1.0/( r_ij[28][27]*(1.0+u_norm_x[28][k]*u_norm_x[27][k]+
                        u_norm_y[28][k]*u_norm_y[27][k]+u_norm_z[28][k]*u_norm_z[27][k]) ) )*
                        ( -f_norm_y[27]*(u_norm_x[28][k]*v_norm_x[27]+u_norm_y[28][k]*v_norm_y[27]+
                        u_norm_z[28][k]*v_norm_z[27])+v_norm_y[27]*(u_norm_x[28][k]*f_norm_x[27]+
                        u_norm_y[28][k]*f_norm_y[27]+u_norm_z[28][k]*f_norm_z[27]) );
                B_plus_y[26]=( 1.0/( r_ij[26][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+      
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*       
                        ( -f_norm_y[26]*(u_norm_x[27][k]*v_norm_x[26]+u_norm_y[27][k]*v_norm_y[26]+  
                        u_norm_z[27][k]*v_norm_z[26])+v_norm_y[26]*(u_norm_x[27][k]*f_norm_x[26]+    
                        u_norm_y[27][k]*f_norm_y[26]+u_norm_z[27][k]*f_norm_z[26]) );       

                f_t_y=( kB*temp/pow(ZI,2.0) )*( a_p_g[27][k]*B_plus_y[27]+
                        a_p_g[26][k]*A_plus_y[26]-a_p_g[26][k]*B_plus_y[26] );
	}
        if(i!=26 && i!=27 && i!=N){
       		A_plus_y[i-1]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( f_norm_y[i]*(v_norm_x[i]*u_norm_x[i-1][k]+v_norm_y[i]*u_norm_y[i-1][k]
                        +v_norm_z[i]*u_norm_z[i-1][k])-v_norm_y[i]*(f_norm_x[i]*u_norm_x[i-1][k]
                        +f_norm_y[i]*u_norm_y[i-1][k]+f_norm_z[i]*u_norm_z[i-1][k]) );
                A_plus_y[i-2]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i-1][k]*u_norm_x[i-2][k]+
                        u_norm_y[i-1][k]*u_norm_y[i-2][k]+u_norm_z[i-1][k]*u_norm_z[i-2][k]) ) )*
                        ( f_norm_y[i-1]*(v_norm_x[i-1]*u_norm_x[i-2][k]+v_norm_y[i-1]*u_norm_y[i-2][k]
                        +v_norm_z[i-1]*u_norm_z[i-2][k])-v_norm_y[i-1]*(f_norm_x[i-1]*u_norm_x[i-2][k]
                        +f_norm_y[i-1]*u_norm_y[i-2][k]+f_norm_z[i-1]*u_norm_z[i-2][k]) );
                B_plus_y[i]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i+1][k]*u_norm_x[i][k]+
                        u_norm_y[i+1][k]*u_norm_y[i][k]+u_norm_z[i+1][k]*u_norm_z[i][k]) ) )*
                        ( -f_norm_y[i]*(u_norm_x[i+1][k]*v_norm_x[i]+u_norm_y[i+1][k]*v_norm_y[i]+
                        u_norm_z[i+1][k]*v_norm_z[i])+v_norm_y[i]*(u_norm_x[i+1][k]*f_norm_x[i]+
                        u_norm_y[i+1][k]*f_norm_y[i]+u_norm_z[i+1][k]*f_norm_z[i]) );
                B_plus_y[i-1]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( -f_norm_y[i-1]*(u_norm_x[i][k]*v_norm_x[i-1]+u_norm_y[i][k]*v_norm_y[i-1]+
                        u_norm_z[i][k]*v_norm_z[i-1])+v_norm_y[i-1]*(u_norm_x[i][k]*f_norm_x[i-1]+
                        u_norm_y[i][k]*f_norm_y[i-1]+u_norm_z[i][k]*f_norm_z[i-1]) );

                f_t_y=( kB*temp/pow(ZI,2.0) )*( a_p_g[i][k]*B_plus_y[i]+
                        a_p_g[i-1][k]*A_plus_y[i-1]-a_p_g[i-1][k]*B_plus_y[i-1]-
                        a_p_g[i-2][k]*A_plus_y[i-2] ); 
	}        
	if(i==N){
                A_plus_y[N-2]=( 1.0/( r_ij[N][N-1]*(1.0+u_norm_x[N-1][k]*u_norm_x[N-2][k]+
                        u_norm_y[N-1][k]*u_norm_y[N-2][k]+u_norm_z[N-1][k]*u_norm_z[N-2][k]) ) )*
                        ( f_norm_y[N-1]*(v_norm_x[N-1]*u_norm_x[N-2][k]+v_norm_y[N-1]*u_norm_y[N-2][k]
                        +v_norm_z[N-1]*u_norm_z[N-2][k])-v_norm_y[N-1]*(f_norm_x[N-1]*u_norm_x[N-2][k]
                        +f_norm_y[N-1]*u_norm_y[N-2][k]+f_norm_z[N-1]*u_norm_z[N-2][k]) );
                B_plus_y[N-1]=( 1.0/( r_ij[N][N-1]*(1.0+u_norm_x[N][k]*u_norm_x[N-1][k]+
                        u_norm_y[N][k]*u_norm_y[N-1][k]+u_norm_z[N][k]*u_norm_z[N-1][k]) ) )*
                        ( -f_norm_y[N-1]*(u_norm_x[N][k]*v_norm_x[N-1]+u_norm_y[N][k]*v_norm_y[N-1]+
                        u_norm_z[N][k]*v_norm_z[N-1])+v_norm_y[N-1]*(u_norm_x[N][k]*f_norm_x[N-1]+
                        u_norm_y[N][k]*f_norm_y[N-1]+u_norm_z[N][k]*f_norm_z[N-1]) );

                f_t_y=( kB*temp/pow(ZI,2.0) )*( 
                        -a_p_g[N-1][k]*B_plus_y[N-1]-
                        a_p_g[N-2][k]*A_plus_y[N-2] );        
	}
}

        free_dvector(A_plus_y,1,N);
        free_dvector(B_plus_y,1,N);

        return f_t_y;
}
double F_t_z_sub(int i, int k, double **r_ij, double **u_norm_x, double **u_norm_y,
        double **u_norm_z, double *f_norm_x, double *f_norm_y, double *f_norm_z,
        double *v_norm_x, double *v_norm_y, double *v_norm_z, double **a_p_g) 
{
        double f_t_z;
        double *A_plus_z,*B_plus_z;

        A_plus_z=dvector(1,N),B_plus_z=dvector(1,N);

if(i<=25){
	if(i==1){
		B_plus_z[1]=( 1.0/( r_ij[1][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*
			( -f_norm_z[1]*(u_norm_x[2][k]*v_norm_x[1]+u_norm_y[2][k]*v_norm_y[1]+
			u_norm_z[2][k]*v_norm_z[1])+v_norm_z[1]*(u_norm_x[2][k]*f_norm_x[1]+
			u_norm_y[2][k]*f_norm_y[1]+u_norm_z[2][k]*f_norm_z[1]) );

                f_t_z=( kB*temp/pow(ZI,2.0) )*a_p_g[1][k]*B_plus_z[1];
        }
        if(i==2){
     		A_plus_z[1]=( 1.0/( r_ij[3][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*
                        ( f_norm_z[2]*(v_norm_x[2]*u_norm_x[1][k]+v_norm_y[2]*u_norm_y[1][k]
                        +v_norm_z[2]*u_norm_z[1][k])-v_norm_z[2]*(f_norm_x[2]*u_norm_x[1][k]
                        +f_norm_y[2]*u_norm_y[1][k]+f_norm_z[2]*u_norm_z[1][k]) );
                B_plus_z[2]=( 1.0/( r_ij[3][2]*(1.0+u_norm_x[3][k]*u_norm_x[2][k]+
                        u_norm_y[3][k]*u_norm_y[2][k]+u_norm_z[3][k]*u_norm_z[2][k]) ) )*
                        ( -f_norm_z[2]*(u_norm_x[3][k]*v_norm_x[2]+u_norm_y[3][k]*v_norm_y[2]+
                        u_norm_z[3][k]*v_norm_z[2])+v_norm_z[2]*(u_norm_x[3][k]*f_norm_x[2]+
                        u_norm_y[3][k]*f_norm_y[2]+u_norm_z[3][k]*f_norm_z[2]) );
                B_plus_z[1]=( 1.0/( r_ij[1][2]*(1.0+u_norm_x[2][k]*u_norm_x[1][k]+      
                        u_norm_y[2][k]*u_norm_y[1][k]+u_norm_z[2][k]*u_norm_z[1][k]) ) )*       
                        ( -f_norm_z[1]*(u_norm_x[2][k]*v_norm_x[1]+u_norm_y[2][k]*v_norm_y[1]+  
                        u_norm_z[2][k]*v_norm_z[1])+v_norm_z[1]*(u_norm_x[2][k]*f_norm_x[1]+    
                        u_norm_y[2][k]*f_norm_y[1]+u_norm_z[2][k]*f_norm_z[1]) );       

                f_t_z=( kB*temp/pow(ZI,2.0) )*( a_p_g[2][k]*B_plus_z[2]+
                        a_p_g[1][k]*A_plus_z[1]-a_p_g[1][k]*B_plus_z[1] );
	}
        if(i!=1 && i!=2 && i!=25){
       		A_plus_z[i-1]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( f_norm_z[i]*(v_norm_x[i]*u_norm_x[i-1][k]+v_norm_y[i]*u_norm_y[i-1][k]
                        +v_norm_z[i]*u_norm_z[i-1][k])-v_norm_z[i]*(f_norm_x[i]*u_norm_x[i-1][k]
                        +f_norm_y[i]*u_norm_y[i-1][k]+f_norm_z[i]*u_norm_z[i-1][k]) );
                A_plus_z[i-2]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i-1][k]*u_norm_x[i-2][k]+
                        u_norm_y[i-1][k]*u_norm_y[i-2][k]+u_norm_z[i-1][k]*u_norm_z[i-2][k]) ) )*
                        ( f_norm_z[i-1]*(v_norm_x[i-1]*u_norm_x[i-2][k]+v_norm_y[i-1]*u_norm_y[i-2][k]
                        +v_norm_z[i-1]*u_norm_z[i-2][k])-v_norm_z[i-1]*(f_norm_x[i-1]*u_norm_x[i-2][k]
                        +f_norm_y[i-1]*u_norm_y[i-2][k]+f_norm_z[i-1]*u_norm_z[i-2][k]) );
                B_plus_z[i]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i+1][k]*u_norm_x[i][k]+
                        u_norm_y[i+1][k]*u_norm_y[i][k]+u_norm_z[i+1][k]*u_norm_z[i][k]) ) )*
                        ( -f_norm_z[i]*(u_norm_x[i+1][k]*v_norm_x[i]+u_norm_y[i+1][k]*v_norm_y[i]+
                        u_norm_z[i+1][k]*v_norm_z[i])+v_norm_z[i]*(u_norm_x[i+1][k]*f_norm_x[i]+
                        u_norm_y[i+1][k]*f_norm_y[i]+u_norm_z[i+1][k]*f_norm_z[i]) );
                B_plus_z[i-1]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( -f_norm_z[i-1]*(u_norm_x[i][k]*v_norm_x[i-1]+u_norm_y[i][k]*v_norm_y[i-1]+
                        u_norm_z[i][k]*v_norm_z[i-1])+v_norm_z[i-1]*(u_norm_x[i][k]*f_norm_x[i-1]+
                        u_norm_y[i][k]*f_norm_y[i-1]+u_norm_z[i][k]*f_norm_z[i-1]) );

                f_t_z=( kB*temp/pow(ZI,2.0) )*( a_p_g[i][k]*B_plus_z[i]+
                        a_p_g[i-1][k]*A_plus_z[i-1]-a_p_g[i-1][k]*B_plus_z[i-1]-
                        a_p_g[i-2][k]*A_plus_z[i-2] ); 
	}
        if(i==25){
                A_plus_z[25-2]=( 1.0/( r_ij[25][25-1]*(1.0+u_norm_x[25-1][k]*u_norm_x[25-2][k]+
                        u_norm_y[25-1][k]*u_norm_y[25-2][k]+u_norm_z[25-1][k]*u_norm_z[25-2][k]) ) )*
                        ( f_norm_z[25-1]*(v_norm_x[25-1]*u_norm_x[25-2][k]+v_norm_y[25-1]*u_norm_y[25-2][k]
                        +v_norm_z[25-1]*u_norm_z[25-2][k])-v_norm_z[25-1]*(f_norm_x[25-1]*u_norm_x[25-2][k]
                        +f_norm_y[25-1]*u_norm_y[25-2][k]+f_norm_z[25-1]*u_norm_z[25-2][k]) );
                B_plus_z[25-1]=( 1.0/( r_ij[25][25-1]*(1.0+u_norm_x[25][k]*u_norm_x[25-1][k]+
                        u_norm_y[25][k]*u_norm_y[25-1][k]+u_norm_z[25][k]*u_norm_z[25-1][k]) ) )*
                        ( -f_norm_z[25-1]*(u_norm_x[25][k]*v_norm_x[25-1]+u_norm_y[25][k]*v_norm_y[25-1]+
                        u_norm_z[25][k]*v_norm_z[25-1])+v_norm_z[25-1]*(u_norm_x[25][k]*f_norm_x[25-1]+
                        u_norm_y[25][k]*f_norm_y[25-1]+u_norm_z[25][k]*f_norm_z[25-1]) );

                f_t_z=( kB*temp/pow(ZI,2.0) )*( 
                        -a_p_g[25-1][k]*B_plus_z[25-1]-
                        a_p_g[25-2][k]*A_plus_z[25-2] );        
	}
}
if(i>25){
	if(i==26){
		B_plus_z[26]=( 1.0/( r_ij[26][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*
			( -f_norm_z[26]*(u_norm_x[27][k]*v_norm_x[26]+u_norm_y[27][k]*v_norm_y[26]+
			u_norm_z[27][k]*v_norm_z[26])+v_norm_z[26]*(u_norm_x[27][k]*f_norm_x[26]+
			u_norm_y[27][k]*f_norm_y[26]+u_norm_z[27][k]*f_norm_z[26]) );

                f_t_z=( kB*temp/pow(ZI,2.0) )*a_p_g[26][k]*B_plus_z[26];
        }
        if(i==27){
     		A_plus_z[26]=( 1.0/( r_ij[28][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*
                        ( f_norm_z[27]*(v_norm_x[27]*u_norm_x[26][k]+v_norm_y[27]*u_norm_y[26][k]
                        +v_norm_z[27]*u_norm_z[26][k])-v_norm_z[27]*(f_norm_x[27]*u_norm_x[26][k]
                        +f_norm_y[27]*u_norm_y[26][k]+f_norm_z[27]*u_norm_z[26][k]) );
                B_plus_z[27]=( 1.0/( r_ij[28][27]*(1.0+u_norm_x[28][k]*u_norm_x[27][k]+
                        u_norm_y[28][k]*u_norm_y[27][k]+u_norm_z[28][k]*u_norm_z[27][k]) ) )*
                        ( -f_norm_z[27]*(u_norm_x[28][k]*v_norm_x[27]+u_norm_y[28][k]*v_norm_y[27]+
                        u_norm_z[28][k]*v_norm_z[27])+v_norm_z[27]*(u_norm_x[28][k]*f_norm_x[27]+
                        u_norm_y[28][k]*f_norm_y[27]+u_norm_z[28][k]*f_norm_z[27]) );
                B_plus_z[26]=( 1.0/( r_ij[26][27]*(1.0+u_norm_x[27][k]*u_norm_x[26][k]+      
                        u_norm_y[27][k]*u_norm_y[26][k]+u_norm_z[27][k]*u_norm_z[26][k]) ) )*       
                        ( -f_norm_z[26]*(u_norm_x[27][k]*v_norm_x[26]+u_norm_y[27][k]*v_norm_y[26]+  
                        u_norm_z[27][k]*v_norm_z[26])+v_norm_z[26]*(u_norm_x[27][k]*f_norm_x[26]+    
                        u_norm_y[27][k]*f_norm_y[26]+u_norm_z[27][k]*f_norm_z[26]) );       

                f_t_z=( kB*temp/pow(ZI,2.0) )*( a_p_g[27][k]*B_plus_z[27]+
                        a_p_g[26][k]*A_plus_z[26]-a_p_g[26][k]*B_plus_z[26] );
	}
        if(i!=26 && i!=27 && i!=N){
       		A_plus_z[i-1]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( f_norm_z[i]*(v_norm_x[i]*u_norm_x[i-1][k]+v_norm_y[i]*u_norm_y[i-1][k]
                        +v_norm_z[i]*u_norm_z[i-1][k])-v_norm_z[i]*(f_norm_x[i]*u_norm_x[i-1][k]
                        +f_norm_y[i]*u_norm_y[i-1][k]+f_norm_z[i]*u_norm_z[i-1][k]) );
                A_plus_z[i-2]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i-1][k]*u_norm_x[i-2][k]+
                        u_norm_y[i-1][k]*u_norm_y[i-2][k]+u_norm_z[i-1][k]*u_norm_z[i-2][k]) ) )*
                        ( f_norm_z[i-1]*(v_norm_x[i-1]*u_norm_x[i-2][k]+v_norm_y[i-1]*u_norm_y[i-2][k]
                        +v_norm_z[i-1]*u_norm_z[i-2][k])-v_norm_z[i-1]*(f_norm_x[i-1]*u_norm_x[i-2][k]
                        +f_norm_y[i-1]*u_norm_y[i-2][k]+f_norm_z[i-1]*u_norm_z[i-2][k]) );
                B_plus_z[i]=( 1.0/( r_ij[i][i+1]*(1.0+u_norm_x[i+1][k]*u_norm_x[i][k]+
                        u_norm_y[i+1][k]*u_norm_y[i][k]+u_norm_z[i+1][k]*u_norm_z[i][k]) ) )*
                        ( -f_norm_z[i]*(u_norm_x[i+1][k]*v_norm_x[i]+u_norm_y[i+1][k]*v_norm_y[i]+
                        u_norm_z[i+1][k]*v_norm_z[i])+v_norm_z[i]*(u_norm_x[i+1][k]*f_norm_x[i]+
                        u_norm_y[i+1][k]*f_norm_y[i]+u_norm_z[i+1][k]*f_norm_z[i]) );
                B_plus_z[i-1]=( 1.0/( r_ij[i][i-1]*(1.0+u_norm_x[i][k]*u_norm_x[i-1][k]+
                        u_norm_y[i][k]*u_norm_y[i-1][k]+u_norm_z[i][k]*u_norm_z[i-1][k]) ) )*
                        ( -f_norm_z[i-1]*(u_norm_x[i][k]*v_norm_x[i-1]+u_norm_y[i][k]*v_norm_y[i-1]+
                        u_norm_z[i][k]*v_norm_z[i-1])+v_norm_z[i-1]*(u_norm_x[i][k]*f_norm_x[i-1]+
                        u_norm_y[i][k]*f_norm_y[i-1]+u_norm_z[i][k]*f_norm_z[i-1]) );

                f_t_z=( kB*temp/pow(ZI,2.0) )*( a_p_g[i][k]*B_plus_z[i]+
                        a_p_g[i-1][k]*A_plus_z[i-1]-a_p_g[i-1][k]*B_plus_z[i-1]-
                        a_p_g[i-2][k]*A_plus_z[i-2] ); 
	}
        if(i==N){
                A_plus_z[N-2]=( 1.0/( r_ij[N][N-1]*(1.0+u_norm_x[N-1][k]*u_norm_x[N-2][k]+
                        u_norm_y[N-1][k]*u_norm_y[N-2][k]+u_norm_z[N-1][k]*u_norm_z[N-2][k]) ) )*
                        ( f_norm_z[N-1]*(v_norm_x[N-1]*u_norm_x[N-2][k]+v_norm_y[N-1]*u_norm_y[N-2][k]
                        +v_norm_z[N-1]*u_norm_z[N-2][k])-v_norm_z[N-1]*(f_norm_x[N-1]*u_norm_x[N-2][k]
                        +f_norm_y[N-1]*u_norm_y[N-2][k]+f_norm_z[N-1]*u_norm_z[N-2][k]) );
                B_plus_z[N-1]=( 1.0/( r_ij[N][N-1]*(1.0+u_norm_x[N][k]*u_norm_x[N-1][k]+
                        u_norm_y[N][k]*u_norm_y[N-1][k]+u_norm_z[N][k]*u_norm_z[N-1][k]) ) )*
                        ( -f_norm_z[N-1]*(u_norm_x[N][k]*v_norm_x[N-1]+u_norm_y[N][k]*v_norm_y[N-1]+
                        u_norm_z[N][k]*v_norm_z[N-1])+v_norm_z[N-1]*(u_norm_x[N][k]*f_norm_x[N-1]+
                        u_norm_y[N][k]*f_norm_y[N-1]+u_norm_z[N][k]*f_norm_z[N-1]) );

                f_t_z=( kB*temp/pow(ZI,2.0) )*( 
                        -a_p_g[N-1][k]*B_plus_z[N-1]-
                        a_p_g[N-2][k]*A_plus_z[N-2] );        
	}
}	
        free_dvector(A_plus_z,1,N);
        free_dvector(B_plus_z,1,N);
	
        return f_t_z;
}

double F_EV_x_sub(int i,int k,double **rx,double **r_ij)
{
        double L,f_EV_x=0.0;
        int l;

        for(l=1;l<=N;l++){
                if(l==i){
                        f_EV_x+=0.0;
                }
                else{
                        if(r_ij[i][l]<CUTOFF){
                                L=-( (6.0*EPSILON*kB*temp*pow(SIGMA,6.0))/pow(r_ij[i][l],7.0) )*
                                ( 1.0 - 2.0*pow(SIGMA,6.0)/pow(r_ij[i][l],6.0) );
                        }
                        else{
                                L=0.0;
                        }

                        f_EV_x+=L*(rx[i][k]-rx[l][k])/r_ij[i][l];
                }
        }

        return f_EV_x;
}

double F_EV_y_sub(int i,int k,double **ry,double **r_ij)
{
        double L,f_EV_y=0.0;
        int l;

        for(l=1;l<=N;l++){
                if(l==i){
                        f_EV_y+=0.0;
                }
                else{
                        if(r_ij[i][l]<CUTOFF){
                                L=-( (6.0*EPSILON*kB*temp*pow(SIGMA,6.0))/pow(r_ij[i][l],7.0) )*
                                ( 1.0 - 2.0*pow(SIGMA,6.0)/pow(r_ij[i][l],6.0) );
                        }
                        else{
                                L=0.0;
                        }

                        f_EV_y+=L*(ry[i][k]-ry[l][k])/r_ij[i][l];
                }
        }

        return f_EV_y;
}

double F_EV_z_sub(int i,int k,double **rz,double **r_ij)
{
        double L,f_EV_z=0.0;
        int l;

        for(l=1;l<=N;l++){
                if(l==i){
                        f_EV_z+=0.0;
                }
                else{
                        if(r_ij[i][l]<CUTOFF){
                                L=-( (6.0*EPSILON*kB*temp*pow(SIGMA,6.0))/pow(r_ij[i][l],7.0) )*
                                ( 1.0 - 2.0*pow(SIGMA,6.0)/pow(r_ij[i][l],6.0) );
                        }
                        else{
                                L=0.0;
                        }

                        f_EV_z+=L*(rz[i][k]-rz[l][k])/r_ij[i][l];
                }
        }

        return f_EV_z;
}

int main(void)
{
	int i,j,k;
	double *bz,*bx,**rx,**ry,**rz,**r_ij;
	double **u_norm_x,**u_norm_y,**u_norm_z;
	double *f_norm_x,*f_norm_y,*f_norm_z,*v_norm_x,*v_norm_y,*v_norm_z;
	double *fx,*fy,*fz,*vx,*vy,*vz;
	double *r_f_x,*r_f_y,*r_f_z;
	double **a_p_g,*beta,arg_beta,**phi,*tau;
	double *F_tot_x,*F_tot_y,*F_tot_z;
	double Wr_sum;/*,Wr1,Wr2,Tw1,Tw2,Tw_sum,Lk1,Lk2;*/
	/*double Wr_tot_sum=0.0,Tw_tot_sum=0.0,Lk_tot_sum=0.0;
	double dens1,dens2;*/
	long idum=(SEED);

	FILE *fptr1; char filename1[]="ml9_twin.out";
	FILE *fptr2; char filename2[]="ml9_twin.dgn";
	/*FILE *fptr2; char filename2[]="ml9_twin.Lk1";	
	FILE *fptr3; char filename3[]="ml9_twin.Lk2";
	FILE *fptr4; char filename4[]="ml9_twin.Lk_tot";
	FILE *fptr5; char filename5[]="ml9_twin.SWr1";
	FILE *fptr6; char filename6[]="ml9_twin.SWr2";
	FILE *fptr7; char filename7[]="ml9_twin.STw1";
 	FILE *fptr8; char filename8[]="ml9_twin.STw2";
	FILE *fptr9; char filename9[]="ml9_twin.3Wr1";
	FILE *fptr10; char filename10[]="ml9_twin.3Wr2";*/

	bz=dvector(1,N-1),bx=dvector(1,N-1);
	rx=dmatrix(1,N+2,1,2),ry=dmatrix(1,N+2,1,2),rz=dmatrix(1,N+2,1,2); 
	r_ij=dmatrix(1,N,1,N); 
	u_norm_x=dmatrix(1,N,1,2),u_norm_y=dmatrix(1,N,1,2),u_norm_z=dmatrix(1,N,1,2);
	f_norm_x=dvector(1,N+2),f_norm_y=dvector(1,N+2),f_norm_z=dvector(1,N+2);
	v_norm_x=dvector(1,N),v_norm_y=dvector(1,N),v_norm_z=dvector(1,N);
	fx=dvector(1,N),fy=dvector(1,N),fz=dvector(1,N);
	vx=dvector(1,N),vy=dvector(1,N),vz=dvector(1,N);
	r_f_x=dvector(1,N+2),r_f_y=dvector(1,N+2),r_f_z=dvector(1,N+2);
	a_p_g=dmatrix(1,N-1,1,2),beta=dvector(1,N-1),phi=dmatrix(1,N,1,2);
	tau=dvector(1,N);
	F_tot_x=dvector(1,N),F_tot_y=dvector(1,N),F_tot_z=dvector(1,N);

/* Initial Configuration */
	
	/* Initial beta[i] */
        for(i=1;i<=N/2;i++)
        {
                beta[i]=-2.0*Pi/(N/2);
        }
	for(i=26;i<=N;i++)
        {
                beta[i]=-2.0*Pi/(N/2);
        }

        /* Initial Positions */
        rx[1][1]=0.0, ry[1][1]=0.0, rz[1][1]=b0;
        for(i=2;i<=N/2;i++)
        {
                bz[i-1]=b0*cos((i-2.0)*beta[i]);
                bx[i-1]=b0*sin((i-2.0)*beta[i]);
                rx[i][1]=rx[i-1][1]+bx[i-1];
                ry[i][1]=0.0;
                rz[i][1]=rz[i-1][1]+bz[i-1];
        }
       	rx[26][1]=0.0, ry[26][1]=0.0, rz[26][1]=-b0;
        for(i=27;i<=N;i++)
        {
                bz[i-1]=-b0*cos((i-27)*beta[i]);
                bx[i-1]=-b0*sin((i-27)*beta[i]);
                rx[i][1]=rx[i-1][1]+bx[i-1];
                ry[i][1]=0.0;
                rz[i][1]=rz[i-1][1]+bz[i-1];
        }
	
	/* Initial Separations */
	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			r_ij[i][j]=sqrt(pow(rx[j][1]-rx[i][1],2.0)+
                                pow(ry[j][1]-ry[i][1],2.0)+
                                pow(rz[j][1]-rz[i][1],2.0));
                }
        }
	
	/* Initial BFCs */
		/* initial u_norm[i] */ 
	for(i=1;i<=24;i++){
		u_norm_x[i][1]=(rx[i+1][1]-rx[i][1])/r_ij[i+1][i];
            	u_norm_y[i][1]=(ry[i+1][1]-ry[i][1])/r_ij[i+1][i];
            	u_norm_z[i][1]=(rz[i+1][1]-rz[i][1])/r_ij[i+1][i];
                	
	}
	u_norm_x[25][1]=(rx[1][1]-rx[25][1])/r_ij[1][25];
        u_norm_y[25][1]=(ry[1][1]-ry[25][1])/r_ij[1][25];
        u_norm_z[25][1]=(rz[1][1]-rz[25][1])/r_ij[1][25];

	for(i=26;i<=N-1;i++){
		u_norm_x[i][1]=(rx[i+1][1]-rx[i][1])/r_ij[i+1][i];
            	u_norm_y[i][1]=(ry[i+1][1]-ry[i][1])/r_ij[i+1][i];
            	u_norm_z[i][1]=(rz[i+1][1]-rz[i][1])/r_ij[i+1][i];
                	
	}
	u_norm_x[50][1]=(rx[26][1]-rx[50][1])/r_ij[26][50];
        u_norm_y[50][1]=(ry[26][1]-ry[50][1])/r_ij[26][50];
        u_norm_z[50][1]=(rz[26][1]-rz[50][1])/r_ij[26][50];

        	/* initial f[i] */
	for(i=1;i<=25;i++){
        	f_norm_x[i]=0.0;     
            	f_norm_y[i]=1.0;
            	f_norm_z[i]=0.0;
	}	
      	for(i=26;i<=N;i++){
        	f_norm_x[i]=0.0;     
            	f_norm_y[i]=1.0;
            	f_norm_z[i]=0.0;
	}
                /* initial v_norm[i] */
	for(i=1;i<=25;i++){
                vx[i]=u_norm_y[i][1]*f_norm_z[i]-u_norm_z[i][1]*f_norm_y[i];
                vy[i]=u_norm_z[i][1]*f_norm_x[i]-u_norm_x[i][1]*f_norm_z[i];
                vz[i]=u_norm_x[i][1]*f_norm_y[i]-u_norm_y[i][1]*f_norm_x[i];

                v_norm_x[i]=vx[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                v_norm_y[i]=vy[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                v_norm_z[i]=vz[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
	}
        for(i=26;i<=N;i++){
                vx[i]=u_norm_y[i][1]*f_norm_z[i]-u_norm_z[i][1]*f_norm_y[i];
                vy[i]=u_norm_z[i][1]*f_norm_x[i]-u_norm_x[i][1]*f_norm_z[i];
                vz[i]=u_norm_x[i][1]*f_norm_y[i]-u_norm_y[i][1]*f_norm_x[i];

                v_norm_x[i]=vx[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                v_norm_y[i]=vy[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                v_norm_z[i]=vz[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
	}

	for(i=1;i<=N;i++){
	        r_f_x[i]=rx[i][1]+0.5*R_hyd*f_norm_x[i];
                r_f_y[i]=ry[i][1]+0.5*R_hyd*f_norm_y[i];
                r_f_z[i]=rz[i][1]+0.5*R_hyd*f_norm_z[i];
	}

	fptr2=fopen(filename2,"w");
	fclose(fptr2);

	/* Write Initial Config. to "*.out" */
	fptr1=fopen(filename1,"w");
       	fprintf(fptr1,"j=0\n"); 
	for(i=50;i>=26;i--){
                fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",rx[i][1]*1e10,
                        ry[i][1]*1e10,rz[i][1]*1e10);
                fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",r_f_x[i]*1e10,
                        r_f_y[i]*1e10,r_f_z[i]*1e10);
        }	
	for(i=1;i<=25;i++){ 
        	fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",rx[i][1]*1e10,
                        ry[i][1]*1e10,rz[i][1]*1e10);
                fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",r_f_x[i]*1e10,
                        r_f_y[i]*1e10,r_f_z[i]*1e10);
        }
        fclose(fptr1);

/* Initial Phi */
	for(i=1;i<=N;i++){
		phi[i][1]=0.0;
	}

/* Initial alpha[i]+gamma[i] (a_p_g[i][1]) */
	for(i=1;i<=24;i++){
                a_p_g[i][1]=asin( (v_norm_x[i]*f_norm_x[i+1]+v_norm_y[i]*f_norm_y[i+1]+
                        v_norm_z[i]*f_norm_z[i+1]-f_norm_x[i]*v_norm_x[i+1]-
                        f_norm_y[i]*v_norm_y[i+1]-f_norm_z[i]*v_norm_z[i+1])/
                        (2.0*cos(beta[i]/2.0)*cos(beta[i]/2.0)) );
	}
	for(i=26;i<=N-1;i++){
                a_p_g[i][1]=asin( (v_norm_x[i]*f_norm_x[i+1]+v_norm_y[i]*f_norm_y[i+1]+
                        v_norm_z[i]*f_norm_z[i+1]-f_norm_x[i]*v_norm_x[i+1]-
                        f_norm_y[i]*v_norm_y[i+1]-f_norm_z[i]*v_norm_z[i+1])/
                        (2.0*cos(beta[i]/2.0)*cos(beta[i]/2.0)) );
	}

/* Initial Lk=Tw+Wr */
			/* Chain 1 */
			/*Wr_sum=0.0;
			for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][1]-
                                rx[i][1],2.0)+pow(ry[k][1]-ry[i][1],2.0)+
                                pow(rz[k][1]-rz[i][1],2.0)),3.0) )*(
                                (rx[k+1][1]-rx[k][1])*(ry[i+1][1]-ry[i][1])*
                                (rz[k][1]-rz[i][1]) - (rx[k+1][1]-rx[k][1])*
                                (rz[i+1][1]-rz[i][1])*(ry[k][1]-ry[i][1]) -
                                (ry[k+1][1]-ry[k][1])*(rx[i+1][1]-rx[i][1])*
                                (rz[k][1]-rz[i][1]) + (ry[k+1][1]-ry[k][1])*
                                (rz[i+1][1]-rz[i][1])*(rx[k][1]-rx[i][1]) +
                                (rz[k+1][1]-rz[k][1])*(rx[i+1][1]-rx[i][1])*
                                (ry[k][1]-ry[i][1]) - (rz[k+1][1]-rz[k][1])*
                                (ry[i+1][1]-ry[i][1])*(rx[k][1]-rx[i][1]) );
                                        }
                                }
                        }
                        Wr1=Wr_sum/(4.0*Pi);*/
			
			/*Wr_sum=0.0;
			fptr5=fopen(filename5,"w");
                        fprintf(fptr5,"t(microseconds)=0.0\n");
                        fclose(fptr5);
			fptr5=fopen(filename5,"a");
                        for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][1]-
                                rx[i][1],2.0)+pow(ry[k][1]-ry[i][1],2.0)+
                                pow(rz[k][1]-rz[i][1],2.0)),3.0) )*(
                                (rx[k+1][1]-rx[k][1])*(ry[i+1][1]-ry[i][1])*
                                (rz[k][1]-rz[i][1]) - (rx[k+1][1]-rx[k][1])*
                                (rz[i+1][1]-rz[i][1])*(ry[k][1]-ry[i][1]) -
                                (ry[k+1][1]-ry[k][1])*(rx[i+1][1]-rx[i][1])*
                                (rz[k][1]-rz[i][1]) + (ry[k+1][1]-ry[k][1])*
                                (rz[i+1][1]-rz[i][1])*(rx[k][1]-rx[i][1]) +
                                (rz[k+1][1]-rz[k][1])*(rx[i+1][1]-rx[i][1])*
                                (ry[k][1]-ry[i][1]) - (rz[k+1][1]-rz[k][1])*
                                (ry[i+1][1]-ry[i][1])*(rx[k][1]-rx[i][1]) );
                                        }
                                }
                        fprintf(fptr5,"Bead=%d\tSWr1=%.5f\n",k,Wr_sum);
			Wr_sum=0.0;
                        }
			fclose(fptr5);   

                        Tw_sum=0.0;
			for(k=1;k<=24;k++)
                        {
                                Tw_sum+=a_p_g[k][1];
                        }
                        Tw1=Tw_sum/(2.0*Pi);

			fptr7=fopen(filename7,"w");
                        fprintf(fptr7,"t(microseconds)=0.0\n");
                        fclose(fptr7);
                        fptr7=fopen(filename7,"a");
			for(k=1;k<=24;k++)
                        {
                                fprintf(fptr7,"Bead=%d\tSTw1=%.5f\n",k,a_p_g[k][1]);
                        }
			fclose(fptr7);

                        Lk1=Tw1+Wr1;
			dens1=( (Lk1+21.4)-21.4 )/21.4;

                        fptr2=fopen(filename2,"w");
                        fprintf(fptr2,"j=0\n");
                        fprintf(fptr2,"Tw1=%.5f\tWr1=%.5f\tLk1=%.5f\tdens1=%.5f\n",
                        Tw1,Wr1,Lk1,dens1);
                        fclose(fptr2);*/
			
			/*Wr_sum=0.0;
                        fptr9=fopen(filename9,"w");
                        fprintf(fptr9,"t(microseconds)=0.0\n");
                        fclose(fptr9);
                        fptr9=fopen(filename9,"a");
                        for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum=0.0;
                                        }
                                        else
                                        {
                                Wr_sum=( 1.0/pow(sqrt(pow(rx[k][1]-
                                rx[i][1],2.0)+pow(ry[k][1]-ry[i][1],2.0)+
                                pow(rz[k][1]-rz[i][1],2.0)),3.0) )*(
                                (rx[k+1][1]-rx[k][1])*(ry[i+1][1]-ry[i][1])*
                                (rz[k][1]-rz[i][1]) - (rx[k+1][1]-rx[k][1])*
                                (rz[i+1][1]-rz[i][1])*(ry[k][1]-ry[i][1]) -
                                (ry[k+1][1]-ry[k][1])*(rx[i+1][1]-rx[i][1])*
                                (rz[k][1]-rz[i][1]) + (ry[k+1][1]-ry[k][1])*
                                (rz[i+1][1]-rz[i][1])*(rx[k][1]-rx[i][1]) +
                                (rz[k+1][1]-rz[k][1])*(rx[i+1][1]-rx[i][1])*
                                (ry[k][1]-ry[i][1]) - (rz[k+1][1]-rz[k][1])*
                                (ry[i+1][1]-ry[i][1])*(rx[k][1]-rx[i][1]) );
                                        }
					fprintf(fptr9,"%d\t%d\t%.5f\n",k,i,Wr_sum);
                                }
				fprintf(fptr9,"\n");
                        }
                        fclose(fptr9);*/

			/* Chain 2 */
			/*Wr_sum=0.0;
                        for(k=26;k<=N-1;k++)
                        {
                                for(i=26;i<=N-1;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][1]-
                                rx[i][1],2.0)+pow(ry[k][1]-ry[i][1],2.0)+
                                pow(rz[k][1]-rz[i][1],2.0)),3.0) )*(
                                (rx[k+1][1]-rx[k][1])*(ry[i+1][1]-ry[i][1])*
                                (rz[k][1]-rz[i][1]) - (rx[k+1][1]-rx[k][1])*
                                (rz[i+1][1]-rz[i][1])*(ry[k][1]-ry[i][1]) -
                                (ry[k+1][1]-ry[k][1])*(rx[i+1][1]-rx[i][1])*
                                (rz[k][1]-rz[i][1]) + (ry[k+1][1]-ry[k][1])*
                                (rz[i+1][1]-rz[i][1])*(rx[k][1]-rx[i][1]) +
                                (rz[k+1][1]-rz[k][1])*(rx[i+1][1]-rx[i][1])*
                                (ry[k][1]-ry[i][1]) - (rz[k+1][1]-rz[k][1])*
                                (ry[i+1][1]-ry[i][1])*(rx[k][1]-rx[i][1]) );
                                        }
                                }
                        }
                        Wr2=Wr_sum/(4.0*Pi);

			Wr_sum=0.0;
                        fptr6=fopen(filename6,"w");
                        fprintf(fptr6,"t(microseconds)=0.0\n");
                        fclose(fptr6);
                        fptr6=fopen(filename6,"a");
                        for(k=26;k<=N-1;k++)
                        {
                                for(i=26;i<=N-1;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][1]-
                                rx[i][1],2.0)+pow(ry[k][1]-ry[i][1],2.0)+
                                pow(rz[k][1]-rz[i][1],2.0)),3.0) )*(
                                (rx[k+1][1]-rx[k][1])*(ry[i+1][1]-ry[i][1])*
                                (rz[k][1]-rz[i][1]) - (rx[k+1][1]-rx[k][1])*
                                (rz[i+1][1]-rz[i][1])*(ry[k][1]-ry[i][1]) -
                                (ry[k+1][1]-ry[k][1])*(rx[i+1][1]-rx[i][1])*
                                (rz[k][1]-rz[i][1]) + (ry[k+1][1]-ry[k][1])*
                                (rz[i+1][1]-rz[i][1])*(rx[k][1]-rx[i][1]) +
                                (rz[k+1][1]-rz[k][1])*(rx[i+1][1]-rx[i][1])*
                                (ry[k][1]-ry[i][1]) - (rz[k+1][1]-rz[k][1])*
                                (ry[i+1][1]-ry[i][1])*(rx[k][1]-rx[i][1]) );
                                        }
                                }
                        fprintf(fptr6,"Bead=%d\tSWr2=%.5f\n",k,Wr_sum);
                        Wr_sum=0.0;
                        }
                        fclose(fptr6);

                        Tw_sum=0.0;
			for(k=26;k<=N-1;k++)
                        {
                                Tw_sum+=a_p_g[k][1];
                        }
                        Tw2=Tw_sum/(2.0*Pi);

			fptr8=fopen(filename8,"w");
                        fprintf(fptr8,"t(microseconds)=0.0\n");
                        fclose(fptr8);
                        fptr8=fopen(filename8,"a");
                        for(k=26;k<=N-1;k++)
                        {
                                fprintf(fptr8,"Bead=%d\tSTw2=%.5f\n",k,a_p_g[k][1]);
                        }
                        fclose(fptr8);

                        Lk2=Tw2+Wr2;
			dens2=( (Lk2+21.4)-21.4 )/21.4;

                        fptr3=fopen(filename3,"w");
                        fprintf(fptr3,"j=0\n");
                        fprintf(fptr3,"Tw2=%.5f\tWr2=%.5f\tLk2=%.5f\tdens2=%.5f\n",
                        Tw2,Wr2,Lk2,dens2);
                        fclose(fptr3);*/

			/*Wr_sum=0.0;
                        fptr10=fopen(filename10,"w");
                        fprintf(fptr10,"t(microseconds)=0.0\n");
                        fclose(fptr10);
                        fptr10=fopen(filename10,"a");
                        for(k=26;k<=49;k++)
                        {
                                for(i=26;i<=49;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum=0.0;
                                        }
                                        else
                                        {
                                Wr_sum=( 1.0/pow(sqrt(pow(rx[k][1]-
                                rx[i][1],2.0)+pow(ry[k][1]-ry[i][1],2.0)+
                                pow(rz[k][1]-rz[i][1],2.0)),3.0) )*(
                                (rx[k+1][1]-rx[k][1])*(ry[i+1][1]-ry[i][1])*
                                (rz[k][1]-rz[i][1]) - (rx[k+1][1]-rx[k][1])*
                                (rz[i+1][1]-rz[i][1])*(ry[k][1]-ry[i][1]) -
                                (ry[k+1][1]-ry[k][1])*(rx[i+1][1]-rx[i][1])*
                                (rz[k][1]-rz[i][1]) + (ry[k+1][1]-ry[k][1])*
                                (rz[i+1][1]-rz[i][1])*(rx[k][1]-rx[i][1]) +
                                (rz[k+1][1]-rz[k][1])*(rx[i+1][1]-rx[i][1])*
                                (ry[k][1]-ry[i][1]) - (rz[k+1][1]-rz[k][1])*
                                (ry[i+1][1]-ry[i][1])*(rx[k][1]-rx[i][1]) );
                                        }
                                        fprintf(fptr10,"%d\t%d\t%.5f\n",k,i,Wr_sum);
                                }
                                fprintf(fptr10,"\n");
                        }
                        fclose(fptr10);

			fptr4=fopen(filename4,"w");
                        fprintf(fptr4,"j=0\n");
                        fprintf(fptr4,"Tw_tot=%.5f\tWr_tot=%.5f\tLk_tot=%.5f\n",
                        Tw1+Tw2,Wr1+Wr2,Lk1+Lk2);
                        fclose(fptr4);*/

			/*Lk_tot_sum=Lk1+Lk2,Tw_tot_sum=Tw1+Tw2,Wr_tot_sum=Wr1+Wr2;*/

/* Initial Torques */
        tau[1]=(kB*temp/(ZI*ZI))*a_p_g[1][1]+(1.0e-20)*exp(-pow(1-1,2.0)/3.56);
	tau[26]=(kB*temp/(ZI*ZI))*a_p_g[26][1]-(1.0e-20)*exp(-pow(26-26,2.0)/3.56);
	for(i=2;i<=24;i++)
        {
                tau[i]=(kB*temp/(ZI*ZI))*(a_p_g[i][1]-a_p_g[i-1][1])+(1.0e-20)*
			exp(-pow(i-1,2.0)/3.56);
        }
        for(i=27;i<=N-1;i++)
        {
                tau[i]=(kB*temp/(ZI*ZI))*(a_p_g[i][1]-a_p_g[i-1][1])-(1.0e-20)*
			exp(-pow(i-26,2.0)/3.56);
        }
	
/* Initial Forces */
fptr2=fopen(filename2,"a");
        for(i=1;i<=N;i++){
                F_tot_x[i]=F_s_x_sub(i,1,r_ij,u_norm_x) + F_b_x_sub(i,1,beta,r_ij,u_norm_x)
                        + F_t_x_sub(i,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g) +
                        F_EV_x_sub(i,1,rx,r_ij);
                F_tot_y[i]=F_s_y_sub(i,1,r_ij,u_norm_y) + F_b_y_sub(i,1,beta,r_ij,u_norm_y)
                        + F_t_y_sub(i,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g) +
                        F_EV_y_sub(i,1,ry,r_ij);
                F_tot_z[i]=F_s_z_sub(i,1,r_ij,u_norm_z) + F_b_z_sub(i,1,beta,r_ij,u_norm_z)
                        + F_t_z_sub(i,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g) +
                        F_EV_z_sub(i,1,rz,r_ij);
        }

fprintf(fptr2,"%d\t%.3f\t%.3f\t%.3f\t%.3f\n",1,F_s_x_sub(1,1,r_ij,u_norm_x)*1e10,F_b_x_sub(1,1,beta,r_ij,u_norm_x)*1e10,
                        F_t_x_sub(1,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g)*1e10,
                        F_EV_x_sub(1,1,rx,r_ij)*1e10);
fprintf(fptr2,"\n%d\t%.3f\t%.3f\t%.3f\t%.3f\n",2,F_s_x_sub(2,1,r_ij,u_norm_x)*1e10,F_b_x_sub(2,1,beta,r_ij,u_norm_x)*1e10,
                        F_t_x_sub(2,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g)*1e10,
                        F_EV_x_sub(2,1,rx,r_ij)*1e10);
fprintf(fptr2,"\n%d\t%.3f\t%.3f\t%.3f\t%.3f\n",25,F_s_x_sub(25,1,r_ij,u_norm_x)*1e10,F_b_x_sub(25,1,beta,r_ij,u_norm_x)*1e10,
                        F_t_x_sub(25,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g)*1e10,
                        F_EV_x_sub(25,1,rx,r_ij)*1e10);
fprintf(fptr2,"%d\t%.3f\t%.3f\t%.3f\t%.3f\n",26,F_s_x_sub(26,1,r_ij,u_norm_x)*1e10,F_b_x_sub(26,1,beta,r_ij,u_norm_x)*1e10,
                        F_t_x_sub(26,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g)*1e10,
                        F_EV_x_sub(26,1,rx,r_ij)*1e10);
fprintf(fptr2,"\n%d\t%.3f\t%.3f\t%.3f\t%.3f\n",27,F_s_x_sub(27,1,r_ij,u_norm_x)*1e10,F_b_x_sub(27,1,beta,r_ij,u_norm_x)*1e10,
                        F_t_x_sub(27,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g)*1e10,
                        F_EV_x_sub(27,1,rx,r_ij)*1e10);
fprintf(fptr2,"\n%d\t%.3f\t%.3f\t%.3f\t%.3f\n",50,F_s_x_sub(50,1,r_ij,u_norm_x)*1e10,F_b_x_sub(50,1,beta,r_ij,u_norm_x)*1e10,
                        F_t_x_sub(50,1,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                        f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g)*1e10,
                        F_EV_x_sub(50,1,rx,r_ij)*1e10);

fclose(fptr2);



/* Time Evolution */
        for(j=1;j<=Nt;j++)
	{
	/* Calculate New Phi */
                for(i=1;i<=24;i++){
                        phi[i][2]=phi[i][1] + (D_rot*dt/(kB*temp))*tau[i] +
                                gasdev(&idum,2.0*D_rot*dt);
                }
		    for(i=26;i<=N-1;i++){
                        phi[i][2]=phi[i][1] + (D_rot*dt/(kB*temp))*tau[i] +
                                gasdev(&idum,2.0*D_rot*dt);
                }

	/* Calculate New Positions */
		rx[1][2]=rx[1][1],ry[1][2]=ry[1][1],rz[1][2]=rz[1][1];
		rx[25][2]=rx[25][1],ry[25][2]=ry[25][1],rz[25][2]=rz[25][1];
		rx[26][2]=rx[26][1],ry[26][2]=ry[26][1],rz[26][2]=rz[26][1];
		rx[50][2]=rx[50][1],ry[50][2]=ry[50][1],rz[50][2]=rz[50][1];
		for(i=2;i<=24;i++){
			rx[i][2]=rx[i][1] + (D_trans*dt/(kB*temp))*F_tot_x[i] +
                                gasdev(&idum,2.0*D_trans*dt);
			ry[i][2]=ry[i][1] + (D_trans*dt/(kB*temp))*F_tot_y[i] +
                                gasdev(&idum,2.0*D_trans*dt);
			rz[i][2]=rz[i][1] + (D_trans*dt/(kB*temp))*F_tot_z[i] +
                                gasdev(&idum,2.0*D_trans*dt);
		}
		for(i=27;i<=N-1;i++){
			rx[i][2]=rx[i][1] + (D_trans*dt/(kB*temp))*F_tot_x[i] +
                                gasdev(&idum,2.0*D_trans*dt);
			ry[i][2]=ry[i][1] + (D_trans*dt/(kB*temp))*F_tot_y[i] +
                                gasdev(&idum,2.0*D_trans*dt);
			rz[i][2]=rz[i][1] + (D_trans*dt/(kB*temp))*F_tot_z[i] +
                                gasdev(&idum,2.0*D_trans*dt);
		}
		
	/* Calculate New Separations */
		for(i=1;i<=N;i++)
                {
                        for(k=1;k<=N;k++)
                        {
                                r_ij[i][k]=sqrt(pow(rx[k][2]-rx[i][2],2.0)+
                                        pow(ry[k][2]-ry[i][2],2.0)+
                                        pow(rz[k][2]-rz[i][2],2.0));
                        }
                }

	/* Calculate New BFCs */
		/* Calculate New u_norm[i] */
		for(i=1;i<=24;i++){
			u_norm_x[i][2]=(rx[i+1][2]-rx[i][2])/r_ij[i+1][i];
                        u_norm_y[i][2]=(ry[i+1][2]-ry[i][2])/r_ij[i+1][i];
                        u_norm_z[i][2]=(rz[i+1][2]-rz[i][2])/r_ij[i+1][i];
		}	
		
		for(i=26;i<=N-1;i++){
			u_norm_x[i][2]=(rx[i+1][2]-rx[i][2])/r_ij[i+1][i];
                        u_norm_y[i][2]=(ry[i+1][2]-ry[i][2])/r_ij[i+1][i];
                        u_norm_z[i][2]=(rz[i+1][2]-rz[i][2])/r_ij[i+1][i];
		}
		u_norm_x[25][2]=u_norm_x[25][1];
                u_norm_y[25][2]=u_norm_y[25][1];
                u_norm_z[25][2]=u_norm_z[25][1];
		u_norm_x[N][2]=u_norm_x[N][1];
                u_norm_y[N][2]=u_norm_y[N][1];
                u_norm_z[N][2]=u_norm_z[N][1];

		/* Calculate New f_norm[i] */
                        /* calculate new f[i] */
                for(i=1;i<=24;i++)
                {
                        fx[i]=f_norm_x[i] + (phi[i][2]-phi[i][1])*v_norm_x[i] -
                                u_norm_x[i][1]*(f_norm_x[i]*u_norm_x[i][2]+f_norm_y[i]*
                                u_norm_y[i][2]+f_norm_z[i]*u_norm_z[i][2]);
                        fy[i]=f_norm_y[i] + (phi[i][2]-phi[i][1])*v_norm_y[i] -
                                u_norm_y[i][1]*(f_norm_x[i]*u_norm_x[i][2]+f_norm_y[i]*
                                u_norm_y[i][2]+f_norm_z[i]*u_norm_z[i][2]);
                        fz[i]=f_norm_z[i] + (phi[i][2]-phi[i][1])*v_norm_z[i] -
                                u_norm_z[i][1]*(f_norm_x[i]*u_norm_x[i][2]+f_norm_y[i]*
                                u_norm_y[i][2]+f_norm_z[i]*u_norm_z[i][2]);
                }
		    for(i=26;i<=N-1;i++)
                {
                        fx[i]=f_norm_x[i] + (phi[i][2]-phi[i][1])*v_norm_x[i] -
                                u_norm_x[i][1]*(f_norm_x[i]*u_norm_x[i][2]+f_norm_y[i]*
                                u_norm_y[i][2]+f_norm_z[i]*u_norm_z[i][2]);
                        fy[i]=f_norm_y[i] + (phi[i][2]-phi[i][1])*v_norm_y[i] -
                                u_norm_y[i][1]*(f_norm_x[i]*u_norm_x[i][2]+f_norm_y[i]*
                                u_norm_y[i][2]+f_norm_z[i]*u_norm_z[i][2]);
                        fz[i]=f_norm_z[i] + (phi[i][2]-phi[i][1])*v_norm_z[i] -
                                u_norm_z[i][1]*(f_norm_x[i]*u_norm_x[i][2]+f_norm_y[i]*
                                u_norm_y[i][2]+f_norm_z[i]*u_norm_z[i][2]);
                }
                        /* "perpendicularize" new f[i] */
                for(i=1;i<=24;i++)
                {
                        fx[i]=fx[i]-u_norm_x[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                        fy[i]=fy[i]-u_norm_y[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                        fz[i]=fz[i]-u_norm_z[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                }
		    for(i=26;i<=N-1;i++)
                {
                        fx[i]=fx[i]-u_norm_x[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                        fy[i]=fy[i]-u_norm_y[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                        fz[i]=fz[i]-u_norm_z[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                }
                        /* normalize new f[i] */
                for(i=1;i<=24;i++)
                {
                        f_norm_x[i]=fx[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                        f_norm_y[i]=fy[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                        f_norm_z[i]=fz[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                }
		    for(i=26;i<=N-1;i++)
                {
                        f_norm_x[i]=fx[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                        f_norm_y[i]=fy[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                        f_norm_z[i]=fz[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                }

		/* Calculate New v_norm[i] */
                for(i=1;i<=24;i++)
                {
                        vx[i]=u_norm_y[i][2]*f_norm_z[i]-u_norm_z[i][2]*f_norm_y[i];
                        vy[i]=u_norm_z[i][2]*f_norm_x[i]-u_norm_x[i][2]*f_norm_z[i];
                        vz[i]=u_norm_x[i][2]*f_norm_y[i]-u_norm_y[i][2]*f_norm_x[i];

                        v_norm_x[i]=vx[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                        v_norm_y[i]=vy[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                        v_norm_z[i]=vz[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                }
		for(i=26;i<=N-1;i++)
                {
                        vx[i]=u_norm_y[i][2]*f_norm_z[i]-u_norm_z[i][2]*f_norm_y[i];
                        vy[i]=u_norm_z[i][2]*f_norm_x[i]-u_norm_x[i][2]*f_norm_z[i];
                        vz[i]=u_norm_x[i][2]*f_norm_y[i]-u_norm_y[i][2]*f_norm_x[i];

                        v_norm_x[i]=vx[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                        v_norm_y[i]=vy[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                        v_norm_z[i]=vz[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                }

	/* Calculate New beta[i] */
                for(i=1;i<=24;i++)
                {
                        arg_beta=(u_norm_x[i][2]*u_norm_x[i+1][2]+u_norm_y[i][2]*
                                u_norm_y[i+1][2]+u_norm_z[i][2]*u_norm_z[i+1][2])/
                                ( sqrt(u_norm_x[i][2]*u_norm_x[i][2]+u_norm_y[i][2]*u_norm_y[i][2]+
                                  u_norm_z[i][2]*u_norm_z[i][2])*sqrt(u_norm_x[i+1][2]*u_norm_x[i+1][2]+
                                  u_norm_y[i+1][2]*u_norm_y[i+1][2]+
                                  u_norm_z[i+1][2]*u_norm_z[i+1][2]) );
                        beta[i]=acos(arg_beta);
                }
		for(i=26;i<=N-1;i++)
                {
                        arg_beta=(u_norm_x[i][2]*u_norm_x[i+1][2]+u_norm_y[i][2]*
                                u_norm_y[i+1][2]+u_norm_z[i][2]*u_norm_z[i+1][2])/
                                ( sqrt(u_norm_x[i][2]*u_norm_x[i][2]+u_norm_y[i][2]*u_norm_y[i][2]+
                                  u_norm_z[i][2]*u_norm_z[i][2])*sqrt(u_norm_x[i+1][2]*u_norm_x[i+1][2]+
                                  u_norm_y[i+1][2]*u_norm_y[i+1][2]+
                                  u_norm_z[i+1][2]*u_norm_z[i+1][2]) );
                        beta[i]=acos(arg_beta);
                }

	/* Calculate New a_p_g[i] */
		for(i=1;i<=24;i++)
                {
                        a_p_g[i][2]=asin( (v_norm_x[i]*f_norm_x[i+1]+v_norm_y[i]*f_norm_y[i+1]+
                        v_norm_z[i]*f_norm_z[i+1]-f_norm_x[i]*v_norm_x[i+1]-
                        f_norm_y[i]*v_norm_y[i+1]-f_norm_z[i]*v_norm_z[i+1])/
                        (2.0*cos(beta[i]/2.0)*cos(beta[i]/2.0)) );
                }
		for(i=26;i<=N-1;i++)
                {
                        a_p_g[i][2]=asin( (v_norm_x[i]*f_norm_x[i+1]+v_norm_y[i]*f_norm_y[i+1]+
                        v_norm_z[i]*f_norm_z[i+1]-f_norm_x[i]*v_norm_x[i+1]-
                        f_norm_y[i]*v_norm_y[i+1]-f_norm_z[i]*v_norm_z[i+1])/
                        (2.0*cos(beta[i]/2.0)*cos(beta[i]/2.0)) );
                }

    	/* Calculate New tau[i] */
                tau[1]=(kB*temp/(ZI*ZI))*a_p_g[1][2]+(1.0e-20)
			*exp(-pow(1-1,2.0)/3.56);
		    tau[26]=(kB*temp/(ZI*ZI))*a_p_g[26][2]-(1.0e-20)
			*exp(-pow(26-26,2.0)/3.56);			
                for(i=2;i<=24;i++)
                {
                        tau[i]=(kB*temp/(ZI*ZI))*(a_p_g[i][2]-a_p_g[i-1][2])+
				(1.0e-20)*exp(-pow(i-1,2.0)/3.56);
                }
		    for(i=27;i<=N-1;i++)
        	    {
                tau[i]=(kB*temp/(ZI*ZI))*(a_p_g[i][2]-a_p_g[i-1][2])-(1.0e-20)*
			exp(-pow(i-26,2.0)/3.56);
		    }
        
			/* Calculate New F_tot[i]: */
                for(i=1;i<=N;i++){
                        F_tot_x[i]=F_s_x_sub(i,2,r_ij,u_norm_x) + F_b_x_sub(i,2,beta,r_ij,u_norm_x)
                                + F_t_x_sub(i,2,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                                f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g) +
                                F_EV_x_sub(i,2,rx,r_ij);
                        F_tot_y[i]=F_s_y_sub(i,2,r_ij,u_norm_y) + F_b_y_sub(i,2,beta,r_ij,u_norm_y)
                                + F_t_y_sub(i,2,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                                f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g) +
                                F_EV_y_sub(i,2,ry,r_ij);
                        F_tot_z[i]=F_s_z_sub(i,2,r_ij,u_norm_z) + F_b_z_sub(i,2,beta,r_ij,u_norm_z)
                                + F_t_z_sub(i,2,r_ij,u_norm_x,u_norm_y,u_norm_z,f_norm_x,f_norm_y,
                                f_norm_z,v_norm_x,v_norm_y,v_norm_z,a_p_g) +
                                F_EV_z_sub(i,2,rz,r_ij);
                }

	/* Calculate New r_f[i] */
                for(i=1;i<=N;i++)
                {
                        r_f_x[i]=rx[i][2]+0.5*R_hyd*f_norm_x[i];
                        r_f_y[i]=ry[i][2]+0.5*R_hyd*f_norm_y[i];
                        r_f_z[i]=rz[i][2]+0.5*R_hyd*f_norm_z[i];
                }
		
			/* Calculate New Lk=Tw+Wr */
                                /* Chain 1 */
                        /*Wr_sum=0.0;
                        for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
                                        }
                                }
                        }
                        Wr1=Wr_sum/(4.0*Pi);

                        Tw_sum=0.0;
                        for(k=1;k<=24;k++)
                        {
                                Tw_sum+=a_p_g[k][2];
                        }
                        Tw1=Tw_sum/(2.0*Pi);

                        Lk1=Tw1+Wr1;
                        dens1=( (Lk1+21.4)-21.4 )/21.4;*/

			 /* Chain 2 */
                        /*Wr_sum=0.0;
                        for(k=26;k<=N-1;k++)
                        {
                                for(i=26;i<=N-1;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
                                        }
                                }
                        }
                        Wr2=Wr_sum/(4.0*Pi);

                        Tw_sum=0.0;
                        for(k=26;k<=N-1;k++)
                        {
                                Tw_sum+=a_p_g[k][2];
                        }
                        Tw2=Tw_sum/(2.0*Pi);

                        Lk2=Tw2+Wr2;
                        dens2=( (Lk2+21.4)-21.4 )/21.4;

			Lk_tot_sum+=Lk1+Lk2,Tw_tot_sum+=Tw1+Tw2,Wr_tot_sum+=Wr1+Wr2;*/

/*printf("j=%d\n",j);*/
                if((j+1000)%1000==0){
			fptr1=fopen(filename1,"a");
                	fprintf(fptr1,"\nj=%d\n",j);
        		for(i=50;i>=26;i--){
                		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",rx[i][2]*1e10,       
                        		ry[i][2]*1e10,rz[i][2]*1e10);   
                		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",r_f_x[i]*1e10,       
                        		r_f_y[i]*1e10,r_f_z[i]*1e10);   
        		}       
        		for(i=1;i<=25;i++){
                		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",rx[i][2]*1e10,
                        		ry[i][2]*1e10,rz[i][2]*1e10);
                		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",r_f_x[i]*1e10,
                        		r_f_y[i]*1e10,r_f_z[i]*1e10);
        		}
			fclose(fptr1);

			/* Calculate New Lk=Tw+Wr */
				/* Chain 1 */
                        /*Wr_sum=0.0;
                        for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
                                        }
                                }
                        }
                        Wr1=Wr_sum/(4.0*Pi);*/

			/*Wr_sum=0.0;
                        fptr5=fopen(filename5,"a");
			fprintf(fptr5,"\n");
                        fprintf(fptr5,"t(microseconds)=%.5f\n",j*5e-6);
			for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
                                        }
                                }
			fprintf(fptr5,"Bead=%d\tSWr1=%.5f\n",k,Wr_sum);
                        Wr_sum=0.0;
                        }
			fclose(fptr5);*/

                        /*Tw_sum=0.0;
                        for(k=1;k<=24;k++)
                        {
                                Tw_sum+=a_p_g[k][2];
                        }
                        Tw1=Tw_sum/(2.0*Pi);

			fptr7=fopen(filename7,"a");
			fprintf(fptr7,"\n");
                        fprintf(fptr7,"t(microseconds)=%.5f\n",j*5e-6);
                        for(k=1;k<=24;k++)
                        {
                                fprintf(fptr7,"Bead=%d\tSTw1=%.5f\n",k,a_p_g[k][2]);
                        }
                        fclose(fptr7);

                        Lk1=Tw1+Wr1;
                        dens1=( (Lk1+21.4)-21.4 )/21.4;
			
			fptr2=fopen(filename2,"a");
                        fprintf(fptr2,"\nj=%d\n",j);
                        fprintf(fptr2,"Tw1=%.5f\tWr1=%.5f\tLk1=%.5f\tdens1=%.5f\n",
                        Tw1,Wr1,Lk1,dens1);
                        fclose(fptr2);*/

			/*Wr_sum=0.0;
                        fptr9=fopen(filename9,"a");
                        fprintf(fptr9,"t(microseconds)=%.5f\n",j*5e-6);
                        for(k=1;k<=24;k++)
                        {
                                for(i=1;i<=24;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum=0.0;
                                        }
                                        else
                                        {
                                Wr_sum=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
 
                                        }
                                        fprintf(fptr9,"%d\t%d\t%.5f\n",k,i,Wr_sum);
                                }
                                fprintf(fptr9,"\n");
                        }
                        fclose(fptr9);*/

				/* Chain 2 */
                        /*Wr_sum=0.0;
                        for(k=26;k<=N-1;k++)
                        {
                                for(i=26;i<=N-1;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
                                        }
                                }
                        }
                        Wr2=Wr_sum/(4.0*Pi);

			Wr_sum=0.0;
                        fptr6=fopen(filename6,"a");
                        fprintf(fptr6,"\n");
                        fprintf(fptr6,"t(microseconds)=%.5f\n",j*5e-6);
			for(k=26;k<=N-1;k++)
                        {
                                for(i=26;i<=N-1;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum+=0.0;
                                        }
                                        else
                                        {
                                Wr_sum+=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );
                                        }
                                }
			fprintf(fptr6,"Bead=%d\tSWr2=%.5f\n",k,Wr_sum);
                        Wr_sum=0.0;
                        }
			fclose(fptr6);

                        Tw_sum=0.0;
                        for(k=26;k<=N-1;k++)
                        {
				Tw_sum+=a_p_g[k][2];
                        }
                        Tw2=Tw_sum/(2.0*Pi);

			fptr8=fopen(filename8,"a");
                        fprintf(fptr8,"\n");
                        fprintf(fptr8,"t(microseconds)=%.5f\n",j*5e-6);
                        for(k=26;k<=N-1;k++)
                        {
                                fprintf(fptr8,"Bead=%d\tSTw2=%.5f\n",k,a_p_g[k][2]);
                        }
                        fclose(fptr8);

                        Lk2=Tw2+Wr2;
                        dens2=( (Lk2+21.4)-21.4 )/21.4;

                        fptr3=fopen(filename3,"a");
                        fprintf(fptr3,"\nj=%d\n",j);
                        fprintf(fptr3,"Tw2=%.5f\tWr2=%.5f\tLk2=%.5f\tdens2=%.5f\n",
                        Tw2,Wr2,Lk2,dens2);
                        fclose(fptr3);*/

			/*Wr_sum=0.0;
                        fptr10=fopen(filename10,"a");
                        fprintf(fptr10,"t(microseconds)=%.5f\n",j*5e-6);
                        for(k=26;k<=49;k++)
                        {
                                for(i=26;i<=49;i++)
                                {
                                        if(i==k)
                                        {
                                                Wr_sum=0.0;
                                        }
                                        else
                                        {
                                Wr_sum=( 1.0/pow(sqrt(pow(rx[k][2]-
                                rx[i][2],2.0)+pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0)),3.0) )*(
                                (rx[k+1][2]-rx[k][2])*(ry[i+1][2]-ry[i][2])*
                                (rz[k][2]-rz[i][2]) - (rx[k+1][2]-rx[k][2])*
                                (rz[i+1][2]-rz[i][2])*(ry[k][2]-ry[i][2]) -
                                (ry[k+1][2]-ry[k][2])*(rx[i+1][2]-rx[i][2])*
                                (rz[k][2]-rz[i][2]) + (ry[k+1][2]-ry[k][2])*
                                (rz[i+1][2]-rz[i][2])*(rx[k][2]-rx[i][2]) +
                                (rz[k+1][2]-rz[k][2])*(rx[i+1][2]-rx[i][2])*
                                (ry[k][2]-ry[i][2]) - (rz[k+1][2]-rz[k][2])*
                                (ry[i+1][2]-ry[i][2])*(rx[k][2]-rx[i][2]) );

                                        }
                                        fprintf(fptr10,"%d\t%d\t%.5f\n",k,i,Wr_sum);
                                }
                                fprintf(fptr10,"\n");
                        }
                        fclose(fptr10);*/

                        /*fptr4=fopen(filename4,"a");
                        fprintf(fptr4,"\nj=%d\n",j);
                        fprintf(fptr4,"Tw_tot=%.5f\tWr_tot=%.5f\tLk_tot=%.5f\n",
                        Tw1+Tw2,Wr1+Wr2,Lk1+Lk2);
                        fclose(fptr4);*/
			
			/*Lk_tot_sum+=Lk1+Lk2,Tw_tot_sum+=Tw1+Tw2,Wr_tot_sum+=Wr1+Wr2;	

			fptr4=fopen(filename4,"w");
			fprintf(fptr4,"\nj=%d\n",j);
			fprintf(fptr4,"Tw_tot_av=%.5f\tWr_tot_av=%.5f\tLk_tot_av=%.5f\n",
			Tw_tot_sum/(j/1000.0+1.0),Wr_tot_sum/(j/1000.0+1.0),Lk_tot_sum/(j/1000.0+1.0));
			fclose(fptr4);*/
		}

	/* r[i][2] -> r[i][1] for next step: */
                for(i=1;i<=N;i++)
                {
                        rx[i][1]=rx[i][2];
                        ry[i][1]=ry[i][2];
                        rz[i][1]=rz[i][2];
                }

        /* u_norm[i][2] -> u_norm[i][1] for next step: */
                for(i=1;i<=24;i++)
                {
                        u_norm_x[i][1]=u_norm_x[i][2];
                        u_norm_y[i][1]=u_norm_y[i][2];
                        u_norm_z[i][1]=u_norm_z[i][2];
                }
                for(i=26;i<=N-1;i++)
                {
                        u_norm_x[i][1]=u_norm_x[i][2];
                        u_norm_y[i][1]=u_norm_y[i][2];
                        u_norm_z[i][1]=u_norm_z[i][2];
                }

	/* phi[i][2] -> phi[i][1] for next step */
                for(i=1;i<=24;i++){
                        phi[i][1]=phi[i][2];
                }
		    for(i=26;i<=N-1;i++){
                        phi[i][1]=phi[i][2];
                }
	}
/*fptr4=fopen(filename4,"w");
fprintf(fptr4,"\n");
fprintf(fptr4,"Tw_tot_av=%.5f\tWr_tot_av=%.5f\tLk_tot_av=%.5f\n",
Tw_tot_sum/20000001.0,Wr_tot_sum/20000001.0,Lk_tot_sum/20000001.0);
fclose(fptr4);*/
}
