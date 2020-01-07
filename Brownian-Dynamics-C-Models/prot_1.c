/* prot_1.c: */

/* Chain Parameters */
#define N       50		
#define b0	31.84e-10	
#define R_hyd	15.92e-10	
#define DELTA	0.008*b0  
#define PSI	sqrt(b0/500.0e-10)		/* bending parameter (in radiants) */ 
#define ZI	sqrt(b0*kB*temp/2.6e-28)	/* torsion parameter (C_t=2.6e-28 J m) */
#define EPSILON 100.0		
#define SIGMA	b0/pow(2.0,1.0/6.0)	
#define CUTOFF	b0	

/* Solvent Parameters */
#define kB	1.3806503e-23	
#define temp	310.0   	
#define ETA	0.01	
#define D_rot   kB*temp/(Pi*ETA*pow(R_hyd,2.0)*b0)
#define D_trans kB*temp/(6.0*Pi*ETA*R_hyd) 

/* Numerical Parameters */
#define dt      5.0e-12		
#define Nt      20000000		
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

        if(i==1){
                f_s_x=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_x[i][k];  
        }
	if(i==N){
                f_s_x=-( kB*temp/pow(DELTA,2.0) )*(r_ij[N][N-1]-b0)*u_norm_x[N-1][k];
        }
	if(i!=1 && i!=N){
		f_s_x=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_x[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_x[i-1][k] );
        }
        return f_s_x;
}
double F_s_y_sub(int i, int k, double **r_ij, double **u_norm_y)
{
        double f_s_y;

	if(i==1){
                f_s_y=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_y[i][k];
        }
        if(i==N){
                f_s_y=-( kB*temp/pow(DELTA,2.0) )*(r_ij[N][N-1]-b0)*u_norm_y[N-1][k];
        }
        if(i!=1 && i!=N){
                f_s_y=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_y[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_y[i-1][k] );
        }
        return f_s_y;
}
double F_s_z_sub(int i, int k, double **r_ij, double **u_norm_z)
{
        double f_s_z;

	if(i==1){
                f_s_z=( kB*temp/pow(DELTA,2.0) )*(r_ij[i][i+1]-b0)*u_norm_z[i][k];
        }
        if(i==N){
                f_s_z=-( kB*temp/pow(DELTA,2.0) )*(r_ij[N][N-1]-b0)*u_norm_z[N-1][k];
        }
        if(i!=1 && i!=N){
                f_s_z=( kB*temp/pow(DELTA,2.0) )*( (r_ij[i][i+1]-b0)*u_norm_z[i][k] -
                (r_ij[i][i-1]-b0)*u_norm_z[i-1][k] );
        }
        return f_s_z;
}

/* Subroutines for Bending Force Components */
double F_b_x_sub(int i, int k, double *beta, double **r_ij, double **u_norm_x)
{
        double f_b_x;
	double *A_x,*B_x; 

	A_x=dvector(1,N),B_x=dvector(1,N);

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
        if(i!=1 && i!=2 && i!=N-1 && i!=N){
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
        if(i==N-1){
                B_x[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_x[i-1][k]-u_norm_x[i][k]*cos(beta[i-1]) );
                A_x[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_x[i][k]-u_norm_x[i-1][k]*cos(beta[i-1]) );
                B_x[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_x[i-2][k]-u_norm_x[i-1][k]*cos(beta[i-2]) );

                f_b_x=( -kB*temp/pow(PSI,2.0) )*( A_x[i-1]+B_x[i]-B_x[i-1] );
        }
        if(i==N){
                B_x[N-1]=( beta[N-2]/(r_ij[N-1][N]*sin(beta[N-2])) )*
                        ( u_norm_x[N-2][k]-u_norm_x[N-1][k]*cos(beta[N-2]) );

                f_b_x=( kB*temp/pow(PSI,2.0) )*B_x[N-1] ;
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
        if(i!=1 && i!=2 && i!=N-1 && i!=N){
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
        if(i==N-1){
                B_y[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_y[i-1][k]-u_norm_y[i][k]*cos(beta[i-1]) );
                A_y[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_y[i][k]-u_norm_y[i-1][k]*cos(beta[i-1]) );
                B_y[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_y[i-2][k]-u_norm_y[i-1][k]*cos(beta[i-2]) );

                f_b_y=( -kB*temp/pow(PSI,2.0) )*( A_y[i-1]+B_y[i]-B_y[i-1] );
        }
        if(i==N){
                B_y[N-1]=( beta[N-2]/(r_ij[N-1][N]*sin(beta[N-2])) )*
                        ( u_norm_y[N-2][k]-u_norm_y[N-1][k]*cos(beta[N-2]) );

                f_b_y=( kB*temp/pow(PSI,2.0) )*B_y[N-1] ;
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
        if(i!=1 && i!=2 && i!=N-1 && i!=N){
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
        if(i==N-1){
                B_z[i]=( beta[i-1]/(r_ij[i][i+1]*sin(beta[i-1])) )*
                        ( u_norm_z[i-1][k]-u_norm_z[i][k]*cos(beta[i-1]) );
                A_z[i-1]=( beta[i-1]/(r_ij[i][i-1]*sin(beta[i-1])) )*
                        ( u_norm_z[i][k]-u_norm_z[i-1][k]*cos(beta[i-1]) );
                B_z[i-1]=( beta[i-2]/(r_ij[i-1][i]*sin(beta[i-2])) )*
                        ( u_norm_z[i-2][k]-u_norm_z[i-1][k]*cos(beta[i-2]) );

                f_b_z=( -kB*temp/pow(PSI,2.0) )*( A_z[i-1]+B_z[i]-B_z[i-1] );
        }
        if(i==N){
                B_z[N-1]=( beta[N-2]/(r_ij[N-1][N]*sin(beta[N-2])) )*
                        ( u_norm_z[N-2][k]-u_norm_z[N-1][k]*cos(beta[N-2]) );

                f_b_z=( kB*temp/pow(PSI,2.0) )*B_z[N-1] ;
        }

        free_dvector(A_z,1,N);
        free_dvector(B_z,1,N);

        return f_b_z;
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
	long idum=(SEED);

	FILE *fptr1; char filename1[]="prot_1.out";

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
        for(i=1;i<=N;i++)
        {
                beta[i]=-2.0*Pi/N;
        }

        /* Initial Positions */
        rx[1][1]=0.0, ry[1][1]=0.0, rz[1][1]=b0;
        for(i=2;i<=N;i++)
        {
                bz[i-1]=b0*cos((i-2.0)*beta[i]);
                bx[i-1]=b0*sin((i-2.0)*beta[i]);
                rx[i][1]=rx[i-1][1]+bx[i-1];
                ry[i][1]=0.0;
                rz[i][1]=rz[i-1][1]+bz[i-1];
        }
	
	/* Initial Separations */
	for(i=1;i<=N-1;i++){
                for(j=i+1;j<=N;j++){
                        r_ij[i][j]=sqrt(pow(rx[j][1]-rx[i][1],2.0)+
                                pow(ry[j][1]-ry[i][1],2.0)+
                                pow(rz[j][1]-rz[i][1],2.0));
                        r_ij[j][i]=r_ij[i][j];
                }
        }
	
	/* Initial BFCs */
		/* initial u_norm[i] */ 
	for(i=1;i<=N-1;i++){
		u_norm_x[i][1]=(rx[i+1][1]-rx[i][1])/r_ij[i+1][i];
            	u_norm_y[i][1]=(ry[i+1][1]-ry[i][1])/r_ij[i+1][i];
            	u_norm_z[i][1]=(rz[i+1][1]-rz[i][1])/r_ij[i+1][i];
                	
	}

        	/* initial f[i] */
	for(i=1;i<=N-1;i++){
        	f_norm_x[i]=0.0;     
            	f_norm_y[i]=1.0;
            	f_norm_z[i]=0.0;
	}	

                /* initial v_norm[i] */
	for(i=1;i<=N-1;i++){
                vx[i]=u_norm_y[i][1]*f_norm_z[i]-u_norm_z[i][1]*f_norm_y[i];
                vy[i]=u_norm_z[i][1]*f_norm_x[i]-u_norm_x[i][1]*f_norm_z[i];
                vz[i]=u_norm_x[i][1]*f_norm_y[i]-u_norm_y[i][1]*f_norm_x[i];

                v_norm_x[i]=vx[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                v_norm_y[i]=vy[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                v_norm_z[i]=vz[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
	}

	for(i=1;i<=N-1;i++){
	        r_f_x[i]=rx[i][1]+0.5*R_hyd*f_norm_x[i];
                r_f_y[i]=ry[i][1]+0.5*R_hyd*f_norm_y[i];
                r_f_z[i]=rz[i][1]+0.5*R_hyd*f_norm_z[i];
	}

	/* Write Initial Config. to "*.out" */
	fptr1=fopen(filename1,"w");
       	fprintf(fptr1,"j=0\n"); 
	for(i=1;i<=N;i++){ 
        	fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",rx[i][1]*1e10,
                        ry[i][1]*1e10,rz[i][1]*1e10);
        }
	for(i=1;i<=N-1;i++){
		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",r_f_x[i]*1e10,
                        r_f_y[i]*1e10,r_f_z[i]*1e10);
	}
        fclose(fptr1);

/* Initial Phi */
	for(i=1;i<=N-2;i++){
		phi[i][1]=0.0;
	}

/* Initial alpha[i]+gamma[i] (a_p_g[i][1]) */
	for(i=1;i<=N-2;i++){
                a_p_g[i][1]=asin( (v_norm_x[i]*f_norm_x[i+1]+v_norm_y[i]*f_norm_y[i+1]+
                        v_norm_z[i]*f_norm_z[i+1]-f_norm_x[i]*v_norm_x[i+1]-
                        f_norm_y[i]*v_norm_y[i+1]-f_norm_z[i]*v_norm_z[i+1])/
                        (2.0*cos(beta[i]/2.0)*cos(beta[i]/2.0)) );
	}

/* Initial Torques */
        tau[1]=(kB*temp/(ZI*ZI))*a_p_g[1][1];
	for(i=2;i<=N-2;i++)
        {
                tau[i]=(kB*temp/(ZI*ZI))*(a_p_g[i][1]-a_p_g[i-1][1]);           
        }
	
/* Initial Forces */
        for(i=1;i<=N;i++){
                F_tot_x[i]=F_s_x_sub(i,1,r_ij,u_norm_x) + F_b_x_sub(i,1,beta,r_ij,u_norm_x) + F_EV_x_sub(i,1,rx,r_ij);
                F_tot_y[i]=F_s_y_sub(i,1,r_ij,u_norm_y) + F_b_y_sub(i,1,beta,r_ij,u_norm_y) + F_EV_y_sub(i,1,ry,r_ij);
                F_tot_z[i]=F_s_z_sub(i,1,r_ij,u_norm_z) + F_b_z_sub(i,1,beta,r_ij,u_norm_z) + F_EV_z_sub(i,1,rz,r_ij);
        }

/* Time Evolution */
        for(j=1;j<=Nt;j++)
	{
	/* Calculate New Phi */
                for(i=1;i<=N-2;i++){
                        phi[i][2]=phi[i][1] + (D_rot*dt/(kB*temp))*tau[i] +
                                gasdev(&idum,2.0*D_rot*dt);
                }

	/* Calculate New Positions */
		for(i=1;i<=N;i++){
			rx[i][2]=rx[i][1] + (D_trans*dt/(kB*temp))*F_tot_x[i] +
                                gasdev(&idum,2.0*D_trans*dt);
			ry[i][2]=ry[i][1] + (D_trans*dt/(kB*temp))*F_tot_y[i] +
                                gasdev(&idum,2.0*D_trans*dt);
			rz[i][2]=rz[i][1] + (D_trans*dt/(kB*temp))*F_tot_z[i] +
                                gasdev(&idum,2.0*D_trans*dt);
		}
		
	/* Calculate New Separations */
	for(i=1;i<=N-1;i++){
                for(k=i+1;k<=N;k++){
                        r_ij[i][k]=sqrt(pow(rx[k][2]-rx[i][2],2.0)+
                                pow(ry[k][2]-ry[i][2],2.0)+
                                pow(rz[k][2]-rz[i][2],2.0));
                        r_ij[k][i]=r_ij[i][k];
                }
        }

	/* Calculate New BFCs */
		/* Calculate New u_norm[i] */
		for(i=1;i<=N-1;i++){
			u_norm_x[i][2]=(rx[i+1][2]-rx[i][2])/r_ij[i+1][i];
                        u_norm_y[i][2]=(ry[i+1][2]-ry[i][2])/r_ij[i+1][i];
                        u_norm_z[i][2]=(rz[i+1][2]-rz[i][2])/r_ij[i+1][i];
		}	

		/* Calculate New f_norm[i] */
                        /* calculate new f[i] */
                for(i=1;i<=N-1;i++)
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
                for(i=1;i<=N-1;i++)
                {
                        fx[i]=fx[i]-u_norm_x[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                        fy[i]=fy[i]-u_norm_y[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                        fz[i]=fz[i]-u_norm_z[i][2]*(u_norm_x[i][2]*fx[i]+u_norm_y[i][2]*
                                fy[i]+u_norm_z[i][2]*fz[i]);
                }
                        /* normalize new f[i] */
                for(i=1;i<=N-1;i++)
                {
                        f_norm_x[i]=fx[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                        f_norm_y[i]=fy[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                        f_norm_z[i]=fz[i]/sqrt(pow(fx[i],2.0)+pow(fy[i],2.0)+pow(fz[i],2.0));
                }

		/* Calculate New v_norm[i] */
                for(i=1;i<=N-1;i++)
                {
                        vx[i]=u_norm_y[i][2]*f_norm_z[i]-u_norm_z[i][2]*f_norm_y[i];
                        vy[i]=u_norm_z[i][2]*f_norm_x[i]-u_norm_x[i][2]*f_norm_z[i];
                        vz[i]=u_norm_x[i][2]*f_norm_y[i]-u_norm_y[i][2]*f_norm_x[i];

                        v_norm_x[i]=vx[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                        v_norm_y[i]=vy[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                        v_norm_z[i]=vz[i]/sqrt(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
                }

	/* Calculate New beta[i] */
                for(i=1;i<=N-2;i++)
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
		for(i=1;i<=N-2;i++)
                {
                        a_p_g[i][2]=asin( (v_norm_x[i]*f_norm_x[i+1]+v_norm_y[i]*f_norm_y[i+1]+
                        v_norm_z[i]*f_norm_z[i+1]-f_norm_x[i]*v_norm_x[i+1]-
                        f_norm_y[i]*v_norm_y[i+1]-f_norm_z[i]*v_norm_z[i+1])/
                        (2.0*cos(beta[i]/2.0)*cos(beta[i]/2.0)) );
                }

    	/* Calculate New tau[i] */
                tau[1]=(kB*temp/(ZI*ZI))*a_p_g[1][2];    
                for(i=2;i<=N-2;i++)
                {
                        tau[i]=(kB*temp/(ZI*ZI))*(a_p_g[i][2]-a_p_g[i-1][2]);  
                }

			/* Calculate New F_tot[i]: */
                for(i=1;i<=N;i++){
                        F_tot_x[i]=F_s_x_sub(i,2,r_ij,u_norm_x) + F_b_x_sub(i,2,beta,r_ij,u_norm_x) + F_EV_x_sub(i,2,rx,r_ij);
                        F_tot_y[i]=F_s_y_sub(i,2,r_ij,u_norm_y) + F_b_y_sub(i,2,beta,r_ij,u_norm_y) + F_EV_y_sub(i,2,ry,r_ij);
                        F_tot_z[i]=F_s_z_sub(i,2,r_ij,u_norm_z) + F_b_z_sub(i,2,beta,r_ij,u_norm_z) + F_EV_z_sub(i,2,rz,r_ij);
                }

	/* Calculate New r_f[i] */
                for(i=1;i<=N-1;i++)
                {
                        r_f_x[i]=rx[i][2]+0.5*R_hyd*f_norm_x[i];
                        r_f_y[i]=ry[i][2]+0.5*R_hyd*f_norm_y[i];
                        r_f_z[i]=rz[i][2]+0.5*R_hyd*f_norm_z[i];
                }
		
                if((j+1000)%1000==0){
			fptr1=fopen(filename1,"a");
                	fprintf(fptr1,"\nj=%d\n",j);
        		for(i=1;i<=N;i++){
                		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",rx[i][2]*1e10,
                        		ry[i][2]*1e10,rz[i][2]*1e10);
        		}
			for(i=1;i<=N-1;i++){
                		fprintf(fptr1,"%.3f\t%.3f\t%.3f\n",r_f_x[i]*1e10,
                        		r_f_y[i]*1e10,r_f_z[i]*1e10);
        		}
			fclose(fptr1);
		}

	/* r[i][2] -> r[i][1] for next step: */
                for(i=1;i<=N;i++)
                {
                        rx[i][1]=rx[i][2];
                        ry[i][1]=ry[i][2];
                        rz[i][1]=rz[i][2];
                }

        /* u_norm[i][2] -> u_norm[i][1] for next step: */
                for(i=1;i<=N-1;i++)
                {
                        u_norm_x[i][1]=u_norm_x[i][2];
                        u_norm_y[i][1]=u_norm_y[i][2];
                        u_norm_z[i][1]=u_norm_z[i][2];
                }

	/* phi[i][2] -> phi[i][1] for next step */
                for(i=1;i<=N-2;i++){
                        phi[i][1]=phi[i][2];
                }
	}
}
