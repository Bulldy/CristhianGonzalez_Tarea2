#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define M 128
#define dx 2.0/256.0
#define dy 2.0/256.0
#define dz 2.0/256.0
#define g 1.4 
#define C 0.5
#define R_air 286.9
#define P_atm 101325.0
#define T0 300.0
#define E_BOOM 10000000000.0

/*---------------------------------------------------------------------------------------*/
/*Units of the problem
 *Constant R only air 286.9 J/(kg*K)
 *Distance are measure in meters
 *Heat capacity ratio gamma = cp/cv = 1.4 for diatomic gas no units
 *Courant number C 0.5 no units
 *Pressure measure in atm 
 *Temperature measure in Kelvins
 *Energy measure in Joules
*/
/*---------------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------------*/
/*Functions to be defined and used*/
void finite_volume(double *wnext,double *wnow, double *wip1, double *wim1, double *wjp1, double *wjm1, double *wkp1, double *wkm1, double dete);
void flux(double *F, double *w);
void glux(double *G, double *w);
void hlux(double *H, double *w);
double delta_time(double ****w);
double pressure(double energy, double rhoo,double u2);
double speed_sound(double pres, double rho);

/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/
/*START OF MAIN*/
int main(){
  // We find density and speed of sound at atmospheric conditions with the ideal gas law
  // We find the energy of atmospheric state
  double rho_atm,e_atm,a_atm;
  rho_atm=P_atm/(R_air*T0);
  e_atm=P_atm/(g-1);
  a_atm=speed_sound(P_atm,rho_atm);
  
  // We turn our initial conditions adimensional
  double Arho_atm,AE_BOOM,Ae_atm;
  Arho_atm=1;
  Ae_atm=e_atm/(rho_atm*a_atm*a_atm);
  AE_BOOM=E_BOOM/(rho_atm*a_atm*a_atm);

  // We initialize our 4d array,indexed by i,j,k for position in space
  // Each i,j,k will hold a 5-component vector
  // We introduce the atmospheric values and the center explosion
  double ****w_now,****w_after;
  w_now=(double ****)malloc(M*sizeof(double ***));
  w_after=(double ****)malloc(M*sizeof(double ***));

  int i,j,k;
  for(i=0;i<M;i++){
    w_now[i]=(double ***)malloc(M*sizeof(double **));
    w_after[i]=(double ***)malloc(M*sizeof(double **));
    for(j=0;j<M;j++){
      w_now[i][j]=(double **)malloc(M*sizeof(double *));
      w_after[i][j]=(double **)malloc(M*sizeof(double *));
      for(k=0;k<M;k++){
	w_now[i][j][k]=(double *)malloc(5*sizeof(double));
	w_after[i][j][k]=(double *)malloc(5*sizeof(double));
	w_now[i][j][k][0]=Arho_atm;
	w_now[i][j][k][1]=0.0;
	w_now[i][j][k][2]=0.0;
	w_now[i][j][k][3]=0.0;
	w_now[i][j][k][4]=e_atm;
	if(i==M/2 && j==M/2 && k==M/2){
	  w_now[i][j][k][4]=AE_BOOM;
	}
      }
    }
  } 
  
  // We advance our system in timesteps of dt calculated at each step
  int Nt,n,l;
  Nt=100;
  double dt,t_total;
  t_total=0;
  for(n=0;n<Nt;n++){
    dt=delta_time(w_now);
    printf("%f y %f\n",dt,t_total);
    t_total=t_total+dt;
    // We fix boundary values at atmospheric conditions
    for(i=0;i<M;i++){
      for(j=0;j<M;j++){
	w_after[i][j][0]=w_now[i][j][0];
	w_after[i][j][M-1]=w_now[i][j][M-1];
	for(k=0;k<M;k++){
	  w_after[0][j][k]=w_now[0][j][k];
	  w_after[M-1][j][k]=w_now[M-1][j][k];
	  w_after[i][0][k]=w_now[i][0][k];
	  w_after[i][M-1][k]=w_now[i][M-1][k];
	}    
      }
    }

    // We calculate the new values inside the cube
    for(i=1;i<M-1;i++){
      for(j=1;j<M-1;j++){
        for(k=1;k<M-1;k++){
          finite_volume(w_after[i][j][k],w_now[i][j][k], w_now[i+1][j][k], w_now[i-1][j][k], w_now[i][j+1][k], w_now[i][j-1][k], w_now[i][j][k+1], w_now[i][j][k-1], dt);
        }
      }
    }

    // We update the system for the next time step
    for(i=0;i<M;i++){
      for(j=0;j<M;j++){
	for(k=0;k<M;k++){
	  for(l=0;l<5;l++){
	    w_now[i][j][k][l]=w_after[i][j][k][l];
	  }
	}
      }
    }
    //printf("%f\n",w_after[94][M/2][M/2][1]);
    // Time step ends
  }
  printf("%f\n",w_after[M/2+1][M/2][M/2][4]);
  return 0;
}

/* END OF MAIN*/
/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------------------*/
/*Function to return the minimum delta time from the Courant condition*/
double delta_time(double ****w){
  int i, j, k;
  double rhoijk,uijk,vijk,wijk,pijk,eijk,aijk,dt_min,dtu,dtv,dtw;
  dt_min=10;
  for(i=0;i<M;i++){
    for(j=0;j<M;j++){
      for(k=0;k<M;k++){
        rhoijk=w[i][j][k][0];
	eijk=w[i][j][k][4];
	uijk=w[i][j][k][1]/rhoijk;
	vijk=w[i][j][k][2]/rhoijk;
	wijk=w[i][j][k][3]/rhoijk;
	pijk=pressure(eijk,rhoijk,uijk*uijk+vijk*vijk+wijk*wijk);
	aijk=speed_sound(pijk,rhoijk);
	dtu=C*dx/(fabs(uijk)+aijk);
	dtv=C*dy/(fabs(vijk)+aijk);
	dtw=C*dz/(fabs(wijk)+aijk);
	if(dtu<dt_min){
	  dt_min=dtu;
	}
	if(dtv<dt_min){
          dt_min=dtv;
	}
	if(dtw<dt_min){
          dt_min=dtu;
	}
      }
    }
  }
  return dt_min;
}

/*---------------------------------------------------------------------------------------*/
/*Function to advance the system in time using a finite volume scheme*/
/*We use a central difference approximation for half values*/
void finite_volume(double *wnext,double *wnow, double *wip1, double *wim1, double *wjp1, double *wjm1, double *wkp1, double *wkm1, double dete){
  int l;
  double *fip12,*fim12,*gjp12,*gjm12,*hkp12,*hkm12;
  double *wip12,*wim12,*wjp12,*wjm12,*wkp12,*wkm12;

  fip12=(double *)malloc(5*sizeof(double));
  fim12=(double *)malloc(5*sizeof(double));
  gjp12=(double *)malloc(5*sizeof(double));
  gjm12=(double *)malloc(5*sizeof(double));
  hkp12=(double *)malloc(5*sizeof(double));
  hkm12=(double *)malloc(5*sizeof(double));

  wip12=(double *)malloc(5*sizeof(double));
  wim12=(double *)malloc(5*sizeof(double));
  wjp12=(double *)malloc(5*sizeof(double));
  wjm12=(double *)malloc(5*sizeof(double));
  wkp12=(double *)malloc(5*sizeof(double));
  wkm12=(double *)malloc(5*sizeof(double));

  for(l=0;l<5;l++){
    wip12[l]=0.5*(wnow[l]+wip1[l]);
    wim12[l]=0.5*(wnow[l]+wim1[l]);
    wjp12[l]=0.5*(wnow[l]+wjp1[l]);
    wjm12[l]=0.5*(wnow[l]+wjm1[l]);
    wkp12[l]=0.5*(wnow[l]+wkp1[l]);
    wkm12[l]=0.5*(wnow[l]+wkm1[l]);
  }

  flux(fip12,wip12);
  flux(fim12,wim12);
  glux(gjp12,wjp12);
  glux(gjm12,wjm12);
  hlux(hkp12,wkp12);
  hlux(hkm12,wkm12);
  
  for(l=0;l<5;l++){
    wnext[l]=wnow[l]-(dete/dx)*(fip12[l]-fim12[l])-(dete/dy)*(gjp12[l]-gjm12[l])-(dete/dz)*(hkp12[l]-hkm12[l]);
  }

  free(wip12);
  free(wim12);
  free(wjp12);
  free(wjm12);
  free(wkp12);
  free(wkm12);
  free(fip12);
  free(fim12);
  free(gjp12);
  free(gjm12);
  free(hkp12);
  free(hkm12);
}

/*---------------------------------------------------------------------------------------*/
/*Function to save the flux vector on F given a certain w vector from the Euler equation */
/* Flux on x direction */
void flux(double *F, double *w){
  F[0]=w[1];
  double p,v2;
  v2=w[1]*w[1]/(w[0]*w[0])+w[2]*w[2]/(w[0]*w[0])+w[3]*w[3]/(w[0]*w[0]);
  p=pressure(w[4],w[0],v2);
  F[1]=w[1]*w[1]/w[0]+p;
  F[2]=w[1]*w[2]/w[0];
  F[3]=w[1]*w[3]/w[0];
  F[4]=(w[4]+p)*(w[1]/w[0]);
}

/*---------------------------------------------------------------------------------------*/
/*Function to save the flux vector on G given a certain w vector from the Euler equation */
/* Flux on y direction */
void glux(double *G, double *w){
  G[0]=w[2];
  G[1]=w[1]*w[2]/w[0];
  double p,v2;
  v2=w[1]*w[1]/(w[0]*w[0])+w[2]*w[2]/(w[0]*w[0])+w[3]*w[3]/(w[0]*w[0]);
  p=pressure(w[4],w[0],v2);
  G[2]=w[2]*w[2]/w[0]+p;
  G[3]=w[2]*w[3]/w[0];
  G[4]=(w[4]+p)*(w[2]/w[0]);
}

/*---------------------------------------------------------------------------------------*/
/*Function to save the flux vector on H given a certain w vector from the Euler equation */
/* Flux on z direction */
void hlux(double *H, double *w){
  H[0]=w[3];
  H[1]=w[1]*w[3]/w[0];
  double p,v2;
  v2=w[1]*w[1]/(w[0]*w[0])+w[2]*w[2]/(w[0]*w[0])+w[3]*w[3]/(w[0]*w[0]);
  p=pressure(w[4],w[0],v2);
  H[2]=w[2]*w[3]/w[0];
  H[3]=w[3]*w[3]/w[0]+p;
  H[4]=(w[4]+p)*(w[3]/w[0]);
}

/*---------------------------------------------------------------------------------------*/
/*Function to return the pressure state given a certain energy, density and velocity*/
double pressure(double energy, double rhoo,double u2){
  double pp;
  pp=(g-1)*(energy-0.5*rhoo*u2);
  return pp;
}

/*---------------------------------------------------------------------------------------*/
/*Function to return the speed of sound given a certain pressure and density*/
double speed_sound(double pres, double rho){
  double aa;
  aa=sqrt(g*pres/rho);
  return aa;
}
