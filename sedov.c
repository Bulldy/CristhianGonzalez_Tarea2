#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define M 128
#define dx 2.0
#define dy 2.0
#define dz 2.0
#define g 1.4 
#define C 0.5
#define R_air 286.9
#define P_atm 101325.0
#define T0 300.0
#define E_BOOM 10000000000.0

/*---------------------------------------------------------------------------------------*/
/*Functions to be defined and used*/
double delta_time(double ****w);
double pressure(double energy, double rhoo,double u2);
double speed_sound(double pres, double rho);

/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/
/*START OF MAIN*/
int main(){
  // We find density and speed of sound at atmospheric conditions with the ideal gas law
  // We find the energy of atmospheric state
  double rho_atm,e_atm;
  rho_atm=P_atm/(R_air*T0);
  e_atm=P_atm/(g-1);

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
	w_now[i][j][k][0]=rho_atm;
	w_now[i][j][k][1]=0.0;
	w_now[i][j][k][2]=0.0;
	w_now[i][j][k][3]=0.0;
	w_now[i][j][k][4]=e_atm;
	if(i==M/2 && j==M/2 && k==M/2){
	  w_now[i][j][k][4]=E_BOOM;
	}
      }
    }
  } 
  
  printf("Inicializo el mierdero\n");
  printf("%f\n",w_now[M/2][M/2][M/2][4]);
  
  // We advance our system in timesteps of dt calculated at each step
  int Nt,n,time,l;
  Nt=1;
  double dt;
  for(n=0;n<Nt;n++){
    dt=delta_time(w_now);
    printf("%f\n",dt);
    
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
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //Aca falta la evolucion temporal 
	  //Deberiamos hacer una funcion que tome como argumentos
	  //los 6 vecinos directos y el dt
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

    // Time step ends
  }

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
  dt_min=1000;
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
