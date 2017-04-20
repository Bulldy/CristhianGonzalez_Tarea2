#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define L 1.0
#define rho0_L 1.0
#define rho0_R 0.1
#define p0_L 1.0
#define p0_R 0.1
#define a_R 0.3741657386773941
#define y 1.4
#define M 1001
#define dx 0.001
#define C 0.5

/*---------------------------------------------------------------------------------------*/
/*Functions to be defined and used*/
void w_init(double **w);
void flux(double *F, double *w);
void laxwendroff_step(double *wnext, double *wi,double *wip1, double *wim1,double dt);
double delta_time(double **w);
double energy(double pres,double rho,double vel);
double pressure(double energy, double ru2);
double speed_sound(double pres, double rho);
void write_final_result(double **ww);

/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/
/*START OF MAIN*/
/*We work with dimensionless variables*/
int main(){

  /*Initialize present and future 2D arrays. The array holds the 3D euler equation vectors*/
  double **w_now, **w_after;
  w_now=(double **)malloc(M*sizeof(double *));
  w_after=(double **)malloc(M*sizeof(double *));
  
  int j;
  for(j=0;j<M;j++){
    w_now[j]=(double *)malloc(3*sizeof(double));
    w_after[j]=(double *)malloc(3*sizeof(double));
  }
  
  w_init(w_now);
  
  int k,Nt,n;
  double dtp,time;
  Nt=1100;
  time=0;
  for(k=0;k<Nt;k++){
    dtp=delta_time(w_now);
    time=time+dtp;
    for(j=1;j<M-1;j++){
      laxwendroff_step(w_after[j],w_now[j],w_now[j+1],w_now[j-1],dtp);
    }
    w_after[0]=w_now[0];
    w_after[M-1]=w_now[M-1];

    for(j=0;j<M;j++){
      for(n=0;n<3;n++){
	w_now[j][n]=w_after[j][n];
      }
    }
  }
  
  write_final_result(w_after);

  FILE *t_file;
  t_file=fopen("time_shock.txt","w");
  fprintf(t_file,"%f \n",time);
  fclose(t_file);

  return 0;
}

/* END OF MAIN*/
/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------------*/
/*Function to initialize our w (3D Euler equation vector) given the initial conditions*/
void w_init(double **w){
  double *rho_init,*p_init,*e_init;
  int j;
  for(j=0;j<M;j++){
    double xx;
    xx=j*dx;
    if(xx<=0.5){
      w[j][0]=rho0_L/rho0_R;
      w[j][1]=0.0;
      w[j][2]=energy((p0_L/(y*p0_R)),(rho0_L/rho0_R),0);
    }
    else{
      w[j][0]=1.0;
      w[j][1]=0.0;
      w[j][2]=energy((1.0/y),1.0,0);
    }
  }
}

/*---------------------------------------------------------------------------------------*/
/*Function to return the minimum delta time from the Courant condition*/
double delta_time(double **w){
  int j;
  double rhoj,uj,pj,ej,aj,dt_min,dt_j;
  dt_min=1000;
  for(j=0;j<M;j++){
    rhoj=w[j][0];
    ej=w[j][2];
    uj=w[j][1]/w[j][0];
    pj=pressure(ej,rhoj*uj*uj);
    aj=speed_sound(pj,rhoj);
    dt_j=C*dx/(fabs(uj)+aj);
    if(dt_j<dt_min){
      dt_min=dt_j;
    }
  }
  return dt_min;
}

/*------------------------------------------------------------------------------------*/
/*Function to advance the system using a Lax-Wendroff step and save on w_next*/
void laxwendroff_step(double *w_next, double *wi,double *wip1, double *wim1,double dt){
  double *fi,*fip1,*fim1,*fip12,*fim12;
  double *wip12,*wim12;

  fi=(double *)malloc(3*sizeof(double));
  fip1=(double *)malloc(3*sizeof(double));
  fim1=(double *)malloc(3*sizeof(double));
  fip12=(double *)malloc(3*sizeof(double));
  fim12=(double *)malloc(3*sizeof(double));

  wip12=(double *)malloc(3*sizeof(double));
  wim12=(double *)malloc(3*sizeof(double));

  flux(fi,wi);
  flux(fip1,wip1);
  flux(fim1,wim1);

  int n;
  for(n=0;n<3;n++){
    wip12[n]=0.5*(wip1[n]+wi[n])-0.5*(dt/dx)*(fip1[n]-fi[n]);
    wim12[n]=0.5*(wi[n]+wim1[n])-0.5*(dt/dx)*(fi[n]-fim1[n]);
  }

  flux(fip12,wip12);
  flux(fim12,wim12);

  for(n=0;n<3;n++){
    w_next[n]=wi[n]-(dt/dx)*(fip12[n]-fim12[n]);
  }
}

/*---------------------------------------------------------------------------------------*/
/*Function to save the flux vector on F given a certain w vector from the Euler equation */
void flux(double *F, double *w){
  F[0]=w[1];
  double rhou2;
  rhou2=w[1]*w[1]/w[0];
  F[1]=rhou2+pressure(w[2],rhou2);
  F[2]=(w[1]/w[0])*(w[2]+pressure(w[2],rhou2));
}

/*---------------------------------------------------------------------------------------*/
/*Function to return the energy state given a certain pressure, density and velocity*/
double energy(double pres,double rho,double vel){
  double ee;
  ee=(pres/(y-1.0))+0.5*rho*vel*vel;
  return ee;
}

/*---------------------------------------------------------------------------------------*/
/*Function to return the pressure state given a certain energy, density and velocity*/
double pressure(double energy, double ru2){
  double pp;
  pp=(y-1)*(energy-0.5*ru2);
  return pp;
}

/*---------------------------------------------------------------------------------------*/
/*Function to return the speed of sound given a certain pressure and density*/ 
double speed_sound(double pres, double rho){
  double aa;
  aa=sqrt(y*pres/rho);
  return aa;
}

/*---------------------------------------------------------------------------------------*/
/*Function to write the files with the final results for density, pressure and velocity*/
void write_final_result(double **ww){
  FILE *rho_file,*pres_file,*vel_file;
  rho_file=fopen("density.txt","w");
  pres_file=fopen("pressure.txt","w");
  vel_file=fopen("velocity.txt","w");
  
  double rhoj,pj,uj;
  int j;
  for(j=0;j<M;j++){
    rhoj=ww[j][0];
    uj=ww[j][1]/ww[j][0];
    pj=pressure(ww[j][2],rhoj*uj*uj);

    fprintf(rho_file,"%f \n",rhoj);
    fprintf(pres_file,"%f \n",pj);
    fprintf(vel_file,"%f \n",uj);
  }
  
  fclose(rho_file);
  fclose(pres_file);
  fclose(vel_file);
}
