#include <cmath>
#include <mkl.h>
#include "base.h"
#include <iostream>
void init(int N, double x[N+1], double xc[N], double U[3][N], struct state left, struct state right){
  double dx = (XMAX - XMIN)/N;
  
  // calculate node points
  for(int i=0; i <= N; i++){
    x[i] = i*dx;
  }
  
  // cell center
  xc[0:N] = 0.5f*(x[1:N] + x[0:N]);

  U[0][0:N/2] = left.rho;
  U[0][N/2:N/2] = right.rho;

  U[1][0:N/2] = left.rho*left.u;
  U[1][N/2:N/2] = right.rho*right.u;

  U[2][0:N/2] = left.p/(GAMMA-1.0f) + 0.5f*left.rho*left.u*left.u;
  U[2][N/2:N/2] = right.p/(GAMMA-1.0f) + 0.5f*right.rho*right.u*right.u;
}

void bc(int N, double U[3][N], double rho[N+2], double u[N+2], double p[N+2], int start, int length){
  if(start + length > N){
    length = N - start;
  }
  rho[start+1:length] = U[0][start:length];
  u[start+1:length] = U[1][start:length]/U[0][start:length];
  p[start+1:length] = GAMMA_M*(U[2][start:length] - 0.5f*(rho[start+1:length]*u[start+1:length]*u[start+1:length]));
  if(start == 0){
    rho[0] = rho[1];
    u[0] = u[1];
    p[0] = p[1];
    
    rho[N+1] = rho[N];
    u[N+1] = u[N];
    p[N+1] = p[N];
  }
}

void reconstruct(int N, double rho[N+2], double u[N+2], double p[N+2], 
		 double rhol[N+1], double rhor[N+1], double ul[N+1], 
		 double ur[N+1], double pl[N+1], double pr[N+1], int start, int length){
  rhol[start:length] = rho[start:length];
  rhor[start:length] = rho[start+1:length];

  ul[start:length] = u[start:length];
  ur[start:length] = u[start+1:length];

  pl[start:length] = p[start:length];
  pr[start:length] = p[start+1:length];
}

__declspec(vector) void roeflux(double rhol, double rhor, double ul, double ur, double pl, double pr,
				double *f_0, double *f_1, double *f_2){
  double el = pl/(GAMMA_M) + 0.5f*rhol*ul*ul;
  double er = pr/(GAMMA_M) + 0.5f*rhor*ur*ur;
  double hl = (el + pl)/rhol;
  double hr = (er + pr)/rhor;

  double sqrtrhol = sqrt(rhol);
  double sqrtrhor = sqrt(rhor);
  double den_inverse = 1/(sqrtrhol + sqrtrhor);
  double uavg = (sqrtrhol*ul + sqrtrhor*ur)*den_inverse;
  double havg = (sqrtrhol*hl + sqrtrhor*hr)*den_inverse;
  double cavg = sqrt(GAMMA_M*(havg - 0.5f*uavg*uavg));
  double cavg_inverse = 1.0f/cavg;

  double d1 = rhor - rhol;
  double d2 = rhor*ur - rhol*ul;
  double d3 = er - el;
  
  double alpha_2 = GAMMA_M*((havg - uavg*uavg)*d1 + uavg*d2 - d3)*cavg_inverse*cavg_inverse;
  double alpha_3 = 0.5f*(d2 + (cavg - uavg)*d1 - cavg*alpha_2)*cavg_inverse;
  double alpha_1 = d1 - alpha_2 - alpha_3;
  
  double lambda_1 =  fabs(uavg - cavg);
  double lambda_2 =  fabs(uavg);
  double lambda_3 =  fabs(uavg + cavg);

  double f1 = lambda_1*alpha_1 + lambda_2*alpha_2 + lambda_3*alpha_3;
  double f2 = lambda_1*alpha_1*(uavg-cavg) + lambda_2*alpha_2*uavg + lambda_3*alpha_3*(uavg+cavg);
  double f3 = lambda_1*alpha_1*(havg-cavg*uavg) + 0.5f*lambda_2*alpha_2*uavg*uavg + lambda_3*alpha_3*(havg+cavg*uavg);

  *f_0 = 0.5f*((rhol*ul + rhor*ur) - f1);
  *f_1 = 0.5f*((rhol*ul*ul + pl + rhor*ur*ur + pr) - f2);
  *f_2 = 0.5f*(ul*hl*rhol + ur*hr*rhor - f3);
}
