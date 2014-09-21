#ifndef _INCL_BASE
#define _INCL_BASE

#define XMAX 10.0f
#define XMIN 0.0f

#define GAMMA 1.4
#define GAMMA_M 0.4

#include <cmath>
#include <mkl.h>
#include "base.h"
#include <iostream>


struct state{
  double rho, u, p;
};

template <class dtype, class xtype>
void init(unsigned int N, xtype x[N+1], xtype xc[N], dtype U[3][N], struct state left, struct state right){
  xtype dx = (XMAX - XMIN)/N;
  
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

template <class dtype>
void bc(unsigned int N, dtype U[3][N], dtype rho[N+2], dtype u[N+2], dtype p[N+2], unsigned int start, unsigned int length){
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

template <class dtype>
void reconstruct(unsigned int N, dtype rho[N+2], dtype u[N+2], dtype p[N+2], 
		 dtype rhol[N+1], dtype rhor[N+1], dtype ul[N+1], 
		 dtype ur[N+1], dtype pl[N+1], dtype pr[N+1], unsigned int start, unsigned int length){
  rhol[start:length] = rho[start:length];
  rhor[start:length] = rho[start+1:length];

  ul[start:length] = u[start:length];
  ur[start:length] = u[start+1:length];

  pl[start:length] = p[start:length];
  pr[start:length] = p[start+1:length];
}

template <class dtype>
__declspec(vector) void roeflux(dtype rhol, dtype rhor, dtype ul, dtype ur, dtype pl, dtype pr,
				dtype *f_0, dtype *f_1, dtype *f_2){
  dtype el = pl/(GAMMA_M) + 0.5f*rhol*ul*ul;
  dtype er = pr/(GAMMA_M) + 0.5f*rhor*ur*ur;
  dtype hl = (el + pl)/rhol;
  dtype hr = (er + pr)/rhor;

  dtype sqrtrhol = sqrt(rhol);
  dtype sqrtrhor = sqrt(rhor);
  dtype den_inverse = 1/(sqrtrhol + sqrtrhor);
  dtype uavg = (sqrtrhol*ul + sqrtrhor*ur)*den_inverse;
  dtype havg = (sqrtrhol*hl + sqrtrhor*hr)*den_inverse;
  dtype cavg = sqrt(GAMMA_M*(havg - 0.5f*uavg*uavg));
  dtype cavg_inverse = 1.0f/cavg;

  dtype d1 = rhor - rhol;
  dtype d2 = rhor*ur - rhol*ul;
  dtype d3 = er - el;
  
  dtype alpha_2 = GAMMA_M*((havg - uavg*uavg)*d1 + uavg*d2 - d3)*cavg_inverse*cavg_inverse;
  dtype alpha_3 = 0.5f*(d2 + (cavg - uavg)*d1 - cavg*alpha_2)*cavg_inverse;
  dtype alpha_1 = d1 - alpha_2 - alpha_3;
  
  dtype lambda_1 =  fabs(uavg - cavg);
  dtype lambda_2 =  fabs(uavg);
  dtype lambda_3 =  fabs(uavg + cavg);

  dtype f1 = lambda_1*alpha_1 + lambda_2*alpha_2 + lambda_3*alpha_3;
  dtype f2 = lambda_1*alpha_1*(uavg-cavg) + lambda_2*alpha_2*uavg + lambda_3*alpha_3*(uavg+cavg);
  dtype f3 = lambda_1*alpha_1*(havg-cavg*uavg) + 0.5f*lambda_2*alpha_2*uavg*uavg + lambda_3*alpha_3*(havg+cavg*uavg);

  *f_0 = 0.5f*((rhol*ul + rhor*ur) - f1);
  *f_1 = 0.5f*((rhol*ul*ul + pl + rhor*ur*ur + pr) - f2);
  *f_2 = 0.5f*(ul*hl*rhol + ur*hr*rhor - f3);
}

#endif
