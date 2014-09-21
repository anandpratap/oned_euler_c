#ifndef _INCL_BASE
#define _INCL_BASE

#define XMAX 1.0f
#define XMIN 0.0f

#define GAMMA 1.4
#define GAMMA_M 0.4

#include <cmath>
#include <mkl.h>
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


#endif
