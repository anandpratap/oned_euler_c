#ifndef _INCL_BC
#define _INCL_BC
#include "base.h"

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

#endif
