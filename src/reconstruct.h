#ifndef _INCL_RECONSTRUCT
#define _INCL_RECONSTRUCT

#include "base.h"

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
#endif
