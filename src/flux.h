#ifndef _INCL_FLUX
#define _INCL_FLUX

#include "base.h"

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
