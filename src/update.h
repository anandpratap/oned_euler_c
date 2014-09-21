#ifndef _INCL_UPDATE
#define _INCL_UPDATE

template <class dtype, class xtype>
  void update(unsigned int N, dtype U[3][N], dtype f_0[N+1], dtype f_1[N+1], dtype f_2[N+1], xtype fac[N], unsigned int start, unsigned int length){     

  if(start + length >= N){
    length = N - start;
  }
  
  U[0][start:length] += -(f_0[start+1:length] - f_0[start:length])*fac[start:length];
  U[1][start:length] += -(f_1[start+1:length] - f_1[start:length])*fac[start:length];
  U[2][start:length] += -(f_2[start+1:length] - f_2[start:length])*fac[start:length];

}
#endif
