#ifndef __UTILS_H
#define __UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include "math.h"

#define GAMMA 1.4
#define GAMMA_M 0.4
#define XMAX 1.0
#define XMIN 0.0

template <class T>
T** allocate_2d_array(int nx, int ny){
    T** A = new T*[nx];
    for(int i(0); i < nx; ++i){
		A[i] = new T[ny];
	}
    return A;
}

template <class T>
void release_2d_array(T** A, int nx, int ny){
    for (int i = 0; i < nx; ++i){
		delete[] A[i];
	}
    delete[] A;
}

struct state{
	double rho, u, p;
};

template <class dtype, class xtype>
	void init(int N, xtype *x, xtype *xc, dtype *U[3], struct state left, struct state right){
	xtype dx = (XMAX - XMIN)/N;
	
	// calculate node points
	for(int i=0; i <= N; i++){
		x[i] = i*dx;
	}
	
	// cell center
	for(int i=0; i < N; i++){
		xc[i] = 0.5f*(x[i+1] + x[i]);
	}

	for(int i=0; i < N; i++){
		if(i < N/2){
			U[0][i] = left.rho;
			U[1][i] = left.rho*left.u;
			U[2][i] = left.p/(GAMMA-1.0f) + 0.5f*left.rho*left.u*left.u;
		}
		else{
			U[0][i] = right.rho;
			U[1][i] = right.rho*right.u;
			U[2][i] = right.p/(GAMMA-1.0f) + 0.5f*right.rho*right.u*right.u;
		}
	}
}

template <class dtype>
void reconstruct(int N, dtype *rho, dtype *u, dtype *p, 
				 dtype *rhol, dtype *rhor, dtype *ul, 
				 dtype *ur, dtype *pl, dtype *pr){
	for(int i=0; i<N+1; i++){
		rhol[i] = rho[i];
		rhor[i] = rho[i+1];
		
		ul[i] = u[i];
		ur[i] = u[i+1];
		
		pl[i] = p[i];
		pr[i] = p[i+1];
	}
}

template <class dtype>
void roeflux(dtype rhol, dtype rhor, dtype ul, dtype ur, dtype pl, dtype pr,
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

template <class dtype>
void bc(int N, dtype *U[3], dtype *rho, dtype *u, dtype *p){
	for(int i=0; i<N; i++){
		rho[i+1] = U[0][i];
		u[i+1] = U[1][i]/U[0][i];
		p[i+1] = GAMMA_M*(U[2][i] - 0.5*(rho[i+1]*u[i+1]*u[i+1]));
	}
	
	rho[0] = rho[1];
    u[0] = u[1];
    p[0] = p[1];
    
    rho[N+1] = rho[N];
    u[N+1] = u[N];
    p[N+1] = p[N];
}

template <class dtype, class xtype>
	void update(int N, dtype *U[3], dtype *f_0, dtype *f_1, dtype *f_2, xtype *fac){     
	
	for(int i=0; i<N; i++){  
		U[0][i] += -(f_0[i+1] - f_0[i])*fac[i];
		U[1][i] += -(f_1[i+1] - f_1[i])*fac[i];
		U[2][i] += -(f_2[i+1] - f_2[i])*fac[i];
	}
	U[0][0] = U[0][1];
	U[1][0] = U[1][1];
	U[2][0] = U[2][1];

	U[0][N-1] = U[0][N-2];
	U[1][N-1] = U[1][N-2];
	U[2][N-1] = U[2][N-2];
}

template <class dtype, class xtype>
	void write_solution(int N, xtype *xc, dtype *U[3]){
	FILE *f;
	f = fopen("output.dat", "w");
	
	fprintf(f, "xc\trho\tu\tp\n");
    for(int i = 1; i < N-1; i++){
		double rho = U[0][i];
		double u = U[1][i]/rho;
		double p = GAMMA_M*(U[2][i] - 0.5f*rho*u*u);
		fprintf(f, "%.10f\t%.10f\t%.10f\t%.10f\n", xc[i], rho, u, p);
	}
	fclose(f);
}

#endif
