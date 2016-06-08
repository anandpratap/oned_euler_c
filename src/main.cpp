#include <stdlib.h>
#include "utils.h"

int main(int argc, char **argv){

	int Ng, N;
	struct state left, right;
	Ng = 21;
	left.rho = 1.0;
	left.u = 0.0;
	left.p = 1.0;
  
	right.rho = 0.125;
	right.u = 0.0;
	right.p = 0.1;
  
	N = Ng - 1;

	double *x = new double[Ng]();
	double *xc = new double[N]();
	double **U = allocate_2d_array<double>(3, N);
	init(N, x, xc, U, left, right);
	double dt, tf, t;
  
	double *rho = new double[N+2]();
	double *u = new double[N+2]();
	double *p = new double[N+2]();
  
	double *rhol = new double[N+1]();
	double *rhor = new double[N+1]();
	double *ul = new double[N+1]();
	double *ur = new double[N+1]();
	double *pl = new double[N+1]();
	double *pr = new double[N+1]();

	double *f_0 = new double[N+1]();
	double *f_1 = new double[N+1]();
	double *f_2 = new double[N+1]();

	double *fac = new double[N+1]();

	dt = 1e-4;
	t = 0.0;
	tf = 0.1;
	for(int i=0; i<N; i++)
		fac[i] = dt/(x[1] - x[0]);


	while(t<tf){
		bc(N, U, rho, u, p);
		reconstruct(N, rho, u, p, rhol, rhor, ul, ur, pl, pr);
		for(int i=0; i<N; i++){
			roeflux(rhol[i], rhor[i], ul[i], ur[i],
					pl[i], pr[i], &f_0[i], &f_1[i], &f_2[i]);
		}
		update(N, U, f_0, f_1, f_2, fac);
		t += dt;
		printf("t = %.4f\n", t);
	}
	write_solution(N, xc, U);
  
	delete[] x;
	delete[] xc;
	delete[] rho;
	delete[] rhol;
	delete[] rhor;
	delete[] u;
	delete[] ul;
	delete[] ur;
	delete[] p;
	delete[] pl;
	delete[] pr;
	delete[] f_0;
	delete[] f_1;
	delete[] f_2;
	delete[] fac;
	release_2d_array(U, 3, N); 
	return 0;
  
}
