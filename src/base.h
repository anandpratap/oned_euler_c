#ifndef _INCL_BASE
#define _INCL_BASE

#define XMAX 10.0f
#define XMIN 0.0f

#define GAMMA 1.4
#define GAMMA_M 0.4

struct state{
  double rho, u, p;
};

void init(int N, double x[N+1], double xc[N], double U[3][N], struct state left, struct state right);
void bc(int N, double U[3][N], double rho[N+2], double u[N+2], double p[N+2], int start, int length);

void reconstruct(int N, double rho[N+2], double u[N+2], double p[N+2], 
		 double rhol[N+1], double rhor[N+1], double ul[N+1], 
		 double ur[N+1], double pl[N+1], double pr[N+1], int start, int length);

__declspec(vector) void roeflux(double rhol, double rhor, double ul, double ur, double pl, double pr,
				double *f_0, double *f_1, double *f_2);

#endif
