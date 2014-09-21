#include <iostream>
#include <fstream>
#include "omp.h"
#include <sys/time.h>

#include "base.h"
#include "io.h"
inline double my_clock(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (1.0e-6*t.tv_usec + t.tv_sec);
}
int main(int argc, char **argv){

  int Ng, N;
  int input;
  input = atoi(argv[1]);
  struct state left, right;
  read_inputs(&left, &right, &Ng);
  N = Ng - 1;

  double x[Ng], xc[N];
  double U[3][N];
  init(N, x, xc, U, left, right);
  double dt, tf, t;
  double rho[N+2], u[N+2], p[N+2];

  double rhol[N+1], rhor[N+1], ul[N+1], ur[N+1], pl[N+1], pr[N+1];

  double f_0[N+1], f_1[N+1], f_2[N+1];
  double res[3][N];
  dt = 0.1*(x[1]-x[0]);
  t = 0.0;
  tf = 2.2;
  int start = 0;
  int length = N+1;
  int remainder;
  int thread_id;
  int thread_count;
  double start_time, end_time;
  start_time = my_clock();
#pragma omp parallel num_threads(input) private(thread_id, thread_count, start, length, remainder, t)
  {
    thread_id = omp_get_thread_num();
    thread_count = omp_get_num_threads();

    for(int i=0; i< 10000; i++){
      length = (N + 1) / thread_count;
      remainder = (N + 1) % thread_count;
      if(thread_id < remainder){
	length += 1;
	start = thread_id*(length);
      }
      else{
	start = thread_id*length + remainder;
      }
#pragma omp barrier
      bc(N, U, rho, u, p, start, length);
#pragma omp barrier

      reconstruct(N, rho, u, p, rhol, rhor, ul, ur, pl, pr, start, length);
      roeflux(rhol[start:length], rhor[start:length], ul[start:length], ur[start:length],
	      pl[start:length], pr[start:length], &f_0[start:length], &f_1[start:length], &f_2[start:length]);
      
#pragma omp barrier
      if(start + length >= N){
	length = N - start;
      }
      U[0][start:length] += -(f_0[start+1:length] - f_0[start:length])/(x[start+1:length] - x[start:length])*dt;
      U[1][start:length] += -(f_1[start+1:length] - f_1[start:length])/(x[start+1:length] - x[start:length])*dt;
      U[2][start:length] += -(f_2[start+1:length] - f_2[start:length])/(x[start+1:length] - x[start:length])*dt;
      t += dt;
    }
  }
  end_time = my_clock();
  printf("time is %lf\n", end_time-start_time);
  write_solution(N, xc, U);
  return 0;
  
}
