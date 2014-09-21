#include <iostream>
#include <fstream>
#include <iomanip>
#include "base.h"
#include "io.h"

void read_inputs(struct state *left, struct state *right, int *N){
  std::fstream inpfile("st.inp", std::ios_base::in);
  inpfile >> *N;
  inpfile >> (*left).rho;  
  inpfile >> (*left).u;
  inpfile >> (*left).p;
  inpfile >> (*right).rho;  
  inpfile >> (*right).u;
  inpfile >> (*right).p;
}


void write_solution(int N, double xc[N], double U[3][N]){
  std::ofstream outfile;
  outfile.open("results.txt");
  outfile << std::fixed << std::setprecision(20);
  for(int i = 0; i < N; i++){
    outfile << xc[i] <<" "<<U[0][i]<<" "<<U[1][i]<<" "<<U[2][i]<<"\n";
  }
  outfile.close();
}
