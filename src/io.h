#ifndef _INCL_IO
#define _INCL_IO

#include "base.h"
#include <iostream>
#include <fstream>
#include <iomanip>

void read_inputs(struct state *left, struct state *right, unsigned int *N){
  std::fstream inpfile("st.inp", std::ios_base::in);
  inpfile >> *N;
  inpfile >> (*left).rho;  
  inpfile >> (*left).u;
  inpfile >> (*left).p;
  inpfile >> (*right).rho;  
  inpfile >> (*right).u;
  inpfile >> (*right).p;
}


template <class dtype, class xtype>
  void write_solution(unsigned int N, xtype xc[N], dtype U[3][N]){
  std::ofstream outfile;
  outfile.open("results.txt");
  outfile << std::fixed << std::setprecision(20);
  for(int i = 0; i < N; i++){
    outfile << xc[i] <<" "<<U[0][i]<<" "<<U[1][i]<<" "<<U[2][i]<<"\n";
  }
  outfile.close();
}
#endif
