#include <iostream>
#include <fstream>
#include <iomanip>
#include "base.h"
#include "io.h"

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
