/*
  Test code for using the template FFT code.

  Author: Tim Molteno, tim@physics.otago.ac.nz

  Copyright Tim Molteno, 2008-2015.
  Licensed under the GPL v3.
*/

#include "fft_complex.h"
#include <Eigen/Dense>
using namespace std;

#define LOG_LENGTH 13
#define N (1<<LOG_LENGTH)

typedef complex<double> fft_complex;

#include <iostream>

int main()
{
  CFFT<N,double> cfft;
  Matrix<complex<double>, N, 1> cdata;
  cdata.setRandom();

  std::cout << "Perform 1000 FFT and INVERSE FFT " << endl;
  std::cout << "FFT Length= " << N << endl;
  
  for (int i=0;i<1000;i++) {
    cfft.fft(cdata);
    cfft.ifft(cdata);
  }
  return 0;
}
