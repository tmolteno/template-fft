/* Test code for the FFT 
 * Author: Tim Molteno, tim@physics.otago.ac.nz
 *
 * Copyright Tim Molteno, 2008-2015.
 * Licensed under the GPL v3.
 * 
 * Perform 1000 FFT followed by Inverse FFT.
 * Print out the calculated number of Megaflops from the formula
 *     Mflops = 5 N log2(N) / t (us)
 * 
 * gcc -O2 -ftree-vectorize fft_noeigen.cpp -lstdc++
 * Real FFT: MFLOPS(N=16384)=1638.4
 * 
 * Real FFT: MFLOPS(N=16384)=1650.19 --march=native --mtune=native
 */
#include "fft_real.h"
#include <iostream>
#include <iomanip>
using namespace std;

#define LOG_LENGTH 14
#define N (1<<LOG_LENGTH)

#define TEST_LOOP 1000
/*
    gcc -O2 -ftree-vectorize fft_tb.cpp -lstdc++
*/


#include <sys/times.h>
#include <unistd.h>
double secnds(void) {
  struct tms buffer;
  
  times(&buffer);
  return 1000.0 * ( (double)(buffer.tms_utime + buffer.tms_stime) ) /
    ( (double) sysconf(_SC_CLK_TCK) );
}

typedef double fft_real;

int main()
{
  FFT<LOG_LENGTH, double> gfft;
  
  unsigned long n = 1<<LOG_LENGTH;
  double* data = new double [2*n];

  {
    FFT<LOG_LENGTH,fft_real> gfft;

    fft_real* data = new fft_real [2*N];
    for (unsigned i=0; i<N; ++i) {
      data[2*i] = 2*i;
      data[2*i+1] = 2*i+1;
    }

    double start_timer, stop_timer;
    start_timer = secnds();

    for (unsigned i=0;i<TEST_LOOP;i++)  {
      gfft.fft(data);
      gfft.ifft(data);
    }

    stop_timer = secnds();
    stop_timer -= start_timer;
    double usec = stop_timer * 1000.0 / (2.0 * TEST_LOOP);
    double Mflops = 5.0*N*LOG_LENGTH / usec;

    cout << "Real FFT: MFLOPS(N=" << N << ")=" << Mflops << endl;

    if (false)  {
      cout<<"-------Real FFT---------" << endl;
      for (unsigned i=0; i<N; ++i)
        cout << setw(10) << setprecision(5) << data[2*i] << "\t"
            << data[2*i+1] << "I" << endl;
    }

    delete [] data;
  }

  return 0;
}
