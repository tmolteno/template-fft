# Template-FFT

Fast Header only C++ FFT and Inverse FFT using templates.

Here is the source code for a complex FFT algorithm that uses template metaprogramming. 
It is quite fast (nearly as fast as fftw) and can be specialized to any type. We will be creating a fixed-point friendly FFT as well in the near future.

## Using the C++ template

Create the file below, and compile it with your favourite C++ compiler. Requires Eigen for handling vectors and matrices. 

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


## Copyright

Author: Tim Molteno, tim@physics.otago.ac.nz

Based on the article "A Simple and Efficient FFT Implementation in C++" by Volodymyr Myrnyy
with an inverse FFT modification.

Licensed under the GPL v3.