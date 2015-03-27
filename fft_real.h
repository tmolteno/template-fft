#ifndef _fft_h_
#define _fft_h_
/*
  A classy FFT and Inverse FFT C++ class library

  Author: Tim Molteno, tim@physics.otago.ac.nz

  Based on the article "A Simple and Efficient FFT Implementation in C++" by Volodymyr Myrnyy
  with just a simple Inverse FFT modification.

  Licensed under the GPL v3.
*/


#include <cmath>

inline void swap(double& a, double& b) {
  double temp=a;
  a = b;
  b = temp;
}

#define TRIG_PRECALC 1

#if TRIG_PRECALC
////// template class SinCosSeries
// common series to compile-time calculation of sine and cosine functions

template<unsigned M, unsigned N, unsigned B, unsigned A>
struct SinCosSeries {
  static double value() {
    return 1-(A*M_PI/B)*(A*M_PI/B)/M/(M+1)
              *SinCosSeries<M+2,N,B,A>::value();
  }
};

template<unsigned N, unsigned B, unsigned A>
struct SinCosSeries<N,N,B,A> {
   static double value() { return 1.; }
};

////// template class Sin
// compile-time calculation of sin(A*M_PI/B) function

template<unsigned B, unsigned A, typename T=double>
struct Sin;

template<unsigned B, unsigned A>
struct Sin<B,A,float> {
   static float value() {
      return (A*M_PI/B)*SinCosSeries<2,24,B,A>::value();
   }
};
template<unsigned B, unsigned A>
struct Sin<B,A,double> {
   static double value() {
      return (A*M_PI/B)*SinCosSeries<2,34,B,A>::value();
   }
};

////// template class Cos
// compile-time calculation of cos(A*M_PI/B) function

template<unsigned B, unsigned A, typename T=double>
struct Cos;

template<unsigned B, unsigned A>
struct Cos<B,A,float> {
   static float value() {
      return SinCosSeries<1,23,B,A>::value();
   }
};
template<unsigned B, unsigned A>
struct Cos<B,A,double> {
   static double value() {
      return SinCosSeries<1,33,B,A>::value();
   }
};

#endif



template<unsigned N, typename T=double>
class DanielsonLanczos
{
  DanielsonLanczos<N/2,T> next;
public:
  void apply(T* data, int iSign) {
    next.apply(data, iSign);
    next.apply(data+N, iSign);

#if TRIG_PRECALC
    T wtemp = iSign*Sin<N,1,T>::value();
    T wpi = -iSign*Sin<N,2,T>::value();
#else
    T wtemp = iSign*std::sin(M_PI/N);
    T wpi = -iSign*std::sin(2*M_PI/N);
#endif
    T wpr = -2.0*wtemp*wtemp;
    T wr = 1.0;
    T wi = 0.0; 

    for (unsigned i=0; i<N; i+=2) {
      int iN = i+N;

      T tempr = data[iN]*wr - data[iN+1]*wi;
      T tempi = data[iN]*wi + data[iN+1]*wr;

      data[iN] = data[i]-tempr;
      data[iN+1] = data[i+1]-tempi;
      data[i] += tempr;
      data[i+1] += tempi;

      wtemp = wr;
      wr += wr*wpr - wi*wpi;
      wi += wi*wpr + wtemp*wpi;
    }
  }
}; 


template<typename T>
class DanielsonLanczos<1,T>
{
public:
  void apply(T* data, int iSign) { }
};


/*!\brief Create a templated FFT/Inverse FFT object

  How to use this FFT
  FFT<LOG_LENGTH, double> gfft;
  
  unsigned long n = 1<<LOG_LENGTH;
  double* data = new double [2*n];

  gfft.fft(data);
  gfft.ifft(data);
*/
template<unsigned P,typename T=double>
class FFT
{
  enum { N = 1<<P };
  DanielsonLanczos<N,T> recursion;

  // reverse-binary reindexing
  void scramble(T* data) {
    int j=1;
    for (int i=1; i<2*N; i+=2) {
      if (j>i) {
        swap(data[j-1], data[i-1]);
        swap(data[j], data[i]);
      }
      int m = N;
      while (m>=2 && j>m) {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
  }

  void rescale(T* data)
  {
    /*  scale all results by 1/n */
    T scale = static_cast<T>(1)/N;
    for (unsigned i=0; i<2*N; i++) {
      data[i] *= scale;
    }
  }

public:
  /*!\brief Replaces data[1..2*N] by its discrete Fourier transform */
  void fft(T* data) {
    scramble(data);
    recursion.apply(data,1);
  }

  /*!\brief Replaces data[1..2*N] by its inverse discrete Fourier transform */
  void ifft(T* data) {
    scramble(data);
    recursion.apply(data,-1);
    rescale(data);
  }
};

#endif /* _fft_h_ */
