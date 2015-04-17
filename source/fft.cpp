//   fft.cpp - impelementation of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

//     modified by Da Cao, Oct.1st/2014.
//     Changed into a FHT method to compute DFT

#include <fstream>
#include <iostream>
#include <math.h>
#include "../include/fft.h"
#include <cmath>
#include <Accelerate/Accelerate.h>
#include "Mixer3D.h"

using namespace std;
const vDSP_Stride Stride = 1;


//   FORWARD FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
bool CFFT::Forward(Complex *const Input, DSPSplitComplex *const Output, const unsigned int N)
{
    //   Check input parameters
    if (!Input || !Output || N < 1 || N & (N - 1))
        return false;

    const int size = sizeof(Input);
    int ln2 = round(log2(size));
    FFTSetup fftSetup = vDSP_create_fftsetup(ln2, FFT_RADIX2);
    
    for (int i = 0; i < size; i++){
        std::cout << i << ": "<< Input[i]<< "\n";
    }
    
    vDSP_ctoz((DSPComplex *)Input, 2*Stride, Output, 1, size/2);
    vDSP_fft_zrip(fftSetup, Output, 1, ln2, FFT_FORWARD);

    return true;
}




//   INVERSE FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
//     Scale  - if to scale result
bool CFFT::Inverse(DSPSplitComplex *const Input, Complex *const Output, const unsigned int N)
{
    //   Check input parameters
    if (!Input || !Output || N < 1 || N & (N - 1))
        return false;
    
    const int size = 2*sizeof(Input);
    int ln2 = round(log2(size));
    FFTSetup fftSetup = vDSP_create_fftsetup(ln2, FFT_RADIX2);
    
    
    float *y = new float[size*2];
    vDSP_fft_zrip(fftSetup, Input, Stride, ln2, FFT_INVERSE);
    vDSP_ztoc(Input, Stride, (DSPComplex*)y, 2, size/2);
    float scale = 0.5/size;
    vDSP_vsmul(y, 1, &scale, y, 1, size);
    
    
    for (int i = 0; i < size; i++) {
        Complex tmp(y[i], 0);
        Output[i] = tmp;
    }
    
    //   Succeeded
    return true;
}






//   Scaling of inverse FFT result
void CFFT::Scale(Complex *const Data, const unsigned int N)
{
    const double Factor = 1. / double(N);
    //   Scale all data entries
    for (unsigned int Position = 0; Position < N; ++Position)
        Data[Position] *= Factor;
}// SIMD?
