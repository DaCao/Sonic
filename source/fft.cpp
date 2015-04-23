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
//    std::cout << "\n";
//    std::cout << "//////////////Da Cao/////////////////////\n";
//    
    
    float InputFloat[N];
    for (int i = 0; i < N; i++){
        InputFloat[i] = Input[i].real();
    }
    
    //   Check input parameters
    if (!Input || !Output || N < 1 || N & (N - 1))
        return false;

    const int size = N;
    int ln2 = log2(N);
    FFTSetup fftSetup = vDSP_create_fftsetup(ln2, FFT_RADIX2);
    
//    std::cout << "forward func, size is: " << N << "   ln2 is: "<< ln2 << "\n";
//    for (int i = 0; i < 10; i++){
//        std::cout << i << "th forward InputFloat: "<< InputFloat[i]<< "\n";
//    }
    
    vDSP_ctoz((DSPComplex *)InputFloat, 2, Output, 1, size/2);
    
//    for (int i = 0; i < 10; i++){
//        std::cout << i << "th forward output before zrip: "<< Output->realp[i] <<"  " << Output->imagp[i]<< "\n";
//    }
    
    vDSP_fft_zrip(fftSetup, Output, 1, ln2, FFT_FORWARD);
    
//    for (int i = 0; i < 10; i++){
//        std::cout << i << "th forward output: "<< Output->realp[i] <<"  " << Output->imagp[i]<< "\n";
//    }

    
    
//    float *y = new float[size];
//    vDSP_fft_zrip(fftSetup, Output, 1, ln2, FFT_INVERSE);
//    vDSP_ztoc(Output, 1, (DSPComplex*)y, 1, size/2);
//    float scale = 0.5/size;
//    vDSP_vsmul(y, 1, &scale, y, 1, size);
    
//    for (int i = 0; i < 10; i++){
//        std::cout << i << "th forward output the inverse back: "<< y[i] <<"\n";
//    }
//    
//    std::cout << "\n";
//    std::cout << "\n";
    
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
    
    const int size = N;
    int ln2 = round(log2(size));
    FFTSetup fftSetup = vDSP_create_fftsetup(ln2, FFT_RADIX2);
    
//    std::cout << "inverse func, size is: " << size << "\n";
//    for (int i = 0; i < 10; i++){
//        std::cout << i << "th inverse input: "<< Input->realp[i] << "  " << Input->imagp[i] << "\n";
//    }
    
    float *y = new float[size];
    vDSP_fft_zrip(fftSetup, Input, Stride, ln2, FFT_INVERSE);
    vDSP_ztoc(Input, Stride, (DSPComplex*)y, 2, size/2);
    float scale = 0.5/size;
    vDSP_vsmul(y, 1, &scale, y, 1, size);
    
//    for (int i = 0; i < 10; i++){
//        std::cout << i << "th y: "<< y[i]<< "\n";
//    }

    
    for (int i = 0; i < size; i++) {
        Complex tmp(y[i], 0);
        Output[i] = tmp;
    }
    
    //   Succeeded
    return true;
}

