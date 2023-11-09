#ifndef FOURIER_H
#define FOURIER_H

#include <complex>
#include <vector>
#include <iostream>
#include <omp.h>

#define complex_vector std::vector<std::complex<double>>
#define complex_matrix std::vector<std::vector<std::complex<double>>>

namespace fourier {
    // FFT Methods
    void fft(complex_vector& v, bool inverse = false);
    void fftRadix2(complex_vector& v, bool inverse = false);

    // Fourier Matrix
    complex_matrix fourierMatrix(int N, bool inverse = false);

    // Applied FFT Methods
    complex_vector convolve(complex_vector f, complex_vector g);

    // Casting
    template <typename T> 
    complex_vector toComplexVector(const std::vector<T> &v);
    std::vector<int> toIntVector(const complex_vector &v);
    std::vector<double> toDoubleVector(const complex_vector &v, int precision = 6);
}

// Printing Functions
template <typename T>
static void printVector(const std::vector<T>& v);
template <typename T>
static void printMatrix(const std::vector<std::vector<T>> &M);

// Auxilliary functions
static void combineHalfVectors(complex_vector &v, const complex_vector &vEven, const complex_vector &vOdd, int n, int inverse);
static int nextPowerOf2(int n);
static bool isPowerOf2(int n);

#endif 