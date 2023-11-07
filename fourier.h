#ifndef FOURIER_H
#define FOURIER_H

#define complex_vector std::vector<std::complex<double>>
#define complex_matrix std::vector<std::vector<std::complex<double>>>

void fft(complex_vector& v, bool inverse = false);
void fftRadix2(complex_vector& v, bool inverse = false);
complex_matrix fourierMatrix(int N, bool inverse = false);
static void combineHalfVectors(complex_vector& v, const complex_vector& vEven, const complex_vector& vOdd, int n, int inverse);
static void printVector(const complex_vector& v);
static void printMatrix(const complex_matrix& M);
static int nextPowerOf2(int n);
static bool isPowerOf2(int n);

#endif 