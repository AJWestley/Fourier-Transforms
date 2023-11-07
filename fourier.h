#ifndef FOURIER_H
#define FOURIER_H

#include <complex>
#include <vector>
#include <iostream>

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> v, bool inverse = false);
std::vector<std::complex<double>> fftRadix2(std::vector<std::complex<double>> v, bool inverse = false);
std::vector<std::vector<std::complex<double>>> fourierMatrix(int N, bool inverse = false);
void printVector(const std::vector<std::complex<double>>& v);
void printMatrix(const std::vector<std::vector<std::complex<double>>>& M);
int nextPowerOf2(int n);
bool isPowerOf2(int n);

#endif 