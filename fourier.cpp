#include <complex>
#include <vector>
#include <iostream>
#include "fourier.h"

const double D_PI = 2 * M_PI;
const std::complex<double> I(0, 1);

int main() {
    std::vector<std::complex<double>> v(4, 1);
    v[1] = 2;
    v[2] = 3;
    v[3] = 4;
    printVector(v);
    fftRadix2(v, false);
    printVector(v);
    fftRadix2(v, true);
    printVector(v);
    return EXIT_SUCCESS;
}

/**
 * @brief Perform a Fast Fourier Transform (FFT) or Inverse FFT (IFFT) on a complex-valued vector.
 *
 * This function computes the FFT or IFFT of a complex-valued vector, using the radix-2 algorithm
 * when the vector size is a power of 2. If the vector size is not a power of 2, it automatically
 * pads the vector with zeros to the nearest power of 2. It can handle both the forward and
 * inverse transformations based on the 'inverse' parameter.
 *
 * @param v The input complex vector for which the FFT or IFFT is computed.
 * @param inverse If set to true, computes the Inverse FFT (IFFT); otherwise, computes the FFT.
 *
 * @return The transformed complex vector (FFT or IFFT result).
 * 
 * @throws std::invalid_argument If the input vector is empty.
 */
void fft(std::vector<std::complex<double>>& v, bool inverse /*= false*/) {
    int n = v.size();

    if (n == 0) {
        throw std::invalid_argument("Input vector is empty.");
    }

    // Ideal case - vector size is a power of 2
    if (isPowerOf2(n)) {
        fftRadix2(v, inverse);
        return;
    }
    
    // Non-ideal case - vector size not a power of 2
    v.resize(nextPowerOf2(n), 0);
    fftRadix2(v, inverse);
}

/**
 * @brief Perform a Fast Fourier Transform (FFT) or Inverse FFT (IFFT) using the radix-2 algorithm.
 *
 * This function computes the FFT or IFFT of a complex-valued vector using the radix-2 algorithm.
 * It can handle both the forward and inverse transformations based on the 'inverse' parameter.
 *
 * @param v The input complex vector for which the FFT or IFFT is computed.
 * @param inverse If set to true, computes the Inverse FFT (IFFT); otherwise, computes the FFT.
 *
 * @return The transformed complex vector (FFT or IFFT result).
 *
 * @throws std::invalid_argument If the length of the input vector 'v' is not a power of 2, an exception is thrown.
 */
void fftRadix2(std::vector<std::complex<double>>& v, bool inverse /*= false*/) {
    // Base case
    int n = v.size();
    if (n == 1) {
        return;
    }
    
    // Error check before run
    if (nextPowerOf2(n) != n) {
        throw std::invalid_argument("Length input vector not a power of 2.");
    }

    // Permute to get odd and even half-vectors
    std::vector<std::complex<double>> vEven(n / 2);
    std::vector<std::complex<double>> vOdd(n / 2);
    for (int i = 0; i < n/2; i++) {
        vEven[i] = v[2*i];
        vOdd[i] = v[2*i+1];
    }

    // Recursively calculate sub-vectors
    fftRadix2(vEven, inverse);
    fftRadix2(vOdd, inverse);

    // Combine to get full vector
    double angle = D_PI / n * (inverse ? 1 : -1);
    std::complex<double> w(1), wn(cos(angle), sin(angle));
    for (int i = 0; i < n/2; i++) {
        v[i] = vEven[i] + w * vOdd[i];
        v[n/2 + i] = vEven[i] - w * vOdd[i];
        if (inverse) {
            v[i] /= 2;
            v[n/2 + i] /= 2;
        }
        w *= wn;
    }
}

/**
 * @brief Generates a Fourier matrix of size N.
 *
 * This function generates a square Fourier matrix of size N, containing complex values.
 *
 * @param N The size of the square Fourier matrix to generate. N must be a positive integer.
 * @param inverse If set to true, computes the inverse Fourier matrix; otherwise, computes the standard Fourier matrix.
 *
 * @return A square Fourier matrix of size N filled with complex values.
 *
 * @throws std::invalid_argument If N is not a positive integer, an exception is thrown.
 */
std::vector<std::vector<std::complex<double>>> fourierMatrix(int N, bool inverse /*= false*/) {
    // N must be positive
    if (N <= 0) {
        throw std::invalid_argument("N must be a positive integer");
    }

    // Initialize Fourier matrix to all 1
    std::vector<std::vector<std::complex<double>>> F(N, std::vector<std::complex<double>>(N, 1.0));

    // Compute fourier coefficient w
    std::complex<double> w = exp((D_PI * I) / (double) N);
    if (inverse) w = pow(w, -1);

    // Populate matrix
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            F[i][j] = std::complex<double>(pow(w, i*j));
        }
    }
    return F;
}

/**
 * @brief Prints a one-dimensional vector of complex numbers to the standard output.
 *
 * This function prints the contents of a one-dimensional vector of complex numbers
 * to the standard output (typically the console). Each element of the vector is
 * separated by a space.
 *
 * @param v The vector of complex numbers to be printed.
 */
void printVector(const std::vector<std::complex<double>>& v) {
    std::cout << "[";
    for (std::complex<double> val : v) {
        std::cout << val << ", ";
    }
    std::cout << "\b\b]" << std::endl;
}

/**
 * @brief Prints a two-dimensional matrix of complex numbers to the standard output.
 *
 * This function prints the contents of a two-dimensional matrix of complex numbers
 * to the standard output (typically the console). Each element of the matrix is
 * separated by a space, and rows are separated by newline characters.
 *
 * @param M The matrix of complex numbers to be printed.
 */
void printMatrix(const std::vector<std::vector<std::complex<double>>>& M) {
    for (std::vector<std::complex<double>> row : M) {
        for (std::complex<double> val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

/**
 * Find the next power of 2 greater than or equal to the given integer.
 *
 * This function calculates the smallest power of 2 that is greater than or equal
 * to the given integer 'n' and returns it.
 *
 * @param n An integer for which the next power of 2 needs to be found.
 * @return The smallest power of 2 greater than or equal to 'n'.
 */
int nextPowerOf2(int n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

/**
 * Find if an integer is a power of 2.
 *
 * This function evaluates whether the given integer is a power of 2.
 *
 * @param n The integer to evaluate.
 * @return True if the integer is a power of 2, false otherwise.
 */
bool isPowerOf2(int n) {
    return ceil(log2(n)) == floor(log2(n));
}