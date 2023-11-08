<h1 align="center">Fourier Transforms</h1>
<p align="center">A C++ library of Discrete Fourier Transform (DFT) methods.</p>

## Overview

This library contains various Discrete Fourier Transform (DFT) methods that can be used in C++. </br>
All the functions live inside the 'fourier' namespace. 

## Custom Datatypes

In the header, there are two 'custom datatypes' that are defined to shorten a few of the function declarations:
 * **complex_vector**: `std::vector<std::complex<double>>`
 * **complex_matrix**: `std::vector<std::vector<std::complex<double>>>`

## Main Methods

### **fft**

Perform a Fast Fourier Transform (FFT) or Inverse FFT (IFFT) on a complex-valued vector. </br>
**fft** can take vectors of any length. If length of the input vector is not a power of 2, it is zero-padded and the transform of the padded vector is returned.

#### Parameters:
 * **v**: The input complex vector for which the FFT or IFFT is computed.
 * **inverse**: If set to true, computes the Inverse FFT (IFFT); otherwise, computes the FFT.

#### Returns:
 * The transformed complex vector (FFT or IFFT result).

#### Throws:
 * **std::invalid_argument**: If the input vector is empty.

#### Example:
```c++
using namespace fourier;
complex_vector v = {4, 0, 0, 0, 0};
fft(v); // v now = {4, 4, 4, 4, 4, 4, 4, 4}
fft(v, true); // v now = {4, 0, 0, 0, 0, 0, 0, 0}
```

### **fftRadix2**

Perform a Fast Fourier Transform (FFT) or Inverse FFT (IFFT) on a complex-valued vector using the radix-2 algorithm. </br>
**fftRadix2** can only take vectors with a length that is a power of 2.

#### Parameters:
 * **v**: The input complex vector for which the FFT or IFFT is computed.
 * **inverse**: If set to true, computes the Inverse FFT (IFFT); otherwise, computes the FFT.

#### Returns:
 * The transformed complex vector (FFT or IFFT result).

#### Throws:
 * **std::invalid_argument**: If the length of the input vector 'v' is not a power of 2, an exception is thrown.

#### Example:
```c++
using namespace fourier;
complex_vector v = {4, 0, 0, 0};
fftRadix2(v); // v now = {4, 4, 4, 4}
fftRadix2(v, true); // v now = {4, 0, 0, 0}
```

### **convolve**

Computes the discrete convolution of two complex-valued input vectors.

#### Parameters:
 * **f**: The first vector to convolve.
 * **g**: The second vector to convolve.

#### Returns:
 * The convolution of input vectors 'f' and 'g'.

#### Example:
```c++
using namespace fourier;
complex_vector F = {1, 2, 3}, G = {4, 5, 6};
complex_vector H = convolve(F, G); 
// H -> {(4,0i), (13,0i), (28,0i), (27,0i), (18,0i), (0,0i), (0,0i), (0,0i)}
```

### **fourierMatrix**

Generates a Fourier matrix of size N.

#### Parameters:
 * **N**: The size of the square Fourier matrix to generate. N must be a positive integer.
 * **inverse**: If set to true, computes the inverse Fourier matrix; otherwise, computes the standard Fourier matrix.

#### Returns:
 * A square Fourier matrix of size N filled with complex values.

#### Throws:
 * **std::invalid_argument**: If N is not a positive integer, an exception is thrown.

#### Example:
```c++
using namespace fourier;
complex_matrix F4;

F4 = fourierMatrix(4);
// F4 ->  {{1,  1,  1,  1},
//         {1,  i, -1, -i},
//         {1, -1,  1, -1},
//         {1, -i, -1,  i}}

F4 = fourierMatrix(4, true);
// F4 ->  {{1,  1,  1,  1},
//         {1, -i, -1,  i},
//         {1, -1,  1, -1},
//         {1,  i, -1, -i}}
```

### **toComplexVector**

Casts a vector of any type to a complex-valued vector.

#### Parameters:
 * **v**: The vector to be cast.

#### Returns:
 * A complex-valued cast of the input vector 'v'.

#### Example: (with convolve)
```c++
using namespace fourier;
std::vector<int> f = {1, 2, 3}, g = {4, 5, 6};
complex_vector F = toComplexVector(f), G = toComplexVector(g);
complex_vector H = convolve(F, G); 
// H -> {(4,0i), (13,0i), (28,0i), (27,0i), (18,0i), (0,0i), (0,0i), (0,0i)}
```

### **toIntVector**

Casts a complex-valued vector to an integer-valued vector.

#### Parameters:
 * **v**: The vector to be cast.

#### Returns:
 * An integer-valued cast of the input vector 'v'.

#### Example: (with convolve)
```c++
using namespace fourier;
std::vector<int> f = {1, 2, 3}, g = {4, 5, 6};
complex_vector F = toComplexVector(f), G = toComplexVector(g);
std::vector<int> h = toIntVector(convolve(F, G)); 
// h -> {4, 13, 28, 27, 18, 0, 0, 0}
```

### **toDoubleVector**

Casts a complex-valued vector to an real-valued vector.

#### Parameters:
 * **v**: The vector to be cast.

#### Returns:
 * A real-valued cast of the input vector 'v'.

#### Example: (with convolve)
```c++
using namespace fourier;
std::vector<int> f = {1, 2, 3}, g = {4, 5, 6};
complex_vector F = toComplexVector(f), G = toComplexVector(g);
std::vector<int> h = toIntVector(convolve(F, G)); 
// h -> {6, 17.65, 33.98, 30.9, 18, 0, 0, 0}
```

## To be added:
 * FFT Shift
 * 2D FFT
 * Mixed-radix FFT
 * Unit Tests
