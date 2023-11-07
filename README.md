<h1 align="center">Fourier Transforms</h1>
<p align="center">A C++ library of discrete Fourier transform methods.</p>

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
std::vector<std::complex<double>> v = {4, 0, 0, 0, 0};
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
std::vector<std::complex<double>> v = {4, 0, 0, 0};
fftRadix2(v); // v now = {4, 4, 4, 4}
fftRadix2(v, true); // v now = {4, 0, 0, 0}
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
std::vector<std::vector<std::complex<double>>> F4;

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

## Auxilliary Methods

### **printVector**

Prints a one-dimensional vector of complex numbers to the standard output.

#### Parameters:
 * **v**: The vector of complex numbers to be printed.

### **printMatrix**

Prints a two-dimensional matrix of complex numbers to the standard output.

#### Parameters:
 * **M**: The matrix of complex numbers to be printed.

### **nextPowerOf2**

Find the next power of 2 greater than or equal to the given integer.

#### Parameters:
 * **n**: An integer for which the next power of 2 needs to be found.

#### Returns:
 * The smallest power of 2 greater than or equal to 'n'.

### isPowerOf2

Find if an integer is a power of 2.

#### Parameters:
 * **n**: The integer to evaluate.

#### Returns:
 * True if the integer is a power of 2, false otherwise.
