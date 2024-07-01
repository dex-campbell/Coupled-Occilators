# Coupled-Occilators

## QR Factorization and Eigenvalue Computation

This repository contains two C programs for calculations of values for a coupled Occilator focusing on matrix operations, specifically QR factorization and eigenvalue calculations.

## Files

- **QR Factorization**
- **Eigenvalue Computation**

## File Description
### QR_Factorisation
- This program performs QR factorization of a square matrix. QR factorization decomposes a matrix into an orthogonal matrix (Q) and an upper triangular matrix (R).

#### KeyFunctions
- print_matrix(int n, double** matrix): Prints the matrix of size n x n.
- qr_factorization(int n, double** A, double** Q, double** R): Performs QR factorization on matrix A and stores the results in matrices Q and R.

### Eigenvalue_Calculation
- This program calculates the eigenvalues of a square matrix and uses them to determine the oscillation frequencies for different mass values.

#### Key Functions
- print_matrix(int n, double** matrix): Prints the matrix of size n x n.
- qr_factorization(int n, double** A, double** Q, double** R): Performs QR factorization on matrix A and stores the results in matrices Q and R.
- main(): Main function to drive the program. Reads input for matrix size, number of different masses, and calculates eigenvalues and oscillation frequencies.

## Dependencies
- Both programs require the math library for mathematical functions. Ensure you link the math library during compilation using the -lm flag.

## Usage
- These programs are designed for educational purposes to demonstrate matrix operations in computational physics. The QR factorization program decomposes a matrix into orthogonal and upper triangular matrices, while the eigenvalue calculation program finds the eigenvalues and corresponding oscillation frequencies for different mass values.

