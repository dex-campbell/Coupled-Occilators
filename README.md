# Coupled-Occilators

# QR Factorization and Eigenvalue Computation

This project implements the QR factorization algorithm to compute the eigenvalues of a square matrix. The program is written in C and provides an interactive console-based interface for users to input the matrix size, matrix values, and the number of iterations for the QR algorithm.

## Features

- **QR Factorization**: Decomposes a given matrix \( A \) into an orthogonal matrix \( Q \) and an upper triangular matrix \( R \).
- **Eigenvalue Computation**: Uses the QR algorithm to iteratively update the matrix \( A \) and compute its eigenvalues.
- **Matrix Printing**: Functions to print matrices in a readable format.

## Getting Started

### Prerequisites

- A C compiler (e.g., `gcc`).

### Compilation

To compile the program, run the following command in your terminal:

```sh
gcc -o qr_factorization Comp_Phys_CW_1.c -lm
