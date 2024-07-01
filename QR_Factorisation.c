#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h> 

void print_matrix(int n, double** matrix) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void qr_factorization(int n, double** A, double** Q, double** R) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
            // Compute the inner product of A's j-th column with Q's i-th column
            double dot_product = 0.0;
            for (int k = 0; k < n; k++) {
                dot_product += A[k][j] * Q[k][i];
            }

            // Subtract the projection of A's j-th column onto Q's i-th column from A's j-th column
            for (int k = 0; k < n; k++) {
                A[k][j] -= dot_product * Q[k][i];
            }

            // Update the corresponding value in R (upper-right triangle)
            R[i][j] = dot_product;
        }

        // Calculate the norm of the j-th column of A
        double norm = 0.0;
        for (int k = 0; k < n; k++) {
            norm += A[k][j] * A[k][j];
        }
        norm = sqrt(norm);

        // Set the corresponding value in R (diagonal)
        R[j][j] = norm;

        // Normalize the j-th column of Q
        for (int k = 0; k < n; k++) {
            Q[k][j] = A[k][j] / norm;
        }
    }
}

int main() {
    int matrix_size;
    bool valid_input = false;

    while (!valid_input) {
        printf("Enter the size of the Square Matrix: ");
        if (scanf("%d", &matrix_size) == 1 && matrix_size > 0) {
            valid_input = true;
        } else {
            printf("Invalid input, please enter a positive integer for the matrix size.\n");
            while (getchar() != '\n');
        }
    }

    // Allocating memory for the arrays and rows in the arrays
    double** A = (double**)malloc(matrix_size * sizeof(double*));
    double** Q = (double**)malloc(matrix_size * sizeof(double*));
    double** R = (double**)malloc(matrix_size * sizeof(double*));

    for (int i = 0; i < matrix_size; i++) {
        A[i] = (double*)malloc(matrix_size * sizeof(double));
        Q[i] = (double*)malloc(matrix_size * sizeof(double));
        R[i] = (double*)malloc(matrix_size * sizeof(double));

        for (int j = 0; j < matrix_size; j++) {
            printf("Please enter the %d row and %d column value: ", i + 1, j + 1);
            if (scanf("%lf", &A[i][j]) != 1) {
                printf("Invalid input\n");
                return 1;
            }
        }
    }

    int iteration;
    valid_input = false;

    while (!valid_input) {
        printf("How many times would you like to iterate it: ");
        if (scanf("%d", &iteration) == 1 && iteration > 0) {
            valid_input = true;
        } else {
            printf("Invalid input, please enter a positive integer for iteration.\n");
            while (getchar() != '\n');
        }
    }

    for (int step = 1; step < iteration; step++) {
        printf("\n Iteration number %d:\n", step);
        qr_factorization(matrix_size, A, Q, R);
        printf("Matrix Q:\n");
        print_matrix(matrix_size, Q);
        printf("Matrix R:\n");
        print_matrix(matrix_size, R);
        
        // Update A for the next step
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                double sum = 0.0;
                for (int k = 0; k < matrix_size; k++) {
                    sum += R[i][k] * Q[k][j];
                }
                A[i][j] = sum;
            }
        }
        printf("A Matrix: \n");
        print_matrix(matrix_size, A);
    }
    printf("Final Matrix: \n");
    print_matrix(matrix_size, A);

    printf("The eigenvalues for this matrix are:\n");
    for (int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++){
            if (i==j){
                printf("%lf\n", A[i][j]);
            }
        }
    }

    // Free memory for the matrices
    for (int i = 0; i < matrix_size; i++) {
        free(A[i]);
        free(Q[i]);
        free(R[i]);
    }
    free(A);
    free(Q);
    free(R);

    return 0;
}
