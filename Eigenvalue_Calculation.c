#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h> 

#define MAX_MATRIX_SIZE 2

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
    // User input for iteration number, range of masses, and spring constant
    int matrix_size = MAX_MATRIX_SIZE;

    double spring_constant;
    bool valid_input = false;

    while (!valid_input) {
        printf("Please enter the value for k : ");
        if (scanf("%lf", &spring_constant) == 1 && spring_constant > 0) {
            valid_input = true;
        } else {
            printf("Invalid input, please try again.\n");
            // Clear the input buffer to avoid an infinite loop on invalid input
            while (getchar() != '\n');
        }
    }

    double mass_low;
    valid_input = false;

    while (!valid_input) {
        printf("Please enter the lower limit for m: ");
        if (scanf("%lf", &mass_low) == 1 && mass_low > 0) {
            valid_input = true;
        } else {
            printf("Invalid input, please try again.\n");
            while (getchar() != '\n');
        }
    }

    double mass_high;
    valid_input = false;

    while (!valid_input) {
        printf("Please enter the upper limit for m: ");
        if (scanf("%lf", &mass_high) == 1 && mass_high > 0 && mass_high > mass_low) {
            valid_input = true;
        } else {
            printf("Invalid input, please try again.\n");
            while (getchar() != '\n');
        }
    }

    int no_diff_mass;
    valid_input = false;

    while (!valid_input) {
        printf("Please enter the number of m you want to test: ");
        if (scanf("%d", &no_diff_mass) == 1 && no_diff_mass > 0) {
            valid_input = true;
        } else {
            printf("Invalid input, please try again.\n");
            while (getchar() != '\n');
        }
    }

    // Calculate mass_diff
    double mass_diff = (mass_high - mass_low) / (no_diff_mass - 1);

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
    

    // Allocate memory for eigenvalues and mass values
    double** eigenvalues = (double**)malloc(no_diff_mass * sizeof(double*));
    double* mass_values = (double*)malloc(no_diff_mass * sizeof(double));

    for (int i = 0; i < no_diff_mass; i++) {
        eigenvalues[i] = (double*)malloc(matrix_size * sizeof(double));
        mass_values[i] = mass_low + i * mass_diff;
    }

    // Define a filename for output
    char filename[] = "mass_vs_Oscillation_frequency.txt";

    FILE* outputFile = fopen(filename, "w");

    if (outputFile == NULL) {
        printf("Error opening the file for writing.\n");
        return 1;
    }

    // Loop through different mass values
    for (int m_index = 0; m_index < no_diff_mass; m_index++) {
        double mass = mass_values[m_index];

        // Allocate and initialize the matrix A based on mass and spring_constant
        double** A = (double**)malloc(matrix_size * sizeof(double*));
        for (int i = 0; i < matrix_size; i++) {
            A[i] = (double*)malloc(matrix_size * sizeof(double));
        }

        // Initialize A based on mass and spring_constant
        A[0][0] = (-2 * spring_constant) / mass;
        A[0][1] = (spring_constant) / mass;
        A[1][0] = (spring_constant) / (10*mass);
        A[1][1] = (-2 * spring_constant) / (10*mass);

        // QR factorization and eigenvalues calculation
        for (int step = 0; step < iteration; step++) {
            double** Q = (double**)malloc(matrix_size * sizeof(double*));
            double** R = (double**)malloc(matrix_size * sizeof(double*));
            for (int i = 0; i < matrix_size; i++) {
                Q[i] = (double*)malloc(matrix_size * sizeof(double));
                R[i] = (double*)malloc(matrix_size * sizeof(double));
            }

            qr_factorization(matrix_size, A, Q, R);

            // Store eigenvalues
            for (int i = 0; i < matrix_size; i++) {
                eigenvalues[m_index][i] = A[i][i];
            }

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

            // Free memory for Q and R
            for (int i = 0; i < matrix_size; i++) {
                free(Q[i]);
                free(R[i]);
            }
            free(Q);
            free(R);
        }

        // Write mass value and eigenvalues to the file
        fprintf(outputFile, "Mass: %lf\n", mass);
        fprintf(outputFile, "Oscillation frequency: ");
        for (int i = 0; i < matrix_size; i++) {
            fprintf(outputFile, "%lf\t", sqrt(-1*eigenvalues[m_index][i]));
        }
        fprintf(outputFile, "\n");

        // Free the memory for matrix A
        for (int i = 0; i < matrix_size; i++) {
            free(A[i]);
        }
        free(A);
    }

    fclose(outputFile);

    // Print or use eigenvalues, oscillation frequency and mass values as needed
    for (int m_index = 0; m_index < no_diff_mass; m_index++) {
        printf("Mass: %lf\n", mass_values[m_index]);
        printf("Eigenvalues: ");
        for (int i = 0; i < matrix_size; i++) {
            printf("%lf\t", eigenvalues[m_index][i]);
        }
        printf("\n");
        printf("Oscillation frequency: ");
        for (int i = 0; i < matrix_size; i++) {
            printf("%lf\t", sqrt(-1*eigenvalues[m_index][i]));
        }
        printf("\n");
    }

    // Free memory for eigenvalues and mass_values
    for (int i = 0; i < no_diff_mass; i++) {
        free(eigenvalues[i]);
    }
    free(eigenvalues);
    free(mass_values);

    return 0;
}
