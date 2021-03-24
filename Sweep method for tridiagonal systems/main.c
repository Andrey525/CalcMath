#include <stdio.h>
#include <stdlib.h>

// #define N 3
#define N 5

void print_matrix(double** matrix, double* vector)
{
    for (int index_row = 0; index_row < N; index_row++) {
        for (int index_col = 0; index_col < N + 1; index_col++) {
            if (index_col != N) {
                printf("%3.0lf  ", matrix[index_row][index_col]);
            } else {
                printf("| %3.0lf  ", vector[index_row]);
            }
        }
        printf("\n");
    }
}

void init_matrix1(double** matrix, double* vector)
{
    matrix[0][0] = 2;
    matrix[0][1] = -1;
    matrix[0][2] = 0;
    matrix[1][0] = 5;
    matrix[1][1] = 4;
    matrix[1][2] = 2;
    matrix[2][0] = 0;
    matrix[2][1] = 1;
    matrix[2][2] = -3;
    vector[0] = 3;
    vector[1] = 6;
    vector[2] = 2;
}

void init_matrix2(double** matrix, double* vector)
{
    matrix[0][0] = 2;
    matrix[0][1] = -1;
    matrix[0][2] = 0;
    matrix[0][3] = 0;
    matrix[0][4] = 0;
    matrix[1][0] = -3;
    matrix[1][1] = 8;
    matrix[1][2] = -1;
    matrix[1][3] = 0;
    matrix[1][4] = 0;
    matrix[2][0] = 0;
    matrix[2][1] = -5;
    matrix[2][2] = 12;
    matrix[2][3] = 2;
    matrix[2][4] = 0;
    matrix[3][0] = 0;
    matrix[3][1] = 0;
    matrix[3][2] = -6;
    matrix[3][3] = 18;
    matrix[3][4] = -4;
    matrix[4][0] = 0;
    matrix[4][1] = 0;
    matrix[4][2] = 0;
    matrix[4][3] = -5;
    matrix[4][4] = 10;
    vector[0] = -25;
    vector[1] = 72;
    vector[2] = -69;
    vector[3] = -156;
    vector[4] = 20;
}

double* sweep_method_for_tridiagonal_systems(double** matrix, double* vector)
{
    double* resolution = malloc(N * sizeof(double));
    double* y = malloc(N * sizeof(double));
    double* alpha = malloc(N * sizeof(double));
    double* beta = malloc(N * sizeof(double));

    /****** Прямая прогонка ******/

    /*** Для первой строки определяем прогоночные коэффициенты ***/
    y[0] = matrix[0][0];
    alpha[0] = -matrix[0][1] / y[0];
    beta[0] = vector[0] / y[0];
    ////////////////////////////

    /*** Для строк от i == 1 до i == N - 2 определяем прогоночные коэффициенты с помощью рекурентных формул ***/
    for (int index_row = 1; index_row < N - 1; index_row++) {
        y[index_row] = matrix[index_row][index_row] + (matrix[index_row][index_row - 1] * alpha[index_row - 1]);
        alpha[index_row] = -matrix[index_row][index_row + 1] / y[index_row];
        beta[index_row] = (vector[index_row] - matrix[index_row][index_row - 1] * beta[index_row - 1]) / y[index_row];
    }
    ////////////////////////////

    /*** Для последней строки определяем прогоночные коэффициенты ***/
    y[N - 1] = matrix[N - 1][N - 1] + matrix[N - 1][N - 2] * alpha[N - 2];
    beta[N - 1] = (vector[N - 1] - matrix[N - 1][N - 2] * beta[N - 2]) / y[N - 1];
    ////////////////////////////

    /*##### Прямая прогонка закончена #####*/

    /****** Обратная прогонка ******/
    resolution[N - 1] = beta[N - 1]; // для последней строки
    for (int index_row = N - 2; index_row >= 0; index_row--) {
        resolution[index_row] = alpha[index_row] * resolution[index_row + 1] + beta[index_row];
    }

    return resolution;

    /*##### Обратная прогонка закончена #####*/
}

void print_resolution(double* resolution)
{
    for (int i = 0; i < N; i++) {
        printf("X[%d] = %.2lf\n", i, resolution[i]);
    }
}

int main()
{
    double *matrix[N], *vector;
    for (int index_row = 0; index_row < N; index_row++) {
        matrix[index_row] = malloc(N * sizeof(double));
        vector = malloc(N * sizeof(double));
    }
    // init_matrix1(matrix, vector);
    init_matrix2(matrix, vector);
    print_matrix(matrix, vector);
    double* resolution;
    resolution = sweep_method_for_tridiagonal_systems(matrix, vector);
    print_resolution(resolution);
    return 0;
}