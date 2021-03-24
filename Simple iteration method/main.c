#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_matrix(double** matrix, double* vector, int n)
{
    for (int index_row = 0; index_row < n; index_row++) {
        for (int index_col = 0; index_col < n + 1; index_col++) {
            if (index_col != n) {
                printf("%.0lf  ", matrix[index_row][index_col]);
            } else {
                printf("|  %.0lf  ", vector[index_row]);
            }
        }
        printf("\n");
    }
}

double* simple_iteration_method(double** matrix, double* vector, int n, double* X)
{
    double e = 0.00001; // точность
    /***** Выводим матрицу *****/
    print_matrix(matrix, vector, n);
    /***** Проверяем на условие сходимости *****/
    for (int index_row = 0; index_row < n; index_row++) {
        double sum = 0;
        for (int index_col = 0; index_col < n; index_col++) {
            if (index_row != index_col) {
                sum += fabs(matrix[index_row][index_col]);
            }
        }
        if (fabs(matrix[index_row][index_row]) > sum) {
            printf("Условие на сходимость строки %d выполнено\n", index_row);
        } else {
            printf("Условие на сходимость строки %d не выполнено\n", index_row);
            exit(-1);
        }
    }
    printf("\nВсе условия на сходимость выполнены\n");

    /***** Метод простой итерации *****/
    X = malloc(sizeof(double) * n);

    double* old_x = malloc(sizeof(double) * n);
    int count_of_iterations = 0;
    while (1) {
        count_of_iterations++;
        int fulfilled_conditions = 0;
        for (int index_row = 0; index_row < n; index_row++) {
            double sum = 0;
            for (int index_col = 0; index_col < n; index_col++) {
                if (index_row != index_col) {
                    sum += matrix[index_row][index_col] * old_x[index_col];
                }
            }
            X[index_row] = (vector[index_row] - sum) / matrix[index_row][index_row];
        }
        printf("#################\n");
        printf("Итерация номер %d\n", count_of_iterations);
        printf("\n");
        for (int index_row = 0; index_row < n; index_row++) {
            printf("X%d = %lf\n", index_row, X[index_row]);
        }
        printf("\n");
        for (int index_row = 0; index_row < n; index_row++) {
            if (fabs(X[index_row] - old_x[index_row]) < e) {
                fulfilled_conditions++;
            }
        }
        printf("%d fulfilled_conditions\n", fulfilled_conditions);
        printf("******************\n");
        if (fulfilled_conditions == n) {
            break;
        } else {
            for (int index_row = 0; index_row < n; index_row++) {
                old_x[index_row] = X[index_row];
            }
        }
    }

    return X;
}

double* method_of_zeidel(double** matrix, double* vector, int n, double* X)
{
    double e = 0.00001; // точность
    /***** Выводим матрицу *****/
    print_matrix(matrix, vector, n);
    /***** Проверяем на условие сходимости *****/
    for (int index_row = 0; index_row < n; index_row++) {
        double sum = 0;
        for (int index_col = 0; index_col < n; index_col++) {
            if (index_row != index_col) {
                sum += fabs(matrix[index_row][index_col]);
            }
        }
        if (fabs(matrix[index_row][index_row]) > sum) {
            printf("Условие на сходимость строки %d выполнено\n", index_row);
        } else {
            printf("Условие на сходимость строки %d не выполнено\n", index_row);
            exit(-1);
        }
    }
    printf("\nВсе условия на сходимость выполнены\n");

    /***** Метод Зейделя *****/
    X = malloc(sizeof(double) * n);

    double* old_x = malloc(sizeof(double) * n);
    int count_of_iterations = 0;
    int fulfilled_conditions = 0;
    while (1) {
        count_of_iterations++;
        fulfilled_conditions = 0;
        for (int index_row = 0; index_row < n; index_row++) {
            double sum = 0;
            for (int index_col = 0; index_col < n; index_col++) {
                if (index_row != index_col) {
                    sum += matrix[index_row][index_col] * old_x[index_col];
                }
            }
            X[index_row] = (vector[index_row] - sum) / matrix[index_row][index_row];
            if (fabs(X[index_row] - old_x[index_row]) < e) {
                fulfilled_conditions++;
            }

            old_x[index_row] = X[index_row];
        }
        if (fulfilled_conditions == n) {
            break;
        }
        printf("#################\n");
        printf("Итерация номер %d\n", count_of_iterations);
        printf("\n");
        for (int index_row = 0; index_row < n; index_row++) {
            printf("X%d = %lf\n", index_row, X[index_row]);
        }
        printf("\n");
        printf("%d fulfilled_conditions\n", fulfilled_conditions);
        printf("******************\n");
    }
    return X;
}

int main()
{

    int n = 3;
    double* matrix[n];
    for (int index_row = 0; index_row < n; index_row++) {
        matrix[index_row] = malloc(n * sizeof(double));
    }
    matrix[0][0] = 10;
    matrix[0][1] = 1;
    matrix[0][2] = 1;
    matrix[1][0] = 2;
    matrix[1][1] = 10;
    matrix[1][2] = 1;
    matrix[2][0] = 2;
    matrix[2][1] = 2;
    matrix[2][2] = 10;
    double vector[3] = { 12, 13, 14 };
    double* X = NULL;

    X = simple_iteration_method(matrix, vector, n, X);
    // X = method_of_zeidel(matrix, vector, n, X);

    for (int index_row = 0; index_row < n; index_row++) {
        printf("X%d = %lf\n", index_row, X[index_row]);
    }

    return 0;
}
