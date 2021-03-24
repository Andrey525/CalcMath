#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void rsly_modified_gauss(double** matrix, double* vector, int n)
{
    int index_max, index_row, index_col, k;
    double element_max, temp;
    /********** Прямой ход **********/
    for (k = 0; k < n; k++) {
        /***** Поиск максимального элемента по абсолютной величине *****/
        index_max = k;
        element_max = fabs(matrix[k][k]);
        for (index_row = k + 1; index_row < n; index_row++) {
            if (fabs(matrix[index_row][k]) > element_max) {
                element_max = fabs(matrix[index_row][k]);
                index_max = index_row;
            }
        }
        /*****/

        /***** Перестановка строк k и index_max *****/
        if (k != index_max) {
            for (index_col = k; index_col < n; index_col++) {
                temp = matrix[k][index_col];
                matrix[k][index_col] = matrix[index_max][index_col];
                matrix[index_max][index_col] = temp;
            }
            temp = vector[k];
            vector[k] = vector[index_max];
            vector[index_max] = temp;
        }
        /*****/

        temp = 1 / matrix[k][k]; // для получения единицы на диагонали
        for (index_row = k; index_row < n; index_row++) {
            matrix[k][index_row] = matrix[k][index_row] * temp; // преобразуем строку (умножаем на обратное числу на диагонали)
        }
        vector[k] = vector[k] * temp;

        for (index_row = k + 1; index_row < n; index_row++) {
            for (index_col = k + 1; index_col < n; index_col++) {
                matrix[index_row][index_col] = matrix[index_row][index_col] - matrix[index_row][k] * matrix[k][index_col];
            }
            vector[index_row] = vector[index_row] - matrix[index_row][k] * vector[k];
        }
    }
    /***********************/

    /********** Обратный ход **********/

    for (index_row = n - 2; index_row >= 0; index_row--) {
        for (index_col = index_row + 1; index_col < n; index_col++) {
            vector[index_row] = vector[index_row] - matrix[index_row][index_col] * vector[index_col];
        }
    }

    /***********************/
}

void rsly_gauss(double** matrix, double* vector, int n)
{
    int index_row, index_col, k;
    double temp;
    /********** Прямой ход **********/
    for (k = 0; k < n; k++) {

        temp = 1 / matrix[k][k]; // для получения единицы на диагонали
        for (index_row = k; index_row < n; index_row++) {
            matrix[k][index_row] = matrix[k][index_row] * temp; // преобразуем строку (умножаем на обратное числу на диагонали)
        }
        vector[k] = vector[k] * temp;

        for (index_row = k + 1; index_row < n; index_row++) {
            for (index_col = k + 1; index_col < n; index_col++) {
                matrix[index_row][index_col] = matrix[index_row][index_col] - matrix[index_row][k] * matrix[k][index_col];
            }
            vector[index_row] = vector[index_row] - matrix[index_row][k] * vector[k];
        }
    }
    /***********************/

    /********** Обратный ход **********/

    for (index_row = n - 2; index_row >= 0; index_row--) {
        for (index_col = index_row + 1; index_col < n; index_col++) {
            vector[index_row] = vector[index_row] - matrix[index_row][index_col] * vector[index_col];
        }
    }

    /***********************/
}

int main()
{
    printf("Введите кол-во переменных: ");
    int n;
    scanf("%d", &n);
    srand(time(NULL));
    double *matrix[n], *matrix2[n];
    for (int index_row = 0; index_row < n; index_row++) {
        matrix[index_row] = malloc(n * sizeof(double));
        matrix2[index_row] = malloc(n * sizeof(double));
    }
    for (int index_row = 0; index_row < n; index_row++) {
        for (int index_col = 0; index_col < n; index_col++) {
            matrix[index_row][index_col] = rand() % 20 - 10;
            matrix2[index_row][index_col] = matrix[index_row][index_col];
        }
    }
    double* vector = malloc(sizeof(double) * n);
    double* vector2 = malloc(sizeof(double) * n);
    for (int index_row = 0; index_row < n; index_row++) {
        vector[index_row] = rand() % 20 - 10;
        vector2[index_row] = vector[index_row];
    }
    // matrix[0][0] = 4;
    // matrix[0][1] = -7;
    // matrix[0][2] = 8;
    // matrix[1][0] = 2;
    // matrix[1][1] = -4;
    // matrix[1][2] = 5;
    // matrix[2][0] = -3;
    // matrix[2][1] = 11;
    // matrix[2][2] = 1;
    // double vector[3] = { -23, -13, 16 };
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
    rsly_gauss(matrix, vector, n);
    rsly_modified_gauss(matrix2, vector2, n);
    printf("Решение, полученное методом Гаусса:\n");
    for (int index_row = 0; index_row < n; index_row++) {
        printf("X%d = %.2lf\n", index_row + 1, vector[index_row]);
    }
    printf("Решение, полученное модиф. методом Гаусса:\n");
    for (int index_row = 0; index_row < n; index_row++) {
    	printf("X%d = %.2lf\n", index_row + 1, vector2[index_row]);
    }

    return 0;
}
