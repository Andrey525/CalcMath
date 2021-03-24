#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double Lagrange(double* x, double* y, int count, double X)
{
    double Y = 0;
    for (int i = 0; i < count; i++) {
        double polynomial_term = 1;
        for (int j = 0; j < count; j++) {
            if (i != j) {
                polynomial_term *= (X - x[j]) / (x[i] - x[j]);
            }
        }
        Y += y[i] * polynomial_term;
    }
    return Y;
}

double Aitken(double* x, double* y, int count, double X)
{
    double P[count][count];
    for (int j = 0; j < count; j++) {
        P[0][j] = y[j];
    }
    for (int i = 1; i < count; i++) {
        for (int j = i; j < count; j++) {
            P[i][j] = (P[i - 1][i - 1] * (x[j] - X) - P[i - 1][j] * (x[i - 1] - X)) / (x[j] - x[i - 1]);
        }
    }
    return P[count - 1][count - 1];
}

void table_of_difference(double** dif, double* x, double* y, int count)
{
    int n = count - 1;
    int row = 0;
    for (int i = 0; i < n; i++) {
        dif[row][i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    }
    row++;
    n--;
    while (row < count - 1) {
        for (int j = 0; j < n; j++) {
            dif[row][j] = (dif[row - 1][j + 1] - dif[row - 1][j]) / (x[row + j + 1] - x[j]);
        }
        row++;
        n--;
    }
}

double Newton(double* x, double* y, int count, double X)
{
    double Y = y[0];
    double* dif[count - 1];
    for (int i = 0; i < count - 1; i++) {
        dif[i] = malloc(sizeof(double) * (count - 1));
    }

    table_of_difference(dif, x, y, count);
    for (int i = 0; i < count - 1; i++) {
        double temp = dif[i][0];
        for (int j = 0; j < i + 1; j++) {
            temp *= (X - x[j]);
        }
        Y += temp;
    }
    return Y;
}

double function(double X)
{
    return (pow(X, 3) - sqrt(pow(X, 2) - 2)) / 6;
}

int main()
{
    double Y, X, Y2, Y3;
    double x[] = { 2, 3, 5, 8, 11, 14, 18, 23, 26, 30 };
    double* y = malloc(sizeof(x));
    for (int i = 0; i < sizeof(x) / sizeof(double); i++) {
        y[i] = function(x[i]);
        printf("y[%d] = %lf, при x[%d] = %lf\n", i, y[i], i, x[i]);
    }
    X = 12;
    Y = Lagrange(x, y, sizeof(x) / sizeof(double), X);
    printf("Ответ, полученный интерполяцией Лагранжа: Y = %.2lf, при X = %.2lf\n", Y, X);
    Y2 = Aitken(x, y, sizeof(x) / sizeof(double), X);
    printf("Ответ, полученный интерполяцией Эйткена: Y2 = %.2lf, при X = %.2lf\n", Y2, X);
    Y3 = Newton(x, y, sizeof(x) / sizeof(double), X);
    printf("Ответ, полученный интерполяцией Ньютона: Y3 = %.2lf, при X = %.2lf\n", Y3, X);
    return 0;
}