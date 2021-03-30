#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double function(double x)
{
    return 8 * x - 3;
}

void linear_approximation(double* x, double* y, int n, double* a, double* b)
{
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x_pow_2 = 0;
    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x_pow_2 += pow(x[i], 2);
    }
    *a = (n * sum_xy - sum_x * sum_y) / (n * sum_x_pow_2 - pow(sum_x, 2));
    *b = (sum_y - *a * sum_x) / n;
}

int main()
{
    srand(time(NULL));
    int n = 10; // количество точек
    double *x, *y;
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        x[i] = i;
        y[i] = function(x[i]) + (rand() % 100 - 50) * 0.007;
    }
    double a, b;
    linear_approximation(x, y, n, &a, &b);
    for (int i = 0; i < n; i++) {
        printf("X[%d] = %lf  Y[%d] = %lf\n", i, x[i], i, y[i]);
    }
    printf("a = %lf  b = %lf\n\n", a, b);
    return 0;
}
