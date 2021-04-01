#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double eps = 0.0001;

double function(double x, double y)
{
    return pow(x, 2) - 2 * y;
}

double dif_Eiler(double a, double b, int n, double y0)
{
    double *y, *x, *func;
    double h;
    y = malloc(sizeof(double) * n);
    x = malloc(sizeof(double) * n);
    func = malloc(sizeof(double) * n);
    h = (b - a) / n; // шаг

    x[0] = a;
    y[0] = y0;
    func[0] = function(x[0], y[0]);
    for (int i = 1; i < n; i++) {
        x[i] = a + h * i;
        y[i] = y[i - 1] + h * func[i - 1];
        func[i] = function(x[i], y[i]);
    }
    // printf("\n\nSolution of the differential equation by the Eiler method:\n\n");
    // for (int i = 0; i < n; i++) {
    //     printf("x[%d] = %.1lf  y[%d] = %lf  func[%d] = %lf  hfunc[%d] = %lf\n", i, x[i], i, y[i], i, func[i], i, h * func[i]);
    // }
    return y[n - 1];
}

void Per_Eiler(double a, double b, int n, double y0)
{
    int n0 = n, k;
    double sq[2], delta = 1;
    for (k = 0; delta > eps; n *= 2, k ^= 1) {

        sq[k] = dif_Eiler(a, b, n, y0);
        if (n > n0) {
            delta = fabs(sq[k] - sq[k ^ 1]);
        }
    }
    printf("Eiler sq[%d] = %lf\n", k, sq[k]);
    // return sq[k ^ 1];
}

double dif_Runge_Kutta(double a, double b, int n, double y0)
{
    double *y, *x, *k1, *k2, *k3, *k4, *delta_y;
    double h;
    h = (b - a) / n; // шаг
    y = malloc(sizeof(double) * n);
    x = malloc(sizeof(double) * n);
    k1 = malloc(sizeof(double) * n);
    k2 = malloc(sizeof(double) * n);
    k3 = malloc(sizeof(double) * n);
    k4 = malloc(sizeof(double) * n);
    delta_y = malloc(sizeof(double) * n);
    x[0] = a;
    y[0] = y0;
    for (int i = 0; i < n; i++) {
        x[i] = a + h * i;
        k1[i] = function(x[i], y[i]);
        k2[i] = function(x[i] + h / 2, y[i] + (h * k1[i]) / 2);
        k3[i] = function(x[i] + h / 2, y[i] + (h * k2[i]) / 2);
        k4[i] = function(x[i] + h, y[i] + h * k3[i]);
        delta_y[i] = (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        y[i + 1] = y[i] + delta_y[i];
    }
    // printf("\n\nSolution of the differential equation by the Runge-Kutta method:\n\n");
    // for (int i = 0; i < n; i++) {
    //     printf("x[%d] = %.1lf  y[%d] = %lf  k1[%d] = %lf  k2[%d] = %lf  k3[%d] = %lf  k4[%d] = %lf  delta_y[%d] = %lf\n", i, x[i], i, y[i], i, k1[i], i, k2[i], i, k3[i], i, k4[i], i, delta_y[i]);
    // }
    return y[n - 1];
}

void Per_Runge_Kutta(double a, double b, int n, double y0)
{
    int n0 = n, k;
    double sq[2], delta = 1;
    for (k = 0; delta > eps; n *= 2, k ^= 1) {

        sq[k] = dif_Runge_Kutta(a, b, n, y0);
        if (n > n0) {
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }
    printf("Runge Kutta sq[%d] = %lf\n", k, sq[k]);
}

int main()
{
    // dif_Eiler(0.0, 1.0, 10, 1);
    Per_Eiler(0.0, 1.0, 10, 1);
    // dif_Runge_Kutta(0.0, 1.0, 10, 1);
    Per_Runge_Kutta(0.0, 1.0, 10, 1);
    return 0;
}
