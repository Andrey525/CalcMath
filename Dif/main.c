#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double function(double x, double y)
{
    return pow(x, 2) - 2 * y;
}

void dif_Eiler(double a, double b, double n, double y0)
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
    for (int i = 1; i <= n; i++) {
        x[i] = a + h * i;
        y[i] = y[i - 1] + h * func[i - 1];
        func[i] = function(x[i], y[i]);
    }
    for (int i = 0; i <= n; i++) {
        printf("x[%d] = %.1lf  y[%d] = %lf  func[%d] = %lf  hfunc[%d] = %lf\n", i, x[i], i, y[i], i, func[i], i, h * func[i]);
    }
}

int main()
{
    dif_Eiler(0.0, 1.0, 10, 1);
    return 0;
}
