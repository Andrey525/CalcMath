#include <math.h>
#include <stdio.h>

const double eps = 1E-7;
int n = 100;

double func(double x)
{
    return x / (pow(sin(2 * x), 3));
}

double medium_rectangle_method_integration(double a, double b)
{
    double s = 0.0;
    double h = (b - a) / n;
    for (int i = 0; i < n; i++) {
        s += func(a + h * (i + 0.5));
    }
    s *= h;
    return s;
}

double Runge_medium_rectangle_method_integration(double a, double b)
{
    int n0 = n, k;
    double sq[2], delta = 1;
    for (k = 0; delta > eps; n *= 2, k ^= 1) {
        double h = (b - a) / n;
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += func(a + h * (i + 0.5));
        }
        sq[k] = s * h;
        if (n > n0) {
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }
    return sq[k ^ 1];
}

int main()
{
    const double a = 0.1;
    const double b = 0.5;
    double S1, S2;
    S1 = medium_rectangle_method_integration(a, b);
    S2 = Runge_medium_rectangle_method_integration(a, b);
    printf("\n\n");
    printf("### medium_rectangle_method_integration ###\nthe integral(S1) = %0.12lf\n\n", S1);
    printf("### Runge_medium_rectangle_method_integration ###\nthe integral(S2) = %0.12lf\n\n", S2);
    return 0;
}
