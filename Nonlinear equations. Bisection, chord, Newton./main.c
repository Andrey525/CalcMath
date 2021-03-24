#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define E 0.0001

double function(double x)
{
    return cos(2 / x) - 2 * sin(1 / x) + 1 / x;
    // return x - cos(x);
}

double first_order_derivative_of_function(double x)
{
    return (2 * sin(2 / x) / (x * x)) + (2 * cos(1 / x) / (x * x)) - 1 / (x * x);
    // return 1 + sin(x);
}

double second_order_derivative_of_function(double x)
{
    return 2 / (x * x * x) - (4 * sin(2 / x) / (x * x * x)) - (4 * cos(2 / x) / (x * x * x * x)) - (4 * cos(1 / x) / (x * x * x)) + (2 * sin(1 / x) / (x * x * x * x));
    // return cos(x);
}

void bisection(double a, double b, double* answer)
{
    printf("Bisection\n");
    double c;
    if (function(a) * function(b) < 0) {
        printf("Условие на сходимость выполнено\n");
        c = (b + a) / 2;
        int iteration = 1;
        printf("iteration\t\ta\t\tb\t\tf(a)\t\tf(b)\n");
        while (fabs(function(c)) > E) {
            c = (b + a) / 2;
            if (function(a) * function(c) < 0) {
                b = c;
            } else {
                a = c;
            }
            printf("%9d\t    %lf\t    %lf\t     %lf\t     %lf\n", iteration, a, b, function(a), function(b));
            iteration++;
            *answer = c;
        }
    } else {
        printf("Условие на сходимость не выполнено\n");
    }
}

void newton(double a, double b, double* answer)
{
    printf("Newton\n");
    int iteration = 1;
    double x;
    int flag = 0;
    if (function(a) * second_order_derivative_of_function(a) > 0) {
        printf("Выполняется условие на сходимость для a\n");
        flag = 1;
        x = a;
    } else if (function(b) * second_order_derivative_of_function(b) > 0) {
        printf("Выполняется условие на сходимость для b\n");
        flag = 1;
        x = b;
    } else {
        printf("Условия на сходимость не выполнены\n");
    }
    if (flag) {
        printf("iteration\t\tx\t\tf(x)\t\tf'(x)\n");
        while (fabs(function(x)) > E) {
            x = x - function(x) / first_order_derivative_of_function(x);
            printf("%9d\t    %lf\t    %lf\t     %lf\n", iteration, x, function(x), first_order_derivative_of_function(x));
            iteration++;
        }
        *answer = x;
    }
}

void chord(double a, double b, double* answer)
{
    printf("Chord\n");
    int iteration = 1;
    double x, temp, previous_x;
    int flag = 0;
    if (function(a) * second_order_derivative_of_function(a) > 0) {
        printf("Выполняется условие на сходимость для a\n");
        flag = 1;
    } else if (function(b) * second_order_derivative_of_function(b) > 0) {
        printf("Выполняется условие на сходимость для b\n");
        flag = 1;
    } else {
        printf("Условия на сходимость не выполнены\n");
    }
    if (flag) {
        x = b;
        previous_x = a;
        printf("iteration\t\tx\t\tprev_x\t\tf(x)\t\tf(prevx)\n");
        while (fabs(x - previous_x) > E) {
            temp = x;
            x = x - function(x) * ((x - previous_x) / (function(x) - function(previous_x)));
            previous_x = temp;
            printf("%9d\t    %lf\t    %lf\t     %lf\t     %lf\n", iteration, x, previous_x, function(x), function(previous_x));
            iteration++;
        }
        *answer = x;
    }
}

int main()
{
    double* answer_bisection = malloc(sizeof(double));
    double* answer_newton = malloc(sizeof(double));
    double* answer_chord = malloc(sizeof(double));
    printf("\n\n\n");
    bisection(1.3, 3.14, answer_bisection);
    printf("\n\n\n");
    newton(1.3, 3.14, answer_newton);
    printf("\n\n\n");
    chord(1.3, 3.14, answer_chord);
    printf("\n\n\n");
    printf("Answer of bisection = %lf\n", *answer_bisection);
    printf("Answer of newton = %lf\n", *answer_newton);
    printf("Answer of chord = %lf\n", *answer_chord);
    return 0;
}
