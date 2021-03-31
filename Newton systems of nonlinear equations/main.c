#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define E 0.0001
int n = 3;

double function(int i, double* x)
{
    if (i == 0) {
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;
    }
    if (i == 1) {
        return 2 * x[0] * x[0] + x[1] * x[1] - 4 * x[2];
    }
    if (i == 2) {
        return 3 * x[0] * x[0] - 4 * x[1] + x[2] * x[2];
    }
    return 0;
}

void dgemm(double** a, double* b, double* c, int n)
{
    for (int i = 0; i < n; i++) {
        c[i] = 0;
        for (int j = 0; j < n; j++) {
            c[i] += a[i][j] * b[j];
        }
    }
}

double sec_opred(double* a)
{
    return a[0] * a[3] - a[2] * a[1];
}
double third_opred(double** a)
{
    return a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2];
}

int Jacobian_calculation(double** a, double* b, int n)
{
    double* nea[n];
    for (int i = 0; i < n; i++) {
        nea[i] = malloc(sizeof(double) * n);
    }
    a[0][0] = 2 * b[0], a[0][1] = 2 * b[0], a[0][2] = 2 * b[0]; // составляем матрицу
    a[1][0] = 4 * b[1], a[1][1] = 2 * b[1], a[1][2] = -4;
    a[2][0] = 6 * b[2], a[2][1] = -4, a[2][2] = 2 * b[2];
    double delta = third_opred(a); // считаем определитель матрицы
    if (delta != 0) {
        double temp[4];
        int num = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    for (int l = 0; l < n; l++) {
                        if (k != i && l != j) {
                            temp[num] = a[k][l], num++;
                        }
                    }
                }
                num = 0;
                nea[i][j] = sec_opred(temp);
            }
        }
        double tmp = 0;
        for (int i = 0; i < n; i++) { // транспонируем
            for (int j = i; j < n; j++) {
                if (i != j) {
                    tmp = nea[i][j];
                    nea[i][j] = nea[j][i];
                    nea[j][i] = tmp;
                }
            }
        }
        nea[0][1] *= -1, nea[1][0] *= -1, nea[1][2] *= -1, nea[2][1] *= -1;
        for (int i = 0; i < n; i++) { // окончательное нахождение Якобиана
            for (int j = 0; j < n; j++) {
                a[i][j] = nea[i][j] / delta;
            }
        }
        return 1;
    }
    return 0;
}
void Newton_SNE(double* X)
{
    double present[n];
    double past[n];
    for (int i = 0; i < n; i++) {
        present[i] = X[i];
        past[i] = present[i];
    }
    double F[n];
    double* a[n];
    double* c = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        a[i] = malloc(sizeof(double) * n);
        for (int j = 0; j < n; j++) {
            a[i][j] = 0;
        }
    }
    int flag = 1;
    do {
        for (int i = 0; i < n; i++) {
            past[i] = present[i];
            F[i] = function(i, past);
        }
        if (Jacobian_calculation(a, past, n)) {
            dgemm(a, F, c, n);
            for (int i = 0; i < n; i++)
                present[i] = past[i] - c[i];
            flag = 0;
            while (fabs(present[flag] - past[flag]) <= E && flag != n) {
                flag++;
            }
            if (flag == n) {
                flag = 0;
            } else {
                flag = 1;
            }
        } else {
            return;
        }
    } while (flag);
    for (int i = 0; i < n; i++) {
        printf("%lf\n", present[i]);
    }
}

int main()
{
    double* x = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        x[i] = 0.5;
    }
    Newton_SNE(x);
    return 0;
}
