#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define STEPS 1000

double f(double x)
{
    return exp(-x*x);
}

double Interpolate(double* X, double* Y, int n, double x)
{
    double w = 0;
    int i,j;

    for(j=0; j<=n; j++)
    {
        double numerator = 1.0;
        double denominator = 1.0;

        for(i=0; i<=n; i++)
        {
            numerator *= x - X[i];

            if(i != j)
                denominator *= X[j] - X[i];
        }

        w += (Y[j]/(x-X[j])) * (numerator/denominator);
    }

    return w;
}

double* ClassicNodes(double x0, double x1, int n)
{
    double* X = malloc(sizeof(double)*(n+1));
    double step = (x1 - x0)/n;
    int i;

    for(i=0; i<=n; i++)
        X[i] = x0 + i*step;

    return X;
}

double* CzebyszewNodes(double x0, double x1, int n)
{
    double* X = malloc(sizeof(double)*(n+1));
    double step = (x1 - x0)/n;
    int m;

    for(m=0; m<=n; m++)
        X[m] = 0.5*((x1 - x0)*cos(3.14*(2*m + 1)/(2*n + 2)) + (x0 + x1));

    return X;
}

double* NodeValues(double* Nodes, int n)
{
    double* Y = malloc(sizeof(double)*(n+1));
    int i;

    for(i=0; i<=n; i++)
        Y[i] = f(Nodes[i]);

    return Y;
}

void SaveComparision(const char* filename, double *X, double* Y, double *X_n, double *Y_n, int n)
{
    FILE* fp = fopen(filename, "w");
    int i;
    double *LagrangeY = malloc(sizeof(double)*STEPS);

    for(i=0; i<STEPS; i++)
    {
        LagrangeY[i] = Interpolate(X_n, Y_n, n, X[i]);

        fprintf(fp, "%g\t%g\t%g\n", X[i], Y[i], LagrangeY[i]);
    }

    free(LagrangeY);
}

int main()
{
    double step = 10.0 / STEPS;
    int i, n;

    double* X = malloc(sizeof(double)*STEPS);
    double* Y = malloc(sizeof(double)*STEPS);

    for(i=0; i<STEPS; i++)
    {
        X[i] = -5 + i*step;
        Y[i] = f(X[i]);
    }

    const char* names[8] = {"data5.dat", "data10.dat", "data15.dat", "data20.dat",
                            "data5c.dat", "data10c.dat", "data15c.dat", "data20c.dat"};
    for(n=5; n<=20; n+=5)
    {
        double* X_n = ClassicNodes(-5, 5, n);
        double* Y_n = NodeValues(X_n, n);
        SaveComparision(names[n/5 - 1], X, Y, X_n, Y_n, n);

        free(X_n);
        free(Y_n);

        X_n = CzebyszewNodes(-5, 5, n);
        Y_n = NodeValues(X_n, n);
        SaveComparision(names[n/5 + 3], X, Y, X_n, Y_n, n);

        free(X_n);
        free(Y_n);
    }

    free(X);
    free(Y);

    return 0;
}