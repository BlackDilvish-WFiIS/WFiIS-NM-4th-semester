#include<stdio.h>
#include<math.h>
#include <stdlib.h>

double f(int m, int k, double x)
{
    return pow(x, m)*sin(k*x);
}

double factorial(int k)
{
    if(k==0 || k==1)
        return 1;
    else
        return (double)k*factorial(k-1);
}

double I(int l, int a, int b, int m, int k)
{
    double s1=0, s2=0;

    for(int i=0; i<l; i++)
    {
        s1 += (pow(-1, i)*pow(0*k, 2*i+m+2))/(pow(k, m+1)*factorial(2*i+1)*(2*i+m+2));
        s2 += (pow(-1, i)*pow(M_PI*k, 2*i+m+2))/(pow(k, m+1)*factorial(2*i+1)*(2*i+m+2));
    }

    return s2 - s1;
}

double S(int n, double a, double b, int m, int k)
{
    double step = (b-a)/(n-1);
    double* nodes = malloc(sizeof(double)*n);

    for(int i=0; i<n; i++)
        nodes[i] = a + i*step;

    double sum = 0;

    for(int i=0; i<n-1; i++)
        sum += f(m, k, nodes[i]) + 4*f(m, k, (nodes[i]+nodes[i+1])/2) + f(m, k, nodes[i+1]);

    return (step/6)*sum;
}

void SaveI(int m, int k, const char* filename, double teo)
{
    FILE* fp = fopen(filename, "w");

    for(int l=1; l<=30; l++)
    {
        double C = I(l, 0, M_PI, m, k);
        fprintf(fp, "%d %g %g\n", l, C, fabs(C - teo));
    }
    fclose(fp);
}

void SaveS(int m, int k, const char* filename, double teo)
{
    FILE* fp = fopen(filename, "w");

    int n[5] = { 11, 21, 51, 101, 201 };

    for(int i=0; i<5; i++)
    {
        double C = S(n[i], 0, M_PI, m, k);
        fprintf(fp, "%d %lf %g\n", n[i], C, fabs(C - teo));
    }
    fclose(fp);
}

int main()
{
    SaveI(0, 1, "Im0k1.dat", 2);
    SaveI(1, 1, "Im1k1.dat", M_PI);
    SaveI(5, 5, "Im5k5.dat", 56.363569);

    SaveS(0, 1, "Sm0k1.dat", 2);
    SaveS(1, 1, "Sm1k1.dat", M_PI);
    SaveS(5, 5, "Sm5k5.dat", 56.363569);
    
    return 0;
}
