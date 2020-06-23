#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define I_TEO 0.2557840245

double U1();
double Q1();

double g(double x, double y)     { return sin(x + y)/log(2 + x + y); }
double f(double x)               { return exp(-x); }

double z(double x, double y)     { return g(x, y) * f(x) * f(y); }
double gCxCy(double x, double y) { return g(x, y)*pow(1 - exp(-1), 2); }

void calc();
void testQ(const char* filename);

int main()
{
    calc("zad1.dat", z, U1);
    calc("zad2.dat", gCxCy, Q1);
    testQ("zad3.dat");

    return 0;
}

double U1()
{
    static long int x1  = 10;
	static long int x2  = 10;
	static long int x3  = 10;	

	int a0 = 1176;
	int a1 = 1476;
	int a2 = 1776;

	long int m = pow(2, 32) - 5;
	long int xi = (a0 * x1 + a1 * x2 + a2 * x3) % m;

    x1 = xi;
    x2 = x1;
	x3 = x2;
	
	return xi / (m + 1.0);
}

double Q1()
{
    return -log(1 - U1()*(1 - exp(-1)));
}

double variant(double* z_x, int n)
{
    double sum1 = 0;
    double sum2 = 0;

    for(int i=1; i<=n; i++)
    {
        sum1 += pow(z_x[i], 2);
        sum2 += z_x[i];
    }

    return (sum1 - pow(sum2, 2)/n)/(n-1);
}

void calc(const char* filename, double(*fun)(double, double), double(*gen)())
{
    FILE* fp = fopen(filename, "w");
    const long int n = 1e5;
    double* z_x = malloc(n * sizeof(double));
    double sum = 0;

    for(int i=1; i<=n; i++)
    {
        double x = gen();
        double y = gen();

        z_x[i] = fun(x, y);
        sum += z_x[i];

        if( i == 10 || i == 100 || i == 1e3 || i == 1e4 || i == 1e5)
        {
            double I = sum / i;

            double deviation = variant(z_x, i) / sqrt(n);

            fprintf(fp, "%d %g %g\n", i, I, deviation);
        }
    }

    fclose(fp);
}

void testQ(const char* filename)
{
    FILE* fp = fopen(filename, "w");
    FILE* fp2 = fopen("chi.dat", "w");
    const long int n = 1e5;
    double n_i[10] = {0};

    for(int i=0; i<n; i++)
    {
        double x = Q1();

        int index = (int)(x * 10);

        n_i[index]++;
    }

    double chi = 0;

    for(int i=1; i<11; i++)
    {
        double P = (1 - exp(-i/10.0))/(1 - exp(-1)) - (1 - exp(-((i-1)/10.0)))/(1 - exp(-1));
        double p = n_i[i-1] / n;

        chi += pow(n_i[i-1] - n * P, 2) / (n * P);

        fprintf(fp, "%d %g %g\n", i, p, P);
    }

    fprintf(fp2, "%g", chi);

    fclose(fp);
    fclose(fp2);
}
