#include<stdio.h>
#include<math.h>

double g1(double x)           { return sin(x);        }
double g2(double x)           { return (x*x)/8;       }
double f(double x)            { return g1(x) - g2(x); }
double derivative_f(double x) { return cos(x)-(x/4);  }

void save_function(const char* filename, double (*fun)(double), double start, double end);
void sieczne(const char* filename, int it_max, double x0, double x1);
void newton(const char* filename, int it_max, double x0);

int main()
{
    save_function("g1.dat", g1, -8.0, 8.0);
    save_function("g2.dat", g2, -8.0, 8.0);
    save_function("f.dat", f, -8.0, 8.0);

    sieczne("sieczne-8.dat", 15, -8.0, -8.1);
    newton("newton-8.dat", 10, -8.0);

    sieczne("sieczne8.dat", 15, 8.0, 8.1);
    newton("newton8.dat", 10, 8.0);

    return 0;
}

void save_function(const char* filename, double (*fun)(double), double start, double end)
{
    FILE* fp = fopen(filename, "w");
    double i;

    for(i=start; i<end; i+=0.1)
        fprintf(fp, "%g %g\n", i, fun(i));

    fclose(fp);
}

void sieczne(const char* filename, int it_max, double x0, double x1)
{
    FILE* fp = fopen(filename, "w");
    int k;
    for(k=1; k<=it_max; k++)
    {
        double x2 = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0));
        fprintf(fp, "%d %g %g %g\n", k, x2, f(x1), f(x0));
        x0 = x1;
        x1 = x2;      
    }
    fclose(fp);
}

void newton(const char* filename, int it_max, double x0)
{
    FILE* fp = fopen(filename, "w");
    int k;
    for(k=1; k<=it_max; k++)
    {
        x0 = x0 - (f(x0)/derivative_f(x0));
        fprintf(fp, "%d %g %g %g\n", k, x0, f(x0), derivative_f(x0));
    }
}