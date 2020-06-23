#include<stdio.h>
#include<math.h>

double f(double x) { return log(pow(x, 5) + 3*x*x + x + 9); }
double g(double x) { return pow(x, 6); }

void findMin(double(*fun)(double), const char* filename, double xa, double xb, double lambda1, double lambda2, double x_teo)
{
    int k = 0;
    FILE *fp = fopen(filename, "w");

    while(fabs(xa - xb) > 10e-6)
    {
        double x1 = xa + lambda1*(xb - xa);
        double x2 = xa + lambda2*(xb - xa);

        if(fun(x1) < fun(x2))
            xb = x2;
        else
            xa = x1;

        double min = (xa + xb)/2;
        fprintf(fp, "%d\t%g\t%g\n", ++k, fabs(min - x_teo), min);
    }

    fclose(fp);
}

int main()
{
    double r = (sqrt(5)-1)/2;
    
    findMin(f, "zlotyF.dat", -0.5, 1.0, r*r, r, -0.1673198);
    findMin(f, "rowne3F.dat", -0.5, 1.0, 1./3, 2./3, -0.1673198);

    findMin(g, "zlotyG.dat", -4.0, 1.0, r*r, r, 0);
    findMin(g, "rowne3G.dat", -4.0, 1.0, 1./3, 2./3, 0);

    return 0;
}