#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double gen1();
double gen2();
double gen3();

void zad1();
void zad2();
void zad3(const char* filename, long int n);

int main()
{
    zad1();
    zad2();

    zad3("g_n2000.dat", 2000);
    zad3("g_n1e4.dat", 1e4);
    zad3("g_n1e7.dat", 1e7);

    return 0;
}

double gen1()
{
    static long int x = 10;
    int a = 17;
    long int m = pow(2, 13) - 1;
    x = (a*x) % m;
    return x / (m + 1.0);
}

double gen2()
{
    static long int x = 10;
    int a = 85;
    long int m = pow(2, 13) - 1;
    x = (a*x) % m;
    return x / (m + 1.0);
}

double gen3()
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

void zad1()
{
    const int N = 2000;
    FILE *fp1 = fopen("gen1.dat", "w");
    FILE *fp2 = fopen("gen2.dat", "w");
    FILE *fp3 = fopen("gen3.dat", "w");

    double* u1 = malloc(sizeof(double)*N);
    double* u2 = malloc(sizeof(double)*N);
    double* u3 = malloc(sizeof(double)*N);

    for(int i=0; i<N; i++)
    {
        u1[i] = gen1();
		u2[i] = gen2();
		u3[i] = gen3();
    }

    for(int i=0; i<N-3; i++)
    {
        fprintf(fp1, "%g %g %g %g\n", u1[i], u1[i+1], u1[i+2], u1[i+3]);
        fprintf(fp2, "%g %g %g %g\n", u2[i], u2[i+1], u2[i+2], u2[i+3]);
        fprintf(fp3, "%g %g %g %g\n", u3[i], u3[i+1], u3[i+2], u3[i+3]);
    }

    fprintf(fp1, "%g %g %g %g\n", u1[N-3], u1[N-2], u1[N-1], 0);
    fprintf(fp1, "%g %g %g %g\n", u1[N-2], u1[N-1], 0, 0);
    fprintf(fp2, "%g %g %g %g\n", u2[N-3], u2[N-2], u2[N-1], 0);
    fprintf(fp2, "%g %g %g %g\n", u2[N-2], u2[N-1], 0, 0);
    fprintf(fp3, "%g %g %g %g\n", u3[N-3], u3[N-2], u3[N-1], 0);
    fprintf(fp3, "%g %g %g %g\n", u3[N-2], u3[N-1], 0, 0);

    free(u1);
    free(u2);
    free(u3);

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
}

void zad2()
{
    FILE* fpSphere = fopen("sphere.dat", "w");
    FILE* fpBall = fopen("ball.dat", "w");

    double ri_s[3];
    double ri_b[3];
    const int N = 2000;

    for(int i=0; i<N; i++)
    {
        double u1 = gen3();
		double u2 = gen3();
		double u3 = gen3();
		double u4 = gen3();

        ri_s[0] = sqrt(-2 * log(1 - u1)) * cos(2 * M_PI * u2);
        ri_s[1] = sqrt(-2 * log(1 - u1)) * sin(2 * M_PI * u2);
        ri_s[2] = sqrt(-2 * log(1 - u3)) * cos(2 * M_PI * u4);

        double norm = sqrt(pow(ri_s[0], 2) + pow(ri_s[1], 2) + pow(ri_s[2], 2));

        ri_s[0] /= norm;
        ri_s[1] /= norm;
        ri_s[2] /= norm;

        fprintf(fpSphere, "%g %g %g\n", ri_s[0],  ri_s[1],  ri_s[2]);

		double si = pow(gen3(), 1/3.0);
		ri_b[0] = si * ri_s[0];
		ri_b[1] = si * ri_s[1];
		ri_b[2] = si * ri_s[2];	

		fprintf(fpBall, "%g %g %g\n", ri_b[0],  ri_b[1],  ri_b[2]);
    }

    fclose(fpSphere);
    fclose(fpBall);
}

void zad3(const char* filename, long int N)
{
    FILE* fp = fopen(filename, "w");

    double ri[3];
    double delta = 1/10.0;
    double nj[11] = {0};

    for (int i = 0; i < N; i++) 
    {
    	double u1 = gen3();
		double u2 = gen3();
		double u3 = gen3();
		double u4 = gen3();

		ri[0] = sqrt(-2 * log(1 - u1)) * cos(2 * M_PI * u2);
		ri[1] = sqrt(-2 * log(1 - u1)) * sin(2 * M_PI * u2);
		ri[2] = sqrt(-2 * log(1 - u3)) * cos(2 * M_PI * u4);

		double norm = sqrt(pow(ri[0], 2) + pow(ri[1], 2) + pow(ri[2], 2));

		ri[0] /= norm; 
		ri[1] /= norm; 
		ri[2] /= norm;

		double si = pow(gen3(), 1.0/3.0);
		ri[0] *= si;
		ri[1] *= si;
		ri[2] *= si;	

		double norm2 = sqrt(pow(ri[0], 2) + pow(ri[1], 2) + pow(ri[2], 2));

		int j = (int)(norm2 / delta) + 1;
		nj[j]++; 
    }

    for(int k = 1; k < 11; k++) 
    {
    	double Rj = delta * k;
        double Rj1 = delta * (k-1);

        double Vj = (4 / 3.0) * M_PI * pow(Rj, 3);
        double Vj1 = (4 / 3.0) * M_PI * pow(Rj1, 3);
        
        double gj = nj[k] / (Vj - Vj1);
        fprintf(fp, "%d %g %g\n", k, nj[k], gj);
    }

    fclose(fp);
}