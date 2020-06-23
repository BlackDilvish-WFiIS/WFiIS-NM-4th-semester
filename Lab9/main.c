#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include</usr/include/gsl/gsl_linalg.h>
#include</usr/include/gsl/gsl_math.h>

double f(double x) { return sin(x); }
double der_f_0(int p) { return !(p%2) ? 0 : pow(-1, (p-1)/2); }
int factorial(int x)
{
    if(x == 0 || x == 1)
        return 1;
    else
        return x*factorial(x-1);
}

double R_NM(double x, double* a, double* b, int N, int M)
{
    double P=0;
    double Q=0;

    for(int i=0; i<=N; i++)
    {
        P += a[i]*pow(x, i);
        Q += b[i]*pow(x, i);
    }

    return P/Q;
}


int main()
{
    int N = 3;
    int M = 3;
    int n = N + M;

    double* C = malloc((n+1)*sizeof(double));

    for(int i=0; i<=n; i++)
        C[i] = der_f_0(i)/factorial(i);

    gsl_matrix* A = gsl_matrix_calloc(M, M);
    gsl_vector* Y = gsl_vector_calloc(M);
    gsl_vector* X = gsl_vector_calloc(M);

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<M; j++)
            gsl_matrix_set(A, i, j, C[N-M+i+j+1]);

        gsl_vector_set(Y, i, -C[N+1+i]);
    }

    gsl_permutation* p = gsl_permutation_calloc(M);
    int signum = 0;
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, Y, X);

    double* b = malloc((M+1)*sizeof(double));

    for(int i=0; i<M-1; i++)
        b[M-i] = gsl_vector_get(X, i);
    b[0] = 1;

    double* a = malloc((N+1)*sizeof(double));

    for(int i=0; i<=N; i++)
    {
        a[i] = 0;
        for(int j=0; j<=i; j++)
            a[i] += C[i-j]*b[j];
    }

    FILE* fp = fopen("MN3.dat", "w");

    for(double x=-6.28; x<=6.28; x += 0.1)
        fprintf(fp, "%lf\t%g\t%g\n", x, f(x), R_NM(x, a, b, N, M));
    
    fclose(fp);

    gsl_matrix_free(A);
    gsl_vector_free(Y);
    gsl_vector_free(X);
    gsl_permutation_free(p);
    free(a);
    free(b);
    free(C); 

    return 0;
}