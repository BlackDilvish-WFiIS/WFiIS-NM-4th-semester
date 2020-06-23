#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 7
#define IT_MAX 12
//#define IT_MAX 30

double** CreateMatrix(int n);
double* CreateVector(int n);
double Norm(double* V);
double Dot(double* V1, double* V2);
double* MultMatrixVector(double** M, double* V);
double** MultMatrixMatrix(double** M1, double** M2);
double** TransposeMatrix(double** M);


int main()
{
    double** A = CreateMatrix(N);
    double** X = CreateMatrix(N);   
    double** W = CreateMatrix(N);

    double* X_old = CreateVector(N);
    double* X_new;

    double lambda = 0;
    FILE* fpl = fopen("lambdy.dat", "w");

    int i,j,k,m;

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
        {
            A[i][j] = (1 + fabs(i+j)) / (1 + fabs(i-j));
            W[i][j] = A[i][j];
        }

    for(k=0; k<N; k++)
    {
        for(i=0; i<N; i++)
            X_old[i] = 1;

        for(i=0; i<N; i++)
            X_old[i] = X_old[i]/Norm(X_old);

        for(m=1; m<=IT_MAX; m++)
        {
            X_new = MultMatrixVector(W, X_old);

            if(fabs(lambda - Dot(X_new, X_old)) < 10e-6) //Warunek STOP
                break;

            lambda = Dot(X_new, X_old);
            fprintf(fpl, "%g ", lambda);

            for(i=0; i<N; i++)
                X_old[i] = X_new[i]/Norm(X_new);

            
        }
        fprintf(fpl, "\n");

        for(i=0; i<N; i++)
            for(j=0; j<N; j++)
                W[i][j] -= lambda * X_old[i] * X_old[j];

        for(i=0; i<N; i++)
            X[i][k] = X_old[i];
    }

    double** D = MultMatrixMatrix(MultMatrixMatrix(TransposeMatrix(X), A), X);


    FILE* fpd = fopen("macierzD.dat", "w");
    FILE* fpx = fopen("macierzX.dat", "w");

    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            fprintf(fpd, "%g ", D[i][j]);
            fprintf(fpx, "%g ", X[i][j]);
        }
        fprintf(fpd, "\n");
        fprintf(fpx, "\n");
    }

    fclose(fpd);
    fclose(fpx);
    fclose(fpl);

    return 0;
}

double** CreateMatrix(int n)
{
    int i,j;
    double** A = malloc(sizeof(double*) * n);

    for(i=0; i<n; i++)
    {
        A[i] = malloc(sizeof(double) * n);
        for(j=0; j<n; j++)
            A[i][j] = 0;
    }

    return A;
}

double* CreateVector(int n)
{
    double* V = malloc(sizeof(double) * n);
    int i;
    for(i=0; i<N; i++)
        V[i] = 0;

    return V;
}

double Norm(double* V)
{
    return sqrt(Dot(V, V));
}

double Dot(double* V1, double* V2)
{
    double dot = 0;
    int i = 0;
    
    for(; i<N; i++)
        dot += V1[i] * V2[i];

    return dot;
}

double* MultMatrixVector(double** M,double* V)
{
    double *W = CreateVector(N);
    int i,j;

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            W[i] += M[i][j] * V[j];

    return W;
}

double** MultMatrixMatrix(double** M1, double** M2)
{
    int i,j,k;
    double** A = CreateMatrix(N);

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            for(k=0; k<N; k++)
                A[i][j] += M1[i][k] * M2[k][j]; 

    return A;
}


double** TransposeMatrix(double** M)
{
    int i,j;
    double** A = CreateMatrix(N);

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            A[i][j] = M[j][i];

    return A;
}