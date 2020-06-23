#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include</usr/include/gsl/gsl_linalg.h>
#include</usr/include/gsl/gsl_math.h>

#define STEPS 1000

double f1(double x) { return 1/(1+x*x); }
double f2(double x) { return cos(2*x); }

void wyznacz_M(double* xw, double* yw, double* m, int n, double alfa, double beta); 
double wyznacz_Sx (double* xw, double* yw, double* m, int n, double x);
void SaveInFile(const char* filename, double* x, double* y, double* s);
double SecondDerivativeF1(double x);

int main()
{
    int n=21;
    double XMIN = -5.0;
    double XMAX = 5.0;
    
    double* x_w = malloc(sizeof(double)*n);
    double* y_w1 = malloc(sizeof(double)*n);
    double* m1 = malloc(sizeof(double)*n);
    
    double* y_w2 = malloc(sizeof(double)*n);
    double* m2 = malloc(sizeof(double)*n);

    double* x_f = malloc(sizeof(double)*STEPS);
    double* y_f1 = malloc(sizeof(double)*STEPS);
    double* s_f1 = malloc(sizeof(double)*STEPS);
    
    double* y_f2 = malloc(sizeof(double)*STEPS);
    double* s_f2 = malloc(sizeof(double)*STEPS);

    double alpha = 0.0;
    double beta = 0.0;

    int n_list[3] = {5, 8, 21};
    const char* f_names[6] = {"f1n5.dat", "f1n8.dat", "f1n21.dat",
                              "f2n5.dat", "f2n8.dat", "f2n21.dat"};

    //Interpolacja
    for(int k=0; k<3; k++)
    {
      n=n_list[k];
      double step = (XMAX - XMIN)/(n-1);
      for(int i=0; i<n; i++)
      {
          x_w[i] = XMIN + i*step;
          y_w1[i] = f1(x_w[i]);
          y_w2[i] = f2(x_w[i]);
      } 
  
      wyznacz_M(x_w, y_w1, m1, n, alpha, beta);
      wyznacz_M(x_w, y_w2, m2, n, alpha, beta);
  
      double step_2 = (XMAX - XMIN)/(STEPS-1);
      for(int i=0; i<STEPS; i++)
      {
          x_f[i] = XMIN + i*step_2;
          y_f1[i] = f1(x_f[i]);
          s_f1[i] = wyznacz_Sx(x_w, y_w1, m1, n, x_f[i]);
          
          y_f2[i] = f2(x_f[i]);
          s_f2[i] = wyznacz_Sx(x_w, y_w2, m2, n, x_f[i]);
      }
  
      SaveInFile(f_names[k], x_f, y_f1, s_f1);
      SaveInFile(f_names[k+3], x_f, y_f2, s_f2);
    
    }
    
    //Porownanie pochodnych
    n=10;
    double step = (XMAX - XMIN)/(n-1);
    for(int i=0; i<n; i++)
    {
        x_w[i] = XMIN + i*step;
        y_w1[i] = f1(x_w[i]);
    } 
    
    wyznacz_M(x_w, y_w1, m1, n, alpha, beta);
    FILE* fp = fopen("Pochodne.dat", "w");
    for(int i=0; i<n; i++)
      fprintf(fp, "%g\t%g\t%g\n", x_w[i], m1[i], SecondDerivativeF1(x_w[i]));
    fclose(fp);
    

    free(x_w);
    free(y_w1);
    free(y_w2);
    free(m1);
    free(m2);
    free(x_f);
    free(y_f1);
    free(y_f2);
    free(s_f1);
    free(s_f2);

    return 0;
}

void wyznacz_M(double* xw, double* yw, double* m, int n, double alfa, double beta)
{
    double h = xw[1] - xw[0];
    double lambda = h / (h+h);
    double mi = 1 - lambda;

    gsl_matrix* A = gsl_matrix_calloc(n, n);
    gsl_vector* m_v = gsl_vector_calloc(n);
    gsl_vector* d = gsl_vector_calloc(n);

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            if((i==0 && i==0) || (i==n-1 && j==n-1))
                gsl_matrix_set(A, i, j, 1.0);
            else if(i==j)
                gsl_matrix_set(A, i, j, 2.0);
            else if(i-1 == j && i != n-1)
                gsl_matrix_set(A, i, j, mi);
            else if(i == j-1 && i!=0)
                gsl_matrix_set(A, i, j, lambda);
            else
                gsl_matrix_set(A, i, j, 0.0);
        }

        if(i>0 && i<n-1)
        {
            double d_val = (6.0/(h+h)) * (((yw[i+1] - yw[i]) / h) - ((yw[i] - yw[i-1]) / h));
            gsl_vector_set(d, i, d_val);
        }       
    }

    gsl_vector_set(d, 0, alfa);
    gsl_vector_set(d, n-1, beta);

    gsl_permutation* p = gsl_permutation_alloc(n);
    int signum = 0;
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, d, m_v);

    for(int i=0; i<n; i++)
        m[i] = gsl_vector_get(m_v, i);

    gsl_matrix_free(A);
    gsl_vector_free(m_v);
    gsl_vector_free(d);
    gsl_permutation_free(p);
}

double wyznacz_Sx (double* xw, double* yw, double* m, int n, double x)
{
    double h = xw[1] - xw[0];
    int i=1;
    for(i=0; i<n; i++)
    {
        if(x>=xw[i-1] && x<=xw[i])
            break;
    }

    double A = ((yw[i] - yw[i-1]) / h) - ((h / 6.0)*(m[i] - m[i-1]));
    double B = yw[i-1] - m[i-1] * (h*h / 6.0);

    double Sx = m[i-1]*(pow(xw[i]-x, 3)/(6*h)) + m[i]*(pow(x-xw[i-1], 3)/(6*h)) + A*(x-xw[i-1]) + B;
    return Sx;
}

void SaveInFile(const char* filename, double* x, double* y, double* s)
{
  FILE* fp = fopen(filename, "w");
  for(int i=0; i<STEPS; i++)
      fprintf(fp, "%g\t%g\t%g\n", x[i], y[i], s[i]);
  fclose(fp);
}

double SecondDerivativeF1(double x)
{
  double dx = 0.01;
  return (f1(x-dx) - 2*f1(x) + f1(x+dx))/(dx*dx);
}