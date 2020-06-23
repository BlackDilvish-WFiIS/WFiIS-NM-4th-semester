#include</usr/include/gsl/gsl_linalg.h>
#include</usr/include/gsl/gsl_math.h>
#include</usr/include/gsl/gsl_eigen.h>
#include<stdio.h>
#include<math.h>

#define L 10
#define n 200
#define N 1
#define DEL_X ((double)L/(n+1.0))

double rho(double x, int alpha)
{
  return 1.0 + 4.0*alpha*x*x;
}

int d_kron(int i, int j)
{
  return i==j ? 1 : 0;
}

void SaveLowestValues(FILE* fp, gsl_matrix* M)
{
  int i, j;
  for(i=0; i<n; i++)
      {
        double x_i = (-L/2.0) + DEL_X*(i+1);
        
        fprintf(fp, "%g ", x_i);
        for(j=0; j<6; j++)
          fprintf(fp, "%g ", gsl_matrix_get(M, i, j));
        fprintf(fp, "\n");
      }
}

int main()
{
  gsl_matrix* A = gsl_matrix_calloc(n,n);
  gsl_matrix* B = gsl_matrix_calloc(n,n);
  gsl_matrix* evec = gsl_matrix_calloc(n,n);
  
  gsl_vector* eval = gsl_vector_calloc(n);
  gsl_eigen_gensymmv_workspace* w = gsl_eigen_gensymmv_alloc(n);
  
  FILE *fp = fopen("najmniejsze_wartosci.dat", "w");
  FILE *f0 = fopen("najnizsze_wartosci_a0.dat", "w");
  FILE *f100 = fopen("najnizsze_wartosci_a100.dat", "w");
  
  int alpha, i, j;
  
  //Iteracja po alfa = [0, 100]
  for(alpha=0; alpha<=100; alpha+=2)
  {
    for(i=0; i<n; i++)
    {
      double x_i = (-L/2.0) + DEL_X*(i+1);
    
      //Uzupelnienie macierzy A i B
      for(j=0; j<n; j++)
      {
        double value_a = (-d_kron(i, j+1) + 2*d_kron(i, j) - d_kron(i, j-1)) / (DEL_X * DEL_X);
        gsl_matrix_set(A, i, j, value_a);
        
        double value_b = (rho(x_i, alpha)/N) * d_kron(i, j);
        gsl_matrix_set(B, i, j, value_b);
      }
    }
    
    //Rozwiazanie problemu oraz sortowanie wektorow i wartosci wlasnych
    gsl_eigen_gensymmv(A, B, eval, evec, w);
    gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    
    //6 najmniejszych wartosci wlasnych
    fprintf(fp, "%d ", alpha);
    for(i=0; i<6; i++)
      fprintf(fp, "%g ", sqrt(gsl_vector_get(eval,i)));
    fprintf(fp, "\n");
    
    //Wektory wlasne dla 6 najnizszych wartosci wlasnych
    if(alpha == 0)
      SaveLowestValues(f0, evec);    
    
    if(alpha == 100)
      SaveLowestValues(f100, evec);
  }
  
  fclose(fp);
  fclose(f0);
  fclose(f100);
  
  free(evec);
  free(w);
  
  free(A);
  free(B);
  free(eval);
  
  return 0;
}
