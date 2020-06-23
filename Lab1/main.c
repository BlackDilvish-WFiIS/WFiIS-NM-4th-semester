#include<stdio.h>
#include<math.h>
#include</usr/include/gsl/gsl_math.h>
#include</usr/include/gsl/gsl_linalg.h>


int main()
{
    int n = 4;

    gsl_matrix *a = gsl_matrix_calloc(n, n);
    gsl_permutation *p = gsl_permutation_calloc(n);
    gsl_matrix *aa = gsl_matrix_calloc(n,n);

    int i,j;
    for(i =0; i<n; i++)
    {
        for(j=0;j<n;j++)
        {
            double value = 1.0/(i + j + 2);
            gsl_matrix_set(a, i, j, value);
	    gsl_matrix_set(aa,i,j,value);
        }
    }

    int signum;
    gsl_linalg_LU_decomp(a, p, &signum);

    FILE *fp;

    fp = fopen("det.dat", "w");
    for(i =0; i<n; i++)
    {
        double value = gsl_matrix_get(a, i, i);
        fprintf(fp, "%g \n", value);
    }

    double wyzn = 1;
    for(i = 0; i<n; i++)
        wyzn *= gsl_matrix_get(a, i, i);
    fprintf(fp, "Wyznacznik: %g\n", wyzn);

    fclose(fp);

    gsl_vector *b = gsl_vector_calloc(n);
    gsl_vector *x = gsl_vector_calloc(n); // wektor rozw
    gsl_matrix *c = gsl_matrix_calloc(n, n);
    int k;

    for(k = 0; k<n; k++)
    {
        for(i=0;i<n;i++)
            gsl_vector_set(b, i, 0.0);

        gsl_vector_set(b, k, 1.0);

        gsl_linalg_LU_solve(a, p, b, x); //rozwiazanie czescowe (jedna z kolumn)

        for(j=0;j<n;j++)
            gsl_matrix_set(c, j, k, gsl_vector_get(x, j)); //uzupelnianie macierzy
    }

    fp = fopen("macierzC.dat", "w");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            fprintf(fp,"%g\t", gsl_matrix_get(c, i, j));
        fprintf(fp,"\n");
    }

   fclose(fp);

   fp = fopen("iloczyn.dat", "w");

   gsl_matrix *d = gsl_matrix_calloc(n, n);

   for(i=0;i<n;i++)
   {
	   for(j=0;j<n;j++)
	   {
		gsl_matrix_set(d,i,j,0);
		for(k=0;k<n;k++)
		{
			double val = gsl_matrix_get(d,i,j) + gsl_matrix_get(aa,i,k) * gsl_matrix_get(c,k,j);
			gsl_matrix_set(d, i, j, val);
		}
		fprintf(fp,"%g ",gsl_matrix_get(d, i, j));
	   }
	   fprintf(fp,"\n");
   }

   fclose(fp);

   fp = fopen("wsk.dat", "w");
    double m1=0,m2=0;

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            double c1 = gsl_matrix_get(aa, i, j);
            double c2 = gsl_matrix_get(c, i, j);

            if(fabs(c1) > m1)
                m1 = fabs(c1);

            if(fabs(c2) > m2)
                m2 = fabs(c2);
        }
    }

   fprintf(fp, "||A|| = %g, ||A-1|| = %g, K = %g", m1, m2, m1*m2);
   fclose(fp);


    return 0;
}
