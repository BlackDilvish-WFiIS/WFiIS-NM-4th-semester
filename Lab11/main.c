#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include</usr/include/gsl/gsl_errno.h>
#include</usr/include/gsl/gsl_fft_complex.h>

#define max(X,Y) ((X)>(Y)?(X):(Y))

double y_0(double omega, int i)
{
    return sin(omega*i) + sin(2*omega*i) + sin(3*omega*i);
}

double noise()
{
    return 2*(rand()/(RAND_MAX + 1.0)) - 1;
}

double y(double omega, int i)
{
    return y_0(omega, i) + noise();
}

int main()
{
    int stride = 1;
    int k = 12;
    double N = pow(2, k);
    double omega = (4*M_PI)/N;

    double* dane = malloc(2*N*sizeof(double));
    double* niezaburzony = malloc(N*sizeof(double));
    double* zaburzony = malloc(N*sizeof(double));

    for(int i=0; i<N; i++)
    {
        dane[2*i] = y(omega, i); //dt=1
        dane[2*i + 1] = 0;
        
        niezaburzony[i] = y_0(omega, i);
        zaburzony[i] = dane[2*i];
    }

    gsl_fft_complex_radix2_forward(dane, stride, N);


    double maxck = dane[0];
    
    FILE* fp = fopen("k12.dat", "w");
    for(int i=0; i<N; i++)
    {
        fprintf(fp, "%d %g %g\n", i, dane[2*i], dane[2*i + 1]);
        maxck = max(maxck, max(dane[2*i], dane[2*i + 1]));
    }
    fclose(fp);
    
    double border = maxck/2;

    fp = fopen("k12cz2.dat", "w");
    for(int i=0; i<N; i++)
    {
        double modulo = sqrt(dane[2*i]*dane[2*i] + dane[2*i + 1]*dane[2*i + 1]);
        fprintf(fp, "%d %g %g\n", i, modulo, border);
        
        if(dane[2*i] < border)
          dane[2*i] = 0;
          
        if(dane[2*i+1] < border)
          dane[2*i+1] = 0;
    }
    fclose(fp);
    
    gsl_fft_complex_radix2_backward(dane, stride, N);
    
    fp = fopen("k12cz3.dat", "w");
    for(int i=0; i<N; i++)
    {
        double unnoised = dane[2*i]/(N/2);
        
        fprintf(fp, "%d %g %g %g\n", i, zaburzony[i], niezaburzony[i], unnoised);
    }
    fclose(fp);
    
    free(dane);
    free(niezaburzony);
    free(zaburzony);

    return 0;
}