#include<stdio.h>

#define xb 2
#define xa 0.5
//#define N 50
#define N 500

int main()
{
  double h = 2*xb / (N-1.0);
  int i, j;

  double a[N+1], d[N+1], c[N+1], x[N+1], rho[N+1];

  //Wyznaczanie elementow macierzy trojdiagonalnej,
  //siatki x oraz ro
  for(i=1;i<N+1;i++) 
  {
	  x[i] = -xb + h * i;
	  d[i] = -2.0 / (h*h);
	  a[i] = c[i] = 1.0/ (h*h);
	  
	  if(x[i] >= -xa && x[i] < 0)
		  rho[i] = 1;
	  else if(x[i] > 0 && x[i] <= xa)
		  rho[i] = -1;
	  else
		  rho[i] = 0;
	  
  }
  d[1] = d[N] = 1;
  c[1] = c[N] = 0;
  
  double u[N+1], l[N+1], y[N+1], v[N+1]; 
  
  //Wyznaczenie elementow macierzy L i U
  u[1] = d[1];
  y[1] = -rho[1];
  for(i = 2 ;i <N+1;i++)
  {
	  l[i] = a[i] / u[i-1];
	  u[i] = d[i] - l[i] * c[i-1];
	  
	  y[i] = -rho[i] - l[i] * y[i-1];
  }
    
  //Wyznaczenie wektora V
  v[N] = y[N]/u[N];
  for(i=N-1; i>=0; i--)
  {
    v[i] = (y[i] - c[i]*v[i+1]) / u[i];
  }

  //Obliczenie wlasciwego wektora V
  double Vrozw[N+1];
  for(int i=1;i<N+1;i++)
  {
    if(x[i]>=-xb && x[i] < -xa)
      Vrozw[i] = x[i]/16 + 1/8.0;
    else if(x[i]>=-xa && x[i] < 0)
      Vrozw[i] = -(x[i] * x[i])/2 - (7.0/16)*x[i];
    else if(x[i]>=0 && x[i] < xa)
      Vrozw[i] = (x[i]*x[i])/2 - (7.0/16)*x[i];
    else
      Vrozw[i] = x[i]/16 - 1/8.0;

  }
    
  FILE* fp = fopen("dane500.dat","w");
  for(i = 1; i < N+1; i++)
    fprintf(fp,"%g\t%g\t%g\n",x[i], v[i], Vrozw[i]);
  
  fclose(fp);
  

  
  return 0;
}
