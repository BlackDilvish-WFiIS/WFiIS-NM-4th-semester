#include<iostream>
#include<cmath>
#include<fstream>

//Wykozystalem uwage praktyczna z poprzedniej wersji zadania na stronie
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))
#define abs(X) ((X)>0? (X):-(X))

#define N 1000
#define m 10

template<typename T>
T* MultMatrixVector(T** M,T* v);

template<typename T>
T Dot(T* v1, T* v2);

template<typename T>
T* VectorsDiff(T* v1, T* v2);

template<typename T>
void NajwiekszySpadek(T defX, const char* fileName);

int main()
{
    NajwiekszySpadek(0.f, "Float0.dat");
    NajwiekszySpadek(1.f, "Float1.dat");
    NajwiekszySpadek(0.0, "Double0.dat");
    NajwiekszySpadek(1.0, "Double1.dat");

    return 0;
}

//Wykozystalem uwage praktyczna z poprzedniej wersji zadania na stronie
template<typename T>
T* MultMatrixVector(T** M,T* v)
{
    T *y = new T[N];

    for(int i=0;i<N;i++)
    {
        int jmin= max(0,i-m);
        int jmax= min(i+m,N-1);

        y[i]=0;

        for(int j=jmin;j<=jmax;j++)
            y[i] += M[i][j] * v[j];
    }

    return y;
}

//Iloczyn skalarny
template<typename T>
T Dot(T* v1, T* v2)
{
    T res = 0;

    for(int i=0; i<N; i++)
        res += v1[i] * v2[i];

    return res;
}

//Roznica wektorow
template<typename T>
T* VectorsDiff(T* v1, T* v2)
{
    T* res = new T[N];

    for(int i=0;i<N;i++)
        res[i] = v1[i] - v2[i];

    return res;
}

template<typename T>
void NajwiekszySpadek(T defX, const char* fileName)
{
    T** A = new T*[N];
    T* b = new T[N];
    T* x = new T[N];

    for(int i=0; i<N; i++)
    {
        A[i] = new T[N];

        for(int j=0;j<N;j++)
            if(abs(i-j) <= m)
                A[i][j] = 1.0 / (1.0 + abs(i-j));
            else
                A[i][j] = 0;

        b[i] = i;
        x[i] = defX;
    }

    int k=0;
    T* r = nullptr;
    std::ofstream f(fileName);

    do{
        delete[] r;
        k++;

        T* Ax = MultMatrixVector(A, x);
        r = VectorsDiff(b, Ax);
        T* Ar = MultMatrixVector(A, r);
        T alfa = Dot(r, r) / Dot(r, Ar);

        for(int i=0; i<N; i++)
            x[i] += alfa * r[i]; 

        f<< k << "\t"<< std::sqrt(Dot(r, r)) << "\t" << alfa << "\t" << std::sqrt(Dot(x, x)) << "\n";

        delete[] Ax;
        delete[] Ar;

    }while(std::sqrt(Dot(r, r)) > 1e-6 && k < 500);

    f.close();

    for(int i=0; i<N; i++)
        delete[] A[i];
    delete[] A;
    delete[] b;
    delete[] x;
}