#include <iostream>
#include <fstream>
#define M_PI       3.14159265358979323846
using namespace std;

int n = 4;
int m = 15;
int NumbPoint(int i, int j, int n) {
    int N;
    if (i == 0 && j == 0) {
        N = 0;
    }
    else {
        N = n * j * (j - 1) / 2 + i + 1;
    }
    return N;
}
int savetostl(double* XX, double* YY, int** T);
int savetostl(double* XX, double* YY, int** T)
{
    int k = 0, ia, ib, ic;
    fstream TRstl;
    TRstl.open("triangulation.stl", ios::out | ios::app);
    TRstl.close();
    TRstl.open("triangulation.stl", ios::out | ios::in);
    TRstl << "solid <Triangulation>\n";
    for (k = 0; k < n * m * m; k++)
    {
        ia = T[k][0];
        ib = T[k][1];
        ic = T[k][2];
        TRstl << "facet normal " << 0.0 << " " << 0.0 << " " << 1.0 << "\n";
        TRstl << "outer loop\n";
        TRstl << "vertex ";
        TRstl << XX[ia] << " " << YY[ia] << " " << 0.0 << "\n";
        TRstl << "vertex ";
        TRstl << XX[ib] << " " << YY[ib] << " " << 0.0 << "\n";
        TRstl << "vertex ";
        TRstl << XX[ic] << " " << YY[ic] << " " << 0.0 << "\n";
        TRstl << "endloop\n";
        TRstl << "endfacet\n";
    }
    TRstl << "endsolid";
    TRstl.close();
    return 0;
}
int surfacesavetostl(double* X, double* Y, double* Z, int** T, int N)
{
    int k = 0, ia, ib, ic;
    double v, x1, y1, z1, x2, y2, z2, nx, ny, nz;
    fstream TRstl;
    remove("triangulation.stl");
    TRstl.open("triangulation.stl", ios::out | ios::app);
    TRstl << "solid <Triangulation>\n";
    for (k = 0; k < N; k++)
    {
        ia = T[k][0];
        ib = T[k][1];
        ic = T[k][2];
        x1 = X[ib] - X[ia];
        y1 = Y[ib] - Y[ia];
        z1 = Z[ib] - Z[ia];
        x2 = X[ic] - X[ia];
        y2 = Y[ic] - Y[ia];
        z2 = Z[ic] - Z[ia];
        nx = (y1 * z2 - y2 * z1);
        ny = (z1 * x2 - x1 * z2);
        nz = (x1 * y2 - x2 * y1);
        v = sqrt(nx * nx + ny * ny + nz * nz);
        nx = nx / v;
        ny = ny / v;
        nz = nz / v;
        TRstl << "facet normal " << nx << " " << ny << " " << nz << "\n";
        TRstl << "outer loop\n";
        TRstl << "vertex ";
        TRstl << X[ia] << " " << Y[ia] << " " << Z[ia] << "\n";
        TRstl << "vertex ";
        TRstl << X[ib] << " " << Y[ib] << " " << Z[ib] << "\n";
        TRstl << "vertex ";
        TRstl << X[ic] << " " << Y[ic] << " " << Z[ic] << "\n";
        TRstl << "endloop\n";
        TRstl << "endfacet\n";
    }
    TRstl << "endsolid";
    TRstl.close();
    return 0;
}
double f(int i, int j) {
    return (double)2* M_PI * i / (j * n);
}
double r(double O) {
    return 5;
}
int main()
{
    int p = 1 + n * m * (m + 1) / 2;
    double* X = new double[p];
    double* Y = new double[p];
    double* Z = new double[p];

    X[0] = 0;
    Y[0] = 0;
    Z[0] = 0;

    for (int j = 1; j <= m; j++) {
        for (int i = 0; i <= n * j - 1; i++) {
            X[NumbPoint(i, j, n)] = j * r(f(i, j)) * cos(f(i, j)) / m;
            Y[NumbPoint(i, j, n)] = j * r(f(i, j)) * sin(f(i, j)) / m;
            Z[NumbPoint(i, j, n)] = f(i,j);
        }
    }
    int** T = new int* [n * m * m];
    for (int i = 0; i < n * m * m; i++) {
        T[i] = new int[3];
    }
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == n - 1) {
            /*T[i][0] = 0;
            T[i][1] = i + 1;
            T[i][2] = 1;
            k++; */
        }
        else {
           /* T[i][0] = 0;
            T[i][1] = i + 1;
            T[i][2] = i + 2;
            k++;*/
        }
    }
    int count = 0;
    for (int s = 1; s <= n; s++) {
        for (int j = 1; j < m; j++) {
            for (int i = (s - 1) * j; i <= s * j; i++) {
                if (i == s * j) {
                   
                    
                        T[k][0] = NumbPoint(i, j, n);
                        T[k][1] = NumbPoint(i + s - 1, j + 1, n);
                        T[k][2] = NumbPoint(i + s, j + 1, n);
                        k++;
                        if (s == n) {
                            cout<< NumbPoint(i, j, n)<<endl;
                            cout << NumbPoint(i + s - 1, j + 1, n) << endl;
                            cout << NumbPoint(i + s, j + 1, n) << endl;
                            cout << "||||||||||||||||||||||||||||"<<endl;
                    }
                    
               
                }
                else {
                    
                        T[k][0] = NumbPoint(i, j, n);
                        T[k][1] = NumbPoint(i + s - 1, j + 1, n);
                        T[k][2] = NumbPoint(i + s, j + 1, n);
                        k++;
                        T[k][0] = NumbPoint(i, j, n);
                        T[k][1] = NumbPoint(i + s, j + 1, n);
                        T[k][2] = NumbPoint(i + 1, j, n);
                        k++;
                        if (s == n) {
                            cout<< NumbPoint(i, j, n)<<endl;
                            cout << NumbPoint(i + s - 1, j + 1, n) << endl;
                            cout << NumbPoint(i + s, j + 1, n) << endl;
                            cout<< NumbPoint(i, j, n)<<endl;
                            cout<< NumbPoint(i + s, j + 1, n)<<endl;
                            cout<< NumbPoint(i + 1, j, n)<<endl;
                            cout << "||||||||||||||||||||||||||||"<<endl;
                        }
                        
                        
                   
                }
            }
        }
    }
    /*int s = n;
    int b = 0;
    for (int j = 1; j < m; j++) {
        for (int i = (s - 1) * j; i <= s * j; i++) {
            if (i == s * j - 1) {
                T[k][0] = NumbPoint(s * j - 1, j, n);
                T[k][1] = NumbPoint(s * (j + 1) - 2, j + 1, n);
                T[k][2] = NumbPoint(s * (j + 1) - 1, j + 1, n);

                k++;
                b++;

                T[k][0] = NumbPoint(s * j - 1, j, n);
                T[k][1] = NumbPoint(s * (j + 1) - 1, j + 1, n);
                T[k][2] = NumbPoint(0, j, n);
                k++;
                b++;

            }
            else {
                if (i == s * j) {
                    T[k][0] = NumbPoint(0, j, n);
                    T[k][1] = NumbPoint(s * (j + 1) - 1, j + 1, n);
                    T[k][2] = NumbPoint(0, j + 1, n);
                    k++;
                    b++;
                }
                else {
                    T[k][0] = NumbPoint(i, j, n);
                    T[k][1] = NumbPoint(i + n - 1, j + 1, n);
                    T[k][2] = NumbPoint(i + n, j + 1, n);
                    k++;
                    b++;

                    T[k][0] = NumbPoint(i, j, n);
                    T[k][1] = NumbPoint(i + n, j + 1, n);
                    T[k][2] = NumbPoint(i + 1, j, n);
                    k++;
                    b++;

                }
            }
        }
    }*/
    surfacesavetostl(X, Y, Z, T, n*m*m-n-14);
    
    cout << "END"<<count;
    delete[] X;
    delete[] Y;
    delete[] Z;
    for (int i = 0; i < n * m * m; i++) {
        delete[] T[i];
    }
    delete[] T;
}
