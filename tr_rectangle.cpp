#include <iostream>
#include <fstream>

using namespace std;
int n = 10;
int m = 10;
int p = (n + 1) * (m + 1);
int N = 2 * n * m;
int R = 25;
int H = 25;

int surfacesavetostl(double* X, double* Y, double* Z, int** T);
int surfacesavetostl(double* X, double* Y, double* Z, int** T)
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
int triang_rectangle(double* X, double* Y,int** T,double a, double b, double c, double d)
{
	int k = 0;
	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			X[k] = a + i * (b - a) / n;
			Y[k] = c + j * (d - c) / m;
			k++;
		}
	}
	k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + j;
			T[k][2] = (m + 1) * i + j + 1;
			k++;
			T[k][0] = (m + 1) * i + j + 1;
			T[k][1] = (m + 1) * (i + 1) + j + 1;
			T[k][2] = (m + 1) * (i + 1) + j;
			k++;
		}
	}
    return 0;
}int main(){	double* Y = new double[p];
	double* Z = new double[p];
	int** T=new int*[N];	for (int i = 0; i < N; i++) {		T[i] = new int[3];	}	triang_rectangle(Y, Z,T, -R, R, 0,H);	double* X = new double[p];	for (int i = 0; i < p; i++) {		X[i] = sqrt(pow(R, 2) - pow(Y[i], 2));	}    surfacesavetostl(X, Y, Z, T);    delete[] Y;    delete[] Z;    delete[] X;    for (int i = 0; i < N; i++) {        delete[] T[i];    }    delete[] T;}