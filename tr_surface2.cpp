#include <iostream>
#include <fstream>
#include <math.h>
#define pi  3.14159265358979323846
using namespace std;




int a = 0;
double b = 2 * pi;
int n = 30;
int m = 10;
int savetostl(double* XX, double* YY, int** T);
int savetostl(double* XX, double* YY, int** T)
{
	int k = 0, ia, ib, ic;
	fstream TRstl;
	TRstl.open("triangulation.stl", ios::out | ios::app);
	TRstl.close();
	TRstl.open("triangulation.stl", ios::out | ios::in);
	TRstl << "solid <Triangulation>\n";
	for (k = 0; k < 2 * n * m; k++)
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
double F(double O) {
	return 0;
}

double Y(double O) {
	return 4;
}

double Gr(double O, double r) {
	return r * Y(O) + (1 - r) * F(O);
}



int main()
{
	double* On = new double[n + 1];
	for (int i = 0; i <= n; i++) {
		On[i] = a + i * (b - a) / n;
	}
	double* r = new double[m + 1];
	for (int i = 0; i <= m; i++) {
		r[i] = (double)i / m;

	}
	pair<double, double>** B = new pair<double, double> *[n + 1];
	for (int i = 0; i < n + 1; i++) {
		B[i] = new pair<double, double>[m + 1];
	}
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < m + 1; j++) {
			B[i][j].first = Gr(On[i], r[j]) * cos(On[i]);
			B[i][j].second = Gr(On[i], r[j]) * sin(On[i]);
			cout << "number" << i << "-" << j << "            " << "x=" << B[i][j].first << "             y=" << B[i][j].second << endl;
		}
	}

	double* X = new double[(n + 1) * (m + 1)];
	double* Y = new double[(n + 1) * (m + 1)];
	double* Z = new double[(n + 1) * (m + 1)];
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < m + 1; j++) {
			X[((m + 1) * i + j)] = B[i][j].first;
			Y[((m + 1) * i + j)] = B[i][j].second;
			Z[((m + 1) * i + j)] = On[i];
		}
	}
	int p = 2 * n * m;
	int** T = new int* [p];
	for (int i = 0; i < p; i++) {
		T[i] = new int[3];
	}
	int k = 0;
	for (int i = 0; i <= n - 1; i++) {
		for (int j = 0; j <= m - 1; j++) {
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + j;
			T[k][2] = (m + 1) * (i + 1) + j + 1;
			k++;
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + j + 1;
			T[k][2] = (m + 1) * (i)+j + 1;
			k++;
		}
	}

	surfacesavetostl(X, Y, Z, T, p);

	delete[] X;
	delete[] Y;
	delete[] Z;
	for (int i = 0; i < p; i++) {
		delete[] T[i];
	}
	delete[] T;
}
