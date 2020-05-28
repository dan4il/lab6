#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double Func(double* x) 
{
	return (-2 * (pow(x[0], 2)) + 20 * x[0] - (pow(x[1], 2)) + 16 * x[1] - 3 * (pow(x[2], 4)) - 48 * (pow(x[2], 3)) - 288 * (pow(x[2], 2)) - 768 * x[2] - 878);
}
//proizvodnaya x1
double Proizvx1(double* x) 
{
	return (-4 * x[0] + 20);
}
//proizvodnaya x2
double Proizvx2(double* x)
{
	return (-2 * x[1] + 16);
}
//proizvodnaya x3
double Proizvx3(double* x) 
{
	return (-12 * (pow(x[2], 3)) - 144 * (pow(x[2], 2)) - 576 * x[2] - 768);
}
//znachen v x+ap
double Tochka(double* x, double a, double* p)
{
	double xn[] = { x[0] + a * p[0], x[1], x[2] + a * p[2] };
	return -(Func(xn));
}
double SkolzOkoshko(double* x, double a, double* p)
{
	double d = 0.025;
	while (!((Tochka(x, a - d, p) >= Tochka(x, a, p)) && (Tochka(x, a, p) <= Tochka(x, a + d, p))))
	{
		if ((Tochka(x, a - d, p) >= Tochka(x, a, p)) && (Tochka(x, a, p) >= Tochka(x, a + d, p)))
			a += d / 2;
		else if (a > -d / 2)
			a -= d / 2;
	}
	return a;
}
//interpolation
double KvadrInt(double* x, double a, double* p)
{
	double d = 0.0000015;
	a = SkolzOkoshko(x, a, p);
	double a1 = a;
	double a2, a3, amin, ast;
	do {
		a2 = a1 + d;
		if (Tochka(x, a1, p) > Tochka(x, a2, p))
			a3 = a1 + 2 * d;
		else
			a3 = a1 - d;
		ast = 0.5 * ((a2 * a2 - a3 * a3) * Tochka(x, a1, p) + (a3 * a3 - a1 * a1) * Tochka(x, a2, p) + (a1 * a1 - a2 * a2) * Tochka(x, a3, p));
		ast = ast / ((a2 - a3) * Tochka(x, a1, p) + (a3 - a1) * Tochka(x, a2, p) + (a1 - a2) * Tochka(x, a3, p));
		if ((Tochka(x, a1, p) <= Tochka(x, a2, p)) && (Tochka(x, a2, p) <= Tochka(x, a3, p)))
			amin = a1;
		else if ((Tochka(x, a2, p) <= Tochka(x, a1, p)) && (Tochka(x, a1, p) <= Tochka(x, a3, p)))
			amin = a2;
		else
			amin = a3;
		if ((a1 <= ast) && (ast <= a3))
		{
			if (Tochka(x, amin, p) < Tochka(x, ast, p))
				a1 = amin;
			else
				a1 = ast;
		}
		else
			a1 = ast;
	} while ((abs(amin - ast)) >= 1e-5);
	return ast;
}
//sopryazhenie gradienti
void SoprGrad(double* x0)
{
	ofstream tr("ans1.dat", ios_base::out);
	double p[] = { Proizvx1(x0),Proizvx2(x0),Proizvx3(x0) };
	double b = 0;
	double a1 = 0.499;
	double a = 0;
	a = KvadrInt(x0, a, p);
	double x[] = { x0[0],x0[1],x0[2] };
	tr << x[0] << "," << x[1] << "," << x[2] << endl;
	double xn[] = { x[0] + a * p[0],x[1] + a * p[1],x[2] + a * p[2] };
	tr << xn[0] << "," << xn[1] << "," << xn[2] << endl;
	double pn[] = { Proizvx1(xn),Proizvx2(xn),Proizvx3(xn) };
	b = (pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]) / (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	pn[0] += b * p[0];
	pn[1] += b * p[1];
	pn[2] += b * p[2];
	a = KvadrInt(xn, a, pn);
	xn[0] += a * pn[0];
	xn[1] += a1 * pn[1];
	xn[2] += a * pn[2];
	tr << xn[0] << "," << xn[1] << "," << xn[2] << endl;
	tr.close();
	cout << "Tochka max funkcii: (" << xn[0] << ", " << xn[1] << ", " << xn[2] << ")" << endl;
	cout << "max funkcii: " << Func(xn) << endl;
}
int main()
{
	double x0[] = { -2, -4,  3 };
	SoprGrad(x0);
	return 0;
}

