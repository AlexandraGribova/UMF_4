#include <iostream>
#include<vector>
#include<fstream>
using namespace std;

int maxiter, N;
double eps, nev;
vector<double> r, z, s, p, boof, boof1, L, U;
vector<double> X, d, di, ggl, ggu, pr;
vector<int> ig, jg;

void input(string file, int N, vector <double>& mas)
{
	ifstream input;
	input.open(file + ".txt");
	for (int i = 0; i < N; i++)
		input >> mas[i];
	input.close();
}
void input_int(string file, int N, vector <int>& mas)
{
	ifstream input;
	input.open(file + ".txt");
	for (int i = 0; i < N; i++)
		input >> mas[i];
	input.close();
}

void r0_s()
{
	int i, j, k;
	int i0, i1;
	for (i = 0; i < N; i++)
	{
		i0 = ig[i]-1;
		i1 = ig[i + 1]-1;
		r[i] = pr[i] - d[i] * X[i];//f-Ax0
		for (k = i0; k < i1; k++)
		{
			j = jg[k]-1;
			r[i] -= ggl[k] * X[j];
			r[j] -= ggu[k] * X[i];
		}
	}
}

void L_1(vector <double> pr, vector <double>& r)
{
	double s;
	int	i0, i1, k;
	for (int i = 0; i < N; i++) {
		s = pr[i];
		i0 = ig[i]-1;
		i1 = ig[i + 1]-1;
		for (k = i0; k < i1; k++)
			s -= L[k] * r[jg[k]-1];
		r[i] = s / di[i];
	}
}

void U_1(vector <double> pr, vector <double>& z)
{
	int	i0, i1, j, k;
	z = pr;
	for (int i = N - 1; i >= 0; i--)
	{
		i0 = ig[i]-1;
		i1 = ig[i + 1]-1;
		z[i] = z[i] / di[i];
		for (k = i0; k < i1; k++)
		{
			j = jg[k]-1;
			z[j] -= U[k] * z[i];
		}
	}
}

void AVec(vector <double> x, vector <double>& y)
{
	for (int i = 0; i < N; i++) {
		y[i] = d[i] * x[i];
		int i0 = ig[i]-1;
		int i1 = ig[i + 1]-1;
		for (int k = i0; k < i1; k++) {
			int j = jg[k]-1;
			y[i] += ggl[k] * x[j];
			y[j] += ggu[k] * x[i];
		}
	}
}

void AVecTrunsp(vector <double> x, vector <double>& y)
{
	for (int i = 0; i < N; i++) {
		y[i] = d[i] * x[i];
		int i0 = ig[i]-1;
		int i1 = ig[i + 1]-1;
		for (int k = i0; k < i1; k++) {
			int j = jg[k]-1;
			y[i] += ggu[k] * x[j];
			y[j] += ggl[k] * x[i];
		}
	}
}

void input()
{
	maxiter = 1500;
	eps = 1e-16;
	ifstream input;
	input.open("input.txt");
	input >> N>>maxiter>>eps;
	input.close();
}

void X_k(double a, vector <double> z)
{
	for (int i = 0; i < N; i++)
		X[i] += a * z[i];
}

void Vec_k(vector <double> &vec, double a, vector <double> p)
{
	for (int i = 0; i < N; i++)
		vec[i] -= a * p[i];
}

void Vec_k_plus (double b, vector <double> some, vector <double> &vec, vector <double> vec_0)
{
	for (int i = 0; i < N; i++)
		vec[i] = some[i] + b * vec_0[i];
}

double scalar_mult(vector <double> v1, vector <double> v2, int size)
{
	double s = 0;
	for (int i = 0; i < size; i++)
		s += v1[i] * v2[i];
	return s;
}

double Norm(vector<double> X)
{
	double norma = 0;
	for (int i = 0; i < N; i++)
	{
		norma += X[i] * X[i];
	}
	norma = sqrt(norma);
	return norma;
}

void BCG()
{
	double nr, nf;
	double a, b;
	vector<double> boof(N), boof1(N);
	vector<double> p_0, r_0, z_0, s_0;
	r0_s();//r=f-A*x
	p = r;
	z = r;
	s = r;
	nr = Norm(r);
	nf = Norm(pr);
	int i;
	for (i = 0; i < maxiter; i++)
	{
		if (nr <= eps * nf)
			break;
		p_0 = p;
		r_0 = r;
		z_0 = z;
		s_0 = s;
		//U_1(z, boof);//boof=(U^-1)*z
		AVec(z, boof);//boof1=A*boof 
		//L_1(boof1, boof);//boof=(L^-1)*boof1
		a = scalar_mult(p, r, N) / scalar_mult(s, boof ,N);//a=(p,r)/(s,boof)  (3.35)
		X_k(a, z);
		Vec_k(r, a, boof);//r=r-a*LAUz
		
		//U_1(s, boof);//boof=(U^-1)*s
		AVecTrunsp(s, boof);//boof1=A^t*boof 
		//L_1(boof1, boof);//boof=(L^-1)*boof1
		Vec_k(p, a, boof);//r=r-a*LA^tUz
		b = scalar_mult(p, r, N) / scalar_mult(p_0, r_0, N);
		Vec_k_plus(b, r, z, z_0);
		Vec_k_plus(b, p, s, s_0);
		nr = Norm(r);
	}
	for (int j = 0; j < N; j++)
		cout << X[j] << endl;
}

void LUS_factorisation()
{
	int i, j, k;
	int i0, i1, ki, kj;
	double suml, sumu, sumd;

	for (i = 0; i < N; i++)
	{
		i0 = ig[i] - 1;
		i1 = ig[i + 1] - 1;
		sumd = 0;
		for (k = i0; k < i1; k++)
		{
			j = jg[k] - 1;
			ki = i0;
			kj = ig[j] - 1;//что такое kj - соответсвующий номер элемента для домножения
			suml = sumu = 0;
			while (ki < k)
				if (jg[kj] == jg[ki]) {
					sumu += U[ki] * L[kj];
					suml += U[kj] * L[ki];
					kj++;
					ki++;
				}
				else
					if (jg[kj] < jg[ki])
						kj++;
					else
						ki++;
			U[k] = (ggu[k] - sumu) / di[j];
			L[k] = (ggl[k] - suml) / di[j];
			sumd += U[k] * L[k];
		}
		di[i] = sqrt(d[i] - sumd);
	}

}

int main()
{
	setlocale(LC_ALL, "Russian");
	input();
	X.resize(N);
	pr.resize(N);
	di.resize(N);
	d.resize(N);
	p.resize(N);
	r.resize(N);
	z.resize(N);
	s.resize(N);
	ig.resize(N + 1);
	input_int("ig", N + 1, ig);
	input("di", N, d);
	input("pr", N, pr);
	int size = ig[N] - 1;
	L.resize(size);
	U.resize(size);
	jg.resize(size);
	ggl.resize(size);
	ggu.resize(size);
	input_int("jg", size, jg);
	input("ggl", size, ggl);
	input("ggu", size, ggu);
	LUS_factorisation();
	BCG();
}