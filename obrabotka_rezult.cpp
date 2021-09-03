#include <iostream>
#include <cmath>
#include <time.h>
#define O1 1
#define O2 2 
#define O3 3
#define g2 0.1 // сигма квадрат
#define dt 0.1 // шаг времени
#define T 10 // интервал времени
#define CN 3 // предназначено для массивов 3 на 3 или 3 на T\dt
#define cn 2 // предназначено для массивов 2 на 2 или 2 на T\dt
#define ZC1 0 // это первая оценка 
#define ZC2 1 //а это вторая оценка 
using namespace std;
   
const float pi = 3.1415926535;// число пи
float RandomFloat();
float normalRand();//нормальное распределение
double y(double);
double e();
void FreedomMem(double **matr, int n);
void Get_matr(double **matr, int n, double **temp_matr, int indRow, int indCol);
double** TransponMtx(double **matr, double **tMatr, int n);// транспонирование матрицы квадратной
double** TransponMtx(double **, double **, int, int);// транспонирование матрицы любой размерности
double Det(double **matr, int n);// вычисление детерминанта
void  obratnya_matr(double, double**, int, double **&);
double **multy(double**, double**, int, int, int, double **&);// перемножение матриц
void jacobi(const unsigned int n, double * const * a, double * d, double * const * v);//нахождение собственных значений и векторов матрицы
int main()
{
	double **z = new double *[T / dt], **C = new double *[T / dt],
		**TC = new double *[CN], **Mat_C = new double *[CN],
		*C_VEC = new double[cn], **C_ZN = new double *[cn],
		**rez = new double *[1];
	double **ans = new double *[CN];
	int j = 0;
	rez[0] = new double[CN];
	C_ZN[0] = new double[cn];
	C_ZN[1] = new double[cn];
	for (double i = 0.0; i < T; i += dt)
	{
		if (j < CN)
		{
			Mat_C[j] = new double[3];
			ans[j] = new double[1];
			TC[j] = new double[T / dt];
		}
		z[j] = new double[1];
		C[j] = new double[3];
		C[j][0] = 1;
		C[j][1] = i;
		C[j][2] = i * i;
		z[j][0] = y(i) + e();
		j++;
	}
	multy(C, TransponMtx(C, TC, T / dt, 3), 3, T / dt, 3, Mat_C);//С*C транспонированная
	FreedomMem(C, T / dt);
	C = new double *[CN];
	for (int i = 0; i < CN; i++)
	{
		C[i] = new double[CN];
		for (int j = 0; j < CN; j++)
			C[i][j] = 0;
	}
	obratnya_matr(Det(Mat_C, CN), Mat_C, CN, C);// обратная матрица от С*C транспонированная размеронсти 3 на 3
	multy(z, TC, 1, T / dt, CN, rez);// C транспонированная * Z
	multy(C, rez, CN, CN, 1, ans);// обратная матрица * Z
	printf("Omnk = %f %f %f\n", ans[0][0], ans[1][0], ans[2][0]);
	FreedomMem(C, CN);
	C = new double *[T / dt];
	j = 0;
	for (double i = 0.0; i < T; i += dt)
	{
		C[j] = new double[cn];
		C[j][0] = TC[ZC1][j];
		C[j][1] = TC[ZC2][j];
		j++;
	}
	FreedomMem(Mat_C, CN);
	FreedomMem(TC, CN);
	TC = new double *[cn];
	Mat_C = new double *[cn];
	Mat_C[0] = new double[cn];
	Mat_C[1] = new double[cn];
	TC[0] = new double[T / dt];
	TC[1] = new double[T / dt];
	multy(C, TransponMtx(C, TC, T / dt, cn), cn, T / dt, cn, Mat_C);//С*C транспонированная
	FreedomMem(C, T / dt);
	C = new double *[cn];
	for (int i = 0; i < cn; i++)
	{
		C[i] = new double[cn];
		for (int j = 0; j < cn; j++)
			C[i][j] = 0;
	}
	obratnya_matr(Det(Mat_C, cn), Mat_C, cn, C);// обратная матрица от С*C транспонированная размеронсти 2 на 2
	cout << endl;
	jacobi(cn, C, C_VEC, C_ZN); // вычисление собственных векторов и значений обратной матрицы
	cout << "Собственные вектора матрицы :" << endl;
	for (int i = 0; i < cn; i++)
	{
		for (j = 0; j < cn; j++)
			cout << C_ZN[i][j] << " ";
		cout << endl;
	}
	cout << endl;
	cout << "Собственные значения матрицы :" << endl;
	for (int i = 0; i < cn; i++)
	{
		cout << C_VEC[i] << " ";
		cout << endl;
	}
	FreedomMem(C, cn);
	FreedomMem(Mat_C, cn);
	FreedomMem(rez, 1);
	FreedomMem(ans, CN);
	FreedomMem(z, CN);
	FreedomMem(TC, cn);
	system("pause");
	return 0;
}
double y(double t)
{
	return O1 + O2 * t + O3 * t*t;
}
double e()
{
	srand(time(NULL));
	return g2 * normalRand();
}
void  obratnya_matr(double det, double** matr, int n, double **& obr_matr)
{

	if (det) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int m = n - 1;
				double **temp_matr = new double *[m];
				for (int k = 0; k < m; k++)
					temp_matr[k] = new double[m];
				Get_matr(matr, n, temp_matr, i, j);
				obr_matr[i][j] = pow(-1.0, i + j + 2) * Det(temp_matr, m) / det;
				FreedomMem(temp_matr, m);
			}
		}
	}
	else
		cout << "Т.к. определитель матрицы = 0,\nто матрица вырожденная и обратной не имеет!!!" << endl;
}

double** TransponMtx(double **matr, double **tmatr, int row, int col)
{

	for (int i = 0; i < col; i++)
	{
		for (int j = 0; j < row; j++)
		{
			tmatr[i][j] = matr[j][i];
		}
	}
	return tmatr;
}
float normalRand()
{
	float z = RandomFloat();
	float f = RandomFloat();

	float r = cos(2 * pi * f) * sqrt(-2 * log(z));

	return r;
}
double** TransponMtx(double **matr, double **tMatr, int n) {//
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			tMatr[j][i] = matr[i][j];
	return tMatr;
}
//Функция освобождения памяти

void FreedomMem(double **matr, int n)
{
	for (int i = 0; i < n; i++)
		delete[] matr[i];
	delete[] matr;
}
void Get_matr(double **matr, int n, double **temp_matr, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matr[ki][kj] = matr[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}
double Det(double **matr, int n)
{
	double temp = 0;   //временная переменная для хранения определителя
	int k = 1;      //степень
	if (n < 1) {
		cout << "Не верный размер матрицы!!!" << endl;
		return 0;
	}
	else if (n == 1)
		temp = matr[0][0];
	else if (n == 2)
		temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
	else {
		for (int i = 0; i < n; i++) {
			int m = n - 1;
			double **temp_matr = new double *[m];
			for (int j = 0; j < m; j++)
				temp_matr[j] = new double[m];
			Get_matr(matr, n, temp_matr, 0, i);
			temp = temp + k * matr[0][i] * Det(temp_matr, m);
			k = -k;
			FreedomMem(temp_matr, m);
		}
	}
	return temp;
}
float RandomFloat()
{
	return (float)rand() / RAND_MAX;
}

double **multy(double** a, double** b, int nrow, int row, int col, double **& c)
{
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < col; j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < row; k++)
			{
				c[i][j] += a[k][i] * b[j][k];

			}
		}
	}
	return c;
}
void jacobi(const unsigned int n, double * const * a, double * d, double * const * v)
{
	if (n == 0) return;
	double * b = new double[n + n];
	double * z = b + n;
	unsigned int i, j;
	for (i = 0; i < n; ++i)
	{
		z[i] = 0.;
		b[i] = d[i] = a[i][i];
		for (j = 0; j < n; ++j)
		{
			v[i][j] = i == j ? 1. : 0.;
		}
	}
	for (i = 0; i < 50; ++i)
	{
		double sm = 0.;
		unsigned int p, q;
		for (p = 0; p < n - 1; ++p)
		{
			for (q = p + 1; q < n; ++q)
			{
				sm += fabs(a[p][q]);
			}
		}
		if (sm == 0) break;
		const double tresh = i < 3 ? 0.2 * sm / (n*n) : 0.;
		for (p = 0; p < n - 1; ++p)
		{
			for (q = p + 1; q < n; ++q)
			{
				const double g = 1e12 * fabs(a[p][q]);
				if (i >= 3 && fabs(d[p]) > g && fabs(d[q]) > g) a[p][q] = 0.;
				else
					if (fabs(a[p][q]) > tresh)
					{
						const double theta = 0.5 * (d[q] - d[p]) / a[p][q];
						double t = 1. / (fabs(theta) + sqrt(1. + theta * theta));
						if (theta < 0) t = -t;
						const double c = 1. / sqrt(1. + t * t);
						const double s = t * c;
						const double tau = s / (1. + c);
						const double h = t * a[p][q];
						z[p] -= h;
						z[q] += h;
						d[p] -= h;
						d[q] += h;
						a[p][q] = 0.;
						for (j = 0; j < p; ++j)
						{
							const double g = a[j][p];
							const double h = a[j][q];
							a[j][p] = g - s * (h + g * tau);
							a[j][q] = h + s * (g - h * tau);
						}
						for (j = p + 1; j < q; ++j)
						{
							const double g = a[p][j];
							const double h = a[j][q];
							a[p][j] = g - s * (h + g * tau);
							a[j][q] = h + s * (g - h * tau);
						}
						for (j = q + 1; j < n; ++j)
						{
							const double g = a[p][j];
							const double h = a[q][j];
							a[p][j] = g - s * (h + g * tau);
							a[q][j] = h + s * (g - h * tau);
						}
						for (j = 0; j < n; ++j)
						{
							const double g = v[j][p];
							const double h = v[j][q];
							v[j][p] = g - s * (h + g * tau);
							v[j][q] = h + s * (g - h * tau);
						}
					}
			}
		}
		for (p = 0; p < n; ++p)
		{
			d[p] = (b[p] += z[p]);
			z[p] = 0.;
		}
	}
	delete[] b;
}