#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#define PI 3.141592653589793
#define N 512

double complex f(double t){
	if ((int)creal(t)%2 == 0){
		return 1.;
	}
	else{
		return -1.;
	}
}


void fourier(double complex *x, double complex *X){
	for (int k = 0; k < N; ++k)
	{
		X[k] = 0;
		for (int n = 0; n < N; ++n)
		{
			X[k] += x[n] * cexp(-2 * I * PI * n * k / N);
		}
	}
}

void fft(double complex *x, double complex *X, int n, int s){
	if (n == 1){
		X[0] = x[0];
	} else{

		fft(x, X, n/2, 2*s);
		fft(x + s, X + n/2, n/2, 2*s);

		for (int k = 0; k < n/2; ++k)
		{
			double complex p = cexp(-2 * I * PI * k / n);
			double complex temp = X[k];

			X[k] = temp + p * X[k+n/2];
			X[k + n/2] = temp - p * X[k+n/2];
		}
	}
}

int reverseBin(int i, int K){
	int j = 0;

	for (int k=K-1; k>-1; k--){
		j += i / (1<<k) * (1<<(K-1-k));
		i %= 1<<k;
	}
	return j;
}


	

void fft_iter(double complex *x, double complex *X){
	
	int logN = 0;
	int n = N;
	while (n >>= 1) ++logN;

	for (int i = 0; i < N; ++i)
	{
		X[i] = x[reverseBin(i,logN)];
	}

	for (int i = 1; i < logN+1; ++i)
	{
		int m = 1<<i;
		double complex p = cexp(-2 * I * PI / m);
		for (int k = 0; k < N; k+=m)
		{
			for (int j = 0; j < m/2; ++j)
			{
				double complex u = X[k+j];
				double complex v = cpow(p, j) * X[k+j+m/2];
				X[k+j] = u + v;
				X[k+j+m/2] = u - v;
			}
		}
	}
}

int main(int argc, char const *argv[])
{
	/* code */
	double complex x[N];
	double complex X[N];

	double t[N];
	double A = 0.;
	double B = 10.;
	double dt = (B-A)/N;

	for (int i = 0; i < N; ++i)
	{

		t[i] = A + dt*i;
		x[i] = f(t[i]);
	}

	//fourier(x, X);
	//fft(x, X, N, 1);
	fft_iter(x,X);

	FILE *fichier = fopen("sol.txt", "w");
	fprintf(fichier,"");
	fclose(fichier);

	fichier = fopen("sol.txt", "a");
	for (int i = 0; i < N; ++i)
	{
		fprintf(fichier, "%f;%f;%f\n", t[i], (double)(x[i]), cabs(X[i]));
	}
	fclose(fichier);

	return 0;
}




