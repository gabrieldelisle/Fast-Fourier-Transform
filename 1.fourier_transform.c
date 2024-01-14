#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.141592653589793
#define f square

double complex square(double t) {
  if ((int)creal(t) % 2 == 0) {
    return 1.;
  } else {
    return -1.;
  }
}

double complex triangle(double t) {
  int t0 = (int)creal(t);
  if (t0 % 2 == 0) {
    return 2 * (t - t0 - 0.5);
  } else {
    return -2 * (t - t0 - 0.5);
  }
}

void fourier(double complex *x, double complex *X, int N) {
  int k, n;
  for (k = 0; k < N; ++k) {
    X[k] = 0;
    for (n = 0; n < N; ++n) {
      X[k] += x[n] * cexp(-2 * I * PI * n * k / N);
    }
  }
}

int main(int argc, char const *argv[]) {
  int N, i;
  /* Find problem size N from command line */
  if (argc < 2) {
    fprintf(stdout, "No size N given\n");
    exit(1);
  }
  N = atoi(argv[1]);

  double complex *x;
  double complex *X;
  X = (double complex *)malloc(2 * N * sizeof(double complex));
  x = (double complex *)malloc(N * sizeof(double complex));

  double t[N];
  double A = 0.;
  double B = 10.;
  double dt = (B - A) / N;

  // generates input vector
  for (i = 0; i < N; ++i) {
    t[i] = A + dt * i;
    x[i] = f(t[i]);
  }

  // fourier transform
  fourier(x, X, N);

  // save output to file
  FILE *fichier = fopen("sol.txt", "w");
  fprintf(fichier, "");
  fclose(fichier);

  fichier = fopen("sol.txt", "a");
  for (i = 0; i < N; ++i) {
    fprintf(fichier, "%f;%f;%f\n", t[i], (double)(x[i]), cabs(X[i]));
  }
  fclose(fichier);

  return 0;
}
