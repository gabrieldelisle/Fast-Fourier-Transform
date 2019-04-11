#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <mpi.h>

#define PI 3.141592653589793
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define f square


double complex square(double t) {
    if ((int) creal(t) % 2 == 0) {
        return 1.;
    } else {
        return -1.;
    }
}

double complex triangle(double t) {
    int t0 = (int) creal(t);
    if (t0 % 2 == 0) {
        return 2*(t-t0-0.5);
    } else {
        return -2*(t-t0-0.5);
    }
}

int bitflip(int nb, int bit_indice) {
    return nb ^ (1 << bit_indice);
}

int reverseBin(int i, int K) {
    int j = 0;

    for (int k = K - 1; k > -1; k--) {
        j += i / (1 << k) * (1 << (K - 1 - k));
        i %= 1 << k;
    }
    return j;
}


int main(int argc, char **argv) {

    int N;
    /* Find problem size N from command line */
    if (argc < 2) {
        fprintf(stdout, "No size N given\n");
        exit(1);
    }
    N = atoi(argv[1]);
    

    /* local variable */
    int p, P, odd_p, notify, global_index, i, Ip, L, R;
    int tag = 0;
    MPI_Status status;


    double complex *x;
    double complex *X;


/* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    if (N < P) {
        fprintf(stdout, "Too few discretization points...\n");
        exit(1);
    }

/* Compute local indices for data distribution */

    int logN = 0;
    int n = N;
    while (n >>= 1) logN++;

    L = N / P;
    R = N % P;
    Ip = (N + P - p - 1) / P;

    X = (double complex *) malloc(2 * Ip * sizeof(double complex));
    x = (double complex *) malloc(Ip * sizeof(double complex));


    double t[N];
    double A = 0.;
    double B = 10.;
    double dt = (B - A) / N;

    for (int i = 0; i < Ip; ++i) {
        global_index = L * p + MIN(p, R) + i;
        t[i] = A + dt*global_index;
        x[i] = f(t[i]);
        X[i] = f(A + dt * reverseBin(global_index,logN));
    }



    int k, j;


    for (i = 1; i < logN+1; ++i){
        int m = 1<<i;
        double complex omega = cexp(-2 * I * PI / m);

        if(m/2>=Ip){
            int other_indice = bitflip(L * p + MIN(p, R) + 0, i-1);
            int p_other = MAX((other_indice) / (L + 1), (other_indice - R) / L);
//            printf("Other processor = %d\r\n", p_other);
//            printf("s = %d \t p =%d \t indice = %d\r\n",i, p,other_indice);

            // might be redundant
            if (p_other != p) {
                if (p > p_other) {
                    MPI_Recv(&X[Ip], Ip, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD, &status);
                    MPI_Send(&X[0], Ip, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD);
                } else {
                    MPI_Send(&X[0], Ip, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD);
                    MPI_Recv(&X[Ip], Ip, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD, &status);
                }

                for (k = 0; k < Ip; k+=m)
                {
                    for (j = 0; j < m/2; ++j)
                    {
//                        printf("s = %d \t p =%d \t indice = %d\r\n",i, p,k+j+Ip);
                        if(p < p_other){
                            double complex u = X[k+j];
                            double complex v = cpow(omega, j) * X[k+j+Ip];
                            X[k+j] = u + v;
                        }else{
                            double complex v = cpow(omega, j) * X[k+j];
                            double complex u = X[k+j+Ip];
                            X[k+j] = u - v;
                        }

                    }
                }
            }
        }
        else{
            for (k = 0; k < Ip; k+=m)
            {
                for (j = 0; j < m/2; ++j)
                {
                    double complex u = X[k+j];
                    double complex v = cpow(omega, j) * X[k+j+m/2];
                    X[k+j] = u + v;
                    X[k+j+m/2] = u - v;
                }
            }
        }


    }

    if(p==0){
        FILE *fichier = fopen("sol.txt", "w");
        fprintf(fichier, "");
        fclose(fichier);
    }



    FILE *fp = fopen("sol.txt", "a");


    if(p!=0){
        MPI_Recv(&notify, 1, MPI_INT,p-1 , tag, MPI_COMM_WORLD, &status);
    }

    /*
     * Write results to the file
     */
    for (i = 0; i < Ip; ++i) {
        fprintf(fp, "%f;%f;%f\n", t[i], (double)(x[i]), cabs(X[i]));
    }

    fclose(fp);
    if(p!=P-1){
        MPI_Send(&notify, 1, MPI_INT, p+1, tag, MPI_COMM_WORLD);
    }

    printf("Process %d over\n", p);
//    /* That's it */
    MPI_Finalize();
    exit(0);
}





