#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <mpi.h>

#define PI 3.141592653589793


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
    int k;
    for (k = K - 1; k > -1; k--) {
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
    int p, P, p_other, global_index, i, j, k;
    int tag = 0;
    MPI_Status status;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    if (N < P*P) {
        fprintf(stdout, "Too few discretization points...\n");
        exit(1);
    }

    /* Compute local indices for data distribution */

    int logN = 0;
    {
        int n = N;
        while (n >>= 1) logN++;
    }
    

    int L = N / P;

    /* Problem definition */
    double complex *x;
    x = (double complex *) malloc(L * sizeof(double complex));

    double t[N];
    double A = 0.;
    double B = 10.;
    double dt = (B - A) / N;

    for (i = 0; i < L; ++i) {
        t[i] = A + dt * (L * p + i);
        x[i] = square(t[i]);
    }



    /* Binary reverse exchange */

    int *origin; //global index from the value in x
    origin = (int *) malloc(L * sizeof(int));
    int *origin2;
    origin2 = (int *) malloc(L * sizeof(int));

    double complex *X;
    X = (double complex *) malloc(2 * L * sizeof(double complex));
    
    int *p_indices; //count values which target is process p
    p_indices = (int *) malloc(P * sizeof(int));
    for (i = 0; i < P; i++) {
        p_indices[i] = 0; //initialize at 0
    }

    //move values and their indice so they are send to the right process
    for (i = 0; i < L; ++i) {
        global_index = L * p + i;
        p_other = reverseBin(global_index, logN)/L;
        X[p_other * (L/P) + p_indices[p_other]] = x[i];
        origin[p_other * (L/P) + p_indices[p_other]] = global_index;
        p_indices[p_other]++;
    }
    free(p_indices);



    MPI_Alltoall(&X[0], L/P, MPI_DOUBLE_COMPLEX, &X[L], L/P, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    MPI_Alltoall(&origin[0], L/P, MPI_INT, &origin2[0], L/P, MPI_INT, MPI_COMM_WORLD);
   
    //reexchange values using their global index
    for (i = 0; i < L; ++i) {
        X[i] = X[L + reverseBin(origin2[i], logN) % L];
    }
    free(origin);
    free(origin2);


    /* iterative fft algorithm */
    double complex u, v, omega;
    int m = 1;
    for (i = 1; i < logN+1; ++i){
        m <<= 1;
        omega = cexp(-2 * I * PI / m);

        // case 1: the values needed are on another processor
        if(m/2>=L){
            p_other = ( (L*p) ^ (m/2) ) / L; 

            if (p > p_other) {
                MPI_Recv(&X[L], L, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD, &status);
                MPI_Send(&X[0], L, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD);
            } else {
                MPI_Send(&X[0], L, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD);
                MPI_Recv(&X[L], L, MPI_DOUBLE_COMPLEX, p_other, tag, MPI_COMM_WORLD, &status);
            }

            for (k = 0; k < L; k+=m)
            {
                for (j = 0; j < L; ++j)
                {
                    if(p < p_other){
                        u = X[k+j];
                        v = cpow(omega, j+(L*p)%(m/2) ) * X[k+j+L];
                        X[k+j] = u + v;
                    }else{
                        v = cpow(omega, j+(L*p)%(m/2) ) * X[k+j];
                        u = X[k+j+L];
                        X[k+j] = u - v;
                    }

                }
            }
        }
        //case 2: the values needed are on the same processor
        else{
            for (k = 0; k < L; k+=m)
            {
                for (j = 0; j < m/2; ++j)
                {
                    u = X[k+j];
                    v = cpow(omega, j) * X[k+j+m/2];
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
    int notify;

    if(p!=0){
        MPI_Recv(&notify, 1, MPI_INT,p-1 , tag, MPI_COMM_WORLD, &status);
    }

    /*
     * Write results to the file
     */
    for (i = 0; i < L; ++i) {
        fprintf(fp, "%f;%f;%f\n", t[i], (double)(x[i]), cabs(X[i]));
    }

    fclose(fp);
    if(p!=P-1){
        MPI_Send(&notify, 1, MPI_INT, p+1, tag, MPI_COMM_WORLD);
    }

//    /* That's it */
    MPI_Finalize();
    exit(0);
}





