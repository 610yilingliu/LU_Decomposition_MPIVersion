#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

#define SIZE 3
// the random number will be generated in range [-RANGE, RANGE]
#define RANGE 3
// divide the generated number, converted it into double. For example, if range = 10 and SCALE = 10, random number will between [-1, 1]
#define SCALE 1
// task id of first task
#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2
#define a(x,y) a[x*SIZE+y] 

clock_t start_t,finish_t;
double total_t;

int nthreads = 1;
// MPI variables
int ntasks = 1;
MPI_Status status;
int taskid;

// arrays for calc
double matrix[SIZE][SIZE];
double L[SIZE][SIZE];
double U[SIZE][SIZE];
double vec[SIZE][1];
double answers[SIZE][1];
double Y[SIZE][1];

// for sub matrix
double *a;
double *mainline;

void exporting(double* arr_2d, int rownum, int colnum, char* fname) {
    // save in csv mode, split by ','
    // 2d array visiting solution: https://stackoverflow.com/questions/16724368/how-to-pass-a-2d-array-by-pointer-in-c
    FILE* fp = NULL;
    fp = fopen(fname, "w");
    double* p = (double*)arr_2d;
    for (int i = 0; i < rownum; i++) {
        for (int j = 0; j < colnum; j++) {
            double ele = p[i * colnum + j];
            fprintf(fp, "%f,", ele);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int find_maxrow(int col) {
    int mx = 0;
    int idx = 0;
    for (int i = col; i < SIZE; i++) {
        double cur = abs(matrix[i][col]);
        if (cur > mx) {
            mx = cur;
            idx = i;
        }
    }
    if (mx == 0) {
        finish_t = clock();
        total_t = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
        printf("Infinite number of answers\n");
        printf("Spent time:%f \n",total_t);
        exit(0);
    }
    return idx;
}

void swap(int a, int b) {
    if (a != b) {
        for (int i = 0; i < SIZE; i++) {
            double tmp = matrix[a][i];
            matrix[a][i] = matrix[b][i];
            matrix[b][i] = tmp;

        }
        double tmp = vec[a][0];
        vec[a][0] = vec[b][0];
        vec[b][0] = tmp;
    }
}

int re_arrange(){
    for(int i = 0; i < SIZE; i++){
        int a = find_maxrow(i);
        swap(a, i);
    }
}
// According to https://zhuanlan.zhihu.com/p/148647732, parrel is unnecessary because T(comp) ~ n and T(comm) ~ p, They are in the same scale.
void vec_generator(){
    for(int i = 0; i < SIZE; i ++){
        vec[i][0] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }
}

// generate random matrix in parallel
// Time Complexity: T(comp) = n/p * n ~ O(n^2), T(comm) = p
void matrix_generator(){
    int numworkers = ntasks - 1;
    int avgrow = SIZE/numworkers;
    int extra = SIZE%numworkers;
    int offset = 0;
    int rows;
    if(taskid == MASTER){
        for(int dest = 1; dest <= numworkers; dest++){
            rows = (dest <= extra) ? avgrow + 1 : avgrow;
            MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&matrix[offset][0], rows * SIZE, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            offset += rows;
        }
        for(int dest = 1; dest <= numworkers; dest++){
            MPI_Recv(&offset, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&matrix[offset][0], rows * SIZE, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &status);
        }
    }
    else{
        MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&matrix[offset][0], rows * SIZE, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        for(int i = offset; i < offset + rows; i++){
            for(int j = 0; j < SIZE; j++){
                // need to set random seed otherwise all rows will be the same
                srand(i * SIZE + j);
                matrix[i][j] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
            }
        }
        MPI_Send(&offset, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&matrix[offset][0], rows * SIZE, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD);
    }
    // sync matrix to all processes
    MPI_Bcast(&matrix[0], SIZE * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void get_lu(){
    for(int i = 0; i < SIZE; i++){
        for(int k = i; k < SIZE; k++){
            double sum = 0;
            double local = 0;
            for(int j = taskid; j < i; j+= ntasks){
                local += L[i][j] * U[j][k];
            }
            MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            U[i][k] = matrix[i][k] - sum;
        }
        for(int k = i; k < SIZE; k++){
            if(i == k){
                L[i][i] = 1;
            }
            else{
                double sum = 0;
                double local = 0;
                for(int j = taskid; j < i; j+= ntasks){
                    local += L[k][j] * U[j][i];
                }
                MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                L[k][i] = (matrix[k][i] - sum)/U[i][i];
            }
        }
    }
}

bool checker(){
    for(int i = 0; i < SIZE; i++){
        double sum = 0;
        double local = 0;
        for(int j = taskid; j < SIZE; j+= ntasks){
            local += answers[j][0] * matrix[i][j];
        }
        MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        // keep 3 digits to compare
        if((int)(sum * 1000 + 0.5)/1000.0 != (int)(vec[i][0] * 1000 + 0.5)/1000.0){
            printf("For line %d, %f != %f\n", i, sum, vec[i][0]);
            return false;
        }
    }
    return true;
}

void count(){
    for(int i = 0; i < SIZE; i++){
        double right = vec[i][0];
        double sum = 0;
        double local = 0;
        for(int j = taskid; j < i; j = ntasks){
            local += L[i][j] * Y[j][0];
        }
        MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        Y[i][0] = (right - sum)/L[i][i];
    }
    for(int i = SIZE - 1; i > -1; i--){
        double right = Y[i][0];
        double sum = 0;
        double local = 0;
        for(int j = SIZE - 1 - taskid; j > i; j-= ntasks){
            local += U[i][j] * answers[j][0];
        }
        MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        answers[i][0] = (right - sum)/U[i][i];
    }
}


int main(int argc, char* argv[]) {
    void *matrix_pointer = matrix;
    void *vector_pointer = vec;
    void *answer_pointer = answers;
    void *L_pointer = L;
    void *U_pointer = U;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
    if(ntasks == 1){
        printf("Two or more tasks are needed\n");
        exit(EXIT_FAILURE);
    }
    matrix_generator();
    if(taskid == MASTER){
        vec_generator();
        exporting(matrix_pointer, SIZE, SIZE, "matrix.csv");
        exporting(vector_pointer, SIZE, 1, "vector.csv");
        re_arrange();
    }
    MPI_Bcast(&matrix[0], SIZE * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&vec[0], SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    get_lu();
    if(taskid == MASTER){
        exporting(L_pointer, SIZE, SIZE, "L.csv");
        exporting(U_pointer, SIZE, SIZE, "U.csv");
    }
    count();
    if(taskid == MASTER){
        exporting(answer_pointer, SIZE, 1, "answer.csv");
    }
    bool res = checker();
    if(taskid == MASTER){
       if(res == true){
           printf("Checked, correct answer");
       }
       else{
           printf("Wrong answer");
       }
    }  
    MPI_Finalize();
    return 0;
}