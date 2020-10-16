#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

#define matrix(i,j) matrix[i * SIZE + j]
#define U(i,j) U[SIZE * i + j]
#define L(i,j) L[SIZE * i +j]

// task id of first task
#define MASTER 0
// We will use MPI_Wtime library instead of time library to avoid error
// Usage of MPI_Wtime https://www.mcs.anl.gov/research/projects/mpi/tutorial/gropp/node139.html
double start_t, finish_t;
double total_t;
// MPI variables
int ntasks = 1;
MPI_Status status;
int taskid;
// size of the matrix, length of the vector
int SIZE;
// the random number will be generated in range [-RANGE, RANGE]
int RANGE;
// divide the generated number, converted it into double. For example, if range = 10 and SCALE = 10, random number will between [-1, 1]
int SCALE;
// SIZE * SIZE
double *matrix;
// SIZE * SIZE
double *L;
// SIZE * SIZE
double *U;
// 1 * SIZE
double *vec;
// 1 * SIZE
double *answers;
// 1 * SIZE
double *Y;
// export: 1, not export: 0
int export_ans;
// error checking, from 
int mpierror;

/*
Debug modules
*/

// Check which process got error. Function from https://github.com/BologneseBandit/HPCassignment2
void mpi_error_check(int mpierror, int taskid){
    if(mpierror != MPI_SUCCESS){
        printf("MPI error at process %d\n", taskid);
        exit(EXIT_FAILURE);
    }
}
// print a matrix on console.
void printmat(double * arr_2d, int rownum ,int colnum){
    for(int i = 0; i < rownum; i++){
        for(int j = 0; j < colnum; j++){
            printf("%f ",arr_2d[i * colnum + j]);
        }
        printf("\n");
    }
}

/*
    MPI I/O
*/

// Export to binary file. 
// Text file is not good for HPC because we cannot locate the index of where a number start. 
// For example, we cannot use multi-processing code to read "1 11111 123" with 3 processes, but in binary file
// it stores like (int)1 (int)11111 (int)123 so each process could read one number.
void exporting(double* arr_2d, int fsize, char* fname) {
    MPI_File wh;
    int ele_per_procs = fsize/ntasks;
    while(ele_per_procs * ntasks < fsize){
        ele_per_procs++;
    }
    // last process treat the rest
    int ele_last_procs = fsize - (ele_per_procs * (ntasks - 1));
    int offset = taskid * ele_per_procs; 

    MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&wh);
    if(taskid == MASTER){
        // the first process should write the size of the matrix
        // will be error if use SIZE instead of &SIZE
        mpierror = MPI_File_write_at(wh, 0, &SIZE, 1, MPI_INT, &status);
        mpi_error_check(mpierror, taskid);
    }
    if(taskid == ntasks - 1){
        // because we already have the dimension at the beginning of the file so the offset should be 1 * sizeof(int) + offset * sizeof(double)
        mpierror = MPI_File_write_at(wh, offset * sizeof(double) + 1 * sizeof(int), &arr_2d[offset], ele_last_procs, MPI_DOUBLE, &status);
    }
    else{
        mpierror = MPI_File_write_at(wh, offset * sizeof(double) + 1 * sizeof(int), &arr_2d[offset], ele_per_procs, MPI_DOUBLE, &status);
    }
    mpi_error_check(mpierror, taskid);
    MPI_File_close(&wh);
}

// another way to do I/O is to use collective method at_all instead of at
// the main reason to get memory error in this section are:
// 1. forgot to use fclose() if the file is already opened by fopen()
// 2. use read_at_all with if..else. All process need to use same variable name for variables put in MPI_File_read_at_all
void read_mat(char* fname, double* m, int fsize){
    MPI_File fh;
    int ele_per_procs = fsize/ntasks;
    MPI_File_open(MPI_COMM_WORLD, fname,MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    while(ele_per_procs * ntasks < fsize){
        ele_per_procs++;
    }
    // last process handle the smallest number of items.
    int ele_last_procs = fsize - (ele_per_procs * (ntasks - 1));
    int offset = taskid * ele_per_procs;
    int ele_in_procs = ele_per_procs;
    // last process
    if(taskid == ntasks - 1) ele_in_procs = ele_last_procs;
    mpierror = MPI_File_read_at_all(fh, offset * sizeof(double) + 1 * sizeof(int), &m[offset], ele_in_procs, MPI_DOUBLE, &status);
    mpi_error_check(mpierror, taskid);
    MPI_File_close(&fh);
}

/*
Find prime element

These three functions designed to find the prime element and swap the current row with the prime row in both matrix A and vector b. 
If we need to use it in the engineering field, this should also generate the transformational matrix P for PAx = Pb, but I will not do that
to make the code simple.
*/
int find_maxrow(int col) {
    int mx = 0;
    int idx = 0;
    for (int i = col; i < SIZE; i++) {
        double cur = abs(matrix(i, col));
        if (cur > mx) {
            mx = cur;
            idx = i;
        }
    }
    if (mx == 0) {
        printf("Infinite number of answers\n");
        return -1;
    }
    return idx;
}

void swap(int a, int b) {
    if (a != b) {
        for (int i = 0; i < SIZE; i++) {
            double tmp = matrix(a, i);
            matrix(a, i) = matrix(b, i);
            matrix(b, i) = tmp;

        }
        double tmp = vec[a];
        vec[a] = vec[b];
        vec[b] = tmp;
    }
}

int re_arrange(){
    for(int i = 0; i < SIZE; i++){
        int a = find_maxrow(i);
        if(a == -1) return -1;
        swap(a, i);
    }
    return 0;
}

/*
    Generate random matrix and vector
*/

// Time Complexity: O(n).Parallel is not needed here
void vec_generator(){
    for(int i = 0; i < SIZE; i ++){
        vec[i]= (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }
}

// generate random matrix in parallel
// Time Complexity: T(comp) = n/p * n ~ O(n^2), T(comm) = p
// Code reference: lecture slides. Use process 0(master process) to allocate tasks to other process and get ther return.
// Finally we should broadcast the matrix in process 0 to all processes because before broadcasting, other process just have the
// part of matrix that they generated themselves.
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
            MPI_Send(&matrix(offset, 0), rows * SIZE, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            offset += rows;
        }
        for(int dest = 1; dest <= numworkers; dest++){
            MPI_Recv(&offset, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&matrix(offset, 0), rows * SIZE, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &status);
        }
    }
    else{
        MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&matrix(offset, 0), rows * SIZE, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        for(int i = offset; i < offset + rows; i++){
            srand(i);
            for(int j = 0; j < SIZE; j++){
                // need to set random seed otherwise all rows will be the same
                matrix(i, j) = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
            }
        }
        MPI_Send(&offset, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&matrix(offset, 0), rows * SIZE, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD);
    }
    // sync matrix to all processes
    MPI_Bcast(&matrix[0], SIZE * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/*
Calculating L and U
*/
// O(n^3) in this step so we have to set a progress presentage here to make sure it is working
void get_lu(){
    int print;
    int print_every;
    int last;
    int cnt = 0;
    if(SIZE > 1000){
        print = 1;
        last = 0;
        print_every = SIZE/100;
    }
    for(int i = 0; i < SIZE; i++){
        if(taskid == MASTER){
            if(print == 1){
                if(i >= last){
                    last += print_every;
                    printf("%d percent of LU finished\r", cnt + 1);
                    fflush(stdout);
                    cnt += 1;
                }
            }
        }
        for(int k = i; k < SIZE; k++){
            double sum = 0;
            double local = 0;
            for(int j = taskid; j < i; j+= ntasks){
                local += L(i, j) * U(j, k);
            }
            MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

            U(i, k) = matrix(i, k) - sum;
        }
        for(int k = i; k < SIZE; k++){
            if(i == k){
                L(i, i) = 1;
            }
            else{
                double sum = 0;
                double local = 0;
                for(int j = taskid; j < i; j+= ntasks){
                    local += L(k, j) * U(j, i);
                }
                MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

                if(U(i, i) == 0){
                    printf("Divide by zero error");
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
                L(k, i) = (matrix(k, i) - sum)/U(i,i);
            }
        }
    }
}

/*
Calculate answer
*/
void count(){
    for(int i = 0; i < SIZE; i++){
        double right = vec[i];
        double sum = 0;
        double local = 0;
        for(int j = taskid; j < i; j += ntasks){
            local += L(i, j) * Y[j];
        }
        MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(L(i, i) == 0){
            printf("Divide by zero error");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        Y[i] = (right - sum)/L(i, i);
    }
    
    for(int i = SIZE - 1; i > -1; i--){
        double right = Y[i];
        double sum = 0;
        double local = 0;
        for(int j = SIZE - 1 - taskid; j > i; j-= ntasks){
            local += U(i, j) * answers[j];
        }
        MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(U(i,i) == 0){
            printf("Divide by zero error");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        answers[i] = (right - sum)/U(i, i);
    }
}

/*
Check Answer
*/
bool checker(){
    for(int i = 0; i < SIZE; i++){
        double sum = 0;
        double local = 0;
        for(int j = taskid; j < SIZE; j+= ntasks){
            local += answers[j] * matrix(i, j);
        }
        MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        // keep 3 digits to compare
        if((int)(sum * 1000 + 0.5)/1000.0 != (int)(vec[i] * 1000 + 0.5)/1000.0){
            printf("For line %d, %f != %f\n", i, sum, vec[i]);
            return false;
        }
    }
    return true;
}

/*
Free the pointer
*/
void freeall(){
    free(matrix);
    free(L);
    free(U);
    free(vec);
    free(answers);
    free(Y);
}


/*
Tool functions. Check if the file is exist or not
*/
int check_file(char* fname){
    FILE *fp = NULL;
    fp = fopen(fname, "r");
    if(fp == NULL) return 0;
    fclose(fp);
    return 1;
}

/*
Get the size of a matrix in a binary file
*/
int get_size(char* fname) {
    FILE *fp = NULL;
    fp = fopen(fname, "rb");
    int num;
    // Reads the dimensions
    fread(&num, sizeof(int), 1, fp);
    fclose(fp);
    return num;
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
    MPI_Barrier(MPI_COMM_WORLD);

    start_t = MPI_Wtime();
    // input error checking section. If error, use MPI_Finalize() to terminate all process
    //read from file
    if(argv[1][0] == 'r' && argc == 2){
        int mat_complete = 1;
        if(taskid == MASTER){
            if(check_file("matrix") == 0 || check_file("vector") == 0) {
                printf("Missing Matrix");
                mat_complete = 0;
        }
        MPI_Bcast(&mat_complete, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(mat_complete == 0){
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }
        // process 0 get size and broadcast to all processes
        if(taskid == MASTER) SIZE = get_size("matrix");
        MPI_Bcast(&SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int rearranged = 1;
        if(taskid == MASTER) printf("Size: %d\n", SIZE);
        // if cannot find the re-arranged matrix and vector, we have to use original matrix to calculate
        if(taskid == MASTER){
            if(check_file("rmatrix") == 0 || check_file("rvector") == 0) rearranged = 0;
            MPI_Bcast(&rearranged, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        if(rearranged == 0){
            matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
            // L and U should be initialized as zero
            L = (double*)calloc(SIZE*SIZE, sizeof(double));
            U = (double*)calloc(SIZE*SIZE, sizeof(double));
            vec = (double*)malloc(SIZE * sizeof(double));
            Y = (double*)calloc(SIZE, sizeof(double));
            answers = (double*)calloc(SIZE, sizeof(double));
            if(taskid == MASTER) printf("Reading matrix and vector\n");
            read_mat("matrix", matrix, SIZE * SIZE);
            if(taskid == MASTER && SIZE < 5){
                printf("Matrix\n");
                printmat(matrix, SIZE, SIZE);
            }
            read_mat("vector", vec, SIZE);
            if(taskid == MASTER && SIZE < 5){
                printf("Vector\n");
                printmat(vec, 1, SIZE);
            }
            re_arrange();
            exporting(matrix, SIZE* SIZE, "rmatrix");
            exporting(vec, SIZE, "rvector");
            get_lu();
            exporting(L, SIZE * SIZE, "L");
            if(taskid == MASTER && SIZE < 5){
                printf("L\n");
                printmat(L, SIZE, SIZE);
            }
            exporting(U, SIZE * SIZE, "U");
            if(taskid == MASTER && SIZE < 5){
                printf("U\n");
                printmat(U, SIZE, SIZE);
            }
        }
        // else we can use L and U to calculate
        else{
            if(taskid == MASTER) printf("Reading L and U\n");
            matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
            L = (double*)malloc(SIZE * SIZE * sizeof(double));
            U = (double*)malloc(SIZE * SIZE * sizeof(double));
            vec = (double*)malloc(SIZE * sizeof(double));
            Y = (double*)calloc(SIZE, sizeof(double));
            answers = (double*)calloc(SIZE, sizeof(double));
            // we need to read the re-arranged matrix and vector instead
            read_mat("rmatrix", matrix, SIZE * SIZE);
            if(taskid == MASTER && SIZE < 5){
                printf("Matrix\n");
                printmat(matrix, SIZE, SIZE);
            }
            read_mat("rvector", vec, SIZE);
            if(taskid == MASTER && SIZE < 5){
                printf("Vector\n");
                printmat(vec, 1, SIZE);
            }
            read_mat("L", L, SIZE * SIZE);
            if(taskid == MASTER && SIZE < 5){
                printf("L\n");
                printmat(L, SIZE, SIZE);
            }
            read_mat("U", U, SIZE * SIZE);
            if(taskid == MASTER && SIZE < 5){
                printf("U\n");
                printmat(U, SIZE, SIZE);
            }
        }
    }
    // generate matrix and calculate
    else if(argv[1][0] == 'g'){
        if(ntasks == 1){
            printf("Two or more tasks are needed\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        if(argc != 6){
            if(taskid == MASTER) printf("You should allocate three variables: size, range, scale, export answer or not(0/1)\n");
            if(taskid == MASTER) printf("Current number of args: %d\n", argc);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        if(argv[2][0] == 0){
            if(taskid == MASTER) printf("Matrix size should greater than 0");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        if(argv[3][0] == 0){
            if(taskid == MASTER) printf("Range should greater than 0");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        if(argv[4][0] == 0){
            if(taskid == MASTER) printf("Scale should greater than 0");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        for(int i = 2; i < argc; i ++){
            for(int j = 0; j < strlen(argv[i]); j++){
                if(argv[i][j] < '0' || argv[i][j] > '9'){
                    if(taskid == MASTER) printf("Invalid argument at arg[%d][%d]", i, j);
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
            }
        }
        if(strlen(argv[5]) > 1||(argv[5][0]!= '0' && argv[5][0]!= '1')){
            if(taskid == MASTER) printf("The 4th argument, export or not should be 0(not export) or 1(export)");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        export_ans = atoi(argv[5]);
        // assign input to variables
        SIZE = atoi(argv[2]);
        RANGE = atoi(argv[3]);
        SCALE = atoi(argv[4]);
        // allocate memory space
        matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
        // L and U should be initialized as zero
        L = (double*)calloc(SIZE*SIZE, sizeof(double));
        U = (double*)calloc(SIZE*SIZE, sizeof(double));
        vec = (double*)malloc(SIZE * sizeof(double));
        Y = (double*)calloc(SIZE, sizeof(double));
        answers = (double*)calloc(SIZE, sizeof(double));

        if(taskid == MASTER) printf("Start Processing......\n");
        vec_generator();
        matrix_generator();
        if(taskid == MASTER && SIZE < 5){
            printf("Matrix\n");
            printmat(matrix, SIZE, SIZE);
        }
        if(taskid == MASTER && SIZE < 5){
            printf("Vector\n");
            printmat(vec, 1, SIZE);
        }
        if(taskid == MASTER) printf("Random vector and Matrix Generated\n");
        if(export_ans == 1){
            exporting(matrix, SIZE * SIZE, "matrix");
            exporting(vec, SIZE, "vector");
        }
        if(taskid == MASTER) printf("Rearranging Matrix...\n");
        re_arrange();
        if(export_ans == 1){
            exporting(matrix, SIZE * SIZE, "rmatrix");
            exporting(vec, SIZE, "rvector");
        }
        if(taskid == MASTER) printf("Calculating LU\n");
        get_lu();
        if(taskid == MASTER && SIZE < 5){
            printf("L\n");
            printmat(L, SIZE, SIZE);
        }
        if(taskid == MASTER && SIZE < 5){
            printf("U\n");
            printmat(U, SIZE, SIZE);
        }
        if(export_ans == 1){
            exporting(L, SIZE * SIZE, "L");
            exporting(U, SIZE * SIZE, "U");
        }

    }
    else{
        if(taskid == MASTER) printf("Invalid Argument\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    if(taskid == MASTER) printf("\nCalculating Result...\n");
    count();
    MPI_Bcast(&answers[0], SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Y[0], SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(export_ans == 1){
        exporting(answers, SIZE, "answers");
    }
    if(taskid == MASTER && SIZE < 5){
        printf("Y\n");
        printmat(Y, 1, SIZE);
    }
    if(taskid == MASTER && SIZE < 5){
        printf("Answers\n");
        printmat(answers, 1, SIZE);
    }
    bool res = checker();
    MPI_Barrier(MPI_COMM_WORLD);
    finish_t = MPI_Wtime();
    if(taskid == MASTER){
       if(res == true){
           printf("\nChecked, correct answer");
       }
       else{
           printf("\nWrong answer");
       }
       total_t = finish_t - start_t;
       printf("\nSpent time: %f seconds\n",total_t);
    }  
    freeall();
    MPI_Finalize();
    return 0;
}