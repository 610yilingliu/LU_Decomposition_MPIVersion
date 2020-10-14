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

int mpierror;

// exit if mpi error
void mpi_error_check(int mpierror){
    if(mpierror != MPI_SUCCESS){
        printf("Error mpi failure\n");
        exit(EXIT_FAILURE);
    }
}


void exporting(double* arr_2d, int fsize, char* fname) {
    // mat_offset: offset in matrix
    // offset: offset in file
    MPI_File wh;
    int ele_per_procs = fsize/ntasks;
    while(ele_per_procs * ntasks < fsize){
        ele_per_procs++;
    }
    // last process treat the rest
    int ele_last_procs = fsize - (ele_per_procs * (ntasks - 1));

    MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&wh);
    if(taskid == MASTER){
        // the first process should write the size of the matrix
        // will be error if use SIZE instead of &SIZE
        mpierror = MPI_File_write_at(wh, 0, &SIZE, 1, MPI_INT, &status);
        mpi_error_check(mpierror);
    }
    if(taskid == ntasks - 1){
        int mat_offset = taskid * ele_per_procs; 
        MPI_Offset offset = mat_offset;
        // because we already have the dimension at the beginning of the file so the offset should be 1 * sizeof(int) + offset * sizeof(double)
        mpierror = MPI_File_write_at_all(wh, offset * sizeof(double) + 1 * sizeof(int), &arr_2d[mat_offset], ele_last_procs, MPI_DOUBLE, &status);
        mpi_error_check(mpierror);
    }
    else{
        int mat_offset = taskid * ele_per_procs;
        MPI_Offset offset = mat_offset;
        mpierror = MPI_File_write_at_all(wh, offset * sizeof(double) + 1 * sizeof(int), &arr_2d[mat_offset], ele_per_procs, MPI_DOUBLE, &status);
        mpi_error_check(mpierror);
    }
    MPI_File_close(&wh);
}

void printmat(double * arr_2d, int rownum ,int colnum){
    for(int i = 0; i < rownum; i++){
        for(int j = 0; j < colnum; j++){
            printf("%f ",arr_2d[i * colnum + j]);
        }
        printf("\n");
    }
}

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


void vec_generator(){
    for(int i = 0; i < SIZE; i ++){
        vec[i]= (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
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


// calculate answer
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

void freeall(){
    free(matrix);
    free(L);
    free(U);
    free(vec);
    free(answers);
    free(Y);
}

void read_mat(char* fname, double* m, int fsize){
    MPI_File fh;
    int ele_per_procs = fsize/ntasks;
    MPI_File_open(MPI_COMM_WORLD, fname,MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    while(ele_per_procs * ntasks < fsize){
        ele_per_procs++;
    }
    // last process handle the smallest number of items.
    int ele_last_procs = fsize - (ele_per_procs * (ntasks - 1));
    // last process
    if(taskid == ntasks - 1){
        int mat_offset = 1 + (taskid * ele_per_procs);
        MPI_Offset offset = mat_offset;
        MPI_File_sync(fh);
        mpierror = MPI_File_read_at_all(fh, offset * sizeof(double), &m[mat_offset], ele_last_procs, MPI_DOUBLE, &status);
        mpi_error_check(mpierror);
    }
    else{
        int mat_offset = 1 + (taskid * ele_per_procs);
        MPI_Offset offset = mat_offset;
        MPI_File_sync(fh);
        mpierror = MPI_File_read_at_all(fh, offset * sizeof(double), &m[mat_offset], ele_per_procs, MPI_DOUBLE, &status);
        mpi_error_check(mpierror);
    }

    MPI_File_close(&fh);
}

int check_file(char* fname){
    return 0;
}

int get_size(char* fname) {
    FILE *fp = NULL;
    fp = fopen(fname, "rb");
    int num;
    // Reads the dimensions
    fread(&num, sizeof(int), 1, fp);
    return num;
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
    if(ntasks == 1){
        printf("Two or more tasks are needed\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    start_t = MPI_Wtime();
        
    // input error checking section. If error, use MPI_Finalize() to terminate all process
    //read from file
    if(argv[1][0] == 'r' && argc == 2){
        if(check_file("matrix") == 0 || check_file("vector") == 0){
            if(taskid == MASTER) printf("Missing Matrix");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
         SIZE = get_size("matrix");
         if(taskid == MASTER) printf("Size: %d\n", SIZE);
        // if cannot find the re-arranged matrix and vector, we have to use original matrix to calculate
        if(check_file("rmatrix") == 0 || check_file("rvector") == 0){
            matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
            // L and U should be initialized as zero
            L = (double*)calloc(SIZE*SIZE, sizeof(double));
            U = (double*)calloc(SIZE*SIZE, sizeof(double));
            vec = (double*)malloc(SIZE * sizeof(double));
            Y = (double*)calloc(SIZE, sizeof(double));
            answers = (double*)calloc(SIZE, sizeof(double));
            if(taskid == MASTER) printf("Reading matrix and vector\n");
            read_mat("matrix", matrix, SIZE * SIZE);
            read_mat("vector", vec, SIZE);
            re_arrange();
            exporting(matrix, SIZE* SIZE, "rmatrix");
            exporting(vec, SIZE, "rvector");
            get_lu();
            exporting(L, SIZE * SIZE, "L");
            exporting(U, SIZE * SIZE, "U");
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
            read_mat("matrix", matrix, SIZE * SIZE);
            read_mat("vector", vec, SIZE);
            read_mat("L", L, SIZE * SIZE);
            read_mat("U", U, SIZE * SIZE);
        }
    }
    // generate matrix and calculate
    else if(argv[1][0] == 'g'){
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