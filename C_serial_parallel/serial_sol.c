#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

// macro to locate the element in matrix
#define matrix(i,j) matrix[i * SIZE + j]
#define U(i,j) U[SIZE * i + j]
#define L(i,j) L[SIZE * i +j]

// size of the matrix, length of the vector
int SIZE;
// the random number will be generated in range [-RANGE, RANGE]
int RANGE;
// divide the generated number, converted it into double. For example, if range = 10 and SCALE = 10, random number will between [-1, 1]
int SCALE;

clock_t start_t,finish_t;

double total_t;
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
int export_ans = 1;


void matrix_generator(){
    for(int i = 0; i < SIZE; i ++){
        for(int j = 0; j < SIZE; j ++){
            matrix(i,j)=(double)(rand() % (RANGE * 2) - RANGE) / SCALE;
        }
    }
}

void vec_generator(){
    for(int i = 0; i < SIZE; i ++){
        vec[i] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }
}

void exporting(double* arr_2d, int rownum, int colnum, char* fname) {
    // save in csv mode, split by ','
    // 2d array visiting solution: https://stackoverflow.com/questions/16724368/how-to-pass-a-2d-array-by-pointer-in-c
    FILE* fp = NULL;
    fp = fopen(fname, "w");
    fprintf(fp, "%d\n", colnum);
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

void printmat(double * arr_2d, int rownum ,int colnum){
    for(int i = 0; i < rownum; i++){
        for(int j = 0; j < colnum; j++){
            printf("%f ",arr_2d[i * colnum + j]);
        }
        printf("\n");
    }
}

// Time complexity(O(n^2))
void count(){
    for(int i = 0; i < SIZE; i++){
        double right = vec[i];
        for(int j = 0; j < i; j++){
            right -= Y[j] * L(i,j);
        }
        // cannot divide by zero, so should do the error checking
        if(L(i, i) == 0){
            // the matrix generated randomly in serial and parallel code are not the same, so we cannot compare the time if they terminate without finish the counting process.
            // so we will not print the time if the code terminate with error.
            printf("Divide by zero error!");
            exit(0);
        }
        Y[i] = right/L(i,i);
    }
    for(int i = SIZE - 1; i > -1; i--){
        double right = Y[i];
        for(int j = SIZE - 1; j > i; j--){
            right -= U(i,j) * answers[j];
        }
        // divide by zero error checking
        if(U(i, i) == 0){
            printf("Divide by zero error!");
            exit(0);
        }
        answers[i] = right/U(i, i);
    }
}

// Time complexity(O(n^3))
void get_lu(){
    int print;
    int print_every;
    int last;
    int cnt = 0;
    if(SIZE > 2000){
        print_every = SIZE/100;
        print = 1;
        last = 0;
    }
    for(int i = 0; i < SIZE; i++){
        if(print == 1){
            if(i >= last){
                last += print_every;
                printf("%d percent of LU finished\r", cnt + 1);
                cnt += 1;
            }
        }
        for(int k = i; k < SIZE; k++){
            double sum = 0;
            for(int j = 0; j < i; j++){
                sum += L(i, j) * U(j, k);
            }
            U(i, k) = matrix(i, k) - sum;
            // printf("U[%d, %d] = %f   ", i, k, U(i, k));
        }
        for(int k = i; k < SIZE; k++){
            if(i == k){
                L(i, i) = 1;
            }
            else{
                double sum = 0;
                for(int j = 0; j < i; j++){
                    sum += L(k, j) * U(j, i);
                }
                // divide by zero error checking
                if(U(i, i) == 0){
                    printf("\nDivide by zero error!");
                    exit(0);
                }
                L(k, i) = (matrix(k, i) - sum)/U(i, i);
            }
        }
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
        printf("\nInvalid Matrix\n");
        exit(0);
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
        swap(a, i);
    }
}

// Time complexity(O(n^2))
bool checker(){
    for(int i = 0; i < SIZE; i++){
        double sum = 0;
        for(int j = 0; j < SIZE; j++){
            sum += answers[j] * matrix(i, j);
        }
        // keep 3 digits to compare
        if((int)(sum * 1000 + 0.5)/1000.0 != (int)(vec[i] * 1000 + 0.5)/1000.0){
            printf("For line %d, %f != %f\n", i, sum, vec[i]);
            return false;
        }
    }
    return true;
}

void freeall(){
    free(matrix);
    free(L);
    free(U);
    free(vec);
    free(answers);
    free(Y);
}

int get_size(char *fname){
    FILE* fp = NULL;
    fp = fopen(fname, "r");
    int colnum = 0;
    char cur[1];
    // stop when finish reading the first line
    while(!feof(fp)){
        cur[0] = fgetc(fp);
        if(cur[0] == '\0' || cur[0] == '\r'||cur[0] == '\n') break;
        colnum = colnum * 10 + atoi(cur);
    }
    fclose(fp);
    return colnum;
}

void read_mat(char *fname, double *m){
    FILE *fp = NULL;
    fp = fopen(fname, "r");
    int linenum = 0;
    int curpos = 0;
    int matpos = 0;
    char buf[30] = "";
    while(!feof(fp)){
        char curchar = fgetc(fp);
        if(curchar == '\n' || curchar == '\r') linenum++;
        // clean up buffer and reset curpos to zero if reach to ','
        else if(curchar == ','){
            m[matpos] = atof(buf);
            char buf[30] = "";
            curpos = 0;
            matpos += 1;
        }
        else{
            if(linenum > 0){
                if(curchar != '\n' && curchar != '\0' && curchar != ' ' && curchar!= '\r'){
                    buf[curpos] = curchar;
                    curpos += 1;
                }
            }
        }
    }
    fclose(fp);
}

// check if target file exists
int check_file(char *fname){
    FILE *fp = NULL;
    fp = fopen(fname, "r");
    if(fp == NULL) return 0;
    fclose(fp);
    return 1;
}

int main(int argc, char* argv[]){
    //read from file
    if(argv[1][0] == 'r' && argc == 2){
        if(check_file("matrix.csv") == 0 || check_file("vector.csv") == 0){
            printf("Missing Matrix");
            exit(0);
        }
         SIZE = get_size("matrix.csv");
         printf("Size: %d\n", SIZE);
        // if cannot find the re-arranged matrix and vector, we have to use original matrix to calculate
        if(check_file("rmatrix.csv") == 0 || check_file("rvector.csv") == 0){
            matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
            // L and U should be initialized as zero
            L = (double*)calloc(SIZE*SIZE, sizeof(double));
            U = (double*)calloc(SIZE*SIZE, sizeof(double));
            vec = (double*)malloc(SIZE * sizeof(double));
            Y = (double*)calloc(SIZE, sizeof(double));
            answers = (double*)calloc(SIZE, sizeof(double));
            printf("Reading matrix and vector\n");
            read_mat("matrix.csv", matrix);
            read_mat("vector.csv", vec);
            re_arrange();
            exporting(matrix, SIZE, SIZE, "rmatrix.csv");
            exporting(vec, 1, SIZE, "rvector.csv");
            get_lu();
            exporting(L, SIZE, SIZE, "L.csv");
            exporting(U, SIZE, SIZE, "U.csv");
        }
        // else we can use L and U to calculate
        else{
            printf("Reading L and U\n");
            matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
            L = (double*)malloc(SIZE * SIZE * sizeof(double));
            U = (double*)malloc(SIZE * SIZE * sizeof(double));
            vec = (double*)malloc(SIZE * sizeof(double));
            Y = (double*)calloc(SIZE, sizeof(double));
            answers = (double*)calloc(SIZE, sizeof(double));
            read_mat("matrix.csv", matrix);
            read_mat("vector.csv", vec);
            read_mat("L.csv", L);
            printmat(L, SIZE, SIZE);
            read_mat("U.csv", U);
        }
    }
    // generate matrix and calculate
    else if(argv[1][0] == 'g'){
        if(argc != 6){
            printf("You should allocate three variables: size, range, scale, export answer or not(0/1)\n");
            printf("Current number of args: %d\n", argc);
            exit(0);
        }
        if(argv[2][0] == 0){
            printf("Matrix size should greater than 0");
            exit(0);
        }
        if(argv[3][0] == 0){
            printf("Range should greater than 0");
            exit(0);
        }
        if(argv[4][0] == 0){
            printf("Scale should greater than 0");
            exit(0);
        }
        for(int i = 2; i < argc; i ++){
            for(int j = 0; j < strlen(argv[i]); j++){
                if(argv[i][j] < '0' || argv[i][j] > '9'){
                    printf("Invalid argument at arg[%d][%d]", i, j);
                    exit(0);
                }
            }
        }
        if(strlen(argv[5]) > 1||(argv[4][0]!= '0' && argv[4][0]!= '1')){
            printf("The 4th argument, export or not should be 0(not export) or 1(export)");
            exit(0);
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

        printf("Start Processing......\n");
        start_t = clock();
        vec_generator();
        matrix_generator();
        // //test input
        // double m[] = {-1.0000, 0.0000, 1.0000, 2.0000,2.0000,1.0000,0.0000,-1.0000,-3.0000};
        // double v[] = {2.0000, 2.0000, 1.0000};
        // matrix = m;
        // vec = v;
        printf("Random vector and Matrix Generated\n");
        if(export_ans == 1){
            exporting(matrix, SIZE, SIZE, "matrix.csv");
            exporting(vec, 1, SIZE, "vector.csv");
        }
        printf("Rearranging Matrix...\n");
        re_arrange();
        if(export_ans == 1){
            exporting(matrix, SIZE, SIZE, "rmatrix.csv");
            exporting(vec, 1, SIZE, "rvector.csv");
        }
        printf("Calculating LU\n");
        get_lu();
        if(export_ans == 1){
            exporting(L, SIZE, SIZE, "L.csv");
            exporting(U, SIZE, SIZE, "U.csv");
        }
    }
    else{
        printf("Invalid Argument\n");
        printf("%s ", argv[0]);
        printf("%s ", argv[1]);
        printf("%s ", argv[2]);
        printf("%s ", argv[3]);
        printf("%s ", argv[4]);
        printf("%s ", argv[5]);
        exit(0);
    }
    printf("\nCalculating Result...\n");
    count();
    printf("Start Checking...\n");
    bool result = checker();
    if(result == true){
        printf("\nAnswer is correct\n");
    }
    else{
        printf("\nNot a correct answer\n");
    }
    if(export_ans == 1){
        exporting(answers, 1, SIZE, "answers.csv");
    }
    finish_t = clock();
    total_t = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
    printf("Spent time: %f seconds\n",total_t);
    freeall();
    return 0;
}