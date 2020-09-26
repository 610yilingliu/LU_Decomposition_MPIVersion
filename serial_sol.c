#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#define SIZE 100
// the random number will be generated in range [-RANGE, RANGE]
#define RANGE 3
// divide the generated number, converted it into double. For example, if range = 10 and SCALE = 10, random number will between [-1, 1]
#define SCALE 1


clock_t start_t,finish_t;
double total_t;
double matrix[SIZE][SIZE];
double L[SIZE][SIZE];
double U[SIZE][SIZE];
double vec[SIZE][1];
double answers[SIZE][1];
double Y[SIZE][1];

void matrix_generator(){
    for(int i = 0; i < SIZE; i ++){
        for(int j = 0; j < SIZE; j ++){
            matrix[i][j] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
        }
    }
}

void vec_generator(){
    for(int i = 0; i < SIZE; i ++){
        vec[i][0] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }
}

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

void count(){

    for(int i = 0; i < SIZE; i++){
        double right = vec[i][0];
        for(int j = 0; j < i; j++){
            right -= L[i][j] * Y[j][0];
        }
        Y[i][0] = right/L[i][i];
    }
    for(int i = SIZE - 1; i > -1; i--){
        double right = Y[i][0];
        for(int j = SIZE - 1; j > i; j--){
            right -= U[i][j] * answers[j][0];
        }
        answers[i][0] = right/U[i][i];
    }
}

void get_lu(){
    for(int i = 0; i < SIZE; i++){
        for(int k = i; k < SIZE; k++){
            double sum = 0;
            for(int j = 0; j < i; j++){
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = matrix[i][k] - sum;
        }
        for(int k = i; k < SIZE; k++){
            if(i == k){
                L[i][i] = 1;
            }
            else{
                double sum = 0;
                for(int j = 0; j < i; j++){
                    sum += L[k][j] * U[j][i];
                }
                L[k][i] = (matrix[k][i] - sum)/U[i][i];
            }
        }
    }
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

bool checker(){
    for(int i = 0; i < SIZE; i++){
        double sum = 0;
        for(int j = 0; j < SIZE; j++){
            sum += answers[j][0] * matrix[i][j];
        }
        // keep 3 digits to compare
        if((int)(sum * 1000 + 0.5)/1000.0 != (int)(vec[i][0] * 1000 + 0.5)/1000.0){
            printf("For line %d, %f != %f\n", i, sum, vec[i][0]);
            return false;
        }
    }
    return true;
}
int main(int argc, char* argv[]){
    printf("Start Processing......\n");
    start_t = clock();
    void *matrix_pointer = matrix;
    void *vector_pointer = vec;
    void *answer_pointer = answers;
    void *L_pointer = L;
    void *U_pointer = U;
    vec_generator();
    matrix_generator(matrix);
    exporting(matrix_pointer, SIZE, SIZE, "matrix.csv");
    exporting(vector_pointer, SIZE, 1, "vector.csv");
    re_arrange();
    // exporting(matrix_pointer, SIZE, SIZE, "arranged_mat.csv");
    // exporting(vector_pointer, SIZE, 1, "arranged_vec.csv");
    get_lu();
    // exporting(L_pointer, SIZE, SIZE, "L.csv");
    // exporting(U_pointer, SIZE, SIZE, "U.csv");
    count();
    exporting(answer_pointer, SIZE, 1, "answer.csv");
    finish_t = clock();
    total_t = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
    printf("Spent time:%f \n",total_t);
    printf("Start Checking...\n");
    bool result = checker();
    if(result == true){
        printf("Answer is correct\n");
    }
    else{
        printf("Not a correct answer\n");
    }
    return 0;
}