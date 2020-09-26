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
// divide the generated number, converted it into double. For example, if range = 10 and SCALE = 10, random number will between [-, 1]
#define SCALE 1
// task id of first task
#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

clock_t start_t,finish_t;
double total_t;

int nthreads = 1;
// MPI variables
int ntasks = 1;
MPI_Status status;
int taskid;

double matrix[SIZE][SIZE];
double L[SIZE][SIZE];
double U[SIZE][SIZE];
double vec[SIZE][1];
double answers[SIZE][1];
double Y[SIZE][1];

// generate matrix with random number. input: 2d matrix pointer, number of thread, is empty or not
extern void matrix_generator(double *);
extern void vec_generator();
// input: number of threads
extern void re_arrange(int);