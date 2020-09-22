#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>
#include <string.h>

// The size of the matrix
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

double matrix[SIZE][SIZE];
double vec[SIZE][1];
double answers[SIZE][1];
// generate matrix with random number. input: number of thread
extern void matrix_generator(int);
extern void vector_generator();
