#include "heads.h"

void line_generator(float *curline, int length){
    for (int i = 0; i < length; i++) {
            &curline[i] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
        }
    }

void matrix_generator(int ntasks){
    if(ntasks == 1){
        for(int i = 0; i < SIZE; i ++){
            line_generator(matrix[i], SIZE);
        }
    }
    else{
        // multiprocessing
        int taskid;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
        for(int i = taskid; i < SIZE; i+= ntasks){
            line_generator(matrix[i], SIZE);
        }
        MPI_Finalize();
    }
}

void vec_generator(){
    for(int i = 0; i < length; i ++){
        vec[i] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }
}