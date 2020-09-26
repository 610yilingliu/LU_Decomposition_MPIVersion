#include "heads.h"

void line_generator(double* m, int curline){
    for (int i = 0; i < SIZE; i++) {
        for(int j = 0; j < SIZE; j ++)
            &m[curline + i][j] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
        }
    }

void matrix_generator(double *m){
    if(ntasks == 1){
        for(int i = 0; i < SIZE; i ++){
            line_generator(m, i);
        }
    }
    else{
        int numworkers = ntasks - 1;
        int avgrow = SIZE/numworkers;
        int extra = SIZE%numworkers;
        int offset = 0;
        int rows;
        if(taskid == MASTER){
            for(int dest = 1; dest < numworkers; dest++){
                rows = (dest <= extra) ? avgrow + 1 : avgrow;
                MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
                MPI_Send(&rows, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
                MPI_Send(&m[offset][0], rows * SIZE, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            }
            for(int dest = 1; dest < numworkers; dest++){
                MPI_Recv(&offset, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(&rows, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(&m[offset][0], rows * SIZE, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &status);
            }
        }
        else{
            MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&m[offset][0], rows * SIZE, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
            for(int i = offset; i < rows; i++){
                line_generator(m, i);
            }
            MPI_Send(&offset, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
            MPI_Send(&m[offset], rows * SIZE, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD);
        }
        // sync to all processes
        MPI_Bcast(&m[0][0], SIZE * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void vec_generator(){
    for(int i = 0; i < length; i ++){
        vec[i][0] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }
}