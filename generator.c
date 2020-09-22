#include "heads.h"

void matrix_generator(int numtasks){
    if numtasks== 1) {
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                matrix[i][j] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
            }
        }
    else{
        int taskid, numworkers, source, dest, mtype, rows, averow, extra, offset
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
        MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        numworkers = numtasks - 1
        if(taskid == MASTER){
            for (int i = 0; i < SIZE; i++) 
                for (int j = 0; j < SIZE; j++) 
                    matrix[i][j] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
                
            mtype = FROM_MASTER;
            averow = SIZE/n_workers;
            extra = SIZE%n_workers;
            for (dest=1; dest<=numworkers; dest++) {
                rows = (dest <= extra) ? averow + 1 : averow;
                MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                MPI_Send(&a[offset][0], rows*NCA, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
                MPI_Send(&b, SIZE * SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
                offset = offset + rows;
                }
            mtype = FROM_WORKER;
            for (i=1; i<=numworkers; i++) {
                source = i;
                MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&c[offset][0], rows*SIZE, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
            }
        }
    }
}

void vector_generator() {
    for (int i = 0; i < SIZE; i++) {
        vec[i][0] = (double)(rand() % (RANGE * 2) - RANGE) / SCALE;
    }