#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
MPI_Status status;
int SIZE = 5;
int MASTER = 0;
int taskid, ntasks, bufsize,mpierror;

int main(int argc, char **argv){
    
    MPI_File fh; //Declaring a File Pointer
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    double *matrix = (double*)malloc(SIZE * SIZE * sizeof(double));
    int i = 0;
    while(i < SIZE * SIZE){
        matrix[i] = i * 0.5;
        i += 1;
    }
    exporting(matrix, SIZE * SIZE, "testdata");
    MPI_Finalize();
return 0;
}

void exporting(double* arr_2d, int fsize, char* fname) {
    // save in csv mode, split by ','
    MPI_File wh;
    int ele_per_procs = fsize/ntasks;
    while(ele_per_procs * ntasks < fsize){
        ele_per_procs++;
    }
    // last process treat the rest
    int ele_last_procs = fsize - (ele_per_procs * (ntasks - 1));

    MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&wh);
    if(taskid == MASTER){
        // the first process only write the size of the matrix
        mpierror = MPI_File_write_at(wh, 0, &SIZE, 1, MPI_INT, &status);
        mpi_error_check(mpierror, taskid);
    }
    if(taskid == ntasks - 1){
        int mat_offset = taskid * ele_per_procs; 
        MPI_Offset offset = mat_offset + 1;
        printf("%d  %d %f\n", taskid, mat_offset, arr_2d[mat_offset]);
        mpierror = MPI_File_write_at_all(wh, offset * sizeof(double), &arr_2d[mat_offset], ele_last_procs, MPI_DOUBLE, &status);
        mpi_error_check(mpierror, taskid);
    }
    else{
        int mat_offset = taskid * ele_per_procs;
        MPI_Offset offset = mat_offset + 1;
        printf("%d  %d %f\n", taskid, mat_offset, arr_2d[mat_offset]);
        mpierror = MPI_File_write_at_all(wh, offset * sizeof(double), &arr_2d[mat_offset], ele_per_procs, MPI_DOUBLE, &status);
        mpi_error_check(mpierror, taskid);
    }
    MPI_File_close(&wh);
}

void mpi_error_check(int mpierror, int taskid){
    if(mpierror != MPI_SUCCESS){
        printf("Error mpi failure\n");
        printf("taskid %d\n", taskid);
        exit(EXIT_FAILURE);
    }
}