#include "heads.h"

int main(int argc, char* argv[]){
    printf("Start Processing......\n");
    start_t = clock();
    if(argc > 2){
        print("Too much Arguments\n");
        exit(EXIT_FAILURE);
    }
    if(argc == 2){
        if(argv[1][0] == '0'){
            printf("Invalid Input, only positive integers are accepted as number of threads\n");
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < strlen(argv[1]); i++){
            if(argv[1][i] < '0' || argv[1][i] > '9'){
                printf("Invalid Input, only integers are accepted as number of threads\n");
                exit(EXIT_FAILURE);
            }
        }
        nthreads = atoi(argv[1]); 
    }
    vec_generator()
    // start multiprocessing
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
    matrix_generator(matrix);
    re_arrange(nthreads);
}