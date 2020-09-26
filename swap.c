#include "heads.h"

int re_arrange(nthreads){
    for(int i = 0; i < SIZE; i++){
        a = find_maxrow(nthreads, i);
        swap(a, i);
    }
}

int find_maxrow(int nthreads, int col) {
    if (nthreads == 1) {
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
    else {
        // initialize two arrays to store max number and its row index
        int* storage = (int*)malloc(nthreads * sizeof(int));
        memset(storage, INT_MIN, nthreads * sizeof(int));
        int* indexes = (int*)malloc(nthreads * sizeof(int));
        for (int i = 0; i < nthreads; i++) {
            indexes[i] = i;
        }
        int jump = nthreads;

        omp_set_num_threads(nthreads);
#pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            for (int i = thread_id + col; i < SIZE; i += jump) {
                double cur = abs(matrix[i][col]);
                if (cur > storage[thread_id]) {
                    indexes[thread_id] = i;
                    storage[thread_id] = cur;
                }
            }
        }
        double mx = 0;
        int idx = 0;
        for (int i = 0; i < nthreads; i++) {
            if (storage[i] > mx) {
                idx = indexes[i];
                mx = storage[i];
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