# LU_Decomposition_MPIVersion

**Author: Yiling Liu**

**Student ID: 22214014**

## Project Scope(Requirement)

1. Random generation of matrix A and vector b, whose elements are double-precision floating-point numbers (same as Project 1);
2. Use a single CPU process to perform LU decomposition. Using the matrices L and U, implement a serial solution to find vector x.
3. Write serial code to save the matrices A, L, U and vector x to files (say, the files name are “mat_A”, “mat_L”, “mat_U” and “vec_x”, respectively). Similarly, write serial code to load/read the files into main memory to the matrices A, L, U and vector x. Make sure the content saved to the files is the same as the content read.
4. Implement a parallel solution using MPI to optimise the serial version done in items 2 and 3 above.
5. A verification function that shows the correctness of your results

## LU Decomposition

### Precondition

1. det(matrix) != 0
2. No pivot number == 0 while eliminating elements

To make it easier in coding, assume `len(matrix[0]) == len(matrix)` in this project.

### Infrastructure
Input: Matrix A

Output: Matrix L and U

n: length of Matrix

**Serial**

```
Initialize: L = 0, U = 0
for i in range(0, n) do
    for k in range(i, n) do
        sum = 0
        for j in range(0, i) do
            sum = sum + L[i][j] * U[i][j]
        end for
        U[i][k] = A[i][k] - sum
    end for
    for k in range(i, n) do
        if i = k then
            L[i][i] = 1
        else
            sum = 0
            for j in range(0, i) do
                sum = sum + L[k][j] * U[j][i]
            end for
            L[k][i] = (A[k][i] - sum)/U[i][i]
        end if
    end for
end for
```
Time Complexity: O(n^3)

**Parallel**

p: number of tasks
```
Initialize: L = 0, U = 0
for i in range(0, n) do
    for k in range(i, n) do
        sum, local = 0
        for j in  range(p, i) , j = j + p, do
            local = local + L[i, j] * U[i, j]
        end for
        MPI.Allreduce(local -> sum)
        U[i, k] = A[i, k] - sum
    end for
    for k in range(i, n) do
        if i = k then
            L[i, i] = 1
        else
            sum, local = 0
            for j in range(p, i) , j = j + p, do
                local = local + L[k, j] * U[j, i]
            end for
            MPI.Allreduce(local -> sum)
            L[k, i] = (A[k, i] - sum)/U[i, i]
        end if
    end for
end for
```
Time Complexity: O(n^3)
1/p time of serial version if the number of process is much smaller than the length of the matrix, but not efficient when a large number of processor is available.For example, if the number of the processor is equal to n, 1, 2... n - 1 processors will be idle in the for loop when i > 0. There is an algorithm that uses O(n^2) processors to achieve O(nlogn) time complexity [here](https://ieeexplore.ieee.org/document/143617)


