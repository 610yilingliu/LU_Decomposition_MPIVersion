# LU_Decomposition_MPIVersion

**Author: Yiling Liu**

**Student ID: 22214014**

## Symbols in this report

procs: number of process in parallel code

T<sub>comm</sub>: Time of communicating

T<sub>comp</sub>: Time of computing


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

 Although the big O time complexity is the same as serial code, but T<sub>comp</sub> of parallel code only spents 1/procs time of serial version if the number of process is much smaller than the length of the matrix, but not efficient when a large number of processor is available.
 
 For example, if the number of the processor is equal to n, 1, 2... n - 1 processors will be idle in the for loop when i > 0. There is an algorithm that uses O(n^2) processors to achieve O(nlogn) time complexity [here](https://ieeexplore.ieee.org/document/143617). I will try this efficient version later if I have time, but this project aims to use MPI API correctly instead of applying better algorithms so I simply leave the simple C version here. 

## Code Structure
To make it easier to compile and debug, I will not use modules in different C files like what I did in project 1. Simply put all code that are available to generate a serial solution in [Serial Code File](./C_serial&parallel/serial_sol.c) and parallel in [Parallel Code File](./C_serial&parallel/mpi_sol.c)

All serial and parallel contains the following modules:

### matrix_generator

**Time complexity: O(n<sup>2</sup>)**

Generate a random `SIZE * SIZE` matrix. The value range depends on static value `RANGE` and `SCALE`, the same as project 1

The difference between serial and parallel code in this module is that parallel code uses master process(process 0) to assign generating task evenly to process 1-n while seral version go through the whole matrix in order. 

When procs==2, T<sub>comp</sub> of parallel version is the same as T<sub>comp</sub> of serial version, but spent extra O(1) time in communication between process 0 and process 1.

### vector_generator

Generate a `1 * SIZE` vector.
It seems unnecessary to do this part in parallel because the time complexity of generate an 1d array is only O(n), so both serial and parallel version simply go through the whole vector to generate it.


### re_arrange, swap, find_maxrow

Find the maximum value in the current column and swap the row which contains this value with the current row

Example:

```
1 2 3
3 1 2
1 3 2
```
will be
```
3 1 2
1 3 2
1 2 3
```
After applying those functions

### get_lu
mentioned in section Infrastructure

### count

calculate the answer based on L matrix and U matrix. The way that the parallel version uses `MPI_Allreduce` is similar with get_lu

### checker

Check if the answer we got is correct or not. Parallel version also uses `MPI_Allreduce` to sum up the local sum in each process.
The final answer we got for checking only keep 3 digits because 1 divide by 3 in computer will be 0.33333...(limited digit) so if we use 3 to multiply this answer, it will not equal to 1.
