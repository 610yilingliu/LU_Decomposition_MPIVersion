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

1. `len(matrix) == len(matrix[0])`
2. det(matrix) != 0
3. No pivot number == 0 while eliminating elements

Condition 1 is already excluded in this project because `len(matrix)` will forced to equal to `len(matrix[0])`.