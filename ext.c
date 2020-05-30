//
// Created by Ce PC on 28/05/2020.
//

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "ext.h"

void printTopo(MPI_Comm com){

    // Check the cartesian:
    int ndims = 0;


    MPI_Cartdim_get(com, &ndims);
    printf("ndim = %d\n", ndims);

    int check_dims[ndims]; // size of each dimension
    int check_periods[ndims]; // whether it's periodic dim or not
    int check_coords[ndims]; // coordonates of the speaking proc
    MPI_Cart_get(com, ndims, check_dims, check_periods, check_coords);


    printf("dims\t prds\t coords\n");
    for (int i = 0 ; i < ndims ; i++) {
        printf("%d\t ", check_dims[i]);
        printf("%d\t ", check_periods[i]);
        printf("%d\n", check_coords[i]);
    }

}

void printCoord(int rank, MPI_Comm com){

    // Whether you want to print coordinate of a specific rank:

    int ndims = 0;
    MPI_Cartdim_get(com, &ndims);

    int check_rank;
    int check_coords[ndims];

    MPI_Cart_coords(com, rank, ndims, check_coords);

    MPI_Cart_rank(com, check_coords, &check_rank); // Just retrieve the rank with the coordinates

    printf("coords = [");
    for (int j = 0; j < ndims; ++j) {
        printf("%d,", check_coords[j]);
    }
    printf("] rank = %d\n", check_rank);


}

int printNeigh(int rank, MPI_Comm com, int dir, int disp){

    // source says: what is the rank of the process which is sending me those data if it calls the same MPI_Cart_shift
    int src_rank, dest_rank;
    MPI_Cart_shift(com, dir, disp, &src_rank, &dest_rank);
    printf("Rank %d; dest: %d; src: %d\n", rank, dest_rank, src_rank);

    return dest_rank;
}

void printMat(int** mat, int row, int col){

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            printf("%d\t", *(*(mat+i)+j));
        }
        printf("\n");
    }
    printf("\n");
}

void initRandMat(int** mat, int row, int col){

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            mat[i][j] = rand()%(RMAX - RMIN) + RMIN;
        }
    }
}

void initZeroMat(int** mat, int row, int col){

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            mat[i][j] = 0;
        }
    }
}

void writeRandMat(char *name, int row, int col) {

    int** mat;
    mat = createMat(row, col);
    initRandMat(mat, row, col);

    writeMat(name, mat, row, col);

    free(mat);

}

void writeMat(char *name, int** mat, int row, int col){
    FILE* file = fopen(name, "w");

    if(file != NULL){

        fprintf(file, "%d %d ", row, col);

        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                fprintf(file,"%d ", mat[i][j]);
            }
        }

        fclose(file);
    } else {
        printf("Error writting file");
    }
}

int** loadMat(char *name) {

    FILE* file = fopen(name, "r");

    if(file != NULL){
        int size[2];

        fscanf(file, "%d %d ", &size[0], &size[1]);

        int** mat = createMat(size[0], size[1]);

        for (int i = 0; i < size[0]; ++i) {
            for (int j = 0; j < size[1]; ++j) {
                fscanf(file, "%d ", &mat[i][j]);
            }

        }
        fclose(file);

        return mat;


    } else {
        printf("Error reading file");
        return NULL;
    }
}

int** loadSubMatFromMat(int** mat, int row, int col, MPI_Comm com) {
    int rank;
    MPI_Comm_rank(com, &rank);
    int coords[NDIM];
    MPI_Cart_coords(com, rank, NDIM, coords);

    int** submat = createMat(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            submat[i][j] = mat[i + row*coords[0] ][j + col*coords[1] ];
        }
    }

    return submat;

}

int** loadSubMatFromFile(char* name, int row, int col, MPI_Comm com) {

    int rank;
    MPI_Comm_rank(com, &rank);
    int coords[NDIM];
    MPI_Cart_coords(com, rank, NDIM, coords);

    FILE* file = fopen(name, "r");

    if(file != NULL){
        int size[2];

        fscanf(file, "%d %d ", &size[0], &size[1]);

        int** submat = createMat(row, col);
        int offset = 0;
        int last = 0;
        int num;
        int x, y = 0;

        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {

                x = i + row*coords[0];
                y = j + col*coords[1];
                last = offset;
                offset = x*size[1] + y + 1; // This +1 makes thing alright but not according to paper calculations...

                for (int k = 0; k < offset - last - 1; ++k) {
                    fscanf(file, "%d ", &num); // skip numbers in between
                }
                fscanf(file, "%d ", &submat[i][j]);

            }

        }

        // printf("debug %d,%d - %d - %s - %d\n", row, col, rank, name, submat[0][0]);
        fclose(file);

        return submat;

    } else {
        printf("Error reading file");
        return NULL;
    }

}

int** createMat(int row, int col){

    /*
    int **arr = (int **)malloc(row * sizeof(int *));
    for (int i=0; i<row; i++)
        arr[i] = (int *)malloc(col * sizeof(int));

    return arr;
    */

    // This allocation below get consecutive memory to simplify the transfer with MPI


    int *data = (int *)malloc(row*col*sizeof(int));
    int **array= (int **)malloc(row*sizeof(int*));
    for (int i=0; i<row; i++)
        array[i] = &(data[col*i]);

    return array;

}

void perfMultiply(int** A, int** B, int** C, int size){
    // Basic implementation that will be optimize later on.

    int i, j, k;
// #pragma omp parallel for
    for (i = 0 ; i < size ; i++) {
        for (k = 0 ; k < size ; k++) {
            for (j = 0 ; j < size ; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

}

int** copyMat(int** in, int row, int col){

    // especially for the back up of pA
    int** out = createMat(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            out[i][j] = in[i][j];
        }
    }

    return out;
}

void freeM(int ** mat, int row){

    // I wanted to free up matrices row by row like a 2D array but it didn't work supposedly because after one free,
    // since the allocation is consecutive in memory, nothing else is needed to be freed up
    // Need to confirm with valgrind or similar
    // if confirmed, int row can be deleted

    if(mat != NULL) {

        free(&(mat[0][0])); // This free the large bundle of allocated memory (see how the memory is allocated above)
        free(mat); // Then we free the array that was referring to starting lines
        mat = NULL;
    }
}

void linesToMat(int*** mat, int row, int col, int subrow, int subcol, MPI_Comm comm){

    int** transformed = createMat(row, col);
    int rank;
    int coord[2];

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {

            coord[0] = i/subrow;
            coord[1] = j/subcol;
            MPI_Cart_rank(comm, coord, &rank);
            // printCoord(rank, comm);
            // printf("debug: [%d,%d] [%d,%d]\n\n", i, j, rank, (i%subrow)*subrow + j%subcol);

            transformed[i][j] = (*mat)[rank][(i%subrow)*subrow + j%subcol];

        }
    }

    freeM(*mat, row);
    *mat = transformed;

}

int equalMat(int** A, int** B, int row, int col){

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (A[i][j] != B[i][j]){
                return 0;
            }
        }
    }
    return 1;
}