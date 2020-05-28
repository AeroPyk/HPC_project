#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "ext.h"


int main(int argc, char *argv[]) {

    srand(time(NULL));

    int wrank, wsize, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);


    // Need a square number of proc and that we can match with square submatrices so far
    double sq = sqrt(wsize);
    if (floor(sq) != sq) {
        MPI_Finalize();
        iELSE printf("A square number of proc is mandatory, found %d", wsize);
        return 1;
    }
    else if(floor(MDIM/sq) != MDIM/sq) {
        MPI_Finalize();
        iELSE printf("With %d proc and size %d the submatrices can't be square", wsize, MDIM);
        return 1;
    }

    // ===============================

    int side = (int) sq; // squarre root of the number of proc
    int subside = MDIM/sq; // size of sub matrix

    int** A = NULL; // matrix A
    int** B = NULL; // matrix B
    int** C = NULL; // matrix C

    int** pA = NULL; // sub matrix A
    int** pAb = NULL; // b for backup like when broadcasting, it is meant to be overwritten
    int** pB = NULL; // sub matrix B
    int** pC = NULL; // sub matrix C


    int periods[NDIM] = {[0 ... NDIM-1] = 1}; // periods link both ends of the mesh (true or false)
    int dims[NDIM] = {0};  // Those values are filled by MPI_Dims_create
    MPI_Comm squareCom; // Com for the square topology

    // Create the topo
    MPI_Dims_create(wsize, NDIM, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NDIM, dims, periods, 1, &squareCom);

    // We create sub com for each line
    MPI_Comm subLines[side];
    MPI_Comm subCol[side];
    int remainLine[2] = {0, 1}; // true / false for dimension
    int remainCol[2] = {1,0};

    for (int i = 0; i < side; ++i) {
        MPI_Cart_sub(squareCom, remainLine, &subLines[i]);
        MPI_Cart_sub(squareCom, remainCol, &subCol[i]);
    }

    // ============================



    // Experiment
    iUNIQUE {
        // writeRandMat("A", MDIM, MDIM); // Create A if needed
        // writeRandMat("B", MDIM, MDIM); // Create B

        A = loadMat("A"); // Load A mostly to show it
        B = loadMat("B");

        iVERBOSE printMat(A, MDIM, MDIM);
        iVERBOSE printMat(B, MDIM, MDIM);


        C = createMat(MDIM, MDIM);
        initZeroMat(C, MDIM, MDIM);

        perfMultiply(A, B, C, MDIM);

        iVERBOSE printMat(C, MDIM, MDIM);

        free(A);
        free(B);
        free(C);

    }

    MPI_Barrier(squareCom); // if A and B are being written, we have to wait for the procedure to end

    pA = loadSubMatFromFile("A", subside, subside, squareCom);
    pAb = copyMat(pA, subside, subside);
    pB = loadSubMatFromFile("B", subside, subside, squareCom);

    iVERBOSE printMat(pA, subside, subside);
    iVERBOSE printMat(pB, subside, subside);

    if(ELSE){

    }

    // Free up

    free(pA);
    free(pB);

    //

    MPI_Finalize(); // Mandatory

    return 0;
}

