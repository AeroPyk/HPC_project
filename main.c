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
        iVERBOSE printf("A square number of proc is mandatory, found %d", wsize);
        return 1;
    }
    else if(floor(MDIM/sq) != MDIM/sq) {
        MPI_Finalize();
        iVERBOSE printf("With %d proc and size %d the submatrices can't be square", wsize, MDIM);
        return 1;
    }

    // ===============================

    int side = (int) sq; // squarre root of the number of proc
    int subside = MDIM/sq; // size of sub matrix

    int** pA; // sub matrix A
    int** pB; // sub matrix B


    int periods[NDIM] = {[0 ... NDIM-1] = 1}; // periods link both ends of the mesh (true or false)
    int dims[NDIM] = {0};  // Those values are filled by MPI_Dims_create
    MPI_Comm squareCom; // Com for the square topology

    // Create the topo
    MPI_Dims_create(wsize, NDIM, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NDIM, dims, periods, 1, &squareCom);

    // Experiment
    iUNIQUE {
        int** A;
        writeMat("A", MDIM, MDIM);
        A = loadMat("A");

        printMat(A, MDIM, MDIM);

        pA = loadSubMatFromMat(A, subside, subside, squareCom);

        printf("\n");

        printMat(pA, subside, subside);

        free(pA);
        free(A);

    }

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


    MPI_Finalize(); // Mandatory

    return 0;
}

