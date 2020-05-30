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
        iUNIQUE printf("A square number of proc is mandatory, found %d", wsize);
        return 1;
    }
    else if(floor(MDIM/sq) != MDIM/sq) {
        MPI_Finalize();
        iUNIQUE printf("With %d proc and size %d the submatrices can't be square", wsize, MDIM);
        return 1;
    }

    // ===============================

    int sideProc = (int) sq; // squarre root of the number of proc
    int subsideMat = MDIM / sideProc; // size of sub matrix

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
    MPI_Comm subLines[sideProc];
    MPI_Comm subCol[sideProc];
    int remainLine[2] = {0, 1}; // true / false for dimension
    int remainCol[2] = {1,0};

    for (int i = 0; i < sideProc; ++i) {
        MPI_Cart_sub(squareCom, remainLine, &subLines[i]);
        MPI_Cart_sub(squareCom, remainCol, &subCol[i]);
    }

    // ============================



    // Experiment

    // Row calculation with one process to compare
    iUNIQUE {
        writeRandMat("A", MDIM, MDIM); // Create A if needed: like when resizing MDIM !!!
        writeRandMat("B", MDIM, MDIM); // Create B

        A = loadMat("A"); // Load A mostly to show it
        B = loadMat("B");

        iVERBOSE printMat(A, MDIM, MDIM);
        iVERBOSE printMat(B, MDIM, MDIM);


        C = createMat(MDIM, MDIM);
        initZeroMat(C, MDIM, MDIM);

        perfMultiply(A, B, C, MDIM);

        iVERBOSE printMat(C, MDIM, MDIM);

        freeM(A, MDIM);
        freeM(B, MDIM);

    }

    MPI_Barrier(squareCom); // if A and B are being written, we have to wait for the procedure to end

    pA = loadSubMatFromFile("A", subsideMat, subsideMat, squareCom);
    pAb = copyMat(pA, subsideMat, subsideMat);
    pB = loadSubMatFromFile("B", subsideMat, subsideMat, squareCom);
    pC = createMat(subsideMat, subsideMat);
    initZeroMat(pC, subsideMat, subsideMat);

    // If you think of it, the broad casting can also be seen as a sequential broadcast of the same column with an offset ... instead of the diagonal + i.

    // Each process gets its rank and coordinates
    int** tmp = NULL;
    int root;
    int rank;
    MPI_Comm_rank(squareCom, &rank);
    int coords[NDIM];
    MPI_Cart_coords(squareCom, rank, NDIM, coords);

    int subrankLine;
    MPI_Comm_rank(subLines[coords[0]], &subrankLine);

    // For B we get the destination before since it doesn't change across the turns
    int src_rank, dest_rank;
    MPI_Cart_shift(subCol[coords[1]], 0, -1, &src_rank, &dest_rank);

    // printf("debug: dest %d , from [%d, %d]; ", dest_rank, coords[0], coords[1]);

    for (int turn = 0; turn < sideProc; ++turn) {


        // ============ A ==============

        root = (turn + coords[0])%sideProc;

        // When you're root you use your backup to calculate
        if(subrankLine == root) {
            tmp = pA;
            pA = pAb;
        }
        // Otherwise the one broadcasted


        MPI_Bcast(&(pA[0][0]), subsideMat * subsideMat, MPI_INT, root, subLines[coords[0]]); // Every process calls it!! That's why the root is precised

        // iVERBOSE printf("from [%d, %d]\n", coords[0], coords[1]);
        // iVERBOSE printMat(pA, subsideMat, subsideMat);
        // iVERBOSE printf("x\n");

        // ============ C ==============
        // iVERBOSE printMat(pB, subsideMat, subsideMat);
        // iVERBOSE printf("=\n");

        perfMultiply(pA, pB, pC, subsideMat);
        // iVERBOSE printMat(pC, subsideMat, subsideMat);

        // ============ B ==============

        MPI_Sendrecv_replace(&(pB[0][0]), subsideMat*subsideMat, MPI_INT, dest_rank, 100, src_rank, MPI_ANY_TAG, subCol[coords[1]], MPI_STATUS_IGNORE);



        // At the end we have to put back the original pA or we're going to lose the allocated space
        if(subrankLine == root) {
            pA = tmp;
        }


        //iUNIQUE printMat(pAb, subsideMat, subsideMat);

    }

    // C result of each process
    // iVERBOSE printf("from [%d, %d]\n", coords[0], coords[1]);
    // iVERBOSE printMat(pC, subsideMat, subsideMat);


    // Gather

    // Caveat: wsize doesnt come from the cartesian
    int** Cb = createMat(wsize, subsideMat*subsideMat);
    MPI_Gather(&(pC[0][0]), subsideMat*subsideMat, MPI_INT, &(Cb[0][0]), subsideMat*subsideMat, MPI_INT, 0, squareCom);

    iUNIQUE {
        linesToMat(&Cb, MDIM, MDIM, subsideMat, subsideMat, squareCom);
        iVERBOSE printMat(Cb, MDIM, MDIM);
        if(!equalMat(C, Cb, MDIM, MDIM)){
            printf("Matrices not equal");
        } else {
            printf("Tutto bene\n");
            writeMat("C", Cb, MDIM, MDIM);
        }
    }

    // Free up

    freeM(pA, subsideMat);
    freeM(pAb, subsideMat);
    freeM(pB, subsideMat);
    freeM(pC, subsideMat);
    freeM(C, MDIM);
    freeM(Cb, MDIM);

    //

    MPI_Finalize(); // Mandatory

    return 0;
}

