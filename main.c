#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>


#include "ext.h"


int main(int argc, char *argv[]) {

    srand((unsigned int) time(NULL));

    int wrank, wsize, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided); // MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
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

    // Timing

    double startNaive = -1, stopNaive = -1, startMulti = -1, stopMulti = -1;

    // ===============================

    int sideProc = (int) sq; // squarre root of the number of proc
    int subsideMat = MDIM / sideProc; // size of sub matrix

    double** A = NULL; // matrix A
    double** B = NULL; // matrix B
    double** C = NULL; // matrix C
    double **Cb = NULL;

    double** pA = NULL; // sub matrix A
    double** pAb = NULL; // b for backup like when broadcasting, it is meant to be overwritten
    double** pB = NULL; // sub matrix B
    double** pC = NULL; // sub matrix C
    double** tmp = NULL;

    // ================= Topology ==================

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

    // Properties for each process
    // Each get its rank inside the square topology and sub line topology and coordinates in square topology
    int sqComRank;
    MPI_Comm_rank(squareCom, &sqComRank);
    int sqComCoords[NDIM];
    MPI_Cart_coords(squareCom, sqComRank, NDIM, sqComCoords);
    int subRankLine;
    MPI_Comm_rank(subLines[sqComCoords[0]], &subRankLine);

    // ===============================================

    // Calculation with one process to compare
    iUNIQUE {
        writeRandMat("A", MDIM, MDIM); // Create A if needed: like when resizing MDIM !!!
        writeRandMat("B", MDIM, MDIM); // Create B

        A = loadMat("A"); // Load A (from file)
        B = loadMat("B"); // Load B

        iVERBOSE printf("A:\n");
        iVERBOSE printMat(A, MDIM, MDIM);
        iVERBOSE printf("B:\n");
        iVERBOSE printMat(B, MDIM, MDIM);


        C = createMat(MDIM, MDIM);
        initZeroMat(C, MDIM, MDIM); // Necessary since we += on C

        startNaive = MPI_Wtime();
        perfMultiply(A, B, C, MDIM);
        stopNaive = MPI_Wtime();

        iVERBOSE printf("C:\n");
        iVERBOSE printMat(C, MDIM, MDIM);

        freeM(&A);
        freeM(&B);

    }

    MPI_Barrier(squareCom); // if A and B are being written, we have to wait for the procedure to end

    // ========== Multi processing ============

    pA = loadSubMatFromFile("A", subsideMat, subsideMat, squareCom);
    pAb = copyMat(pA, subsideMat, subsideMat);
    pB = loadSubMatFromFile("B", subsideMat, subsideMat, squareCom);
    pC = createMat(subsideMat, subsideMat);
    initZeroMat(pC, subsideMat, subsideMat);

    // If you think of it, the broad casting can also be seen as a sequential broadcast of the same column with an offset ... instead of the diagonal + i.

    int root; // Will hold the root rank for each broadcast

    // For B we get the destination before since it doesn't change across the turns
    int src_rank, dest_rank; // will recv pB from src_rank and send it to dest_rank
    MPI_Cart_shift(subCol[sqComCoords[1]], 0, -1, &src_rank, &dest_rank);

    // printf("debug: dest %d , from [%d, %d]; ", dest_rank, sqComCoords[0], sqComCoords[1]);

    startMulti = MPI_Wtime();
    for (int turn = 0; turn < sideProc; ++turn) {

        // ============ A ==============

        root = (turn + sqComCoords[0]) % sideProc;

        // When you're root you use your backup to calculate
        if(subRankLine == root) {
            tmp = pA;
            pA = pAb;
        }
        // Otherwise the one broadcasted

        MPI_Bcast(&(pA[0][0]), subsideMat * subsideMat, MPI_DOUBLE, root, subLines[sqComCoords[0]]); // Every process calls it!! That's why the root is precised

        // iVERBOSE printf("from [%d, %d]\n", sqComCoords[0], sqComCoords[1]);
        // iVERBOSE printMat(pA, subsideMat, subsideMat);
        // iVERBOSE printf("x\n");

        // ============ C ==============
        // iVERBOSE printMat(pB, subsideMat, subsideMat);
        // iVERBOSE printf("=\n");

        perfMultiply(pA, pB, pC, subsideMat);
        // From class about intel matrix multiplication: not tested yet since I cannot run this library at the time on my computer.
        // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, subsideMat, subsideMat, subsideMat, 1.0, &(A[0][0]), subsideMat, &(B[0][0]), subsideMat, 1.0, &(C[0][0]), subsideMat);

        // iVERBOSE printMat(pC, subsideMat, subsideMat);

        // ============ B ==============

        MPI_Sendrecv_replace(&(pB[0][0]), subsideMat*subsideMat, MPI_DOUBLE, dest_rank, 100, src_rank, MPI_ANY_TAG, subCol[sqComCoords[1]], MPI_STATUS_IGNORE);

        // At the end we have to put back the original pA or we're going to lose the allocated space
        if(subRankLine == root) {
            pA = tmp;
            tmp = NULL; // tmp is no more useful at this point and we put it back at NULL
        }

        //iUNIQUE printMat(pAb, subsideMat, subsideMat);

    }
    stopMulti = MPI_Wtime();

    // C result of each process
    // iVERBOSE printf("from [%d, %d]\n", sqComCoords[0], sqComCoords[1]);
    // iVERBOSE printMat(pC, subsideMat, subsideMat);

    // ================== Gather ====================

    // We split here to avoid all process but one creating Cb
    iUNIQUE {
        // Caveat: wsize doesnt come from the cartesian
        Cb = createMat(wsize, subsideMat * subsideMat);
        MPI_Gather(&(pC[0][0]), subsideMat * subsideMat, MPI_DOUBLE, &(Cb[0][0]), subsideMat * subsideMat, MPI_DOUBLE, uRANK, squareCom);
    }
    iELSE MPI_Gather(&(pC[0][0]), subsideMat * subsideMat, MPI_DOUBLE, NULL, subsideMat * subsideMat, MPI_DOUBLE, uRANK, squareCom);

    // Once the data gathered:
    iUNIQUE {
        linesToMat(&Cb, MDIM, MDIM, subsideMat, subsideMat, squareCom);
        iVERBOSE printMat(Cb, MDIM, MDIM);
        if(!equalMat(C, Cb, MDIM, MDIM)){
            printf("Matrices not equal");
        } else {
            printf("Va tutto bene\nNaive %lfs\nMulti %lfs", stopNaive-startNaive, stopMulti-startMulti);
            writeMat("C", Cb, MDIM, MDIM);
        }
    }

    // Free up

    freeM(&pA);
    freeM(&pAb);
    freeM(&pB);
    freeM(&pC);
    freeM(&A);
    freeM(&B);
    freeM(&C);
    freeM(&Cb);
    freeM(&tmp);

    //

    MPI_Finalize(); // Mandatory

    return 0;
}

