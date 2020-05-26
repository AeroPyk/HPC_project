#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define NDIM 2
#define MDIM 4
#define UNIQUE wrank == 0 // Only wrank proc print
#define VERBOSE 1 && UNIQUE

void printTopo(MPI_Comm com);
void printCoord(int rank, MPI_Comm com);
int printNeigh(int rank, MPI_Comm com, int dir, int disp);
void printMat(int[][MDIM]);
void initMat(int[][MDIM]);

int main(int argc, char *argv[]) {

    srand(time(NULL));
    int A[MDIM][MDIM];
    int B[MDIM][MDIM];
    initMat(A);
    initMat(B);

    int wrank, wsize, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    if(!VERBOSE) {
        printMat(A);
        printf("\n");
        printMat(B);
    }

    // Need a square number of proc
    double sq = sqrt(wsize);
    if (floor(sq) != sq) {
        MPI_Finalize();
        if(wrank == 0) printf("A square number of proc is mandatory, found %d", wsize);
        return 1;
    }
    int side = (int) sq;

    int periods[NDIM] = {[0 ... NDIM-1] = 1}; // periods link both ends of the mesh (true or false)
    int dims[NDIM] = {0};  // Those values are filled by MPI_Dims_create
    MPI_Comm squareCom;

    MPI_Dims_create(wsize, NDIM, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NDIM, dims, periods, 1, &squareCom);

    if(!VERBOSE) {
        printTopo(squareCom);
    }

    MPI_Comm subLines[side];
    MPI_Comm subCol[side];
    int remainLine[2] = {0, 1}; // true / false for dimension
    int remainCol[2] = {1,0};

    // We create sub com for each line
    for (int i = 0; i < side; ++i) {
        MPI_Cart_sub(squareCom, remainLine, &subLines[i]);
        MPI_Cart_sub(squareCom, remainCol, &subCol[i]);
    }



    if(VERBOSE) {
        printNeigh(0, subLines[0], 0, 1);
        printTopo(subCol[0]);
    }


    MPI_Finalize();

    return 0;
}

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

void printMat(int mat[MDIM][MDIM]){

    for (int i = 0; i < MDIM; ++i) {
        for (int j = 0; j < MDIM; ++j) {
            printf("%d\t", mat[i][j]);
        }
        printf("\n");
    }
}

void initMat(int mat[MDIM][MDIM]){

    for (int i = 0; i < MDIM; ++i) {
        for (int j = 0; j < MDIM; ++j) {
            mat[i][j] = rand();
        }
    }
}
