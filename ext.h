//
// Created by Ce PC on 28/05/2020.
//

#ifndef PROJECTMM_EXT_H
#define PROJECTMM_EXT_H

#define NDIM 2 // For topology dimension should not be changed
#define MDIM 4 // For matrix size (MDIM * MDIM) // Nb proc <= MDIM*MDIM // Tested with MDIM = 128 and 64 proc on my laptop

// MIN and MAX value to write the matrices
#define RMAX 50
#define RMIN -50

// macro to manage whether one or several process do something
#define uRANK 0
#define UNIQUE (wrank == uRANK) // Only wrank proc print
#define ELSE 1 && !UNIQUE
#define VERBOSE 0 // 0: no peculiar message; 1: show the matrices

#define iUNIQUE if(UNIQUE) // only one process
#define iELSE if(ELSE) // All the other process
#define iVERBOSE if(VERBOSE)

// debug topo
void printTopo(MPI_Comm com);
void printCoord(int rank, MPI_Comm com);
int printNeigh(int rank, MPI_Comm com, int dir, int disp);


// Matrices

// Just shows the matrice
void printMat(double** mat, int row, int col);
// Initialize the matrice in argument with random values between RMAX and RMIN (double)
void initRandMat(double**, int, int);
// Initialize with 0s
void initZeroMat(double** mat, int row, int col);
// Write in the file {name} a random matrices of row*col size.
void writeRandMat(char *name, int row, int col);
// Write {mat} in the file {name}: format = "row col value val val val val ..."
void writeMat(char *name, double** mat, int row, int col);
// From {name} generate a matrix (memory allocation included)
double** loadMat(char *name);

// Load partial matrices (submatrices)
// Either from the full matrix
double** loadSubMatFromMat(double** mat, int row, int col, MPI_Comm com);
// Or directly from the file
double** loadSubMatFromFile(char* name, int row, int col, MPI_Comm com);

// Allocate memory in a special way to have consecutive memory for the matrix
double** createMat(int row, int col);

// Multiply within a process
void perfMultiply(double** A, double** B, double** C, int size);

// Copy a matrix (includes memory allocation) (back up of pA)
double** copyMat(double** in, int row, int col);

// free up the matrix and assign a NULL value to avoid crash if freeing up twice or more
void freeM(double ***pmat);

// When gathering the data from every processes (all sub matrices), we need to reformat since one line correspond to one sub matrix
// /!\ this replace definiely the input matrice (the previous one is freed up)
void linesToMat(double*** mat, int row, int col, int subrow, int subcol, MPI_Comm comm);

// Check whether both matrices are identical (+/-1e(-8))
int equalMat(double** A, double** B, int row, int col);

// Generate a random double between RMAX and RMIN
double drand ();

#endif //PROJECTMM_EXT_H
