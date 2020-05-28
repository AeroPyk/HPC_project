//
// Created by Ce PC on 28/05/2020.
//

#ifndef PROJECTMM_EXT_H
#define PROJECTMM_EXT_H

#define NDIM 2 // For dimension
#define MDIM 2 // For matrix size (MDIM * MDIM)

#define UNIQUE (wrank == 2) // Only wrank proc print
#define ELSE 1 && !UNIQUE
#define VERBOSE 1

#define iUNIQUE if(UNIQUE)
#define iELSE if(ELSE)
#define iVERBOSE if(VERBOSE)

void printTopo(MPI_Comm com);
void printCoord(int rank, MPI_Comm com);
int printNeigh(int rank, MPI_Comm com, int dir, int disp);

void printMat(int** mat, int row, int col);
void initRandMat(int**, int, int);
void initZeroMat(int** mat, int row, int col);
void writeRandMat(char *name, int row, int col);
void writeMat(char *name, int** mat, int row, int col);
int** loadMat(char *name);

int** loadSubMatFromMat(int** mat, int row, int col, MPI_Comm com);
int** loadSubMatFromFile(char* name, int row, int col, MPI_Comm com);

int** createMat(int row, int col);

void perfMultiply(int** A, int** B, int** C, int size);
int** copyMat(int** in, int row, int col);

#endif //PROJECTMM_EXT_H
