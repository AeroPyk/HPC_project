//
// Created by Ce PC on 28/05/2020.
//

#ifndef PROJECTMM_EXT_H
#define PROJECTMM_EXT_H

#define NDIM 2 // For dimension
#define MDIM 8 // For matrix size (MDIM * MDIM)
#define UNIQUE (wrank == 2) // Only wrank proc print
#define VERBOSE 1 && !UNIQUE

#define iUNIQUE if(UNIQUE)
#define iVERBOSE if(VERBOSE)

void printTopo(MPI_Comm com);
void printCoord(int rank, MPI_Comm com);
int printNeigh(int rank, MPI_Comm com, int dir, int disp);

void printMat(int** mat, int row, int col);
void initMat(int**);
void writeMat(char name[], int row, int col);
int** loadMat(char *name);

int** loadSubMatFromMat(int** mat, int row, int col, MPI_Comm com);
int** loadSubMatFromFile(char* name, int row, int col, MPI_Comm com);

int** createMat(int row, int col);

#endif //PROJECTMM_EXT_H
