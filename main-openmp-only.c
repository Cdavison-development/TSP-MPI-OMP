//#include "ompcInsertion.c"
//#include "ompfInsertion.c"
//#include "ompnAddition.c"

#include <stdio.h>
#include <stdlib.h>
//#include "utils.h"
#include <omp.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "coordReader.h"

#ifndef OMPCINSERTION_H_
#define OMPCINSERTION_H_
    int *cheapestInsertion(double **distanceMatrix, int numOfCoords, int startVertex);
#endif
#ifndef OMPFINSERTION_H_
#define OMPFINSERTION_H_
    int *farthestInsertion(double **distanceMatrix, int numOfCoords, int startVertex);
#endif
#ifndef OMPNADDITION_H_
#define OOMPNADDITION_H_
    int *nearestInsertion(double **distanceMatrix, int numOfCoords, int startVertex);
#endif



 double euclideanDistance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return sqrt(dx * dx + dy * dy);
}

/* double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return matrix;
}  


double calculateTotalDistance(int *tour, double **distanceMatrix, int tourLength) {
    double totalDistance = 0.0;
    for (int i = 0; i < tourLength - 1; i++) {
        totalDistance += distanceMatrix[tour[i]][tour[i + 1]];
    }
    totalDistance += distanceMatrix[tour[tourLength - 1]][tour[0]]; // Complete the loop
    printf("Total distance: %f\n", totalDistance);
    return totalDistance;
}


/* void findBestTour(double **distanceMatrix, int numOfCoords, int alg, char *ciFilename, char *fiFilename, char *niFilename) {
    double globalBestCost = DBL_MAX;
    int *globalBestTour = malloc((numOfCoords + 1) * sizeof(int));
    

    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        int *currentTour;
        double currentCost;

        switch (alg) {
            case 0:
                currentTour = cheapestInsertion(distanceMatrix, numOfCoords, startVertex);
                break;
            case 1:
                currentTour = farthestInsertion(distanceMatrix, numOfCoords, startVertex);
                break;
            case 2:
                currentTour = nearestInsertion(distanceMatrix, numOfCoords, startVertex);
                break;
        }

        currentCost = calculateTotalDistance(currentTour, distanceMatrix, numOfCoords);

        if (currentCost < globalBestCost) {
            globalBestCost = currentCost;
            memcpy(globalBestTour, currentTour, (numOfCoords + 1) * sizeof(int));
        }

        free(currentTour);
    }

    printf("Best tour for algorithm %d: ", alg);
    for (int i = 0; i < numOfCoords + 1; i++) {
        printf("%d ", globalBestTour[i]);
    }
    printf("with a cost of: %f\n", globalBestCost);

    switch (alg) {
        case 0:
            writeTourToFile(globalBestTour, numOfCoords + 1, ciFilename);
            break;
        case 1:
            writeTourToFile(globalBestTour, numOfCoords + 1, fiFilename);
            break;
        case 2:
            writeTourToFile(globalBestTour, numOfCoords + 1, niFilename);
            break;
    }

    free(globalBestTour);
} */

typedef int* (*InsertionFunction)(double **, int, int);
void findBestTour(double **distanceMatrix, int numOfCoords, int hrs,char *filenames[]) {
    InsertionFunction insertionFunctions[3] = {cheapestInsertion, farthestInsertion, nearestInsertion};
    
    double globalBestCost = DBL_MAX;
    int *globalBestTour = malloc((numOfCoords + 1) * sizeof(int));
    

    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        int *currentTour = insertionFunctions[hrs](distanceMatrix, numOfCoords, startVertex);
        double currentCost = calculateTotalDistance(currentTour, distanceMatrix, numOfCoords);

        if (currentCost < globalBestCost) {
            globalBestCost = currentCost;
            memcpy(globalBestTour, currentTour, (numOfCoords + 1) * sizeof(int));
        }

        free(currentTour);
    }

    char *filename = filenames[hrs];
    writeTourToFile(globalBestTour, numOfCoords + 1, filename);

    free(globalBestTour);
}
int main(int argc, char *argv[]) {
    
    if (argc != 5) {
        printf("Usage: %s <coord_file> <cheapest output file> <furthest output file> <nearest output file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);
    char *filenames[3] = {argv[2], argv[3], argv[4]};

    for (int alg = 0; alg < 3; alg++) {
        findBestTour(distanceMatrix, numOfCoords, alg, filenames);
    }

    return 0;
}

    // ... Free any other allocated memory and clean up ...

   