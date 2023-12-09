#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords();
void *writeTourToFile();
#endif


int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);

double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix[i][j] = distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return matrix;
}


void findFurthestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *furthestIndex, double *furthestDistance) {
    *furthestDistance = DBL_MIN; // Start with a negative value
    *furthestIndex = -1;

    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode) {
            double dist = distanceMatrix[currentNode][i];
            if (dist > *furthestDistance) {
                *furthestDistance = dist;
                *furthestIndex = i;
            }
        }
    }
}

    

int *furthestInsertion(double **distanceMatrix, int numOfCoords) {
    // Allocate memory for tour and unvisited nodes
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    int *unvisited = malloc(numOfCoords * sizeof(int));
    for (int i = 1; i < numOfCoords; i++) {
        unvisited[i - 1] = i;
    }
    int unvisitedCount = numOfCoords - 1;

    // Start with the first node
    tour[0] = 0;
    int tourSize = 1;

    // Find the furthest neighbor to 0 and add to the tour
    int furthestIndex;
    double furthestDistance;

    // Find the furthest neighbor from node 0
    findFurthestNeighbor(distanceMatrix, numOfCoords, 0, &furthestIndex, &furthestDistance);

    tour[1] = furthestIndex;
    tour[2] = 0; // Complete the initial loop
    tourSize = 3;

    // Remove the furthest neighbor from the unvisited list
     for (int i = 0; i < unvisitedCount; i++) {
        if (unvisited[i] == furthestIndex) {
            unvisited[i] = unvisited[unvisitedCount - 1];
            unvisitedCount--;
            break;
        }
    } 

    // Furthest insertion algorithm
    while (tourSize < numOfCoords + 1) {
        double MaxDist = DBL_MIN;
        int MaxDistIndex = -1, InsertPosition = -1;

    
    double maxNodeDist = -1;
        for (int idx = 0; idx < unvisitedCount; idx++) {
        int node = unvisited[idx];
            for (int j = 0; j < tourSize; j++) {
            double nodeDist = distanceMatrix[tour[j]][node];
                if (nodeDist > maxNodeDist) {
                    maxNodeDist = nodeDist;
                    MaxDistIndex = node;
                }
            }
        }

    double minInsertDist = DBL_MAX;
        for (int j = 0; j < tourSize - 1; j++) {
            double insertDist = distanceMatrix[tour[j]][MaxDistIndex] + distanceMatrix[MaxDistIndex][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
            if (insertDist < minInsertDist) {
                minInsertDist = insertDist;
                InsertPosition = j + 1;
            }
        }
        printf("Selected Node: %d, Insert Position: %d\n", MaxDistIndex, InsertPosition);
        

        if (MaxDistIndex != -1) {
    printf("Tour before insertion: ");
    for (int i = 0; i < tourSize; i++) {
        printf("%d ", tour[i]);
    }
    printf("\n");

    // Shift elements to make space for the new node
    for (int i = tourSize; i > InsertPosition; i--) {
        tour[i] = tour[i - 1];
    }

    // Insert the furthest node at the calculated position
    tour[InsertPosition] = MaxDistIndex;

    printf("Tour after insertion: ");
    for (int i = 0; i < tourSize + 1; i++) {
        printf("%d ", tour[i]);
    }
    printf("\n");

    printf("Unvisited list before update: ");
    for (int i = 0; i < unvisitedCount; i++) {
        printf("%d ", unvisited[i]);
    }
    printf("\n");

    // Remove the inserted node from the unvisited list
    for (int i = 0; i < unvisitedCount; i++) {
        if (unvisited[i] == MaxDistIndex) {
            unvisited[i] = unvisited[unvisitedCount - 1];
            unvisitedCount--;
            break;
        }
    }

    printf("Unvisited list after update: ");
    for (int i = 0; i < unvisitedCount; i++) {
        printf("%d ", unvisited[i]);
    }
    printf("\n");

    // Increment the tour size
    tourSize++;
    printf("Tour size updated to: %d\n", tourSize);
        }
    }
free(unvisited);
return tour;
}


int main(int argc, char *argv[]) {

     // Check if the correct number of arguments are passed
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <coordinate filename> <output filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // argv[1] is the coordinate file name
    // argv[2] is the output file name
    char *inputFilename = argv[1];
    char *outputFilename = argv[2];

    int numOfCoords = readNumOfCoords(inputFilename);

    if (numOfCoords == -1) {
        fprintf(stderr, "Error reading number of coordinates from %s.\n", inputFilename);
        return EXIT_FAILURE;
    }

    double **coords = readCoords(inputFilename, numOfCoords);

    if (coords == NULL) {
        fprintf(stderr, "Error reading coordinates from %s.\n", inputFilename);
        return EXIT_FAILURE;
    }

    // Calculate the distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);
    if (distanceMatrix == NULL) {
        fprintf(stderr, "Error calculating distance matrix.\n");
        return EXIT_FAILURE;
    }

    // Find the cheapest tour
    int *tour = furthestInsertion(distanceMatrix, numOfCoords);
    if (tour == NULL) {
        fprintf(stderr, "Error calculating cheapest insertion tour.\n");
        return EXIT_FAILURE;
    }

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFilename); // Note: numOfCoords + 1 to include the return to the starting node

    // Free the memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);
    free(tour);

    return EXIT_SUCCESS;
}