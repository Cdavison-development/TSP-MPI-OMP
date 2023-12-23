#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>
#include <float.h>
#include <time.h>

#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords();
void *writeTourToFile();
#endif


    //find nearest neighbor for a node 
    void findNearestNeighbor(double *matrix1D, int numOfCoords, int currentNode, int *nearestIndex, double *nearestDistance) {
    //global values to write local values to outside the parallelised area
    //initialise global values with max possible value and invalid index    
    double globalNearestDistance = DBL_MAX;
    int globalNearestIndex = -1;

    //parallel region
    #pragma omp parallel
    {
        //initialise local variables
        double localNearestDistance = DBL_MAX;
        int localNearestIndex = -1;
        //parallel loop
        #pragma omp for
        for (int i = 0; i < numOfCoords; i++) {
            if (i != currentNode) {             
                double dist = matrix1D[currentNode * numOfCoords + i];  //calculate distance from i to node using 1D matrix
                if (dist < localNearestDistance) {
                    localNearestDistance = dist;
                    localNearestIndex = i;
                }
            }
        }
        //update global values
        #pragma omp critical
        {
            if (localNearestDistance < globalNearestDistance) {
                globalNearestDistance = localNearestDistance;
                globalNearestIndex = localNearestIndex;
            }
        }
    }
    *nearestDistance = globalNearestDistance;
    *nearestIndex = globalNearestIndex;
}
 

int *cheapestInsertion(double *matrix1D, int numOfCoords, int startVertex) {
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    int *visited = malloc(numOfCoords * sizeof(int));
    double (*matrix2D)[numOfCoords] = (double (*)[numOfCoords])matrix1D;
    
    printf("Starting cheapestInsertion2\n");
    // Initialize visited array
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0;
    }
    // intialise startvertex
    
    tour[0] = startVertex;
    visited[startVertex] = 1;
    int tourSize = 1; 

    // Find the nearest neighbor to startvertex and add to the tour
    double nearestDistance;
    int nearestIndex;
    findNearestNeighbor(matrix1D, numOfCoords, startVertex, &nearestIndex, &nearestDistance); 
   
     
    //complete initial loop
    tour[1] = nearestIndex; 
    tour[2] = startVertex; 
    visited[nearestIndex] = 1;
    tourSize = 3; 

    //allocate global variables
    double globalminCost;
    int globalminCostIndex, globalinsertPosition;

    //allocate memory for arrays
    int numThreads = omp_get_max_threads();
    double *localMinCosts = malloc(numThreads * sizeof(double));
    int *localMinCostIndex = malloc(numThreads * sizeof(int));
    int *localInsertPositions = malloc(numThreads * sizeof(int));
    //loop until tour is complete
    while (tourSize < numOfCoords + 1) {
        globalminCost = DBL_MAX;
        globalminCostIndex = -1;
        globalinsertPosition = -1;
        //find cheapest node in parallel section
        #pragma omp parallel
        {
            //initialize local variables with threads 
            int threadID = omp_get_thread_num();
            localMinCosts[threadID] = DBL_MAX;
            localMinCostIndex[threadID] = -1;
            localInsertPositions[threadID] = -1;
        //loop over visited array 
        #pragma omp for
            for (int i = 0; i < numOfCoords; i++) {
                if (!visited[i]) {  // Check if i is unvisited
                    for (int j = 0; j < tourSize - 1; j++) {
                        double cost = matrix1D[tour[j] * numOfCoords + i] + matrix1D[i * numOfCoords + tour[j + 1]] - matrix1D[tour[j] * numOfCoords + tour[j + 1]];
                        if (cost < localMinCosts[threadID]) {
                            localMinCosts[threadID] = cost;
                            localMinCostIndex[threadID] = i;
                            localInsertPositions[threadID] = j + 1;
                        }
                    }
                }
            }
        //update global values
        #pragma omp critical
        {
           if (localMinCosts[threadID] < globalminCost) {
                    globalminCost = localMinCosts[threadID];
                    globalminCostIndex = localMinCostIndex[threadID];
                    globalinsertPosition = localInsertPositions[threadID];
                } 
            }
        }
        //remove visited node
        if (globalminCostIndex != -1) {
            for (int i = tourSize; i > globalinsertPosition; i--) {
                tour[i] = tour[i - 1];
            }
            tour[globalinsertPosition] = globalminCostIndex;
            visited[globalminCostIndex] = 1;  // Mark this node as visited
            tourSize++;
        }
    }
    
    free(visited);
    free(localMinCosts);
    free(localMinCostIndex);
    free(localInsertPositions);
    return tour;

    //free(tour);
}

 


/* 
int main(int argc, char *argv[]) {
    // Ensure correct usage
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    char *outputFile = argv[2];


    
    // Read coordinates
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    double *matrix1D = malloc(numOfCoords * numOfCoords * sizeof(double));
    // Generate distance matrix
    GenerateDistanceMatrix(matrix1D, coords, numOfCoords);

    // Apply the cheapest insertion algorithm
    int *tour = cheapestInsertion2(matrix1D, numOfCoords);

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);

    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
    }
    
    
    free(coords);
    free(matrix1D);
    free(tour);

    return 0;
}   */