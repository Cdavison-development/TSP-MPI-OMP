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
 

//Find furthest node 
void findFurthestNeighbor(double *matrix1D, int numOfCoords, int currentNode, int *furthestIndex, double *furthestDistance) {
    double maxDistance = DBL_MIN;
    int maxIndex = -1;
    //loop over all nodes to find furthest neighbor
    #pragma omp parallel for reduction(max:maxDistance)
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode) {
            double dist = matrix1D[currentNode * numOfCoords + i];
            if (dist > maxDistance) {
                maxDistance = dist;
                //update furthest neighbor index
                #pragma omp critical
                {
                    if (dist > *furthestDistance) {
                        *furthestDistance = dist;
                        *furthestIndex = i;
                    }
                }
            }
        }
    }
}

int *farthestInsertion(double *matrix1D, int numOfCoords, int startVertex) {
    // Allocate memory for tour and unvisited nodes
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    int *visited = malloc(numOfCoords * sizeof(int));
    printf("Starting furthestInsertion2\n");
    // Initialize visited array
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0;
    }
    
    // Start start vertex
    tour[0] = startVertex;
    visited[startVertex] = 1;
    int tourSize = 1; 

    // Find the furthest neighbor to 0 and add to the tour
    int furthestIndex;
    double furthestDistance;

    // Find the furthest neighbor from node 0
    findFurthestNeighbor(matrix1D, numOfCoords, startVertex, &furthestIndex, &furthestDistance);

    tour[1] = furthestIndex; 
    tour[2] = startVertex; 
    visited[furthestIndex] = 1;
    tourSize = 3; 

    double MaxDist;
    int MaxDistIndex,InsertPosition;

    int numThreads = omp_get_max_threads();
    double *localMaxDist = malloc(numThreads * sizeof(double));
    int *localMaxDistIndex = malloc(numThreads * sizeof(int));
    int *localInsertPositions = malloc(numThreads * sizeof(int));
    double *localMinInsertDist = malloc(numThreads * sizeof(double));
    // Furthest insertion algorithm

    while (tourSize < numOfCoords + 1) {
        MaxDist = DBL_MIN;
        MaxDistIndex = -1, InsertPosition = -1;
        double minInsertDist = DBL_MAX;
        
        //find furthest node 
        #pragma omp parallel
        {
         
            int threadID = omp_get_thread_num();
            localMaxDist[threadID] = DBL_MIN;
            localMaxDistIndex[threadID] = -1;
            localInsertPositions[threadID] = -1;

            #pragma omp for
            for (int idx = 0; idx < numOfCoords; idx++) {
                if (!visited[idx]){
                for (int j = 0; j < tourSize; j++) {
                    double nodeDist = matrix1D[tour[j] * numOfCoords + idx];
                    if (nodeDist > localMaxDist[threadID]) {
                        localMaxDist[threadID] = nodeDist;
                        localMaxDistIndex[threadID] = idx;
                    }
                }
            }
            }
            #pragma omp critical
           {
                if (localMaxDist[threadID] > MaxDist) {
                    MaxDist = localMaxDist[threadID];
                    MaxDistIndex =  localMaxDistIndex[threadID];
                }
            }
        }

        //double minInsertDist = DBL_MAX;
        #pragma omp parallel
        {
            int threadID = omp_get_thread_num();
            localMinInsertDist[threadID] = DBL_MAX;
            localInsertPositions[threadID] = -1;

            for (int j = 0; j < tourSize - 1; j++) {
                double insertDist = matrix1D[tour[j] * numOfCoords + MaxDistIndex] + matrix1D[MaxDistIndex * numOfCoords + tour[j + 1]] - matrix1D[tour[j] * numOfCoords + tour[j + 1]];
                if (insertDist < localMinInsertDist[threadID]) {
                    localMinInsertDist[threadID] = insertDist;
                    localInsertPositions[threadID] = j + 1;
                }
            }

            #pragma omp critical
            {
                if (localMinInsertDist[threadID] < minInsertDist) {
                    minInsertDist = localMinInsertDist[threadID];
                    InsertPosition = localInsertPositions[threadID];
                }
            }
        }

        if (MaxDistIndex != -1) {
            // Shift elements to make space for the new node
            for (int i = tourSize; i > InsertPosition; i--) {
                tour[i] = tour[i - 1];
            }

            // Insert the furthest node at the calculated position
            tour[InsertPosition] = MaxDistIndex;
            visited[MaxDistIndex] = 1;  // Mark this node as visited

            // Increment the tour size
            tourSize++;
            //printf("Tour size updated to: %d\n", tourSize);
        }
    }
    
    free(localMaxDist);
    free(localMaxDistIndex);
    free(localInsertPositions);
    free(localMinInsertDist);
    free(visited);
    return tour;

    //free(tour);
}
 
