//import libraries and coordReader Files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>
#include <float.h>
#include <time.h>
#include "utils.h"
#include "coordReader.h"


//Provided functions
//int readNumOfCoords(char *filename);
//double **readCoords(char *filename, int numOfCoords);
//void *writeTourToFile(int *tour, int tourLength, char *filename);

//Function to calculate the Euclidean distance between two points
  double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

//Function to generate a distance matrix given a set of coordinates
double **generateDistanceMatrix(double **coords, int numOfCoords) {
    //Dynamically allocate memory for 2-D array
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }
    //generate empty distance matrix of size numOfCoords using parallel for loop  to distribute loop iterations within the team of threads
    #pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix[i][j] = distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return matrix;
}  

 //Find the nearest neighbor for a given node in distance matrix.
void findFurthestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *furthestIndex, double *furthestDistance) {
    double maxDistance = DBL_MIN; //Initialize with the smallest possible double value
    int maxIndex = -1;             //Initialize with an invalid index


    //Parallel loop over all nodes to find the furthest neighbor
    #pragma omp parallel for reduction(max:maxDistance)
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode) {                             //Skip the currentNode 
            double dist = distanceMatrix[currentNode][i];   //Get the distance from currentNode to node i
            if (dist > maxDistance) {
                maxDistance = dist;
                //Critical section to update the furthest neighbor's index in a thread-safe manner
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

//Cheapest insertion algorithm for TSP problem
int *furthestInsertion(double **distanceMatrix, int numOfCoords,int start_vertex) {
    // Allocate memory for tour and unvisited nodes
    int *tour = malloc((numOfCoords + 1) * sizeof(int));    //Allocate tour memory
    int *unvisited = malloc(numOfCoords * sizeof(int));     //Allocate unisited memory for list of unvisited nodes
    // Initialize unvisited nodes starting from 1
    for (int i = 1; i < numOfCoords; i++) {
        unvisited[i - 1] = i;
    }
    int unvisitedCount = numOfCoords;

    // Start with the first node
    //int startVertex = 4;
    tour[0] = start_vertex;
    int tourSize = 1;

    // Find the furthest neighbor to startVertex
    int furthestIndex;
    double furthestDistance;

    // Find the furthest neighbor from startVertex
    findFurthestNeighbor(distanceMatrix, numOfCoords, start_vertex, &furthestIndex, &furthestDistance);

    // Add the furthest neighbor to the tour and close the loop
    tour[1] = furthestIndex;
    tour[2] = start_vertex; // Close the loop by coming back to the start vertex
    tourSize = 3;

    /*  int newUnvisitedCount = 0;
    for (int i = 0; i < unvisitedCount; i++) {
        if (unvisited[i] != startVertex && unvisited[i] != furthestIndex) {
            unvisited[newUnvisitedCount++] = unvisited[i];
        }
    }
        unvisitedCount = newUnvisitedCount;  */
    for (int i = 0; i < unvisitedCount; i++) {
    if (unvisited[i] == start_vertex) {
        for (int j = i; j < unvisitedCount - 1; j++) {
            unvisited[j] = unvisited[j + 1];
        }
        unvisitedCount--;
        break;
    }
}

// Remove the furthest neighbor from the unvisited list
/* for (int i = 0; i < unvisitedCount; i++) {
    if (unvisited[i] == furthestIndex) {
        unvisited[i] = unvisited[unvisitedCount - 1];
        unvisitedCount--;
        break;
    }
}  */
    printf("Initial tour: ");
    for (int i = 0; i < tourSize; i++) {
    printf("%d ", tour[i]);
    }
    printf("\n");  

    //Prepare for the parallelized calculation of the cheapest insertion
    //Allocating memory for relevant array variables

    double MaxDist;
    int MaxDistIndex,InsertPosition;
    int numThreads = omp_get_max_threads();
    double *localMaxDist = malloc(numThreads * sizeof(double));
    int *localMaxDistIndex = malloc(numThreads * sizeof(int));
    int *localInsertPositions = malloc(numThreads * sizeof(int));
    double *localMinInsertDist = malloc(numThreads * sizeof(double));

    //Tour construction until all nodes are visited
    while (tourSize < numOfCoords + 1) {
        MaxDist = DBL_MIN;
        MaxDistIndex = -1, InsertPosition = -1;
        double minInsertDist = DBL_MAX;
        
        //printf("Finding furthest node to insert...\n");
        //Parallel region to find the furthest node to insert
        #pragma omp parallel
        {
            //Assign threads to each local variable
            int threadID = omp_get_thread_num();
            localMaxDist[threadID] = DBL_MIN;
            localMaxDistIndex[threadID] = -1;
            localInsertPositions[threadID] = -1;
            //Parallel loop over unvisited nodes
            #pragma omp for
            //Find the node farthest from any node in the tour and its insertion position
            for (int idx = 0; idx < unvisitedCount; idx++) {
                int node = unvisited[idx];
                //printf("Checking node %d\n", unvisited[idx]);
                for (int j = 0; j < tourSize; j++) {
                    double nodeDist = distanceMatrix[tour[j]][node];    //Calculate distance between current node and node being considered for insertion
                    if (nodeDist > localMaxDist[threadID]) {            //Update Max distance and max distance index
                        localMaxDist[threadID] = nodeDist;
                        localMaxDistIndex[threadID] = node;
                    }
                }
            }
            //critical section to update the global maximum distance and insertion position
            #pragma omp critical
           {
                if (localMaxDist[threadID] > MaxDist) {
                    MaxDist = localMaxDist[threadID];
                    MaxDistIndex =  localMaxDistIndex[threadID];
                }
            }
        }
        //printf("MaxDistIndex: %d\n", MaxDistIndex);
        //double minInsertDist = DBL_MAX;
        #pragma omp parallel
        {
            //Obtain the thread ID for managing local variables specific to each thread
            //Initialize local minimum insertion distance and position for this thread
            int threadID = omp_get_thread_num();
            localMinInsertDist[threadID] = DBL_MAX;
            localInsertPositions[threadID] = -1;
            
            
            //loop to find best insertion point
            for (int j = 0; j < tourSize - 1; j++) {
                //calculate additional cost of inserting furthest node
                double insertDist = distanceMatrix[tour[j]][MaxDistIndex] + distanceMatrix[MaxDistIndex][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                //Update local minimum insertion distance and position
                if (insertDist < localMinInsertDist[threadID]) {
                    localMinInsertDist[threadID] = insertDist;
                    localInsertPositions[threadID] = j + 1;
                }
            }
            //critical section to update global minimum insertion distance and position
           #pragma omp critical
           {
                if (localMinInsertDist[threadID] < minInsertDist) {
                    minInsertDist = localMinInsertDist[threadID];
                    InsertPosition = localInsertPositions[threadID];
                }
           }
        }
        //printf("InsertPosition: %d\n", InsertPosition);
        if (MaxDistIndex != -1) {
            // Shift elements to make space for the new node
            for (int i = tourSize; i > InsertPosition; i--) {
                tour[i] = tour[i - 1];
            }
            // Insert the furthest node at the calculated position
            tour[InsertPosition] = MaxDistIndex;

            // Remove the inserted node from the unvisited list
            for (int i = 0; i < unvisitedCount; i++) {
                if (unvisited[i] == MaxDistIndex) {
                    unvisited[i] = unvisited[unvisitedCount - 1];
                    unvisitedCount--;
                    break;
                }
            }
            // Increment the tour size
            tourSize++;

            printf("final tour: ");
            for (int i = 0; i < tourSize; i++) {
            printf("%d ", tour[i]);
            }
            printf("\n");
    
            }
        }
    //Free allocated memory
    free(localMaxDist);
    free(localMaxDistIndex);
    free(localInsertPositions);
    free(localMinInsertDist);
    free(unvisited);
    return tour;

    } 



int main(int argc, char *argv[]) {

     // Check if the correct number of arguments are passed
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <coordinate filename> <output filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    //Initialising arguments
    char *inputFilename = argv[1];
    char *outputFile = argv[2];

    int numOfCoords = readNumOfCoords(inputFilename);
    //test case in case of break
    if (numOfCoords == -1) {
        fprintf(stderr, "Error reading number of coordinates from %s.\n", inputFilename);
        return EXIT_FAILURE;
    }

    double **coords = readCoords(inputFilename, numOfCoords);
    //test case in case of break
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
    writeTourToFile(tour, numOfCoords + 1, outputFile); // Note: numOfCoords + 1 to include the return to the starting node
    
    
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