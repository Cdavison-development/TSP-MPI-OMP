//libaries and header files 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#ifndef coordReader_H_
#define coordReader_H_
    int readNumOfCoords(char *filename);
    double **readCoords(char *filename, int numOfCoords);           
    void writeTourToFile(int *tour, int tourLength, char *filename);
#endif

/* Checks is vertex is in tour by iterating over the tour and comparing each vertex
if the vertex is in the tour, returns 1, else returns 0 */
int IsVertexInTour(int vertex, int *tour, int tourSize) {
    for (int i = 0; i < tourSize; i++) {       
        if (tour[i] == vertex) {           
            return 1;
        }
    }    
    return 0;
} 

//Structure to group certain global variables relevant to the NearestNotInTour function
  typedef struct {
    int globalNearestVertex;
    int globalClosestTourVertex;
    double globalMinDistance;
} NearestVertexResult;

/* 
    iterates over NumOfCoords and checks if the vertex is in the tour using the IsVertexInTour function
    if it is not in the tour, it finds the distance to each vertex and if the distance is less than the current min distance
    it updates the min distance and the nearest vertex, eventually finding the nearest vertex not in the tour
     */
 NearestVertexResult NearestNotInTour(double *matrix1D, int *tour, int numOfCoords, int tourSize) {
    //initialise the result struct
    NearestVertexResult result;
    //initialise the global variables to be used in the parallel section
    result.globalNearestVertex = -1;
    result.globalClosestTourVertex = -1;
    result.globalMinDistance = DBL_MAX;

    //parallel section to find the nearest vertex not in the tour
    #pragma omp parallel
    {
        //initialise the local variables to be used within each thread to find the nearest vertex not currently in the tour
        int localNearestVertex = -1;
        int localClosestTourVertex = -1;
        double localMinDistance = DBL_MAX;  
        //parallel for loop to find closest vertex not currently in tour
        #pragma omp for
        for (int i = 0; i < numOfCoords; i++) {
            if (!IsVertexInTour(i, tour, tourSize)) { //params are vertex, tour, tourSize                
                for (int j = 0; j < tourSize; j++) { // Iterate over current tour size
                    double distance = matrix1D[i * numOfCoords + tour[j]]; // Finding Distance using the 1D matrix defined in main-mpi, tour[j] used to get index number              
                    if (distance < localMinDistance) {
                        //update the local min distance, local nearest vertex and local closest tour vertex
                        localMinDistance = distance;  
                        localNearestVertex = i;
                        localClosestTourVertex = j;                     
                    }
                }
            }
        }
        //critical section to update the global min distance, global nearest vertex and global closest tour vertex
        #pragma omp critical
        {
            if (localMinDistance < result.globalMinDistance) {
                result.globalMinDistance = localMinDistance;
                result.globalNearestVertex = localNearestVertex;
                result.globalClosestTourVertex = localClosestTourVertex;
            }
        }
    }
    //return the result struct
    return result;
}  
  
/* 
    iterates over the tour and finds the optimal insertion position for the vertex
    calculates the price of insertion at either side of the closest tour vertex, found by the NearestNotInTour function
    
    as the insertion is rather cheap,parallelisation is not necessary and causes an increase in processing time 
 */
void OptimalPositionInsertion(double *matrix1D, int vertex, int *tour, int size, int closestTourVertex, int numOfCoords) {
    //initialise variables
    double minIncrease = DBL_MAX;
    int optimalPosition = -1;

    /* 
    handling of an edge case in which the checking of the starting and ending vertex is made due to a situation
    where the starting and ending vertex were percieved as neighbours and were being compared as so. 
    */
    if (closestTourVertex == 0) {   
        int lastVertexAtEnd = tour[size - 2]; //finds the vertex to compare V[max] to
        int firstVertexAfterStart = tour[1]; // finds the vertex to compare V[0] to

        //manual calculation of the distance between v[max] and v[max - 1]
        double increaseEnd = matrix1D[lastVertexAtEnd * numOfCoords + vertex]
                            + matrix1D[vertex * numOfCoords + tour[0]] 
                            - matrix1D[lastVertexAtEnd * numOfCoords + tour[0]];
        
        //manual calculation of the distance between v[0] and v[1]
        double increaseStart = matrix1D[tour[0] * numOfCoords + vertex] 
                     + matrix1D[vertex * numOfCoords + firstVertexAfterStart] 
                     - matrix1D[tour[0] * numOfCoords + firstVertexAfterStart];

        // Compare and set the optimal position
        if (increaseEnd < minIncrease) {
            minIncrease = increaseEnd;
            optimalPosition = size - 1;
        }
        if (increaseStart < minIncrease) {
            minIncrease = increaseStart;
            optimalPosition = 0;
        }
    } else {
        //iterates over the two potential insertion points, these being the positions before and after the closest tour vertex
        for (int i = 0; i < 2; i++) {
            int pos = (closestTourVertex + i) % size; // ensures the cyclical nature of the tour by wrapping around to the start of the tour if it is greater than tour size
            int lastVertex = tour[pos - 1]; // last vertex is equal to the vertex before the current position 
            int nextVertex = tour[pos]; // next vertex is equal to the current position

            //calculates the increase in distance by inserting the vertex at the current position using the 1D array defined in main-mpi
            double increase = matrix1D[lastVertex * numOfCoords + vertex] 
                + matrix1D[vertex * numOfCoords + nextVertex] 
                - matrix1D[lastVertex * numOfCoords + nextVertex];

            //updates values if the increase is less than the current min increase
            if (increase < minIncrease) {
                minIncrease = increase;
                optimalPosition = pos;
            }
        }
    }

    //inserts the vertex at the optimal position
    if (optimalPosition >= 0) {
        for (int j = size; j > optimalPosition; j--) {
            tour[j] = tour[j - 1];
        }
        tour[optimalPosition] = vertex;
        
    }
} 


//nearest insertion function makes use of the NearestNotInTour and OptimalPositionInsertion functions to find the optimal tour
int *nearestInsertion(double *matrix1D, int numOfCoords, int startVertex) {
    //allocate memory for all arrays
    double (*matrix2D)[numOfCoords] = (double (*)[numOfCoords])matrix1D;
    int *visited = (int *)calloc(numOfCoords, sizeof(int));
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));
    printf("Starting nearestInsertion\n");
    //initalise variables
    double minDistance = DBL_MAX;
    int firstNearestVertex = -1;

   //finds the first nearest vertex to the start vertex
    for (int vertex = 0; vertex < numOfCoords; vertex++) {
        if (vertex != startVertex) {
            if (matrix1D[startVertex * numOfCoords + vertex] < minDistance) {
                minDistance = matrix1D[startVertex * numOfCoords + vertex];
                firstNearestVertex = vertex;
            }
        }
    }

    //creates the initial loop
    visited[startVertex] = 1;
    visited[firstNearestVertex] = 1;
    tour[0] = startVertex;
    tour[1] = firstNearestVertex;
    tour[2] = startVertex; 
    int tour_size = 3;


    /* while the tour is incomplete, complete the full nearest insertion process of identifying the nearest vertex not in the initial tour,
    and finding the optimal position to insert it. 
    */
    while (tour_size < numOfCoords + 1) {        
        NearestVertexResult nearestResult = NearestNotInTour(matrix1D, tour, numOfCoords,tour_size);
        int nextVertex = nearestResult.globalNearestVertex;
        int closestTourVertex = nearestResult.globalClosestTourVertex;
        double distance = matrix1D[nextVertex * numOfCoords + tour[closestTourVertex]];        
        OptimalPositionInsertion(matrix1D, nextVertex, tour, tour_size, closestTourVertex, numOfCoords);       
        visited[nextVertex] = 1;       
        tour_size++;
    }
    //return tour and free allocated memory
    free(visited);
    return tour;
} 


/* int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    char *outputFile = argv[2];

    double start, end;
    start = omp_get_wtime(); 

    
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    int *tour = nearestinsertion(distanceMatrix, numOfCoords);
    
    writeTourToFile(tour, numOfCoords + 1, outputFile);
    end = omp_get_wtime(); 
    printf("Work took %f seconds\n", end - start);
    
    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);
    free(tour);

    return 0;
}
 */