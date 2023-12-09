
//import libraries and coordReader Files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords();
void *writeTourToFile();
#endif
//Provided functions
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);


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

    //generate empty distance matrix of size numOfCoords
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix[i][j] = distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    printf("Distance Matrix:\n");
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    return matrix;
}

    //Find the nearest neighbor for a given node in distance matrix.
    void findNearestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *nearestIndex, double *nearestDistance) {
    double nearestDist = DBL_MAX; //Initialise with max possible value 
    int nearestIdx = -1;          //Initialise with invalid index

    //Iterate over all coordinates to find nearest neighbor
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode) {                             //Skip current node
            double dist = distanceMatrix[currentNode][i];   //Get distance from current node to i
            if (dist < nearestDist) {                       //Update nearest neighbor if a closer node is found
                nearestDist = dist;
                nearestIdx = i;
            }
        }
    }

    *nearestDistance = nearestDist; // Update the value pointed by nearestDistance
    *nearestIndex = nearestIdx;     // Update the value pointed by nearestIndex
}
 


//Cheapest insertion algorithm for TSP problem
int *cheapestInsertion(double **distanceMatrix, int numOfCoords) {

    int *tour = malloc((numOfCoords + 1) * sizeof(int)); //Allocate tour memory
    int *unvisited = malloc(numOfCoords * sizeof(int));  //Allocate unisited memory for list of unvisited nodes

    // Initialize unvisited nodes starting from 1
    for (int i = 1; i < numOfCoords; i++) {
        unvisited[i - 1] = i;
    }
    int unvisitedCount = numOfCoords - 1;

    // Start with vertex 0
    tour[0] = 0;
    int tourSize = 1; 

    // Find the nearest neighbor to 0 and add to the tour
    double nearestDistance;
    int nearestIndex;
    findNearestNeighbor(distanceMatrix, numOfCoords, 0, &nearestIndex, &nearestDistance); 
   
    tour[1] = nearestIndex; //Add nearest neighbor to tour
    tour[2] = 0;            //Complete the inital loop back to the starting node
    tourSize = 3; 

    //NOTE: Following section is included twice as the tour initialises first node twice, i.e. instead of 0,4,0 returns 0,4,4,0
    //This fix removes the first duplicate.

    // Remove the nearest neighbor from the unvisited list
    for (int i = 0; i < unvisitedCount; i++) {
        if (unvisited[i] == nearestIndex) {
            unvisited[i] = unvisited[unvisitedCount - 1]; //Replace found node
            unvisitedCount--; //Decrease count of Unvisited nodes
            break;
        }
    }
    
    //Tour construction until all nodes are visited
    while (tourSize < numOfCoords + 1) {
        double minCost = DBL_MAX;
        int minCostIndex = -1, insertPosition = -1;

        //Find cheapest node and insert into current tour
        for (int idx = 0; idx < unvisitedCount; idx++) {
            int i = unvisited[idx];
            for (int j = 0; j < tourSize - 1; j++) {
                // Calculate cost of inserting node between tour[j] and tour[j + 1] 
                double cost = distanceMatrix[tour[j]][i] + distanceMatrix[i][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                if (cost < minCost) {   //Update minimum cost and corresponding node and insertion position
                    minCost = cost;
                    minCostIndex = i;
                    insertPosition = j + 1;
                }
            }
        }

        //Insert the cheapest node into the tour
        if (minCostIndex != -1) {
            for (int i = tourSize; i > insertPosition; i--) {
                tour[i] = tour[i - 1]; //Shift elements to make space for new node
            }
            tour[insertPosition] = minCostIndex; //Insert node

            // Remove inserted node from list of unvisited nodes
            for (int i = 0; i < unvisitedCount; i++) {
                if (unvisited[i] == minCostIndex) {
                    unvisited[i] = unvisited[unvisitedCount - 1];
                    unvisitedCount--;
                    break;
                }
            }

            tourSize++; //Increment the size of the tour
        }
    }
    
    free(unvisited); //Free memory allocated for unvisited
    return tour;
}

int main(int argc, char *argv[]) {
    // Ensure correct usage
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    //Initialising arguments
    char *inputFile = argv[1];
    char *outputFile = argv[2];


    
    //Read coordinates
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    
    //Generate distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    //Apply the cheapest insertion algorithm
    int *tour = cheapestInsertion(distanceMatrix, numOfCoords);

    //Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);
    
    //Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
       
    free(coords);
    free(distanceMatrix);
    free(tour);

    return 0;
}  