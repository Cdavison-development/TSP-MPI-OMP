// libraries and header files

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <float.h>
#include <string.h>


#ifndef OMPCINSERTION_H_
#define OMPCINSERTION_H_
    int *cheapestInsertion(double *matrix1D, int numOfCoords, int startVertex);
#endif
#ifndef OMPFINSERTION_H_
#define OMPFINSERTION_H_
    int *farthestInsertion(double *matrix1D, int numOfCoords, int startVertex);
#endif
#ifndef OMPNADDITION_H_
#define OOMPNADDITION_H_
    int *nearestInsertion(double *matrix1D, int numOfCoords, int startVertex);
#endif
#ifndef coordReader_H_
#define coordReader_H_
    int readNumOfCoords(char *filename);
    double **readCoords(char *filename, int numOfCoords);           
    void writeTourToFile(int *tour, int tourLength, char *filename);
#endif

//calculates the euclidean distance between two points
double EuclideanDistance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return sqrt(dx * dx + dy * dy);
}
//generates the distance matrix using a 1D array to ensure contiguous memory
void GenerateDistanceMatrix(double *matrix1D,double **coords, int numOfCoords) {
    double (*matrix2D)[numOfCoords] = (double (*)[numOfCoords])matrix1D;
    //initialize vertices where i == j to 0
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            matrix2D[i][j] = 0.0; 
        }
    }
    //calculate the distance between each pair of vertices
    
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix2D[i][j] = EuclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix2D[j][i] = matrix2D[i][j]; // Use symmetry, avoid redundant calculation
        }
    }

}  
double calculateTotalDistance(int *tour, double *matrix1D, int tourLength,  int numOfCoords) {
    double totalDistance = 0.0; //initialize distance variable
    //loop through the all vertexes in the tour other than the last vertex
    
    for (int i = 0; i < tourLength - 1; i++) {
        totalDistance += matrix1D[tour[i] * numOfCoords + tour[i + 1]]; //find value of i in tour and add to total distance
    }
   totalDistance += matrix1D[tour[tourLength - 1] * numOfCoords + tour[0]]; //add the final vertex in the tour to the starting vertex
    return totalDistance; //return the total distance
} 

//pointer to function for the insertion algorithms
typedef int* (*InsertionFunction)(double *, int, int);

void FindBestTour(double *matrix1D, int numOfCoords, int hrs,char *filenames[]) {
    //array of function points for each insertion algorithm
    InsertionFunction insertionFunctions[3] = {cheapestInsertion, farthestInsertion, nearestInsertion};
    //initialise the global variables for the best tour and its cost
    double globalBestCost = DBL_MAX;
    //allocate memory for the global best tour
    int *globalBestTour = malloc((numOfCoords + 1) * sizeof(int));
    
    //iterate over the number of coordinates and find the best tour for each insertion algorithm
    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        //calculate current tour using the current insertion algorithm
        int *currentTour = insertionFunctions[hrs](matrix1D, numOfCoords, startVertex);
        //calculate the distance of the tour using the calculateTotalDistance function
        double currentCost = calculateTotalDistance(currentTour, matrix1D, numOfCoords+1, numOfCoords);

        //update global values if current values have lower costs
        if (currentCost < globalBestCost) {
            globalBestCost = currentCost;
            memcpy(globalBestTour, currentTour, (numOfCoords + 1) * sizeof(int));
        }
        //free allocated memory for current tour
        free(currentTour);
    }
    //write the best tour to the output file
    char *filename = filenames[hrs];
    writeTourToFile(globalBestTour, numOfCoords + 1, filename);
    //free allocated memory for global best tour
    free(globalBestTour);
}


int main(int argc, char *argv[]) {
    
    //check if the correct number of arguments are used
    if (argc != 5) {
        printf("USE : %s <coord_file> <cheapest insertion output file> <furthest insertion output file> <nearest insertion output file>\n", argv[0]);
        return 1; //exit if incorrect amount of arguments are used
    }

    //Read num of coords and coords from the input file
    char *inputFile = argv[1];    
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    //generates dfistance matrix using the coords
    double *matrix1D = malloc(numOfCoords * numOfCoords * sizeof(double));
     GenerateDistanceMatrix(matrix1D, coords, numOfCoords);

    //set up filenames for output files
    char *filenames[3] = {argv[2], argv[3], argv[4]};

    //find the best tour for each insertion algorithm
    for (int hrs = 0; hrs < 3; hrs++) {
        FindBestTour(matrix1D, numOfCoords, hrs, filenames);
    }

    return 0; //exit program
}

    

   