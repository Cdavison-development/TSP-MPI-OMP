//#include "ompcInsertion.c"
//#include "ompfInsertion.c"
//#include "ompnAddition.c"
#include <mpi.h>
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



 double euclideanDistance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return sqrt(dx * dx + dy * dy);
}

/* double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double *matrix1D = malloc(numOfCoords * numOfCoords * sizeof(double));

    double (*matrix2D)[numOfCoords] = (double (*)[numOfCoords])matrix1D;
    
    for (int i = 0; i < numOfCoords; i++) {
        matrix2D[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix2D[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix2D[j][i] = matrix2D[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return (double **)matrix2D;
}  */

void generateDistanceMatrix(double *matrix1D,double **coords, int numOfCoords) {
     printf("generateDistanceMatrix called with numOfCoords = %d\n", numOfCoords);
    double (*matrix2D)[numOfCoords] = (double (*)[numOfCoords])matrix1D;

    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            matrix2D[i][j] = 0.0; // Initialize with -1 for non-diagonal elements
        }
    }
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix2D[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix2D[j][i] = matrix2D[i][j]; // Use symmetry, avoid redundant calculation

        //Print statement for debugging: show each distance calculated
           //printf("Distance between %d and %d: %f\n", i, j, matrix2D[i][j]);
        }
    }
   printf("Completed distance matrix:\n");
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            printf("%f ", matrix2D[i][j]);
        }
        printf("\n");
    }
   
} 

double calculateTotalDistance(int *tour, double *matrix1D, int tourLength,  int numOfCoords) {
    double totalDistance = 0.0;
    for (int i = 0; i < tourLength - 1; i++) {
        totalDistance += matrix1D[tour[i] * numOfCoords + tour[i + 1]];
    }
   totalDistance += matrix1D[tour[tourLength - 1] * numOfCoords + tour[0]]; // Complete the loop
    printf("Total distance: %f\n", totalDistance);
    return totalDistance;
}


/* typedef int* (*InsertionFunction)(double *, int, int);
void findBestTour(double *matrix1D, int numOfCoords, int hrs,char *filenames[], int startVertex) {
    InsertionFunction insertionFunctions[3] = {cheapestInsertion, farthestInsertion, nearestInsertion};
    
    
    double globalBestCost = DBL_MAX;
    int *globalBestTour = malloc((numOfCoords + 1) * sizeof(int));
    int tourLength = numOfCoords + 1; 
    printf("findBestTour called with hrs = %d\n", hrs);
    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        int *currentTour = insertionFunctions[hrs](matrix1D, numOfCoords, startVertex);
        double currentCost = calculateTotalDistance(currentTour, matrix1D, tourLength, numOfCoords);


        if (currentCost < globalBestCost) {
            globalBestCost = currentCost;
            memcpy(globalBestTour, currentTour, (numOfCoords + 1) * sizeof(int));
        }

        free(currentTour);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char *filename = filenames[hrs];
    if (rank == 0) {
        char *filename = filenames[hrs];
        writeTourToFile(globalBestTour, numOfCoords + 1, filename);
    }

    free(globalBestTour);
}


int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Usage: %s <coord_file> <cheapest output file> <furthest output file> <nearest output file>\n", argv[0]);
        return 1;
    }
    
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double start_time = MPI_Wtime();
    char *inputFile = argv[1];
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    double *matrix1D = malloc(numOfCoords * numOfCoords * sizeof(double));
    int verticesPerProcess = numOfCoords / size;
    int startVertex = rank * verticesPerProcess;
    int endVertex = (rank == size - 1) ? numOfCoords : (rank + 1) * verticesPerProcess;

    if (matrix1D == NULL) {
        fprintf(stderr, "Failed to allocate memory for matrix1D\n");
        exit(EXIT_FAILURE);
    }

    if (rank == 0) {
        // Generate the distance matrix only in the root process
        generateDistanceMatrix(matrix1D, coords, numOfCoords);
    }

    // Broadcast the distance matrix to all processes
    int bcast_result = MPI_Bcast(matrix1D, numOfCoords * numOfCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (bcast_result != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Bcast failed.\n");
        MPI_Abort(MPI_COMM_WORLD, bcast_result);
        // Exit or handle the error as necessary
    }   
    char *filenames[3] = {argv[2], argv[3], argv[4]};

    
    for (int v = startVertex; v < endVertex; v++) {
        for (int alg = 0; alg < 3; alg++) {
            findBestTour(matrix1D, numOfCoords, alg, filenames, v);  // Pass the starting vertex
        }
    }
    double end_time = MPI_Wtime();
    double duration = end_time - start_time;

    if (rank == 0) {
        // Print the duration (only in the root process)
        printf("Time taken: %f seconds\n", duration);
    }
    
    free(matrix1D); // Free the memory

    // Free other dynamically allocated memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
    }
    free(coords);

    MPI_Finalize(); // Finalize MPI
    return 0;
} 
*/
typedef struct {
    double cost;
    int *tour;
} TourResult;
typedef int* (*InsertionFunction)(double *, int, int);
TourResult findBestTour(double *matrix1D, int numOfCoords, int alg, int startVertex) {
    InsertionFunction insertionFunctions[3] = {cheapestInsertion, farthestInsertion, nearestInsertion};

    double globalBestCost = DBL_MAX;
    int *globalBestTour = malloc((numOfCoords + 1) * sizeof(int));
    if (globalBestTour == NULL) {
        fprintf(stderr, "Failed to allocate memory for globalBestTour\n");
        // Handle allocation failure
    } else {
        // Optionally initialize and print initial values
    for (int i = 0; i <= numOfCoords; i++) {
        globalBestTour[i] = -1;  // Initialize with a placeholder value
       // printf("%d ", globalBestTour[i]);  // Print initial value
    }
    printf("\n");
}

    TourResult result = { .cost = DBL_MAX, .tour = NULL };
    int *currentTour = insertionFunctions[alg](matrix1D, numOfCoords, startVertex);
    if (currentTour == NULL) {
        fprintf(stderr, "Algorithm %d failed at start vertex %d\n", alg, startVertex);
        // Handle the error, e.g., skip this iteration or abort
    }
    double currentCost = calculateTotalDistance(currentTour, matrix1D, numOfCoords + 1, numOfCoords);
    if (currentCost < globalBestCost) {
        printf("currentCost: %f\n", currentCost);
        globalBestCost = currentCost;
        memcpy(globalBestTour, currentTour, (numOfCoords + 1) * sizeof(int));
        printf("New global best tour found with cost %f: ", globalBestCost);
        for (int i = 0; i <= numOfCoords; i++) {
            printf("%d ", globalBestTour[i]);
        }
        printf("\n");
    }
    free(currentTour);
    result.cost = globalBestCost;
    result.tour = globalBestTour;
    return result;
}

typedef struct {
    double cost;
    int rank;
} CostRank;

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Usage: %s <coord_file> <cheapest output file> <furthest output file> <nearest output file>\n", argv[0]);
        return 1;
    }

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double start_time = MPI_Wtime();

    char *inputFile = argv[1];
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);

    double *matrix1D = malloc(numOfCoords * numOfCoords * sizeof(double));
    if (matrix1D == NULL) {
        fprintf(stderr, "Failed to allocate memory for matrix1D\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (rank == 0) {
        generateDistanceMatrix(matrix1D, coords, numOfCoords);
    }

    MPI_Bcast(matrix1D, numOfCoords * numOfCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    char *filenames[3] = {argv[2], argv[3], argv[4]};
    TourResult localBestTours[3];
    for (int i = 0; i < 3; i++) {
        localBestTours[i].cost = DBL_MAX;
        localBestTours[i].tour = malloc((numOfCoords + 1) * sizeof(int));
    }

    int verticesPerProcess = numOfCoords / size;
    int remainder = numOfCoords % size;
    int startVertex = rank * verticesPerProcess + (rank < remainder ? rank : remainder);
    int endVertex = startVertex + verticesPerProcess + (rank < remainder ? 1 : 0);

    if (startVertex >= numOfCoords) {
        startVertex = endVertex = numOfCoords;
    }

    for (int v = startVertex; v < endVertex; v++) {
        for (int alg = 0; alg < 3; alg++) {
            TourResult currentResult = findBestTour(matrix1D, numOfCoords, alg, v);
            if (currentResult.cost < localBestTours[alg].cost) {
                memcpy(localBestTours[alg].tour, currentResult.tour, (numOfCoords + 1) * sizeof(int));
                localBestTours[alg].cost = currentResult.cost;
            }
            free(currentResult.tour);
        }
    }

    CostRank localCostRank[3], globalCostRank[3];
    for (int alg = 0; alg < 3; alg++) {
        localCostRank[alg].cost = localBestTours[alg].cost;
        localCostRank[alg].rank = rank;
    }

    MPI_Allreduce(localCostRank, globalCostRank, 3, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    int *allTours = NULL;
    if (rank == 0) {
        allTours = malloc(numOfCoords * size * 3 * sizeof(int));
    }

    for (int alg = 0; alg < 3; alg++) {
        MPI_Gather(localBestTours[alg].tour, numOfCoords, MPI_INT, allTours + alg * numOfCoords * size, numOfCoords, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* if (rank == 0) {
        for (int alg = 0; alg < 3; alg++) {
            printf("Global best tour for algorithm %d with cost %f: ", alg, globalCostRank[alg].cost);
            int *finalBestTour = allTours + alg * numOfCoords * size + globalCostRank[alg].rank * numOfCoords;
            for (int i = 0; i < numOfCoords + 1; i++) { // Add +1 to include the return to the initial city
                printf("%d ", finalBestTour[i % numOfCoords]); // Use modulo to loop back to the first city
            }
            printf("\n");
            writeTourToFile(finalBestTour, numOfCoords + 1, filenames[alg]); // Pass numOfCoords + 1 to include the full loop
        }
        free(allTours);
    } */
    if (rank == 0) {
    for (int alg = 0; alg < 3; alg++) {
        printf("Global best tour for algorithm %d with cost %f: ", alg, globalCostRank[alg].cost);
        int *finalBestTour = allTours + alg * numOfCoords * size + globalCostRank[alg].rank * numOfCoords;
        for (int i = 0; i < numOfCoords + 1; i++) { // Add +1 to include the return to the initial city
            printf("%d ", finalBestTour[i % numOfCoords]); // Use modulo to loop back to the first city
        }
        // Create a temporary array for the complete tour
        int *completeTour = malloc((numOfCoords + 1) * sizeof(int));
        memcpy(completeTour, finalBestTour, numOfCoords * sizeof(int));
        completeTour[numOfCoords] = completeTour[0]; // Add the starting city at the end

        // Now write using the complete tour
        writeTourToFile(completeTour, numOfCoords + 1, filenames[alg]);

        // Free the temporary array
        free(completeTour);

        // Continue with the rest of the code...
    }
    free(allTours);
}
    double end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Time taken: %f seconds\n", end_time - start_time);
    }

    free(matrix1D);
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
    }
    free(coords);

    MPI_Finalize();
    return 0;
}  
    

   