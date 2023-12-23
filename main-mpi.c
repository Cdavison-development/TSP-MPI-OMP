//libraries and header files
#include <mpi.h>
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
#define OMPNADDITION_H_
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
    #pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix2D[i][j] = EuclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix2D[j][i] = matrix2D[i][j]; // Use symmetry, avoid redundant calculation
        }
    }

} 


//define a struct to store the tour result as well as its cost and path
typedef struct {
    double cost;
    int *tour;
} TourResult;

typedef struct {
    double cost;
    int rank;
} CostRank;

//pointer to function for the insertion heuristic
typedef int* (*InsertionFunction)(double *, int, int);
   void findBestToursForEachHeuristic(double *matrix1D, int numOfCoords, int rank, int size, char *filenames[3]) {
    InsertionFunction insertionFunctions[3] = {cheapestInsertion, farthestInsertion, nearestInsertion};
    TourResult localBestTours[3];
    CostRank localCostRank[3], globalCostRank[3];
    int *allTours = NULL;  // Buffer to gather all tours
    if (rank == 0) {
        allTours = malloc((numOfCoords + 1) * size * 3 * sizeof(int));
        if (allTours == NULL) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
    //initalise array to store each heuristics tour
    for (int i = 0; i < 3; i++) {
        localBestTours[i].cost = DBL_MAX;
        localBestTours[i].tour = malloc((numOfCoords + 1) * sizeof(int));
        if (localBestTours[i].tour == NULL) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

   
    //find the first vertex index in nuum of coords
    int startVertex = rank * (numOfCoords / size) + (rank < (numOfCoords % size) ? rank : (numOfCoords % size));
    //find the final vertex index in coords
    int endVertex = startVertex + (numOfCoords / size) + (rank < (numOfCoords % size) ? 1 : 0);
   
    if (startVertex >= numOfCoords) {
        startVertex = endVertex = numOfCoords;
    }

    #pragma omp parallel
    {
    //iterate over each start vertex     
    for (int i = startVertex; i < endVertex; i++) {
        #pragma omp for 
        for (int hrs = 0; hrs < 3; hrs++) {
            //initialise current tour for each insertion heuristic
            int *currentTour = insertionFunctions[hrs](matrix1D, numOfCoords, i);
            if (currentTour == NULL) {
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }                    
            double currentCost = 0.0;
            //calculate distance of current tour
            #pragma omp parallel for reduction(+:currentCost)
            for (int j = 0; j < numOfCoords; j++) {
                double segmentCost = matrix1D[currentTour[j] * numOfCoords + currentTour[(j + 1) % numOfCoords]];
                currentCost += segmentCost;           
            }
            //update tour if new tour cost is cheaper than current tour cost in thread safe environment
            #pragma omp critical
            {
                if (currentCost < localBestTours[hrs].cost) {
                    memcpy(localBestTours[hrs].tour, currentTour, (numOfCoords + 1) * sizeof(int));
                    localBestTours[hrs].cost = currentCost;
                } 
          
            }
            free(currentTour);
        }
            
        }
    }


    // Store local best tour costs and ranks from the struct for each heuristic
    for (int hrs = 0; hrs < 3; hrs++) {
        localCostRank[hrs].cost = localBestTours[hrs].cost;
        localCostRank[hrs].rank = rank;
        //gathers the values of local best tour from all processes
        MPI_Gather(localBestTours[hrs].tour, numOfCoords + 1, MPI_INT, allTours + hrs * (numOfCoords + 1) * size, numOfCoords + 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    //aggregates the total results from every process and finds the best tours across all heuristics and processes
    MPI_Allreduce(localCostRank, globalCostRank, 3, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    //writes each tour to file 
    if(rank == 0){
        for (int hrs = 0; hrs < 3; hrs++) {
            printf("Global best cost for heuristic %d is %f found by process %d.\n", hrs, globalCostRank[hrs].cost, globalCostRank[hrs].rank);
            printf("Global best tour for heuristic %d: ", hrs);
            for (int j = 0; j < numOfCoords + 1; j++) {
                printf("%d ", allTours[hrs * (numOfCoords + 1) * size + globalCostRank[hrs].rank * (numOfCoords + 1) + j]);
                // Call the function to write the tour to file
                char* outputFilename = filenames[hrs];
                writeTourToFile(&allTours[hrs * (numOfCoords + 1) * size + globalCostRank[hrs].rank * (numOfCoords + 1)], numOfCoords + 1, outputFilename);
            }       
        }
   }
    //free allTours on rank 0
    if (rank == 0) {
        if (allTours != NULL) {
            free(allTours);
        } else {
            fprintf(stderr, "all Tours null %d", rank);
        }
    }
    

    // Free localBestTours in all processes
    for (int i = 0; i < 3; i++) {
        if (localBestTours[i].tour != NULL) {
            free(localBestTours[i].tour);
        } else {
            fprintf(stderr, "local best tours null %d\n", rank);
        }
    }

}
   
int main(int argc, char *argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check the correct number of arguments
    if (argc != 5) {
        if (rank == 0) {
            printf("USE: %s <coord_file> <cheapest insertion output file> <furthest insertion output file> <nearest insertion output file>\n", argv[0]);
        }
        
        return 1;
    }

    

    // Read the number of coordinates and the coordinates themselves
    char *inputFile = argv[1];
    int numOfCoords = readNumOfCoords(inputFile);  // Implement this function as needed
    double **coords = readCoords(inputFile, numOfCoords);  // Implement this function as needed
    double startTime, endTime;
    if (rank == 0) {
        startTime = MPI_Wtime();
    }
    // Allocate memory for the distance matrix
    double *matrix1D = malloc(numOfCoords * numOfCoords * sizeof(double));
    if (matrix1D == NULL) {
        fprintf(stderr, "Failed to allocate memory\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Generate and broadcast the distance matrix (only in process 0)
    if (rank == 0) {
        GenerateDistanceMatrix(matrix1D, coords, numOfCoords);  // Implement this function as needed
    }
    //Bcast the distance matrix from root to all other processes
    MPI_Bcast(matrix1D, numOfCoords * numOfCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Store the output file names
    char *filenames[3] = {argv[2], argv[3], argv[4]};

    // Call the function to find and process the best tours
    findBestToursForEachHeuristic(matrix1D, numOfCoords, rank, size, filenames);

     if (rank == 0) {
        endTime = MPI_Wtime();
        printf("Total execution time: %f seconds\n", endTime - startTime);
    }

    // Clean up and finalize MPI
    free(matrix1D);
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
    }
    free(coords);

    MPI_Finalize();
    return 0;
}
