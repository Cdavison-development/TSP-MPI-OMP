#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "utils.h"
#include "coordReader.h"
#include <omp.h>
/* #ifndef COORDREADER_H
#define COORDREADER_H
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);
#endif  */



/* double euclideanDistance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return sqrt(dx * dx + dy * dy);
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
            matrix[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return matrix;
} 
 */
int IsVertexInTour(int vertex, int *tour, int tourSize) {
    for (int i = 0; i < tourSize; i++) {       
        if (tour[i] == vertex) {           
            return 1;
        }
    }    
    return 0;
} 


  typedef struct {
    int global_nearest_vertex;
    int global_closest_tour_vertex;
    double global_min_distance;
    //int local_nearest_vertex;
    //int local_closest_tour_vertex;
    //double local_min_distance;
} NearestVertexResult;


 NearestVertexResult NearestNotInTour(double **distanceMatrix, int *tour, int numOfCoords, int tourSize) {
    NearestVertexResult result;
    result.global_nearest_vertex = -1;
    result.global_closest_tour_vertex = -1;
    result.global_min_distance = DBL_MAX;

    #pragma omp parallel
    {
        int local_nearest_vertex = -1;
        int local_closest_tour_vertex = -1;
        double local_min_distance = DBL_MAX;  

        #pragma omp for
        for (int i = 0; i < numOfCoords; i++) {
            if (!IsVertexInTour(i, tour, tourSize)) {
                //printf("Vertex %d is not in the tour.\n", i);
                for (int j = 0; j < tourSize; j++) { // Iterate over current tour size
                    double distance = distanceMatrix[i][tour[j]]; // Use tour[j] to get the vertex number
                  //printf("Checking distance from vertex %d (not in tour) to vertex %d (in tour): %f\n", i, tour[j], distance);
                    if (distance < local_min_distance) {
                        local_min_distance = distance;
                        local_nearest_vertex = i;
                        local_closest_tour_vertex = j; // j is the index in the tour
                    //printf("New nearest vertex: %d with distance: %f, Closest tour vertex (index): %d\n", i, distance, j);
                    }
                }
            }
        }
        #pragma omp critical
        {
            if (local_min_distance < result.global_min_distance) {
                result.global_min_distance = local_min_distance;
                result.global_nearest_vertex = local_nearest_vertex;
                result.global_closest_tour_vertex = local_closest_tour_vertex;
            }
        }
    }
    return result;
}  
  

void OptimalPositionInsertion(double **distanceMatrix, int vertex, int *tour, int size, int closest_tour_vertex, int numOfCoords) {
    double minIncrease = DBL_MAX;
    int OptimalPosition = -1;

    //printf("Closest tour vertex index: %d, Vertex number: %d\n", closest_tour_vertex, tour[closest_tour_vertex]);

    if (closest_tour_vertex == 0) {
        // Special handling when the closest tour vertex is the start of the tour
        int lastVertexAtEnd = tour[size - 2]; // Second-to-last vertex (339)
        int firstVertexAfterStart = tour[1]; // First vertex after start (248)

        // Calculate increase for inserting between last vertex and start
        double increaseEnd = distanceMatrix[lastVertexAtEnd][vertex] + distanceMatrix[vertex][tour[0]] - distanceMatrix[lastVertexAtEnd][tour[0]];
        //printf("Checking position between vertices %d and %d, Increase: %f\n", lastVertexAtEnd, tour[0], increaseEnd);

        // Calculate increase for inserting between start and first vertex after start
        double increaseStart = distanceMatrix[tour[0]][vertex] + distanceMatrix[vertex][firstVertexAfterStart] - distanceMatrix[tour[0]][firstVertexAfterStart];
        //printf("Checking position between vertices %d and %d, Increase: %f\n", tour[0], firstVertexAfterStart, increaseStart);

        // Compare and set the optimal position
        if (increaseEnd < minIncrease) {
            minIncrease = increaseEnd;
            OptimalPosition = size - 1;
        }
        if (increaseStart < minIncrease) {
            minIncrease = increaseStart;
            OptimalPosition = 0;
        }
    } else {
        // Regular case
        for (int i = 0; i < 2; i++) {
            int pos = (closest_tour_vertex + i) % size;
            int lastVertex = (pos == 0) ? tour[size - 2] : tour[pos - 1];
            int nextVertex = tour[pos];

            double increase = distanceMatrix[lastVertex][vertex] + distanceMatrix[vertex][nextVertex] - distanceMatrix[lastVertex][nextVertex];
            //printf("Checking position between vertices %d and %d, Increase: %f\n", lastVertex, nextVertex, increase);

            if (increase < minIncrease) {
                minIncrease = increase;
                OptimalPosition = pos;
            }
        }
    }

    // Insert vertex at the optimal position
    if (OptimalPosition >= 0) {
        for (int j = size; j > OptimalPosition; j--) {
            tour[j] = tour[j - 1];
        }
        tour[OptimalPosition] = vertex;
        //printf("Inserted vertex %d at position %d\n", vertex, OptimalPosition);
    } else {
        //printf("No optimal position found for vertex %d\n", vertex);
    }
} 



int *nearestInsertion(double **distanceMatrix, int numOfCoords, int start_vertex) {
    int *visited = (int *)calloc(numOfCoords, sizeof(int));
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));

    double min_distance = DBL_MAX;
    int first_nearest_vertex = -1;

   
    for (int vertex = 0; vertex < numOfCoords; vertex++) {
        if (vertex != start_vertex && distanceMatrix[start_vertex][vertex] < min_distance) {
            min_distance = distanceMatrix[start_vertex][vertex];
            first_nearest_vertex = vertex;
        }
    }

    
    visited[start_vertex] = 1;
    visited[first_nearest_vertex] = 1;
    tour[0] = start_vertex;
    tour[1] = first_nearest_vertex;
    tour[2] = start_vertex; // Completing the initial loop
    int tour_size = 3;

    /* printf("Initial tour: ");
    for (int i = 0; i < tour_size; i++) {
    printf("%d ", tour[i]);
}
printf("\n"); */

    //printf("Starting vertex: %d\n", start_vertex);
    //printf("First nearest vertex to start: %d\n", first_nearest_vertex);
    
    while (tour_size < numOfCoords + 1) {
        //printf("Iteration: %d\n", k);
        //printf("numOfCoords: %d\n",numOfCoords);
       // printf("tour_size: %d\n",tour_size);
        NearestVertexResult nearestResult = NearestNotInTour(distanceMatrix, tour, numOfCoords,tour_size);
        int next_vertex = nearestResult.global_nearest_vertex;
        int closest_tour_vertex = nearestResult.global_closest_tour_vertex;
        double distance = distanceMatrix[next_vertex][tour[closest_tour_vertex]];
        
         //printf("Next nearest vertex not in tour: %d (Distance to closest tour vertex %d: %f)\n", 
                //next_vertex, tour[closest_tour_vertex], distance);
        

        //printf("Next nearest vertex not in tour: %d\n", next_vertex);
        OptimalPositionInsertion(distanceMatrix, next_vertex, tour, tour_size, closest_tour_vertex, numOfCoords);
        //printf("Inserted vertex %d at position %d in the tour\n", next_vertex, tour_size);
        visited[next_vertex] = 1;
        //printf("visited[next_vertex]: %d\n",visited[next_vertex]);
        tour_size++;
        //printf("tour_size: %d\n",tour_size);

    }

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