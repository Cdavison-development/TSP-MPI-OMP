#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <string.h>
#ifndef COORDREADER_H
#define COORDREADER_H
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);
#endif 



double euclideanDistance(double x1, double y1, double x2, double y2) {
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
            matrix[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
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
int IsVertexInTour(int vertex, int *tour, int numOfCoords) {
    for (int i = 0; i < numOfCoords; i++) {
        if (tour[i] == vertex) {
            return 1;
        }
    }
    return 0;
}

typedef struct {

    double min_distance;
    double distance;
    int nearest_vertex;

}   NearestVertex;








typedef struct {
    int nearest_vertex;
    int closest_tour_vertex;
} NearestVertexResult;



NearestVertexResult NearestNotInTour(double **distanceMatrix, int *tour, int numOfCoords) {
    NearestVertexResult result;
    result.nearest_vertex = -1;
    result.closest_tour_vertex = -1;
    double min_distance = DBL_MAX;

    for (int i = 0; i < numOfCoords; i++) {
        if (!IsVertexInTour(i, tour, numOfCoords)) {
            for (int j = 0; j < numOfCoords; j++) {
                if (IsVertexInTour(j, tour, numOfCoords)) {
                    double distance = distanceMatrix[i][j];
                    if (distance < min_distance) {
                        min_distance = distance;
                        result.nearest_vertex = i;
                        // Find the index of vertex 'j' in the tour array
                        for (int k = 0; k < numOfCoords; k++) {
                            if (tour[k] == j) {
                                result.closest_tour_vertex = k;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}
/*
void OptimalPositionInsertion(double **distanceMatrix, int vertex, int *tour, int size, int closest_tour_vertex, int numOfCoords) {
    double minIncrease = DBL_MAX;
    int OptimalPosition = -1;

    // Define the two positions to check: before and after the closest tour vertex
    int positions[2];

    // Check if the closest tour vertex is the first in the tour
    positions[0] = (closest_tour_vertex == 0) ? size - 1 : closest_tour_vertex - 1;
    positions[1] = closest_tour_vertex;

    printf("Inserting vertex: %d\n", vertex);
    printf("Closest tour vertex index: %d, Vertex number: %d\n", closest_tour_vertex, tour[closest_tour_vertex]);

    for (int i = 0; i < 2; i++) {
        int pos = positions[i];
        int nextVertex = (pos == size - 1) ? tour[0] : tour[pos + 1];
        //int nextVertex = tour[(pos + 1) % size];
        int lastVertex = tour[pos];

        double increase = distanceMatrix[lastVertex][vertex] + distanceMatrix[vertex][nextVertex] - distanceMatrix[lastVertex][nextVertex];
        printf("Checking position between vertices %d and %d (index %d), Increase: %f\n", lastVertex, nextVertex, pos, increase);

        // Check if the increase is less than or equal to the minimum found so far
        // If equal, prefer the position before the closest tour vertex
        if (increase < minIncrease || (increase == minIncrease && i == 0)) {
            minIncrease = increase;
            OptimalPosition = pos + (i == 0 ? 0 : 1);
            printf("Current optimal position: %d (between index %d and %d), Increase: %f\n", OptimalPosition, pos, (pos + 1) % size, minIncrease);
        }
    }

    // Insert vertex
    if (OptimalPosition >= 0) {
        for (int j = size; j > OptimalPosition; j--) {
            tour[j] = tour[j - 1];
        }
        tour[OptimalPosition] = vertex;
        printf("Inserted vertex %d at position %d (index %d)\n", vertex, OptimalPosition, OptimalPosition);
    }else {
        printf("No optimal position found for vertex %d\n", vertex);
    }

}
*/
void OptimalPositionInsertion(double **distanceMatrix, int vertex, int *tour, int size, int closest_tour_vertex, int numOfCoords) {
    double minIncrease = DBL_MAX;
    int OptimalPosition = -1;

    // Positions to check: just before and after the closest tour vertex
    int positions[2] = {
        closest_tour_vertex, // Position before the closest tour vertex
        (closest_tour_vertex + 1) % size // Position after the closest tour vertex (with wrap-around)
    };

    printf("Inserting vertex: %d\n", vertex);
    printf("Closest tour vertex index: %d, Vertex number: %d\n", closest_tour_vertex, tour[closest_tour_vertex]);

    for (int i = 0; i < 2; i++) {
        int pos = positions[i];
        int nextVertex = (pos == size) ? tour[0] : tour[pos];
        int lastVertex = (pos == 0) ? tour[size - 1] : tour[pos - 1];

        double increase = distanceMatrix[lastVertex][vertex] + distanceMatrix[vertex][nextVertex] - distanceMatrix[lastVertex][nextVertex];
        printf("Checking position between vertices %d and %d, Increase: %f\n", lastVertex, nextVertex, increase);

        if (increase < minIncrease) {
            minIncrease = increase;
            OptimalPosition = pos;
            printf("Current optimal position: %d, Increase: %f\n", OptimalPosition, minIncrease);
        }
    }

    // Adjust the position for insertion if needed
    if (OptimalPosition == size) {
        OptimalPosition = 0;
    }

    // Insert vertex
    if (OptimalPosition >= 0) {
        for (int j = size; j > OptimalPosition; j--) {
            tour[j] = tour[j - 1];
        }
        tour[OptimalPosition] = vertex;
        printf("Inserted vertex %d at position %d\n", vertex, OptimalPosition);
    } else {
        printf("No optimal position found for vertex %d\n", vertex);
    }
}

// Update nearestinsertion function and main accordingly to use the new NearestNotInTour return type


int *nearestinsertion(double **distanceMatrix, int numOfCoords) {
    int *visited = (int *)calloc(numOfCoords, sizeof(int));
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));

    int start_vertex = 1; // Start from vertex 0
    double min_distance = DBL_MAX;
    int first_nearest_vertex = -1;

    // Find the nearest vertex to the start vertex
    for (int vertex = 0; vertex < numOfCoords; vertex++) {
        if (vertex != start_vertex && distanceMatrix[start_vertex][vertex] < min_distance) {
            min_distance = distanceMatrix[start_vertex][vertex];
            first_nearest_vertex = vertex;
        }
    }

    // Initialize the tour with the start vertex and its nearest vertex
    visited[start_vertex] = 1;
    visited[first_nearest_vertex] = 1;
    tour[0] = start_vertex;
    tour[1] = first_nearest_vertex;
    tour[2] = start_vertex; // Completing the initial loop
    int tour_size = 3;

    printf("Initial tour: ");
    for (int i = 0; i < tour_size; i++) {
    printf("%d ", tour[i]);
}
printf("\n");

    printf("Starting vertex: %d\n", start_vertex);
    printf("First nearest vertex to start: %d\n", first_nearest_vertex);
    
    while (tour_size < numOfCoords + 1) {
        //printf("Iteration: %d\n", k);
        //printf("numOfCoords: %d\n",numOfCoords);
       // printf("tour_size: %d\n",tour_size);
        NearestVertexResult nearestResult = NearestNotInTour(distanceMatrix, tour, numOfCoords);
        int next_vertex = nearestResult.nearest_vertex;
        int closest_tour_vertex = nearestResult.closest_tour_vertex;
        double distance = distanceMatrix[next_vertex][tour[closest_tour_vertex]];
        
         printf("Next nearest vertex not in tour: %d (Distance to closest tour vertex %d: %f)\n", 
                next_vertex, tour[closest_tour_vertex], distance);
        

        //printf("Next nearest vertex not in tour: %d\n", next_vertex);
        OptimalPositionInsertion(distanceMatrix, next_vertex, tour, tour_size, closest_tour_vertex, numOfCoords);
        //printf("Inserted vertex %d at position %d in the tour\n", next_vertex, tour_size);
        visited[next_vertex] = 1;
        //printf("visited[next_vertex]: %d\n",visited[next_vertex]);
        tour_size++;
        //printf("tour_size: %d\n",tour_size);


        printf("Current Tour: ");
        for (int i = 0; i < tour_size; i++) {
            printf("%d ", tour[i]);
        }
        printf("\n");
    }



    printf("Final Tour: ");
    for (int i = 0; i < tour_size; i++) {
        printf("%d ", tour[i]);
    }
    printf("\n");

    free(visited);
    return tour;
} 


int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    char *outputFile = argv[2];

    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    int *tour = nearestinsertion(distanceMatrix, numOfCoords);
    
    writeTourToFile(tour, numOfCoords + 1, outputFile);
    //register int i asm("esp"); //add this line         
    //printf("$esp = %#010x\n", i); // and this line will print esp

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
