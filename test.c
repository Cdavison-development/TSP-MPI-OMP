#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
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

int find_nearest_vertex(double **adj_matrix,  int current_vertex, int numOfCoords, int *visited) {
    double min_distance = DBL_MAX;
    int nearest_vertex = -1;

    for (int i = 0; i < numOfCoords; i++) {
        if (!visited[i] && i != current_vertex) {
            double dist = adj_matrix[current_vertex][i];
            if (dist < min_distance) {
                nearest_vertex = i;
                min_distance = dist;
            }
        }
    }
    printf("Nearest vertex found: %d, Distance: %f\n", nearest_vertex, min_distance);
    return nearest_vertex;
}


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

int *nearest_insertion_matrix(double **adj_matrix, int numOfCoords) {
    int *visited = (int *)calloc(numOfCoords, sizeof(int));
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));
    
    // Start with the first vertex (v[0])
    int start_vertex = 0;
    double min_distance = DBL_MAX;
    int first_nearest_vertex = -1;

// Find the nearest vertex to v[0]
    for (int j = 0; j < numOfCoords; j++) {
        if (j != start_vertex && adj_matrix[start_vertex][j] < min_distance) {
            first_nearest_vertex = j;
            min_distance = adj_matrix[start_vertex][j];
        }
    }

    visited[start_vertex] = 1;
    visited[first_nearest_vertex] = 1;
    tour[0] = start_vertex;
    tour[1] = first_nearest_vertex;
    tour[2] = start_vertex; // Completing the initial loop
    int tour_size = 3;

    printf("Initial Tour: ");
    for (int j = 0; j < tour_size; j++) {
        printf("%d ", tour[j]);
        }
    printf("\n");

    while (tour_size < numOfCoords + 1) {
        int last_vertex_in_tour = tour[tour_size - 1];
        int next_vertex = find_nearest_vertex(adj_matrix, last_vertex_in_tour, numOfCoords, visited);
        printf("Next Vertex to insert: %d\n", next_vertex);

        double min_insert_cost = DBL_MAX;
        int optimal_position = 0;

        // Find the best position to insert the next vertex
        for (int i = 0; i < tour_size; i++) {
            int next_index = (i + 1) % tour_size;
            double cost_with_next_vertex = adj_matrix[tour[i]][next_vertex] + adj_matrix[next_vertex][tour[next_index]];
            //printf("cost_with_next_vertex: %d\n", cost_with_next_vertex);
            double original_cost = adj_matrix[tour[i]][tour[next_index]];
            //printf("original_cost: %d\n", original_cost);
            double current_cost = cost_with_next_vertex - original_cost;
            //printf("current_cost: %d\n", current_cost);

             printf("Checking insertion between %d and %d, Cost: %f\n", tour[i], tour[next_index], current_cost);

            if (current_cost < min_insert_cost) {
                min_insert_cost = current_cost;
                optimal_position = i + 1;
            }
        }

        printf("Inserting %d between %d and %d with cost %f at position %d\n", next_vertex, tour[optimal_position - 1], tour[optimal_position], min_insert_cost, optimal_position);
        // Shift elements to make space for the new vertex
        for (int j = tour_size; j >= optimal_position; j--) {
            tour[j] = tour[j - 1];
        }

        // Insert the next vertex at the optimal position
        tour[optimal_position] = next_vertex;
        visited[next_vertex] = 1;
        tour_size++;

        printf("Updated Tour: ");
        for (int j = 0; j < tour_size; j++) {
            printf("%d ", tour[j]);
        }
        printf("\n");
    }

    // Add the start vertex to the end to complete the tour
    
    printf("Final Tour: ");
    for (int j = 0; j < tour_size; j++) {
        printf("%d ", tour[j]);
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

    int *tour = nearest_insertion_matrix(distanceMatrix, numOfCoords);
    
    writeTourToFile(tour, numOfCoords, outputFile);

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