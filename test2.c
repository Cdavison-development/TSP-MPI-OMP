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



double distance(double **adj_matrix, int from_vertex, int to_vertex) {
    return adj_matrix[from_vertex][to_vertex];
}



double tour_length_matrix(double **adj_matrix, int *tour_order, int num_vertices) {
    double length = 0.0;
    int i;

    for (i = 0; i < num_vertices - 1; i++) {
        length += distance(adj_matrix, tour_order[i], tour_order[i + 1]);
    }

    // Add distance from last to first vertex to complete the tour
    length += distance(adj_matrix, tour_order[num_vertices - 1], tour_order[0]);
    return length;
}

// Function to generate a distance matrix from coordinates
double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            double dx = coords[i][0] - coords[j][0];
            double dy = coords[i][1] - coords[j][1];
            matrix[i][j] = matrix[j][i] = sqrt(dx * dx + dy * dy); // Calculate Euclidean distance
        }
    }
    return matrix;
}

int find_nearest_vertex(double **adj_matrix, int *tour, int tour_size, int num_vertices, int *visited) {
    double min_distance = DBL_MAX;
    int nearest_vertex = -1;
    int i;

    for (i = 0; i < num_vertices; i++) {
        if (!visited[i]) {
            double dist = distance(adj_matrix, tour[tour_size - 1], i);
            if (dist < min_distance) {
                nearest_vertex = i;
                min_distance = dist;
            }
        }
    }

    return nearest_vertex;
}

int *nearest_insertion_matrix(double **adj_matrix, int num_vertices) {
    int *visited = (int *)calloc(num_vertices, sizeof(int));
    int *tour = (int *)malloc(num_vertices * sizeof(int));
    int tour_size = 0;
    int i, j, start_vertex;

    srand(time(NULL)); // Seed for random number generation
    start_vertex = rand() % num_vertices;
    tour[tour_size++] = start_vertex;
    visited[start_vertex] = 1;

    while (tour_size < num_vertices) {
        int next_vertex = find_nearest_vertex(adj_matrix, tour, tour_size, num_vertices, visited);
        double min_insert_cost = DBL_MAX;
        int optimal_position = 0;

        for (i = 0; i < tour_size; i++) {
            int next_index = (i + 1) % tour_size;
            double current_cost = distance(adj_matrix, tour[i], next_vertex) + distance(adj_matrix, next_vertex, tour[next_index]);

            if (current_cost < min_insert_cost) {
                min_insert_cost = current_cost;
                optimal_position = i;
            }
        }

        // Shift elements to make space for the new vertex
        for (j = tour_size; j > optimal_position + 1; j--) {
            tour[j] = tour[j - 1];
        }

        // Insert the next vertex
        tour[optimal_position + 1] = next_vertex;
        visited[next_vertex] = 1;
        tour_size++;
    }

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