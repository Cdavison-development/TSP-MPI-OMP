#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coordReader.h"

    double euclideanDistance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return sqrt(dx * dx + dy * dy);
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