#ifndef COORDREADER_H
#define COORDREADER_H

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);
#endif // COORDREADER_H