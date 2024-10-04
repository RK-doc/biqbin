#ifndef ReadAndOutput_h
#define ReadAndOutput_h

#include <vector>
#include "Structures.h"

// note: 
// optimal = 1: optimal solution
// optimal = 0: approximate solution

void read_params(Params &params);
int read_input(char *argv[], int &n, int &m, std::vector<Vertex> &vertices, std::vector<Edge> &edges);
void print_solution(int n, int m, int sol[], int val, int counterBaBNodes, double timeTotal, double density, bool optimal);
void outputFile(char *argv[], int n, int m, int sol[], int val, int counterBaBNodes, double timeTotal, double density, bool optimal);
void outputFile_after_start(char *argv[], int n, int m, int sol[], int val,  int counterBaBNodes, double timeTotal, double density, bool begin);

#endif /* ReadAndOutput_h */
