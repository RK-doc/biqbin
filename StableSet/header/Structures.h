/*
 * Structures.h
 *
 *  Created on: 30.08.2017
 *      Author: mesiebenhofe
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <map>

struct Params{
	//program
	int DFS = 1;
	int minSizeGraphSizeForFullEnum = 1;
	int subproblemWithVertexInSetFirst = 1;
	int strategyBranchVar = 1; //1: random, 2: min, 3:max, 4:closest to mean, 5: closest to 0.5, 6: median
	int lowerBoundForEveryNth = 1;
	int lowerBoundInRoot = 4;
	int maxExecutionTime = 60*60*1; //1 hour
	std::map<int,int> numberExactSubgraphs; //(size, numberOfGraphs)
	//bundle
	double mL;
	double epsStop;
	int Lmax;
	double muSSFactor;
	double muNSFactor;
	double muMin;
	double mr;

	//output
	int printlevel = 1;
	int outputEdgesSubproblem = 0;
	int outputVerticesSubproblem = 0;
	int outputVertexLocal = 0;
	int outputVertexGlobal = 1;
	int logfile = 0;
};


struct Vertex{
	int gNumber; //global number
	int lNumber; //local number
	Vertex(int g, int l): gNumber(g), lNumber(l){}
};

struct Edge{
	int vertex1;
	int vertex2;
	Edge(int i,int j): vertex1(i), vertex2(j){
	}
};


class Structures {
public:
	Structures();
	virtual ~Structures();
};

#endif /* STRUCTURES_H_ */
