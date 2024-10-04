/*
 * Graph.h
 *
 *  Created on: 08.08.2017
 *      Author: mesiebenhofe
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include "Structures.h"


extern Params params;

class Graph {

protected:
	int m_n;
	int m_m;
	std::vector<Edge> m_edges;
	std::vector<Vertex> m_vertices;




public:
	Graph();
	virtual ~Graph();

	Graph(/*int n, int m,*/ std::vector<Vertex>& vertices,std::vector<Edge>& edges);
	int fixVertex(int vertexNumber, int inSet); //1 = in Set, 0 = not in Set

	std::vector<int> getMaximumStableSet();
	bool isStableSet(std::vector<int>::iterator i, std::vector<int>::iterator j);

	void printInfo();
	void printRudyFormatToFile(std::string filename);

	//getter
	int n(){ return m_n;}
	int m(){ return m_m;}
	std::vector<Edge>& edges(){ return m_edges;}
	std::vector<Vertex>& vertices(){ return m_vertices;}



private:
	//help variables for getBiggestStableSet
	bool m_stableSetFound;
	int* m_adjacency;

	inline void setLocalNumberVertices();
	inline void setLocalNumberEdges();
	inline void deleteEdges(int vertexNumber);
	inline void deleteVertex(int vertexNumber);

	//Graph(const Graph&);
	//Graph& operator=(const Graph&);
};

#endif /* GRAPH_H_ */
