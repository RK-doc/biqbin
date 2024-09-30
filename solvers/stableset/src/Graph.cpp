/*
 * Graph.cpp
 *
 *  Created on: 08.08.2017
 *      Author: mesiebenhofe
 */

#include "Graph.h"
#include "CombAndPerm.h"

#include <functional>
#include <fstream>      // std::ofstream

Graph::Graph(){

}

Graph::Graph(std::vector<Vertex>& vertices, std::vector<Edge>& edges){
	m_vertices = vertices;
	m_n = m_vertices.size();
	m_edges = edges;
	m_m = m_edges.size();

}

int Graph::fixVertex(int vertexNumber, int inSet){
	if(inSet == 0){
		//1) delete edges:
		deleteEdges(vertexNumber);

		//2) delete vertex
		deleteVertex(vertexNumber);

		//3) set local number
		setLocalNumberEdges();
		setLocalNumberVertices();
	}

	if(inSet == 1){
		//1) delete vertex
		deleteVertex(vertexNumber);

		//2) delete all neighbors j from vertex and their edges
		std::set<int> neighbors;
		for (std::vector<Edge>::iterator it = m_edges.begin() ; it != m_edges.end(); ++it){
			if(vertexNumber == it->vertex1){
				neighbors.insert(it->vertex2);
			} else if(vertexNumber == it->vertex2){
				neighbors.insert(it->vertex1);
			}
		}
		for(auto& j: neighbors){
			deleteVertex(j);
			deleteEdges(j);
		}

		//3) set local number
		setLocalNumberEdges();
		setLocalNumberVertices();
	}

	m_m = m_edges.size();
	m_n = m_vertices.size();


	return 0;
}

Graph::~Graph() {
}


std::vector<int> Graph::getMaximumStableSet(){
	//generate adjacency matrix
	std::vector<int> stableSet;
	if(m_n == 0){ return stableSet;}

	int adjacency[m_n][m_n];// = {0};  //Initialize with 0
	for(int i = 0; i<m_n; i++){
		for(int j = 0; j<m_n; j++){
			adjacency[i][j] = 0;
		}
	}
	m_adjacency = &adjacency[0][0];
	for(std::vector<Edge>::iterator ite = m_edges.begin() ; ite != m_edges.end(); ++ite){
		adjacency[ite->vertex1-1][ite->vertex2-1] = 1;
		adjacency[ite->vertex2-1][ite->vertex1-1] = 1;
	}


	for(int i = 0; i < m_n; i++){
		stableSet.push_back(i+1);
	}
	m_stableSetFound = false;

	for(int k = m_n; k > 1; k--){

		//for all S with |S| = k  (k = 1 --> stableSet)
		//check if stableSet
		//for_each_combination(stableSet.begin(),stableSet.begin()+k, stableSet.end(), isStableSet); //stops, if stableSet found
		for_each_combination(stableSet.begin(),stableSet.begin()+k, stableSet.end(),
						std::bind(&Graph::isStableSet, this, std::placeholders::_1, std::placeholders::_2));

		if(m_stableSetFound){
#ifdef DEBUG
			if(params.printlevel>1){
				std::cout << "stable set found on subgraph with size: " << k << std::endl;
			}
#endif

			//stableSet.erase(stableSet.begin()+k, stableSet.end());
			stableSet.resize(k);
			return stableSet; //return the first k elements
		}

	}

	stableSet.resize(1);
	return stableSet;
}

bool Graph::isStableSet(std::vector<int>::iterator start, std::vector<int>::iterator end){
	for (std::vector<int>::iterator itv1 = start ; itv1 != end; ++itv1){
		for (std::vector<int>::iterator itv2 = itv1+1 ; itv2 != end; ++itv2)
			if(/*m_adjacency[*itv1][*itv2]*/ m_adjacency[m_n*((*itv1)-1)+((*itv2)-1)] == 1){
				return false; //get next combination
			}
	}
	m_stableSetFound = true;
	return true; //stop iterator for combinations
}



void Graph::printInfo(){
	std::cout << "===================================================" << std::endl;
	std::cout << "Graph" <<  std::endl;
	std::cout << "number of vertices: " << m_n << std::endl;
	std::cout << "number of edges: " << m_m << std::endl;
	std::cout << "edges: "  << std::endl;
	for (auto it : m_edges){
		std::cout <<  "(" << it.vertex1 << "," << it.vertex2 << ")" << std::endl;
	}
	std::cout << "vertices:" << std::endl;
	std::cout <<  "(g,l)" << std::endl;
	for (auto it : m_vertices){
		std::cout <<  "(" << it.gNumber << "," << it.lNumber << ")\n";
	}
	std::cout << "\n";
}


inline void Graph::setLocalNumberVertices(){
	int i = 1;
	for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ++it){
		it->lNumber = i++;
	}
}

inline void Graph::setLocalNumberEdges(){
	int i = 1;
	for (std::vector<Vertex>::iterator itv = m_vertices.begin() ; itv != m_vertices.end(); ++itv){
		if(itv->lNumber != i){ //if local number of vertex changed
			for (std::vector<Edge>::iterator ite = m_edges.begin() ; ite != m_edges.end(); ++ite){
				if(ite->vertex1 == itv->lNumber){
					ite->vertex1 = i;
				}
				if(ite->vertex2 == itv->lNumber){
					ite->vertex2 = i;
				}
			}
		}
		i++;
	}
}

inline void Graph::deleteEdges(int vertexNumber){
	for (std::vector<Edge>::iterator it = m_edges.begin() ; it != m_edges.end(); ){
	    if(it->vertex1==vertexNumber || it->vertex2==vertexNumber){
	    	it = m_edges.erase(it);
	    } else {
	    	++it;
	    }
	}
	m_m = m_edges.size();
}

inline void Graph::deleteVertex(int vertexNumber){
	for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ){
		if(it->lNumber == vertexNumber){
			it = m_vertices.erase(it);
		} else {
			++it;
		}
	}
	m_n = m_vertices.size();
}

void Graph::printRudyFormatToFile(std::string filename){
	std::ofstream rudyFile;
	rudyFile.open(filename);

	rudyFile << m_n << " " << m_m << std::endl;
	for(auto i: m_edges){
		rudyFile << i.vertex1 << " " << i.vertex2 << " " << 1 << std::endl;
	}

	rudyFile.close();
}