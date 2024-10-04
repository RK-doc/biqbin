/*
 * Subproblem.cpp
 *
 *  Created on: 09.08.2017
 *      Author: mesiebenhofe
 */

#include "Subproblem.h"
#include <stdlib.h> //rand()

// constrcut original problem from graph
Subproblem::Subproblem(Graph& g)
{
	m_depth = 0;
	m_graph = Graph(g.vertices(), g.edges());
	m_upperBound = m_graph.n();
	m_constValue = 0;
	m_nextBranchVar = 0;	
}

// construct subproblem from received info
Subproblem::Subproblem(std::vector<int>& global_vertices, std::vector<Edge>& edges_original_graph, int depth, double upper_bound, int constValue, std::set<int>& stableSet)
{
	 
	std::vector<Vertex> vertices; 
	std::vector<Edge> edges;

	// set local numbers
	for(int i = 1; i <= global_vertices.size(); ++i)
	{
		vertices.emplace_back(global_vertices.at(i-1),i);
	}
	
	// determine edges in subgraph --> need to add LOCAL numbers of edges!!!!
	int i = 0;
	int edge_localVertex1;
	int edge_localVertex2;
	
	for (std::vector<Edge>::iterator it = edges_original_graph.begin() ; it != edges_original_graph.end(); ++it)
	{
		for (std::vector<Vertex>::iterator itv = vertices.begin() ; itv != vertices.end(); ++itv)
		{
			if ((itv->gNumber == it->vertex1) || (itv->gNumber == it->vertex2))
				++i;
		}

		if (i == 2) // both
		{
			for(int j = 0; j < global_vertices.size(); ++j)
			{
				if (global_vertices.at(j) == it->vertex1)
					edge_localVertex1 = j+1;

				else if (global_vertices.at(j) == it->vertex2)
					edge_localVertex2 = j+1;

				
			}
			edges.emplace_back(edge_localVertex1,edge_localVertex2);
		}

		i = 0;
	}	
	

	m_graph = Graph(vertices, edges);
	m_upperBound = upper_bound;
	m_depth = depth;
	m_nextBranchVar = 0;
	m_constValue = constValue;
	m_stableSet = stableSet;

}
// construct problem from branching
Subproblem::Subproblem(Subproblem& parent, int branchVar, int inSet){

	m_depth = parent.m_depth+1;
	m_graph = Graph(parent.m_graph.vertices(),parent.m_graph.edges());
	m_nextBranchVar = 0;
	m_upperBound = parent.m_upperBound;
	m_constValue = parent.m_constValue + inSet;
	m_stableSet = parent.m_stableSet;
	if(inSet == 1){
		m_stableSet.insert(m_graph.vertices()[branchVar-1].gNumber);
	}

	m_graph.fixVertex(branchVar, inSet);

}

Subproblem::~Subproblem() {

}

int Subproblem::computeUpperBound(){
	double upperBound = m_upperBound;
	int status = g_ubndComp->getUpperBound(&m_graph,upperBound, m_nextBranchVar, m_constValue);
	//update upper bound
	if(status == -1){
		//if bound computation was not successful:
		//upperBound remains unchanged
		//random branchvar (between 1 and n)
		if(m_graph.n()!=0){
			m_nextBranchVar = std::rand() % m_graph.n() + 1;
		} else {
			m_nextBranchVar = 0;
		}
		std::cerr << "Problem computing upper bound" << std::endl;
		counterBoundComputationFailed++;
		return status;
	}
	if(upperBound + m_constValue < m_upperBound){ //if upperBound lower than inherited upperBound -> update
		m_upperBound = upperBound + m_constValue;
	}

	return status;
}

int Subproblem::computeLowerBound(int* globalSol, int n){;
	int status = 0;
	int sizeGraph = m_graph.n();
	int sol[sizeGraph];
	int lowerBound = g_lowerBound - m_constValue;

	if(sizeGraph > params.minSizeGraphSizeForFullEnum){
		status = g_lbndComp->getLowerBound(&m_graph,lowerBound, sol); //updates lowerBound and sol if better lower bound found
		lowerBound += m_constValue;
	}
	else{
		std::vector<int> stable = m_graph.getMaximumStableSet();

		//get solution from stableSet
		for(int i = 0; i<sizeGraph; i++){
			sol[i] = 0;
		}
		for(auto it: stable){
			sol[it-1] = 1;
		}

		lowerBound = stable.size() + m_constValue;
	}

	if(status != 0){
		//if bound computation was not successful --> do nothing
		std::cerr << "Problem computing lower bound" << std::endl;
		counterBoundComputationFailed++;
		return status;
	}

	if(g_lowerBound < lowerBound){
		//std::cout << "update global solution\n";
		g_lowerBound = lowerBound;
		for(int i = 0; i<n; i++){
			globalSol[i] = 0;
		}
		for (auto it : m_stableSet){
			globalSol[it-1] = 1;
		}
		std::vector<Vertex> vertices = m_graph.vertices();
		for(int i = 0; i<sizeGraph; i++){
			if(sol[i]==1){
				globalSol[vertices.at(i).gNumber - 1] = 1;
			}
		}
	}
	return status;
}

void Subproblem::printInfo1(int count){
	std::cout << "===================================================" << std::endl;
	std::cout << "Node " << count << " depth: "<< m_depth <<std::endl;
	std::cout << "upper bound: " << m_upperBound << std::endl;
	std::cout << "constant term: " << m_constValue << std::endl;
	std::cout << "Stable Set: " ;
	for (auto it : m_stableSet){
		std::cout << it << " ";
	}
	std::cout << std::endl;

	if(params.outputVerticesSubproblem==1 || params.outputEdgesSubproblem==1){
		std::cout << "-----------------------" << std::endl;
		std::cout << "number of vertices: " << m_graph.n() << std::endl;
		std::cout << "number of edges: " << m_graph.m() << std::endl;
		if(params.outputEdgesSubproblem==1){
			std::cout << "edges: "  << std::endl;
			std::vector<Vertex> vertices = m_graph.vertices();
			if(params.outputVertexGlobal==1 && params.outputVertexLocal==1){
				std::cout << "global/local" << std::endl;
				for(auto it : m_graph.edges()){
					std::cout <<"[" << "(" << vertices.at(it.vertex1 - 1).gNumber << "," << vertices.at(it.vertex2 - 1).gNumber << ")"; //edges with global number
					std::cout << "," <<  "(" << it.vertex1 << "," << it.vertex2 << ")" << "]"; //edges with local number
					std::cout << std::endl;
				}
			}
			else if(params.outputVertexLocal==1){
				for(auto it : m_graph.edges()){
					std::cout <<  "(" << it.vertex1 << "," << it.vertex2 << ")" ; //edges with local number
					std::cout << std::endl;
				}
			}
			else if(params.outputVertexGlobal==1){
				for(auto it : m_graph.edges()){
					std::cout << "(" << vertices.at(it.vertex1 - 1).gNumber << "," << vertices.at(it.vertex2 - 1).gNumber << ")"; //edges with global number
					std::cout << std::endl;
				}
			}
		}
		if(params.outputVerticesSubproblem==1){
			std::cout << "vertices:" << std::endl;
			if(params.outputVertexGlobal==1 && params.outputVertexLocal==1){
				std::cout <<  "(g,l)" << std::endl;
				for(auto it : m_graph.vertices()){
					std::cout <<  "[" << it.gNumber << "," << it.lNumber << "]" << std::endl;
				}
			}
			else if(params.outputVertexLocal==1){
				for (auto it : m_graph.vertices()){
					std::cout << it.lNumber << " ";
				}
				std::cout << std::endl;
			}
			else if(params.outputVertexGlobal==1){
				for(auto it : m_graph.vertices()){
					std::cout << it.gNumber << " ";
				}
			}
		}
	}
}

void Subproblem::printInfo2(){
	std::cout << "-----------------------" << std::endl;
	std::cout << "new upper bound for this subproblem: " << m_upperBound << std::endl;
	std::cout << "best lower bound: " << g_lowerBound << std::endl;
	//std::cout << "next branching variable: " << m_nextBranchVar << std::endl;
}

void Subproblem::printInfo3(){
	std::cout << "new best lower bound: " << g_lowerBound << std::endl;
}

void Subproblem::printInfo4(){

	if(params.outputVertexGlobal==1 && params.outputVertexLocal==1){
		std::cout << "next branching variable: " << "[" << m_graph.vertices().at(m_nextBranchVar - 1).gNumber << "," << m_nextBranchVar  << "]" << std::endl;
	}
	else if(params.outputVertexGlobal==1){
		std::cout << "next branching variable: "  << m_graph.vertices().at(m_nextBranchVar - 1).gNumber << std::endl;
	}
	else if(params.outputVertexLocal==1){
		std::cout << "next branching variable: " << m_nextBranchVar << std::endl;
	}
	else{
		std::cout << "next branching variable: "  << m_graph.vertices().at(m_nextBranchVar - 1).gNumber << std::endl;
	}
}
