#include <iostream>
#include <stdlib.h>
#include "Graph.h"
#include "CombAndPerm.h"
#include <functional>       // for std::bind


Graph::Graph(){}
Graph::~Graph(){}

Graph::Graph(std::vector<Vertex> &vertices, std::vector<Edge> &edges) : m_n(vertices.size()), m_m(edges.size()), m_vertices(vertices), m_edges(edges) {}


void Graph::fixVertex(int vertexNumber, int inSet){
    if(inSet == 0){
        // 1) delete edges:
        deleteEdges(vertexNumber);
        
        // 2) delete vertex
        deleteVertex(vertexNumber);
    
        // 3) set local number
        setLocalNumberEdges();
        setLocalNumberVertices();
    }
    if(inSet == 1){
        // 1) delete vertex
        deleteVertex(vertexNumber);
        
        // 2) delete all neighbors j from vertex and their edges
        
        // find neightbors
        std::set<int> neighbors;
        for (std::vector<Edge>::iterator it = m_edges.begin(); it != m_edges.end(); ++it){
            if(vertexNumber == it->vertex1){
                neighbors.insert(it->vertex2);
            } else if(vertexNumber == it->vertex2){
                neighbors.insert(it->vertex1);
            }
        }
        
        for(auto j : neighbors){
            deleteVertex(j);
            deleteEdges(j);
        }
        
        // 3) set local number
        setLocalNumberEdges();
        setLocalNumberVertices();
    }
    
    // compute new sizes of vertices and edges
    m_m = m_edges.size();
    m_n = m_vertices.size();
    
}

std::vector<int> Graph::getMaximumStableSet(double upperBound, int lowerBound){
    if(lowerBound < 1)
        lowerBound = 1;     //lowerBound needs to be 1

    int upperBoundInt = m_n;

    if(upperBound < m_n)
        upperBoundInt = static_cast<int>(upperBound);   //floor
    
    std::vector<int> stableSet;

    if(m_n == 0) {return stableSet;}
    
    // m_adjacency matrix in row format
    m_adjacency = (int*)malloc( (m_n*m_n) * sizeof(int) );

    for(std::vector<Edge>::iterator ite = m_edges.begin() ; ite != m_edges.end(); ++ite){
        m_adjacency[m_n * (ite->vertex1-1) + (ite->vertex2-1)] = 1;
        m_adjacency[m_n * (ite->vertex2-1) + (ite->vertex1-1)] = 1;
    }

    //stableSet = {1,...,n}
    for(int i = 0; i < m_n; i++){
        stableSet.push_back(i+1);
    }
    m_stableSetFound = false;
    
    //for(int k = m_n; k > 1; k--){
    for(int k = upperBoundInt; k > lowerBound; k--){
        // from k = upperBound because it can not be larger and to lowerbound + 1
        // because we want to get a better lower bound
        // for all S with |S| = k  (k = 1 --> stableSet)
        // check if stableSet
        for_each_combination(stableSet.begin(),stableSet.begin()+k, stableSet.end(),
                             std::bind(&Graph::isStableSet, this, std::placeholders::_1, std::placeholders::_2));
        
        if(m_stableSetFound){            
            //first k entries in stableSet are a maximum stable set
            stableSet.resize(k);
            free(m_adjacency);
            return stableSet; //return the first k elements
        }
        
    }
    
    //no stable set of size > lowerbound found --> return any vertex as stable set
    stableSet.resize(1);

    free(m_adjacency);
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


void Graph::setLocalNumberVertices(){
    // renumber vertices from 1 to n
    int i = 1;
    for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ++it){
        it->lNumber = i++;
    }
    m_n = m_vertices.size();
}

void Graph::setLocalNumberEdges(){
    // renumber the edges
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
    m_m = m_edges.size();
}

inline void Graph::deleteEdges(int vertexNumber){
    for (std::vector<Edge>::iterator it = m_edges.begin() ; it != m_edges.end(); ){
        if( (it->vertex1 == vertexNumber) || (it->vertex2==vertexNumber) ){
            it = m_edges.erase(it);
        }
        else {
            ++it;
        }
    }
    m_m = m_edges.size();
}

inline void Graph::deleteVertex(int vertexNumber){
    for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ){
        if(it->lNumber == vertexNumber){
            it = m_vertices.erase(it);
            break;
        }
        else {
            ++it;
        }
    }
    m_n = m_vertices.size();
}
