#ifndef Graph_h
#define Graph_h

#include <vector>
#include <set>
#include "Structures.h"

extern Params params;


class Graph {
    
protected:
    // number of vertices
    int m_n;        

    // number of edges
    int m_m;

    // vertices and edges
    std::vector<Vertex> m_vertices;
    std::vector<Edge> m_edges;
    
public:
    Graph();
    virtual ~Graph();
 
    Graph(std::vector<Vertex> &vertices, std::vector<Edge> &edges);
    
    // getters
    int n() const { return m_n; }
    int m() const { return m_m; }
    std::vector<Edge>& edges(){ return m_edges; }
    std::vector<Vertex>& vertices(){ return m_vertices; }
    
    // fix vertex: 1 = in Set, 0 = not in Set
    void fixVertex(int vertexNumber, int inSet);
    
    // for FULL enumeration on small graphs
    // input arguments: size of maximum stable set can be between 
    // minimum stable set found and upper bound
    std::vector<int> getMaximumStableSet(double upperBound, int lowerBound);
    
    // made public only for for_each_combination (from ComAndPerm.h)
    // function to check if given set is stable
    bool isStableSet(std::vector<int>::iterator i, std::vector<int>::iterator j);
    
private:
    // help variables for getMaximumStableSet

    // is stable set found
    bool m_stableSetFound;

    // adjacency matrix
    int* m_adjacency;
    
    // functions to renumber vertices and edges
    void setLocalNumberVertices();
    void setLocalNumberEdges();
    inline void deleteEdges(int vertexNumber);
    inline void deleteVertex(int vertexNumber);
    
    
};

#endif /* Graph_h */
