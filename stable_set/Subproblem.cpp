#include "Subproblem.h"
#include <iostream>


Subproblem::~Subproblem() {}
Subproblem::Subproblem(Graph &g) : m_graph(g) {

    m_upperBound = g.n();
    m_constValue = 0;           // number of already fixed vertices
    m_nextBranchVar = 0;        // determined after the computation of upper bound

}

Subproblem::Subproblem(Subproblem &parent, int branchVar, int inSet) : m_graph(parent.m_graph), m_stableSet(parent.m_stableSet) {
    
    m_nextBranchVar = 0;
    m_constValue = parent.m_constValue + inSet;     // add fixed vertex
   
    if(inSet == 1){
        m_stableSet.insert(m_graph.vertices()[branchVar-1].gNumber);
    }
    m_graph.fixVertex(branchVar, inSet);
    
    // inherited upper bound vs. size of the graph    m_upperBound = (parent.m_upperBound < m_graph.n()) ? parent.m_upperBound : m_graph.n(); 

}


// constructor for Subproblem when received over message (MPI) --> to send and receive subproblem we only need global numbers of vetices of subgraph
// used for PARALLEL BRANCH & BOUND
Subproblem::Subproblem(std::vector<int> &global_vertices, std::vector<Edge> &edges_original_graph, double upper_bound, int constValue, std::set<int> &stableSet){
    
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    
    // set local numbers
    for(int i = 1; i <= global_vertices.size(); ++i)
    {
        vertices.emplace_back(global_vertices[i-1],i);
    }
    
    // determine edges in subgraph --> need to add LOCAL numbers of vertices!!
    int i = 0;
    int edge_localVertex1;
    int edge_localVertex2;
    
    for (std::vector<Edge>::iterator it = edges_original_graph.begin(); it != edges_original_graph.end(); ++it)
    {
        for (std::vector<Vertex>::iterator itv = vertices.begin(); itv != vertices.end(); ++itv)
        {
            if ( itv->gNumber == it->vertex1 ){
            	edge_localVertex1 = itv->lNumber;
                ++i;
            }
            if ( itv->gNumber == it->vertex2 ){
            	edge_localVertex2 = itv->lNumber;
                ++i;
            }
            
            if ( i == 2 ) {	// new edge discovered
            	edges.emplace_back(edge_localVertex1,edge_localVertex2);
            	break;		// edge found
            }
        }       
        i = 0;		// reset i
    }	
    
    m_graph = Graph(vertices, edges);
    m_upperBound = upper_bound;
    m_nextBranchVar = 0;
    m_constValue = constValue;
    m_stableSet = stableSet; 
}

// compute upper bound
void Subproblem::computeUpperBound(){
    
    // upper bound and next branching variable are changed in getUpperBound routine!!
    // temp var -> upper bound without taking into account fixed variables
    double upperBound; 

    g_bndComp->getUpperBound(&m_graph, upperBound, m_nextBranchVar);

/*  TODO: what if upper bound is not computed successfully  
    if(status == -1){
        // if bound computation was not successful:
        // upperBound remains unchanged
        // random branchvar (between 1 and n)
        if(m_graph.n() != 0){
            m_nextBranchVar = std::rand() % m_graph.n() + 1;
        }
        else {
            m_nextBranchVar = 0;
        }
        std::cerr << "Problem computing upper bound" << std::endl;
        return status;
    }
 */   

    // if upperBound lower than inherited upperBound -> update
    if(upperBound + m_constValue < m_upperBound){
        m_upperBound = upperBound + m_constValue;
    }   
}


void Subproblem::computeLowerBound(int* globalSol, int n){
    
    int sizeGraph = m_graph.n();
    int *sol = new int[sizeGraph];					// int sol[sizeGraph] is variable length
    int lowerBound = g_lowerBound - m_constValue;   // local lower bound
    
    //updates lowerBound and sol if better lower bound found
    g_bndComp->getLowerBound(&m_graph,lowerBound, sol);
    
    lowerBound += m_constValue;

    if(g_lowerBound < lowerBound){
        
        g_lowerBound = lowerBound;
        for(int i = 0; i<n; i++){
            globalSol[i] = 0;
        }
        
        // already fixed vertices in the stable set (from before) --> global numbers!!
        for (auto it : m_stableSet){
            globalSol[it-1] = 1;
        }
        
        // add new vertices
        std::vector<Vertex> vertices = m_graph.vertices();
        for(int i = 0; i < sizeGraph; i++){
            if( sol[i]==1 ){
                globalSol[vertices.at(i).gNumber - 1] = 1;
            }
        }
    }
    delete[] sol;
}


void Subproblem::doFullEnumeration(int *globalSol, int n){
    
    // use local lower bound: g_lowerBound - m_constValue
    std::vector<int> stable = m_graph.getMaximumStableSet(m_upperBound, g_lowerBound - m_constValue);
    
    // this would be the new global lower bound
   int lowerBound = stable.size() + m_constValue;
    
    // update if better lower bound found
    if(g_lowerBound < lowerBound){
        
        g_lowerBound = lowerBound;
        for(int i = 0; i<n; i++){
            globalSol[i] = 0;
        }
        
        // already fixed vertices in the stable set (from before) --> global numbers!!
        for (auto it : m_stableSet){
            globalSol[it-1] = 1;
        }
        
        // new vertices in the stable set
        std::vector<Vertex> vertices = m_graph.vertices();

        for(auto index : stable){
            globalSol[vertices[index - 1].gNumber - 1] = 1;
        }
    }
}
