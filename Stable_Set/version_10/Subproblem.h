#ifndef Subproblem_h
#define Subproblem_h

#include "Graph.h"
#include "BoundComputationTheta.h"

extern BoundComputationTheta* g_bndComp;
extern int g_lowerBound;                    // global lower bound


class Subproblem{
protected:
    Graph m_graph;
    double m_upperBound;                    // local upper bound
    int m_constValue;                       // add to objective function (number of already fixed vertices)
    std::set<int> m_stableSet;              // already fixed vertices in the stable set
    int m_nextBranchVar;
    
    
public:
    Subproblem(){}
    virtual ~Subproblem();
    
    // initial problem with whole graph
    Subproblem(Graph &g);   

    // construct subproblem from parent by branching on variable
    Subproblem(Subproblem& parent, int branchVar, int inSet);
    
    // constructor for Subproblem when received over message (MPI)
    Subproblem(std::vector<int>& vertices, std::vector<Edge>& edges, double upper_bound, int constValue, std::set<int> &stableSet);
    
    // getters
    Graph& graph(){ return m_graph; }
    double upperBound() const { return m_upperBound; }
    int nextBranchVar() const { return m_nextBranchVar; }
    int getConstValue(){ return m_constValue; }
    std::set<int>& stableSet(){ return m_stableSet; }
    
    // methods for computing upper and lower bounds
    void computeUpperBound();                        
    void computeLowerBound(int* globalSol, int n);  
    void doFullEnumeration(int* globalSol, int n);

    
};

#endif /* Subproblem_h */
