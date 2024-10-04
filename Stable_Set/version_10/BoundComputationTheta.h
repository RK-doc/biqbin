#ifndef BoundComputationTheta_h
#define BoundComputationTheta_h

#include <vector>
#include <armadillo>
#include "Graph.h"
#include "Heuristic.h"

using namespace arma;

// declarations for Boundary Point Method and Interior Point Method
void bpm(std::vector<Edge> &edge_list, int n, double &uBound, mat &X);
void ipm_Lovasz_theta(std::vector<Edge> &edge_list, int n, double &uBound, mat &X);


class BoundComputationTheta {
public:

    BoundComputationTheta();
    virtual ~BoundComputationTheta();

    void getUpperBound(Graph *g, double &upperBound, int &branchVar);
    void getLowerBound(Graph *g, int &lowerBound, int *sol);
    
protected:
    // vector x from theta --> is used to determine the next branching variable
    std::vector<float> m_vector;
    
    // boolean to indicate whether x was set by the upper bound computation
    bool m_vectorSet;
    
    // boolean to indicate if we are in the root node
    // set to true in root node and set to false after first lower bound computation
    bool m_rootNode;
   
    // set vector x from theta matrix X for lower bound computation
    void setVector(mat &X, double theta, int n);
    
    // set next branching varaible
    void setBranchVar(int &branchVar, int n);
    
    // Burer-Monteiro heuristic: combinations: 11, 21, 25
    std::vector<int> getStableSetWithAO(Graph* g, int rank, int numResets, int maxSeconds); 
};


#endif /* BoundComputationTheta_h */
