#include "BoundComputationTheta.h"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
      
using namespace arma; 

BoundComputationTheta::BoundComputationTheta(){
    // first computation: node in branch and bound tree is the root
    m_rootNode = true;
    m_vectorSet = false;
    
    //set seed --> for setting branching variable randomly if params
    std::srand(static_cast<unsigned int>(std::time(0)));
}
 
BoundComputationTheta::~BoundComputationTheta() {}

// computation of upper bound --> THETA FUNCTION computed with Interior or Boundary Point Method
// depending on the number of constraints and density of the graph
void BoundComputationTheta::getUpperBound(Graph *g, double &upperBound, int &branchVar){ //, int constValueUpperBound){

    m_vectorSet = false;
    int n = g->n();
    
    // case: empty graph
    if( n == 0 ) {
        branchVar = 0;
        upperBound = 0;
        return;
    }
    
    // compute density of the graph
    double density = 2 * static_cast<double>(g->m()) / (n*(n-1));
    
    // decide which SDP solver to use:
    // interior-point method (sparse graps) or boundary Point Method (dense graphs)
    if ( (n <= 100) || ((density < 0.1) && (g->m() <= 15000)) ){
    
    	// std::cout << "Using Interior Point Method.\n";
        
        // ***************** INTERIOR-POINT METHOD *****************
        mat X = 1.0/n * eye<mat>(n,n); // strictly feasible solution for IPM

        ipm_Lovasz_theta(g->edges(), n, upperBound, X);

        setVector(X, upperBound, n);
        setBranchVar(branchVar, n);
	}   
    else // use Boundary Point Method
    {
        // std::cout << "Using Boundary Point Method.\n";
        
        // ***************** BOUNDARY POINT METHOD *****************
        mat X = zeros<mat>(n,n);
        
        bpm(g->edges(), n, upperBound, X);

        setVector(X, upperBound, n);
        setBranchVar(branchVar, n);
    }
}

// LOWER BOUND COMPUTATION --> local lower bound for the reduced graph!
// p.computeLowerBound() routine checks if better found and set it to global variables
void BoundComputationTheta::getLowerBound(Graph *g, int &lowerBound, int *sol){
    
    int n = g->n();

    if( n == 0 ){
        return;
    }
 
    // try to find a better lower bound with Burer-Monteiro heuristic
    std::vector<int> newStableSet;
    
    if(m_rootNode){
        newStableSet = getStableSetWithAO(g, 2, 5, 50);     // take longer time for root node: BIGGER GRAPHS
        m_rootNode = false;
    }
    else if(n >= 200){
        newStableSet = getStableSetWithAO(g, 2, 1, 5);
    }
    else{
        newStableSet = getStableSetWithAO(g, 2, 1, 1);
    }
    
    
    // update lower bound
    if( lowerBound < newStableSet.size() ){

        lowerBound = newStableSet.size();
        
        // update solution vector
        for(int i = 0; i < n; i++){
            sol[i] = 0;
        }
        for(int i = 0; i < newStableSet.size(); ++i){
            sol[newStableSet[i] - 1] = 1;
        }
    }
    
}

// Set m_vector --> from Interior or Boundary Point Method solution X
void BoundComputationTheta::setVector(mat &X, double theta, int n){
    
    // vector: theta * diag(X)
    m_vector.clear();
    m_vector = conv_to<std::vector<float>>::from( theta*(X.diag()) );
    m_vectorSet = true;
    
}

void BoundComputationTheta::setBranchVar(int& branchVar, int n){
    // set branching variable
    switch (params.strategyBranchVar) {
        case 1: {
            // set random branching variable
            branchVar = std::rand() % n + 1;  // random branchvar (between 1 and n)
            break;
        }    

        case 2: {
            // position of min in vector
            int index = 0;
            int i = 0;
            double min = m_vector[0];
            for(std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); ++it){
                if(*it < min){
                    min = *it;
                    index = i;
                }
                ++i;
            }
            branchVar = index + 1;
            break;
        }

        case 3: {
            // position of max in vector
            int index = 0;
            int i = 0;
            double max = m_vector[0];

            for(std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); ++it){
                if(*it > max){
                    max = *it;
                    index = i;
                }
                ++i;
            }
            branchVar = index + 1;
            break;
        }

        case 4: {
            // closest to 0.5
            int index = 0;
            int i = 0;
            double min = abs(m_vector[0]-0.5);

            for (std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); ++it){
                if(abs(*it-0.5) < min){
                    min = abs(*it-0.5);
                    index = i;
                }
                ++i;
            }
            
            branchVar = index + 1;
            break;
        }

        default: // set to random
            branchVar = std::rand() % n + 1;
    }

}

// Burer-Monteiro heuristic
std::vector<int> BoundComputationTheta::getStableSetWithAO(Graph* g, int rank, int numResets, int maxSeconds){
    //rank has to be 1 or 2!
    Heuristic h;
    return h.getBestStableSetFound(g->n(), g->edges(), rank, numResets, maxSeconds); //11, 21, 25
}
