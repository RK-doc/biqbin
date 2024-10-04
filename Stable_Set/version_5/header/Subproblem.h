/*
 * Subproblem.h
 *
 *  Created on: 09.08.2017
 *      Author: mesiebenhofe
 */

#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_

#include "Graph.h"
#include "LowerBoundComputation.h"
#include "UpperBoundComputation.h"

extern LowerBoundComputation* g_lbndComp;
extern UpperBoundComputation* g_ubndComp;
extern int g_lowerBound;
extern Params params;
extern int counterBoundComputationFailed;

class Subproblem {
protected:
	Graph m_graph;
	double m_upperBound; //local upper bound
	int m_depth; //depth of subproblem in BaB tree
	int m_constValue; //add to objective function
	std::set<int> m_stableSet;
	int m_nextBranchVar;
	
	


public:
	//Timotej
	// added branching variable vector and in/out vector --> needs to be public
	//std::vector<int> encode_subproblem;						// contains variables on which we branched
	//std::vector<int> encode_in_or_out; 	// contains 0 and 1 (in or out)

	Subproblem(Graph& g);

	// Timotej
	Subproblem(std::vector<int>& vertices, std::vector<Edge>& edges, int depth, double upper_bound, int constValue, std::set<int>& stableSet);

	Subproblem(Subproblem& parent, int branchVar, int inSet);
	int computeUpperBound(); //returns 0 if successful
	int computeLowerBound(int* globalSol, int n); //returns 0 if successful
	void printInfo1(int count);
	void printInfo2();
	void printInfo3();
	void printInfo4();

	virtual ~Subproblem();
	Subproblem(){};

	Graph graph(){return m_graph;};
	double upperBound(){return m_upperBound;};
	int nextBranchVar(){return m_nextBranchVar;};
	int depth(){return m_depth;}
	int get_constValue(){return m_constValue;};
	std::set<int> stableSet(){return m_stableSet;};
/*
	//Timotej added settter for upperBound
	void set_upperBound(double upper_bound)
	{
		m_upperBound = upper_bound;
	}
*/
private:
	//Subproblem();
	//Subproblem(const Subproblem&);
	//Subproblem& operator=(const Subproblem&);
};

#endif /* SUBPROBLEM_H_ */
