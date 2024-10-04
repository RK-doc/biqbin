/*
 * Oracle_SS.cpp
 *
 *  Created on: 26.09.2017
 *      Author: mesiebenhofe
 */

#include "Oracle_SS.h"

#include <limits> //std::numeric_limits<double>::min()

//#include <time.h>

using namespace mosek::fusion;
using namespace monty;

Oracle_SS::Oracle_SS(Graph* graph, std::vector<std::vector<int>> exactSubgraphs, std::vector<int> tI,
		 std::vector<int> bI, std::vector<std::vector<int> > yINeeded, std::vector<std::vector<std::vector<int>>> DI) {
	n = graph->n();
	edges = graph->edges();

	m_exactSubgraphs = exactSubgraphs;
	m_tI = tI;
	m_bI = bI;
	m_yINeeded = yINeeded;
	m_DI = DI;


	///constraints type 1: X_{i,i} = x_{i}
	for(int i = 1; i<n+1; i++){
		///X_{0,i} = X_{i,0} = X_{i,i} <==>
		///-0.5 X_{0,i} - 0.5 X_{i,0} + X_{i,i} = 0 <==>
		/// <A,X> = 0,  a_{0,i} = -0.5, a_{i,0} = -0.5, a_{i,i} = 1
		auto rows = new_array_ptr<int,1>({0,i,i});
		auto cols = new_array_ptr<int,1>({i,0,i});
		auto values = new_array_ptr<double,1>({-0.5,-0.5,1});
		//auto A = Matrix::sparse(n+1,n+1,rows,cols,values);
		staticConstraints.push_back(Matrix::sparse(n+1,n+1,rows,cols,values));
		//M->constraint(Expr::dot(A,X),Domain::equalsTo(0.0));
	}

	///constraints type 2: X_{i,j} = 0  \forall (i,j) \in E <==>
	///0.5 X_{i,j} + 0.5 X_{j,i} = 0
	for (std::vector<Edge>::iterator it = edges.begin() ; it != edges.end(); it++ ){
		auto rows = new_array_ptr<int,1>({it->vertex1,it->vertex2});
		auto cols = new_array_ptr<int,1>({it->vertex2, it->vertex1});
		auto values = new_array_ptr<double,1>({0.5,0.5});
		//auto A = Matrix::sparse(n+1,n+1,rows,cols,values);
		staticConstraints.push_back(Matrix::sparse(n+1,n+1,rows,cols,values));
		//M->constraint(Expr::dot(A,X),Domain::equalsTo(0.0));
	}

	//indices to get submatrix of X
	//remove first column and row
	//store only the lower triangular matrix
	int numberOfEntries = (n*(n+1))/2;
	rowsSubX = new_array_ptr<int,1>(numberOfEntries);
	colsSubX = new_array_ptr<int,1>(numberOfEntries);
	int indexCounter = 0;
	for(int i = 0; i < n; i++){
		for(int j = 0; j <= i; j++){
			(*rowsSubX)[indexCounter] = i+1;
			(*colsSubX)[indexCounter] = j+1;
			indexCounter++;
		}
	}

}

Oracle_SS::~Oracle_SS() {

}

int Oracle_SS::oracle(std::vector<double> y, double& fmax, double& h, std::vector<double>& g, std::vector<double>& XMatrix){


	int sumBI = 0;
	for(std::vector<int>::iterator it = m_bI.begin(); it!= m_bI.end(); it++){
		sumBI += *it;
	}
	if(sumBI != y.size()){
		std::cerr << "[ERROR] Oracle: size of y and sum(bI) do not match! " << std::endl;
		return -1;
	}


	/************************************************************************************
	 *
	 *	y (sum(bI)x1 array) ... point, at which the oracle should be evaluated
	 *
	 *	exactSubgraphs (qx1 cell) ... list of Subgraphs
	 *		exactSubgraphs_{i} = i-th subgraph to be exact
	 *
	 *	n ... size of graph
	 *
	 *	bI (qx1 array) ... #equalities to hold for the exact subgraph constraint
	 *
	 *	DI (qx1 cell)
	 *		DI_{i} (tI(i)xbI(i) Matrix) ... needed in the max term of f
	 *
	 *	yINeeded (qx1 cell)
	 *		yINeeded_{i} (kI*(kI+1)/2 x 1 array)
	 *
	 ************************************************************************************/


	/**
	 *
	 * fmax(yJ) = sum_{I \in J}  max_{1 \leq i \leq t_i} [D_I(y_I)]_{i}
	 *
	 */

	int yIndex = 0;
	double fmaxTemp;
	double rowSum;
	fmax = 0;

	for(unsigned int indexSubgraph = 0; indexSubgraph < m_exactSubgraphs.size(); indexSubgraph++){ //for I \in J
		//DI(y)
		//DI ...  (tI)x(bI)
		fmaxTemp = std::numeric_limits<double>::min(); //max_{1\leq i \leq t_I} [DI(yI)]_i
		for(int i = 0; i < m_tI[indexSubgraph]; i++){
			//compute DI(y)_i
			rowSum = 0;
			for(int j = 0; j < m_bI[indexSubgraph]; j++){
				rowSum += m_DI[indexSubgraph][i][j] * y[yIndex + j];
			}
			if(rowSum > fmaxTemp){
				fmaxTemp = rowSum;
			}
		}

		fmax += fmaxTemp;

		yIndex += m_bI[indexSubgraph];

	}

#ifdef DEBUG
	if(params.printlevel > 5){
		std::cout << "----------------------------------" << std::endl;
		std::cout << "Oracle: " << std::endl;
		std::cout << "fmax = " << fmax << std::endl;
	}
#endif


	/**
	 *   h(yJ)= max <eye(n,n) - sum_{I \in J} PI'(MI(yI)), X>
	 *   		s.t diag(X) = x
	 *   			X_ij = 0 \forall (i,j) \in E
	 *   			X - xx' psd
	 */




	/// create matrix (nxn) C = eye(n,n) - sum_{I \in J} PI'(MI(yI))

	auto arrayForMatrix = new_array_ptr<double,1>(n*n);

	// I_n
	for(int i = 0; i < n; i++){
		(*arrayForMatrix)[i*n+i] = 1;
	}
	// - sum_{I \in J} PI'(MI(yI))
	// PI'(MI(yI)) ... Matrix (nxn) aus Vektor der Größe bI
	yIndex = 0;

	for(unsigned int indexSubgraph = 0; indexSubgraph < m_exactSubgraphs.size(); indexSubgraph++){ // \forall I \in J
		//some helper variables
		int sizeOfSubgraph = m_exactSubgraphs[indexSubgraph].size();
		int indexYNeeded = 0;
		int counterY = 0;
		bool diagEntry = true;

		for(int i = 0; i < sizeOfSubgraph; i++){
			diagEntry = true;					// first entry in each row is the diag entry (add just once)

			for(int j = i; j < sizeOfSubgraph; j++){
				if(m_yINeeded[indexSubgraph][indexYNeeded] == 1){
					(*arrayForMatrix)[(m_exactSubgraphs[indexSubgraph][i]-1)*n + (m_exactSubgraphs[indexSubgraph][j]-1)] -= y[yIndex + counterY];

					if(diagEntry){
						diagEntry = false;
					} else {
						(*arrayForMatrix)[(m_exactSubgraphs[indexSubgraph][j]-1)*n + (m_exactSubgraphs[indexSubgraph][i]-1)] -= y[yIndex + counterY];
					}

					counterY++;
				}
				indexYNeeded++;
			}
		}

		yIndex += m_bI[indexSubgraph];
	}


	auto C = Matrix::dense(n,n,arrayForMatrix);

#ifdef DEBUG
	if(params.printlevel > 6){
		std::cout << "Matrix C =  " << C->toString() << std::endl;
	}
#endif

	///create Mosek model
	Model::t M = new Model(); auto _M = finally([&]() { M->dispose(); });

	///var X ... (n+1)x(n+1) matrix, psd
	Variable::t X = M->variable(Domain::inPSDCone(n+1));

	///objective function:
	//max <C,X_{2,n+1}{2,n+1}>   (X in first row/column is x resp. x')
	M->objective(ObjectiveSense::Maximize, Expr::dot(C, X->slice(new_array_ptr<int,1>({1,1}), new_array_ptr<int,1>({n+1,n+1}) ) ));

	///add constraints:
	///constraint 1: X_{1,1} = 1
	M->constraint(X->index(0,0)->asExpr(), Domain::equalsTo(1.0));

	///add constraints
	///type 1: X_{i,i} = x_{i} and
	///type 2: X_{i,j} = 0  \forall (i,j) \in E
	for(int i = 0; i < staticConstraints.size(); i++){
		M->constraint(Expr::dot(staticConstraints[i],X),Domain::equalsTo(0.0));
	}

	///solve SDP
	//clock_t time = clock();
	try{
		M->solve();
	}
	catch(...){
		std::cerr << "mosek failed to solve SOC problem" << std::endl;
		return -1;
	}
	//time = clock() - time;

	//std::cout << "time to solve oracle in mosek = " <<  ((float)time)/CLOCKS_PER_SEC << std::endl;


	///check solution status
	if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal || M->getPrimalSolutionStatus() == SolutionStatus::Optimal){
		///set value of h(yJ)
		h = M->primalObjValue();

		if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal){
			counterMosekSolutionNearOptimal++;
		}

#ifdef DEBUG
		if(params.printlevel>5){
			if(M->getPrimalSolutionStatus() == SolutionStatus::Optimal){
				std::cout << "Problem solution status optimal" << std::endl;
			}
			else if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal){
				std::cout << "Problem solution status near optimal" << std::endl;
			}
			std::cout << "h = " << h << std::endl;
		}
#endif
	} else{
		counterBoundComputationFailed++;
		std::cerr << "Problem solution not optimal (Oracle)" << std::endl;
		return -1;
	}



	/**
	 *
	 * gI = -MI'PI(X*)
	 * gJ = (gI)_{I \in J}
	 *
	 */
	g.clear();
	bool diagEntry = true;
	for(unsigned int indexSubgraph = 0; indexSubgraph < m_exactSubgraphs.size(); indexSubgraph++){
		int sizeOfSubgraph = m_exactSubgraphs[indexSubgraph].size();
		int indexYNeeded = 0;
		for(int i = 0; i < sizeOfSubgraph; i++){
			diagEntry = true; //first entry is diagentry
			for(int j = i; j < sizeOfSubgraph; j++){
				if(m_yINeeded[indexSubgraph][indexYNeeded] == 1){
					if(diagEntry){
						g.push_back( - (*(X->index(m_exactSubgraphs[indexSubgraph][i],m_exactSubgraphs[indexSubgraph][j])->level()))[0]); //exclude 1st row/column (x resp. x')
						diagEntry = false;
					}
					else{
						g.push_back( -2* (*(X->index(m_exactSubgraphs[indexSubgraph][i],m_exactSubgraphs[indexSubgraph][j])->level()))[0]); //exclude 1st row/column (x resp. x')
					}
				}
				indexYNeeded++;
			}
		}
	}

#ifdef DEBUG
	if(params.printlevel > 5){
		std::cout << "g = (";
		for(auto i: g){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;
	}
#endif

	XMatrix.clear();
	auto arr = *(Var::reshape(X->pick(rowsSubX,colsSubX),(n*(n+1))/2,1)->level());
	for(int i = 0; i < (n*(n+1))/2; i++){
		XMatrix.push_back(arr[i]);
	}


#ifdef DEBUG
	if(params.printlevel > 5){
		std::cout << "X = (";
		for(auto i: XMatrix){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;
	}
#endif


}
