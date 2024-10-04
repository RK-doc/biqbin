/*
 * bdl_master_dual_soc_rotated.cpp
 *
 *  Created on: 04.10.2017
 *      Author: mesiebenhofe
 */

#include "bdl_master_dual_soc_rotated.h"
#include <time.h>	//* clock_t, clock, CLOCKS_PER_SEC */

using namespace mosek::fusion;
using namespace monty;

bdl_master_dual_soc_rotated::bdl_master_dual_soc_rotated(std::vector<int>& tI, std::vector<int>& bI,
		std::vector<std::vector<std::vector<int>>>& DI) {

	m_tI = tI;
	m_bI = bI;
	m_DI = DI;

#ifdef DEBUG
	if(params.printlevel > 5){
		std::cout << "----------------------------------" << std::endl;
		std::cout << "create dual_soc_rotated: " << std::endl;
		std::cout << "create constant constraint matrix" << std::endl;
	}
#endif

	///generate constant constraint matrix
	//n = 1 + 2q + sum(bI)
	//m = sum(tI) + 2q + sum(bI)
	//# entries = sum(bI_i * tI_i) + sum(bI) +				//thrid block
	//				+ q +									//second block
	//				+ sum(tI)								//first block

	///get number of matrix entries:
	m_sum_bI = 0;
	m_sum_tI = 0;
	int sum_tIbI = 0;

	//TODO: if programmed correct not needed
	if(tI.size() != bI.size()){
		std::cerr << "error: dimension of bI and tI not equal!" << std::endl;
	}
	if(tI.size() != DI.size()){
		std::cerr << "error: dimension of DI and tI not equal!" << std::endl;
	}

	m_q = tI.size(); //number of exact subgraphs

	//calculate sum(bI), sum(tI), <bI,tI>
	for(int i = 0; i < m_q; i++){
		m_sum_bI += bI[i];
		m_sum_tI += tI[i];
		sum_tIbI += bI[i] * tI[i];
	}

	int numberOfMatrixEntries = sum_tIbI + m_sum_bI + m_q + m_sum_tI;

	///set entries of matrix

	auto rows = new_array_ptr<int,1>(numberOfMatrixEntries);
	auto cols = new_array_ptr<int,1>(numberOfMatrixEntries);
	auto values = new_array_ptr<double,1>(numberOfMatrixEntries);

	int index_arrayForMatrix = 0;

	//constraint type: X1I = 1/2 \forall I \in J
	int row = m_q + 1;
	int col = m_sum_tI;
	for(int i = 0; i < m_q; i++){
		(*rows)[index_arrayForMatrix] = row;
		(*cols)[index_arrayForMatrix] = col;
		(*values)[index_arrayForMatrix] = 1;

		row++;
		col += bI[i] + 2;
		index_arrayForMatrix++;
	}

	//constraint type: \beta^I \in \Delta_t_I
	row = 1;
	col = 0;
	for(int i = 0; i<m_q; i++){
		for(int j = 0; j < tI[i]; j++){
			(*rows)[index_arrayForMatrix] = row;
			(*cols)[index_arrayForMatrix] = col;
			(*values)[index_arrayForMatrix] = 1;

			col++;
			index_arrayForMatrix++;
		}
		row++;
	}


	//constraint type: xI = GI(\alpha) + DI'(\betaI)
	row = 1 + 2*m_q;
	col = 0;
	int col2 = m_sum_tI + 2;

	for(int i = 0; i < m_q; i++){
		for(int j = 0; j < bI[i]; j++){
			for(int k = 0; k < tI[i]; k++){
				(*rows)[index_arrayForMatrix] = row;
				(*cols)[index_arrayForMatrix] = col;
				(*values)[index_arrayForMatrix] = DI[i][k][j]; //DI'_i{j,k} = DI_i{k,j}

				index_arrayForMatrix++;
				col++;
			}
			(*rows)[index_arrayForMatrix] = row;
			(*cols)[index_arrayForMatrix] = col2;
			(*values)[index_arrayForMatrix] = -1;

			index_arrayForMatrix++;
			col2++;
			row++;
			col -= tI[i];

		}
		col += tI[i];
		col2 += 2;
	}

	constConstraintMatrix = Matrix::sparse(1+2*m_q+m_sum_bI, m_sum_tI+2*m_q+m_sum_bI, rows, cols, values);

#ifdef DEBUG
	if(params.printlevel > 6){
		std::cout << "const matrix: " << std::endl;
		std::cout << Matrix::dense(constConstraintMatrix)->toString() << std::endl;
	}
#endif

	///generate left side of constraint (b)

	b = new_array_ptr<double,1>(1 + 2*m_q + m_sum_bI);

	for(int i = 0; i < m_q+1; i++){
		(*b)[i] = 1;
	}

	for(int i = 0; i < m_q; i++){
		(*b)[m_q+1+i] = 0.5;
	}

	for(int i = 0; i < m_sum_bI; i++){
		(*b)[1 + 2*m_q + i] = 0;       //Todo: maybe this is not necessary because entries are 0 in general
	}


#ifdef DEBUG
	if(params.printlevel > 6){
		std::cout << "b = ( " ;
		for(int i = 0; i < 1 + 2*m_q + m_sum_bI; i++){
			std::cout << (*b)[i] << ", ";
		}
		std::cout << ")" << std::endl;
	}
#endif


}

bdl_master_dual_soc_rotated::~bdl_master_dual_soc_rotated() {

}

int bdl_master_dual_soc_rotated::get_new_trial_point(std::vector<double>& yCenter, std::vector<std::vector<double>>& G,
														std::vector<double>& e, double tau, std::vector<double>& yTrial,
														std::vector<double>& alpha, double& timeQP){


#ifdef DEBUG
	if(params.printlevel > 5){
		std::cout << "---------------------------------" << std::endl;
		std::cout << "get new trial point" << std::endl;

		std::cout << "yCenter = (";
		for(auto i: yCenter){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;

		std::cout << "e = (";
		for(auto i: e){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;

		std::cout << "tau = " << tau << std::endl;

	}
#endif

	int r = G.size(); //number of bundle elements

	///generate submatrix from matrix A (for constraint Ax = b)
	Matrix::t subMatrixOfA;
	int numberOfMatrixEntries = r*(m_sum_bI + 1);
	auto rows = new_array_ptr<int,1>(numberOfMatrixEntries);
	auto cols = new_array_ptr<int,1>(numberOfMatrixEntries);
	auto values = new_array_ptr<double,1>(numberOfMatrixEntries);

	int index_arrayForMatrix = 0;


	//constraint type sum(\alpha) = 1
	for(int i = 0; i < r; i++){
		(*rows)[index_arrayForMatrix] = 0;
		(*cols)[index_arrayForMatrix] = i;
		(*values)[index_arrayForMatrix] = 1;
		index_arrayForMatrix++;
	}

	//constraint type xI = GI(\alpha) + DI'(\betaI)
	/*
	int rowindexMatrix = 1 + 2*m_q;
	int rowindexG = 0;
	for(int i = 0; i < r; i++){ //for each element in bundle
		for(int j = 0; j < m_q; j++){ //for each subgraph
			for(int k = 0; k < m_bI[j]; k++){
				(*rows)[index_arrayForMatrix] = rowindexMatrix;
				(*cols)[index_arrayForMatrix] = i; //[0,r-1]
				(*values)[index_arrayForMatrix] = G[i][rowindexG];

				rowindexMatrix++;
				rowindexG++;
				index_arrayForMatrix++;
			}

		}
		rowindexG = 0;
		rowindexMatrix -= m_sum_bI;
	}
	*/
	//easier?

	int rowIndexMatrix = 1 + 2*m_q;
	for(int j = 0; j < r; j++){
		for(int i = 0; i < m_sum_bI; i++){
			(*rows)[index_arrayForMatrix] = rowIndexMatrix + i;
			(*cols)[index_arrayForMatrix] = j;
			(*values)[index_arrayForMatrix] = G[j][i];

			index_arrayForMatrix++;
		}
		rowIndexMatrix = 1 + 2*m_q;
	}


	subMatrixOfA = Matrix::sparse(1+2*m_q+m_sum_bI, r, rows, cols, values);

#ifdef DEBUG
	if(params.printlevel > 6){
		std::cout << "submatrix of A: " << std::endl;
		std::cout << Matrix::dense(subMatrixOfA)->toString() << std::endl;
	}
#endif

	///generate matrix C
	Matrix::t C;
	rows = new_array_ptr<int,1>(r + m_sum_tI + m_q);
	cols = new_array_ptr<int,1>(r + m_sum_tI + m_q);
	values = new_array_ptr<double,1>(r + m_sum_tI + m_q);

	int colIndexMatrix = 0;
	index_arrayForMatrix = 0;

	double entry;


	//sum(e_j,\alpha_j)
	for(int i = 0; i < r; i++){
		(*rows)[index_arrayForMatrix] = 0;
		(*cols)[index_arrayForMatrix] = colIndexMatrix;
		(*values)[index_arrayForMatrix] = e[i];

		colIndexMatrix++;
		index_arrayForMatrix++;
	}

	//- sum<DI,betaI>
	int yIndex = 0;
	for(int i = 0; i < m_q; i++){			//for each subgraph
		for(int j = 0; j < m_tI[i]; j++){
			(*rows)[index_arrayForMatrix] = 0;
			(*cols)[index_arrayForMatrix] = colIndexMatrix;

			//-DIi(yCenter_Ii)
			entry = 0;
			for(int k = 0; k < m_bI[i]; k++){
				entry -= m_DI[i][j][k] * yCenter[yIndex + k];
			}

			(*values)[index_arrayForMatrix] = entry;

			index_arrayForMatrix++;
			colIndexMatrix++;
		}
		yIndex += m_bI[i];
	}

	// sum tau/2 x2I
	for(int i = 0; i < m_q; i++){
		colIndexMatrix++;
		(*rows)[index_arrayForMatrix] = 0;
		(*cols)[index_arrayForMatrix] = colIndexMatrix;
		(*values)[index_arrayForMatrix] = tau/2;

		index_arrayForMatrix++;
		colIndexMatrix += 1 + m_bI[i];
	}

	C = Matrix::sparse(1, r + m_sum_tI + 2*m_q + m_sum_bI,rows,cols,values);

#ifdef DEBUG
	if(params.printlevel > 6){
		std::cout << "matrix C: " << std::endl;
		std::cout << Matrix::dense(C)->toString() << std::endl;
	}
#endif

	///create mosek model
	Model::t M = new Model(); auto _M = finally([&]() { M->dispose(); });
	Variable::t X = M->variable(r + m_sum_tI + 2*m_q + m_sum_bI);

	M->objective(ObjectiveSense::Minimize, Expr::mul(C, X));

	M->constraint(X->slice(0,r+m_sum_tI), Domain::greaterThan(0.0)); //\alpha \geq 0, \beta \geq 0

	//constraint xI1 \cdot xI2 \geq ||xI||^2 (q constraints)
	int startIndex = r + m_sum_tI;
	int endIndex = startIndex + 2 + m_bI[0];
	for(int i = 0; i < m_q - 1; i++){
		M->constraint(X->slice(startIndex,endIndex), Domain::inRotatedQCone());
		startIndex = endIndex;
		endIndex += 2 + m_bI[i+1];
	}
	//last pair of variables (cannot access m_bI[q])
	M->constraint(X->slice(startIndex,endIndex), Domain::inRotatedQCone());

	//constraint Ax = b
	Expression::t expr1 = Expr::mul(subMatrixOfA, X->slice(0,r));
	Expression::t expr2 = Expr::mul(constConstraintMatrix, X->slice(r, r+m_sum_tI + 2*m_q + m_sum_bI));
	M->constraint(Expr::add(expr1, expr2), Domain::equalsTo(b));



	clock_t time;
	///solve
	try{
		time = clock();
		M->solve();
		time = clock() - time;
		timeQP = ((double)time)/CLOCKS_PER_SEC;
	}
	catch(...){
		std::cerr << "mosek failed to solve SOC problem" << std::endl;
		return -1;
	}


	///check solution status
	if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal || M->getPrimalSolutionStatus() == SolutionStatus::Optimal){

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
			std::cout << "objective function value = " << M->primalObjValue() << std::endl;
		}
#endif
	} else{
		counterBoundComputationFailed++;
		std::cerr << "Problem solution not optimal (Oracle)" << std::endl;
		return -1;
	}


	///calculate yTrial (store alpha and beta)
	yTrial.clear();
	std::vector<double> beta;
	alpha.clear();
	for(int i = 0; i < r; i++){
		alpha.push_back((*(X->index(i)->level()))[0]);	//solution at index (i,0)
	}
	for(int i = 0; i < m_sum_tI; i++){
		beta.push_back((*(X->index(r+i)->level()))[0]); //solution at index (r+i,0)
	}

	/*
	std::vector<double> xAusgabe;
	for(int i = 0; i < r + m_sum_tI + 2*m_q + m_sum_bI; i++){
		xAusgabe.push_back((*(X->index(i)->level()))[0]);
	}
	std::cout << "x = [";
	for(auto i : xAusgabe){
		std::cout << i << ",";
	}
	std::cout << "]" << std::endl;
	*/

	//yI_i = yCenterI_i - tau (sum_{j=1}^{r} GI_{i,j}alpha_j + sum_{j=1}^{tI} DI_{j,i}beta_j)
	yIndex = 0;
	int betaIndex = 0;
	double tempValue;
	for(int g = 0; g < m_q; g++){			//for each subgraph
		for(int i = 0; i < m_bI[g]; i++){	//for each index i of yI
			tempValue = 0;
			for(int j = 0; j < r; j++){
				tempValue += G[j][yIndex + i] * alpha[j];
			}
			for(int j = 0; j < m_tI[g]; j++){
				tempValue += m_DI[g][j][i] * beta[betaIndex + j];
			}

			yTrial.push_back(yCenter[yIndex + i] - tau*tempValue);

		}
		yIndex += m_bI[g];
		betaIndex += m_tI[g];
	}


#ifdef DEBUG
if(params.printlevel > 5){
	std::cout << "new trial point = (";
	for(auto i: yTrial){
		std::cout << i << ", ";
	}
	std::cout <<  ")" << std::endl;

	std::cout << "alpha = (";
	for(auto i: alpha){
		std::cout << i << ", ";
	}
	std::cout << ")" << std::endl;

	std::cout << "beta = (";
	for(auto i: beta){
		std::cout << i << ", ";
	}
	std::cout << ")" << std::endl;

}
#endif

	return 0;


}
