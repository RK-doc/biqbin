/*
 * BoundComputationBundle.cpp
 *
 *  Created on: 20.09.2017
 *      Author: mesiebenhofe
 */

#include "BoundComputationBundle.h"
#include <math.h>       /* exp */

using namespace mosek::fusion;
using namespace monty;

BoundComputationBundle::BoundComputationBundle() {

}

BoundComputationBundle::~BoundComputationBundle() {

}

int BoundComputationBundle::getUpperBound(Graph* graph, double& upperBound, int& branchVar, int constValueUpperBound){
	int status = BoundComputationTheta::getUpperBound(graph, upperBound, branchVar, constValueUpperBound);
	if(status != -1){
		if(upperBound + OPTEPS < g_lowerBound - constValueUpperBound + 1){
			//prune -> don't need to start bundle!
			//variables are already set
#ifdef DEBUG
			if(params.printlevel > 2){
				std::cout << "theta function found upper bound -> prune" << std::endl;
			}
#endif
			return 0;
		}
	}

	m_vectorSet = false;

	int n = graph->n();

	if(n==0){
		branchVar = 0;
		return 0;
	}

	///variables for bundle
	std::vector<std::vector<int>> subgraphs;
	std::vector<int> tI;
	std::vector<int> bI;
	std::vector<std::vector<int> > yINeeded;
	std::vector< std::vector<std::vector<int>> > DI;




	/// create adjacency from list of edges
	int adjacency[n][n];
	for(int i = 0; i<n; i++){
		for(int j = 0; j<n; j++){
			adjacency[i][j] = 0;
		}
	}
	std::vector<Edge> edges = graph->edges();
	for(std::vector<Edge>::iterator ite = edges.begin() ; ite != edges.end(); ++ite){
		adjacency[ite->vertex1-1][ite->vertex2-1] = 1;
		adjacency[ite->vertex2-1][ite->vertex1-1] = 1;
	}



	///create subgraphs
	//random device for std::shuffle


	std::random_device rd;
	std::mt19937 g1(rd());

	for (std::map<int,int>::iterator it= params.numberExactSubgraphs.begin(); it!=params.numberExactSubgraphs.end(); ++it){
		int sizeOfSubgraph = it->first;
//		int sizeOfSubgraph = 5; //NEW

#ifdef DEBUG
		//print size of subgraphs
		if(params.printlevel>5){
			std::cout << "------------" << std::endl;
			std::cout << "size of subgraphs: " << sizeOfSubgraph << std::endl;
		}
#endif

		if(sizeOfSubgraph < n){
			//create vector of indices
			std::vector<int>indicesGraph;
			for(int i = 0; i < n; i++){
				indicesGraph.push_back(i+1);
			}

#ifdef DEBUG
			//print number of subgraphs for this subgraph-size
			if(params.printlevel>3){
				std::cout << "add " << it->second << " subgraph constraints" << std::endl;
			}
#endif
			//add subgraph of size (param in paramfile) (param in paramfile) times
			// <==> add subgraph constraints it->second times
			for(unsigned iteratorSubgraphs = 0; iteratorSubgraphs < it->second; iteratorSubgraphs++){
//			for(int i = 0; i < n - 4; i++){	//NEW
				//shuffle indices
				std::shuffle(indicesGraph.begin(), indicesGraph.end(), g1);

				//subgraph = first sizeSubgraph entries in indices
				std::vector<int> subgraph (indicesGraph.begin(),indicesGraph.begin() + sizeOfSubgraph) ;
//				int verticesSubgraph[] = {i+1,i+2,i+3,i+4,i+5};	  //NEW
//				std::vector<int> subgraph(verticesSubgraph, verticesSubgraph+sizeof(verticesSubgraph)/sizeof(int)); //NEW


				subgraphs.push_back(subgraph);

#ifdef DEBUG
				if(params.printlevel>5){
					std::cout << "Subgraph: ";
					for(auto s: subgraph){
						std::cout << s << " ";
					}
					std::cout << std::endl;
				}
#endif
				//get all stable sets
				std::vector<std::vector<int>> stableSets = getStableSets(subgraph, &adjacency[0][0], n);

#ifdef DEBUG
				if(params.printlevel > 6){
					int counter = 1;
					std::cout << "Stable Sets:" << std::endl;
					for(auto ss : stableSets){
						std::cout << counter++ << ")  ";
						for(auto s : ss){
							std::cout << s << " ";
						}
						std::cout << std::endl;
					}
				}
#endif

				///store number of stable sets in tI
				tI.push_back(stableSets.size());


#ifdef DEBUG
				if(params.printlevel>5){
					std::cout << "tI{i} = " << tI.back() << std::endl;
				}
#endif


				///get all stable set matrices
				std::vector<Matrix::t> stableSetMatrices = getAllStableSetMatrices(stableSets, sizeOfSubgraph);

#ifdef DEBUG
				if(params.printlevel > 6){
					std::cout << "Stable Set Matrices: " << std::endl;
					for(unsigned l = 0; l < stableSetMatrices.size(); l++){
						std::cout << l+1 << ")" << std::endl;
						for(int i = 0; i < stableSetMatrices[l]->numRows(); i++){
							for(int j = 0; j < stableSetMatrices[l]->numColumns(); j++){
								std::cout << stableSetMatrices[l]->get(i,j) << " ";
							}
							std::cout << std::endl;
						}
						std::cout << std::endl;
					}
				}
#endif


				///get b and yNeeded
				int b;
				std::vector<int> yNeeded;
				getbAndyNeeded(b, yNeeded, &adjacency[0][0], subgraph, n);

				///store b and yNeeded in bI and yINeeded
				bI.push_back(b);
				yINeeded.push_back(yNeeded);

#ifdef DEBUG
				if(params.printlevel > 5){
					std::cout << "b = " << bI.back() << std::endl;
					std::cout << "yNeeded = ";
					for(auto i : yINeeded.back()){
						std::cout << i << " ";
					}
					std::cout << std::endl;
				}
#endif


				///get D
				std::vector<std::vector<int> > D = getMatrixD(b, sizeOfSubgraph, yNeeded, stableSetMatrices);
				DI.push_back(D);

#ifdef DEBUG
				if(params.printlevel > 6){
					std::cout << "D: " << std::endl;
					for(auto i: DI.back()){
						for(auto elem: i){
							std::cout << elem << " ";
						}
						std::cout << std::endl;
					}


					//ACHTUNG: Matrix::t nicht mit get!
					//TODO: merken, dass wenn nicht nxn --> Index anders berechnen
				}
#endif

			}
		} //endif(sizeSubgraph < n)
	} //endfor(subgraphs of different size)


	//create Oracle
	Oracle_SS oracle = Oracle_SS(graph, subgraphs, tI, bI, yINeeded, DI);

	//create bdl_master_dual_soc_rotated
	bdl_master_dual_soc_rotated soc_rotated = bdl_master_dual_soc_rotated(tI, bI, DI);

	//initialize y
	std::vector<double> y; //size: sum(bI)
	//set y to (1,2,1,2,...)
	/*
	int counter = 0;
	for(std::vector<int>::iterator it = bI.begin() ; it != bI.end(); ++it){
			for(int i = 0; i < *it; i++){
				y.push_back(counter++%2 + 1);
			}
	}
	*/


	//print y:
	//std::cout << "y = (";
	//for(auto i : y){
	//	std::cout << i << " ";
	//}
	//std::cout << ")" << std::endl;

	//double fmax;
	//double h;
	//std::vector<double> g;

	//call oracle
	//oracle.oracle(y, fmax, h, g);

	//call next trial point
	/*
	if(subgraphs.size()!=0){
		std::vector<double> yTrial;

		std::vector<std::vector<double> > G;
		G.push_back(g);

		std::vector<double> e; //size of e = r
		e.push_back(0);
		double tau = 1/0.5;		//tau = 1/mue

		soc_rotated.get_new_trial_point(y, h, G, e, tau, yTrial);
	}
	*/

	double fOpt = n;
	std::vector<double> yOpt;
	std::vector<double> gOpt;
	std::vector<double> X;

	/*
	//TODO: delete
	std::cout << "bI = (";
	for(int i = 0; i < bI.size(); i++){
		std::cout << bI[i] << " ";
	}
	std::cout <<")" <<std::endl;
	//------
	 */

	if(subgraphs.size()!=0){
		int sum_bI = 0;
		for(auto i : bI){
			sum_bI += i;
		}

		for(int i = 0; i < sum_bI; i++){
			y.push_back(0);
		}

		//std::cout << "sumBI = " << sum_bI << std::endl;
		//std::cout << "n = " << n << std::endl;
		//std::cout << "(sum_bi)^1/3: " << std::pow(sum_bI, ((double)1)/3) << std::endl;

		//double muStart = ((double)12 * std::pow(sum_bI,((double)1)/3))/(double)n; //old mu
		double density = (double)graph->m()*2/(n*(n-1)); //density = #edges/#possibleEdges = m/(nChoose2)
		double muStart = 7.92*(265+sum_bI)*(3.96+exp(-10*n))*(0.04+exp(-8.25*density))/(std::pow(n,2));

		if(muStart == 0){
			muStart = 1/10000000; //TODO: prüfe, wie damit umgehen, dass mü 0 wird
		}

		status = bdl_dgr::bundle_for_exact_subgraph_constraints(oracle, soc_rotated, y, muStart, "SS", g_lowerBound - constValueUpperBound,fOpt, yOpt, gOpt, X);
		if(status == -1){
			return -1;
		}

		setVector(X, n);
		setBranchVar(branchVar, n);
	}
	else{
		branchVar = std::rand() % n + 1;  ///< random branchvar (between 1 and n)
	}


	upperBound = fOpt;
	return 0;

}

std::vector<Matrix::t> BoundComputationBundle::getAllStableSetMatrices(std::vector<std::vector<int> >& stableSets, int sizeOfSubgraph){
	std::vector<Matrix::t> stableSetMatrices;
	for(std::vector<std::vector<int> >::iterator it1 = stableSets.begin() ; it1 != stableSets.end(); ++it1){
		std::vector<int> stableSet = *it1;
		//sparse matrix with sizeof(stableSet)^2 entries
		int numberMatrixEntries = stableSet.size()*stableSet.size();
		auto rows = new_array_ptr<int,1>(numberMatrixEntries);
		auto cols = new_array_ptr<int,1>(numberMatrixEntries);
		auto values = new_array_ptr<double,1>(numberMatrixEntries);

		int i = 0;
		for(std::vector<int>::iterator it2 = stableSet.begin() ; it2 != stableSet.end(); ++it2){
			//S_{j,j} = 1 \forall j \in stableSet
			(*rows)[i] = *it2 - 1;
			(*cols)[i] = *it2 - 1;
			(*values)[i] = 1;
			i++;
			//S_{i,j} = 1 \forall i,j \in stableSet, i \neq j
			for(std::vector<int>::iterator it3 = it2 + 1; it3 != stableSet.end(); ++it3){
				(*rows)[i] = *it2 - 1;
				(*cols)[i] = *it3 - 1;
				(*values)[i] = 1;
				i++;
				(*rows)[i] = *it3 - 1;
				(*cols)[i] = *it2 - 1;
				(*values)[i] = 1;
				i++;
			}
		}
		stableSetMatrices.push_back(Matrix::sparse(sizeOfSubgraph,sizeOfSubgraph,rows,cols,values));
	}
	return stableSetMatrices;
}

void BoundComputationBundle::getbAndyNeeded(int& b, std::vector<int>& yNeeded, int* adjacency, std::vector<int>& subgraph, int n){
	int sizeOfSubgraph = subgraph.size();
	b = 0; /// b = numbers of nonzeros in yNeeded
	for(int i = 0; i < sizeOfSubgraph; i++){
		for(int j = i; j < sizeOfSubgraph; j++){
			if(adjacency[(subgraph[i]-1)*n + (subgraph[j]-1)]==1){
				yNeeded.push_back(0);
			}
			else{
				yNeeded.push_back(1);
				b++;
			}
		}
	}
}

std::vector<std::vector<int>> BoundComputationBundle::getMatrixD(int b, int sizeSubgraph, std::vector<int>& yNeeded, std::vector<Matrix::t>& stableSetMatrices){

	std::vector<std::vector<int> > D;
	std::vector<int> rowTemp;

	int t = stableSetMatrices.size();
	//auto arrayForMatrix = new_array_ptr<double,1>(t*b); //t*b entries
	//int arrayIndex = 0;

	bool diagEntry = true;

	for(int l = 0; l < t; l++){ // for all stable set matrices

		rowTemp.clear();
		int yIndex = 0;

		for(int i = 0; i < sizeSubgraph; i++){
			diagEntry = true;
			//(*arrayForMatrix)[arrayIndex] = - stableSetMatrices[l]->get(i,i); //index {i,i} factor = 1
			for(int j = i; j < sizeSubgraph; j++){
				if(diagEntry){
					rowTemp.push_back(stableSetMatrices[l]->get(i,i));
					diagEntry = false;
				} else if(yNeeded[yIndex] == 1){
					rowTemp.push_back(2*(stableSetMatrices[l]->get(i,j)));
				}

				//if(yNeeded[counter] == 1){
					//(*arrayForMatrix)[arrayIndex] += 2 * stableSetMatrices[l]->get(i,j);
					//arrayIndex++;
				//}
				yIndex++;
			}
		}
		D.push_back(rowTemp);
	}



	//DEBUG OUTPUT
	//std::cout << "array for Matrix: " << std::endl;
	//for(int i = 0; i < t*b; i++){
	//	std::cout << (*arrayForMatrix)[i] << std::endl;
	//}

	//Matrix::t m = Matrix::dense(t,b,arrayForMatrix);
	//std::cout << "Matrix: " << std::endl;
	//std::cout << m->toString() << std::endl;

	//get liefert falsche werte, wenn matrix nicht nxn sondern mxn
	//siehe test1.cpp
	/*
	int row = m->numRows();
	int col = m->numColumns();
	for(int i = 0; i < row; i++){
		for(int j = i; j < col+i; j++){
			std::cout << m->get(i,j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/

	return D;
}

//TODO: set Vector, set BranchVar
void BoundComputationBundle::setVector(std::vector<double> X, int n){
	///set vector (if no errors occurred)
	//X is Lower Triangular Matrix (need entries 0,2,5,9,14,...)
	m_vector.clear();
	int index = 0;
	for(int i = 0; i<n; i++){
		m_vector.push_back(X[index]);
		index += i+2;
		//m_vector.push_back((*(X->index(0,i+1)->level()))[0]);
	}
	m_vectorSet = true;
}

void BoundComputationBundle::searchForGoodSubgraphs(std::vector<double>& X, std::vector<double>& B, std::vector<int>& bestFoundSubgraph){




}
