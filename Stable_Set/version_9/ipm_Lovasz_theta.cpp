/***************************************************************************************************
* Interior point method for computing Lovasz theta number
* C++ serial version
* 
* solves: max <J,X> s.t. A_G*X(:) = b; X psd
*       (means X(i,j) = 0 for all [ij] in E(G), tr(X)=1) 
*
* data:
* the edge list of the graph is given in e
*
*
* Code written by Timotej Hrga (timotej.hrga@gmail.com)
* Based on MATLAB code from paper
* Helmberg, Rendl, Vanderbei, Wolkowisz, An Interior-Point Method for Semidefinite Programming
****************************************************************************************************
*/

#include <iostream>
#include <armadillo>
#include <vector>
#include "Structures.h"

using namespace std;
using namespace arma;


/* ------------------------- MAIN ROUTINE ------------------------- */
// NOTE:  armadillo starts indexing at 0!

void ipm_Lovasz_theta(std::vector<Edge> &edge_list, int n, double &uBound, mat &X){

	/* ------------------------- Data ------------------------- */
    // number of edges
    int m1 = edge_list.size();
    
    // edgelist taken apart
    uvec I(m1,fill::zeros);
    uvec J(m1,fill::zeros);
    
    for (int i = 0; i < m1; ++i){
        I(i) = edge_list[i].vertex1;
        J(i) = edge_list[i].vertex2;
    }
    
	// number of constraints is |E| + 1
	int m = m1 + 1;
                 
	// right hand side
	vec b(m,fill::zeros);
	b(m-1) = 1;           // armadillo starts indexing at 0!


	/* ------------------------ INITIALIZE: Set starting matrices and starting values ------------------------ */
	vec y(m,fill::zeros);
	y(m-1) = n + 1;
	mat Z = (n+1) * eye<mat>(n,n) - ones<mat>(n,n);
	
	// initial primal and dual objective
	double phi = dot(b,y);
	double psi = sum(sum(X));

	// initial complementarity
	double mu = dot(Z,X) / (2.0*n);

	// default tolerance
	int digits = 6;
	
	// iteration count
	int iter = 1;     

       
 	/* -------------------------- Loop -----------------------------*/


	// while duality gap too large
	while ( (phi - psi) > ( (1 > std::abs(phi)) ? 1 : std::abs(phi) ) * pow(10,-digits) ) {
	
    
    	mat Zi = inv_sympd(Z);                    // explicitly compute inv(Z)
    	Zi = ( Zi + Zi.t() ) / 2.0;

    	// construct matrix M for system My = rhs
    	mat M(m,m,fill::zeros);
    
    	// (m,m) element
    	M(m-1,m-1) = dot(Zi,X);
    
   		// the rest of the matrix M
   		// last row and last column
    	M( m-1, span(0,m1-1) ) = sum( Zi.cols(I-1) % X.cols(J-1) + Zi.cols(J-1) % X.cols(I-1) );
    	M( span(0,m1-1), m-1 ) = M( m-1, span(0,m1-1) ).t();

    	// submatrix (1:m1) x (1:m1)
    	M( span(0,m1-1), span(0,m1-1) ) = Zi(I-1,J-1) % X(J-1,I-1) + Zi(J-1,I-1) % X(I-1,J-1) + Zi(I-1,I-1) % X(J-1,J-1) + Zi(J-1,J-1) % X(I-1,I-1);


    	
   		// determine right hand side of linear system to get y
    	// we need to form A(.)
    
    	vec rhs(m,fill::zeros);                      // initialize
    	rhs(m-1) = trace(Zi);                  		 // last component is trace 
    	vec tmp = vectorise(Zi);                     // make tmp a long vector

    	rhs.subvec(0,m1-1) = tmp( (I-1)*n + J - 1) * 2.0; 
   	

   		// get dy
    	vec dy = solve(M, (mu*rhs - b), solve_opts::fast);   			 // solve for dy
   
    	
    	// now compute A^t(y)
        // back substitute for dZ = A^T(y)
    	vec Aty_vec(n*n, fill::zeros);

    	Aty_vec( (I-1)*n + J - 1) = dy(span(0,m1-1));

    	mat Aty = reshape(Aty_vec,n,n);
	    
	    mat dZ = Aty + Aty.t() + dy(m-1)*eye<mat>(n,n);
                      
    	mat dX = mu * Zi - X - Zi*dZ*X;         	// back substitute for dX 
    	dX = (dX + dX.t()) / 2.0;


    	// line search on primal: X  = X + alpha_p * dX  psd matrix
    	double alpha_p = 1;

    	mat R;
    	
    	bool isPosDef = chol(R, X + alpha_p*dX);		//  test if positive definite

    	
    	while (!isPosDef){
    		alpha_p = 0.8 * alpha_p;
    		isPosDef = chol(R, X + alpha_p*dX);
    	}
    
    	if (alpha_p < 1){
       		alpha_p = 0.95 * alpha_p;        			// stay away from boundary 
    	}

    	// line search on dual
    	double alpha_d = 1;

    	isPosDef = chol(R, Z + alpha_d*dZ);		//  test if positive definite
    
    	while (!isPosDef){
    		alpha_d = 0.8 * alpha_d;
    		isPosDef = chol(R, Z + alpha_d*dZ);
    	}

    	if (alpha_d < 1){
       		alpha_d = 0.95 * alpha_d;        			// stay away from boundary 
    	}

    
    	// update
    	X = X + alpha_p * dX;
    	y = y + alpha_d * dy;
    	Z = Z + alpha_d * dZ;
    	mu = dot(Z,X) / (2.0*n);
    
    	if ( (alpha_p + alpha_d) > 1.8 ){
    		mu = 0.5 * mu;								// speed up for long steps
    	} 			
    
    	// objective values
    	phi = dot(b,y);
    	psi = sum(sum(X));
    

    	// start new iteration
	    ++iter;    
	                  

	} // end while loop
	
	/* ----------------------------------------------------------------------*/
    
	uBound = phi;
}
