/***************************************************************************************************
 * Boundary point method, based on augmented Lagrangian for computing Lovasz theta number
 * C++ serial version
 * 
 * solves: max <J,X>  
 *		   subject to:	A*X = b   (means X(i,j) = 0 for all [ij] in E(G), tr(X)=1) 
 *						 X psd
 *						 X >= 0 (nonnegativity constraint)
 *
 * Data:
 * the edge list of the graph is given in std::vector<Edge> e and the size n of the (sub)graph
 *
 * 
 * NOTE: From this variant of SDP program for Lovasz theta function we do not obtain vector x
 * which is further used in deciding what the next branching variable is.
 * 
 * We get vector x by x = theta * diag(X), where X is optimal solution of above SDP.
 *
 *
 * For details see:
 * J. Povh, F. Rendl and A. Wiegele. A boundary point method to solve semidefinite programs.
 * Computing, 78(3):277-286, November 2006
 *
 * Code written by Timotej Hrga (timotej.hrga@gmail.com)
 * Based on MATLAB code theta_bp.m available at https://www.aau.at/mathematik/publikationen/software/
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

void bpm(std::vector<Edge> &edge_list, int n, double &uBound, mat &X)
{
    
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
    
    // all ones matrix (cost function)
	mat L(n,n,fill::ones); 
    
    // right hand side
    vec b(m,fill::zeros);
    b(m-1) = 1;
    
    
    /* ------------------- Auxiliary data ---------------------------*/
    // A*A', which is diagonal
	vec DA(m,fill::ones);
	DA *= 2;
	DA(m-1) = n;
	
	// identity matrix
	mat Id = eye(n,n); 
    
    
    /* ------------------------ INITIALIZE: Set starting matrices and starting values ------------------------ */
    mat Z = zeros<mat>(n,n);
    mat S = zeros<mat>(n,n);
	vec y = zeros<vec>(m);
	vec Aty_vec(n*n, fill::zeros);
	mat Aty;
    
    // sigma is tuned for random graphs
    double sigma;
    if (n <= 250){
        sigma = 0.1/n; // before 20.0
    }
    else {
        sigma = 0.05/n; // before 10.0 
    }
    
    // default tolerance and max-iter
    double tol = 1e-5;
    int max_iter = 5000;
    
    // iteration count
    int iter = 1;
   
    /* -------------------------- Outer loop -----------------------------*/
    
    double err;
    
    // variable to hold approximal value of SDP if algorithm reaches max_iter
    double approx_value;
    
    
    // MAIN LOOP
    while (iter <= max_iter) {
        
        mat temp = Z + S + L + X/sigma;
        
	    vec rhs(m,fill::zeros); 	// initialize
	    rhs(m-1) = trace(temp); 	// last component is trace
	    vec tmp = vectorise(temp); 	// make tmp a long vector
	    
	    rhs.subvec(0,m1-1) = tmp( (I-1)*n + J - 1) * 2.0;

	    rhs = rhs - b/sigma; 		// the final right hand side
        
        // now compute y
	    y = rhs/DA;   				// solve sigma*AA^Ty = sigma*A(L+Z) + A(X)-b by element division

	    // now compute A^t(y)
	 	Aty_vec((I-1)*n + J-1) = y.subvec(0,m1-1);
  
	    Aty = reshape(Aty_vec,n,n);
	    
	    Aty = Aty + Aty.t() + y(m-1)*Id;
        
        
        // compute S from A'y - Z - C - X/sigma
    	S = max(Aty - Z - L - X/sigma,zeros<mat>(n,n));

	    // compute W(y)
	    mat W = Aty - S - L - X/sigma; 
        
        // now compute projections to get Wp and Wn 
	    mat Wp(n,n,fill::zeros);
	    mat Wn(n,n,fill::zeros);

	    vec lam(n,fill::zeros);
	    mat eigvec(n,n,fill::zeros);
	    eig_sym(lam, eigvec, W);

	    // find indeces of positive elemets
	    uvec Ip = find(lam > 0);
	    int j = Ip.n_elem;
	    
	    if (j > n/2.0){           
	    	mat eigvec_poz(n,j,fill::zeros);
	    	for (int r = 0; r<j; r++){
		    	int ic = Ip(r); 
			    eigvec_poz.col(r) = eigvec.col(ic)*sqrt(lam(ic)); 
	    	}

		    Wp = eigvec_poz * eigvec_poz.t();     	// the projection
	    	Wp = (Wp+Wp.t())/2;   						// should be symmetric
	    	Wn = W - Wp;
	    }
	    else {
		    uvec In = find(lam < 0);
		    int jj = In.n_elem;
		    mat eigvec_poz(n,jj,fill::zeros);

	    for (int r = 0; r < jj; r++){    
		    int ic = In(r);
			    eigvec_poz.col(r) = eigvec.col(ic)*sqrt(-lam(ic)); 
		    }
	  
	    	Wn = -eigvec_poz * eigvec_poz.t();   	// the projection
	    	Wn = (Wn+Wn.t())/2;   						// should be symmetric
	    	Wp = W - Wn; 
	    }
        
        // determine V = X and Z
        Z = Wp;
        X = -sigma * Wn;
        
        // compute inner error, even though we iterate only once anyway
	    temp = X; 
	    rhs.zeros(m,1);                     
	    rhs(m-1) = trace(temp);               	// last component is trace
    
	    tmp.zeros(m,1);
	    tmp = vectorise(temp);                  	// make tmp a long vector

	    rhs.subvec(0,m1-1) = tmp( (I-1)*n + J - 1) * 2.0; 
	
	    vec res_p = rhs - b;                     // A(X)-b

	    double err_p = norm(res_p);              // inner error
        
        
        // compute outer error (relative dual error)
	    mat res_d = Z + S - Aty + L;
	    double err_d = norm(res_d,"fro"); 

	    // compute primal and dual objective
	    double primal = sum(sum(X));
	    double dual = dot(b,y);

        // update sigma
        if (err_d > 5000 * err_p) {
            sigma *= 1.005;
        }
        
        iter++;
        
        err = max( err_p/2, err_d/(n+1) );
        
        approx_value = (primal > dual) ? primal : dual;
        
        if (err < tol) {
            uBound = primal;
            break;
        }
        
    } // end while loop
    
    uBound = approx_value;  // if max iterations reached

	// safe upper bound
	Z = Aty - L - S;
	vec lam(n,fill::zeros);
	mat eigvec(n,n,fill::zeros);
	eig_sym(lam, eigvec, Z);

	if ( lam(0) < 0 )
         uBound = y(m-1) - lam(0);
	
    
}
