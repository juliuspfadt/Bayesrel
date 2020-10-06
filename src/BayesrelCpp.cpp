

#include <RcppArmadillo.h>
extern "C" {
#include "declarations.h"
}
#include "customsdp.h"


//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::dvec csdpArma(
              int n_p,
              int nconstraints_p,
              int nblocks_p,
              const arma::ivec& blocktypes_p,
              const arma::ivec& blocksizes_p,
              const Rcpp::List& C_p,
              const Rcpp::List& A_p,
              const arma::dvec& b_p,
              const arma::cube& car)
{

    struct blockmatrix C;
    struct constraintmatrix *constraints;
    struct blockmatrix X, Z;
    double *y, *b;
    double pobj, dobj;
    int status;
    arma::dvec out(car.n_slices);


    /*
     * setup C
     */
    C = blkmatrix_R2csdpArma(C_p);

    /*
     * setup constraints
     */
    constraints = constraints_R2csdpArma(A_p);

    /*
     * Allocate storage for RHS
     */
    b = double_vector_R2csdpArma(nconstraints_p,b_p);

    /*
     * Create an initial solution. This allocates space for X, y, and Z,
     * and sets initial values
     */
    initsoln(n_p,nconstraints_p,C,b,constraints,&X,&y,&Z);

    /*
     * Solve the problem
     */
    status = custom_sdpCpp(n_p,nconstraints_p,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj, car, out);

    /*
     * Grab the results
     */

    free_prob(n_p,nconstraints_p,C,b,constraints,X,y,Z);

    return out;
}
