

#include <RcppArmadillo.h>
extern "C" {
#include "declarations.h"
}
#include "customsdp.h"


//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double csdpArma(
              int n_p,
              int nconstraints_p,
              int nblocks_p,
              const arma::ivec& blocktypes_p,
              const arma::ivec& blocksizes_p,
              const Rcpp::List& C_p,
              const Rcpp::List& A_p,
              const arma::dvec& b_p)
{

    struct blockmatrix C;
    struct constraintmatrix *constraints;
    struct blockmatrix X, Z;
    double *y, *b;
    double pobj, dobj, ret;
    int status;


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
    status = custom_sdpCpp(n_p,nconstraints_p,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj);

    /*
     * Grab the results
     */

    /*
     * Grab X
     */
    Rcpp::List X_p = blkmatrix_csdp2RArma(X);

    /*
     * Grab Z
     */
    Rcpp::List Z_p = blkmatrix_csdp2RArma(Z);

    /* Copy y */
    arma::dvec y_p = double_vector_csdp2RArma(nconstraints_p, y);

    free_prob(n_p,nconstraints_p,C,b,constraints,X,y,Z);

    y_p(0) = 0;
    ret = arma::accu(y_p);
    return ret;
}
