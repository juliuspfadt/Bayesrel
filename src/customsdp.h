

#ifndef customsdp_h
#define customsdp_h
#include <RcppArmadillo.h>


arma::ivec int_vector_csdp2RArma(int n, int *y);

arma::dvec double_vector_csdp2RArma(int n, double *y);

int * int_vector_R2csdpArma(int n, const arma::ivec& y);

double * double_vector_R2csdpArma(int n, const arma::dvec& y);

struct blockmatrix blkmatrix_R2csdpArma(const Rcpp::List& X);

Rcpp::List blkmatrix_csdp2RArma(const blockmatrix& X);

struct constraintmatrix *constraints_R2csdpArma(const Rcpp::List& A);

int custom_sdpCpp(
    int n,
    int k,
    const blockmatrix& C,
    double *a,
    struct constraintmatrix *constraints,
    double constant_offset,
    struct blockmatrix *pX,
    double **py,
    struct blockmatrix *pZ,
    double *ppobj,
    double *pdobj);

void printBlockMat(const blockmatrix& C);

void printConstMat(const constraintmatrix& S, int k);

/*
 A block record describes an individual block within a matrix.
 */

struct blockrec2 {
    union blockdatarec data;
    enum blockcat blockcategory;
#ifndef NOSHORTS
    unsigned short blocksize;
#else
    int blocksize;
#endif
};

/*
 A block matrix contains an entire matrix in block diagonal form.
 */

struct blockmatrix2 {
    int nblocks;
    struct blockrec2 *blocks;
};
#endif /* customsdp_h */
