//implementation of SimplicialCholeskyhandler
#include <Eigen/Sparse>
#include <vector>

#include "cusp_eigen.h"

Eigen::SparseMatrix<double, Eigen::RowMajor, int> convert( const cusp::coo_matrix<int, double, cusp::host_memory>& cm)
{
    typedef Eigen::Triplet<double> T;
    //first generate a vector of Eigen Triplets
    //thrust::host_vector<T> triplets( cm.num_entries);
    std::vector<T> triplets( cm.num_entries);
    for( unsigned i=0; i<cm.num_entries; i++)
        triplets[i] = T( cm.row_indices[i], cm.column_indices[i], cm.values[i]);
    //now construct the Eigen matrix from triplets (will even sort and reduce triplets)
    Eigen::SparseMatrix<double, Eigen::RowMajor, int> em( cm.num_rows, cm.num_cols);
    em.setFromTriplets( triplets.begin(), triplets.end());
    return em;
}

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int> EMatrix;
typedef Eigen::SimplicialLDLT<EMatrix> SimplicialCholeskyImpl;

namespace dg
{
SimplicialCholesky::SimplicialCholesky(): pImpl( new SimplicialCholeskyImpl) {}
SimplicialCholesky::SimplicialCholesky(const HMatrix& matrix): 
    pImpl( new SimplicialCholeskyImpl( convert(matrix))) {}
bool SimplicialCholesky::compute( const HMatrix& matrix) 
{
    pImpl->compute( convert( matrix));
    if( pImpl->info() != Succeeded ) return false;
    return true;
}
bool SimplicialCholesky::solve( double *x, const double* b, unsigned N) 
{
    if( x == b)
    {
        Eigen::Map< VectorXd> xmap( x, N);
        xmap = pImpl->solve( xmap);
    }
    else
    {
        Eigen::Map< VectorXd> xmap( x, N);
        Eigen::Map< const VectorXd> bmap( b, N);
        xmap = pImpl->solve( xmap);
    }
    if( pImpl->info() != Succeeded) return false;
    return true;
}
} //namespace dg



