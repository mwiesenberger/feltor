#pragma once
#include "projection.h"
#include "gridX.h"

/*!@file

  @brief Creation of projection matrices for X point topology
 */
namespace dg{
namespace create{
///@addtogroup interpolation
///@{

/**
 * @brief Create a projection between two grids with X-point topology
 *
 * The projection is the one between the grids without X-point topology
 * <tt>g_new.grid()</tt> and <tt>g_old.grid()</tt>
 * @param g_new The new (coarse) grid
 * @param g_old The old (fine) grid
 * @copydoc hide_method
 * @return Projection matrix
 * @sa dg::create::projection(const aRealTopology<real_type,Nd>&,const aRealTopology<real_type,Nd>&,std::string)
 */
template<class real_type>
dg::SparseMatrix< int, real_type, thrust::host_vector> projection( const RealGridX1d<real_type>& g_new, const RealGridX1d<real_type>& g_old,std::string method = "dg") {
    return projection(g_new.grid(), g_old.grid(),method);
}

///@copydoc projection(const RealGridX1d<real_type>&,const RealGridX1d<real_type>&,std::string)
template<class real_type>
dg::SparseMatrix< int, real_type, thrust::host_vector> projection( const aRealTopologyX2d<real_type>& g_new, const aRealTopologyX2d<real_type>& g_old,std::string method = "dg") {
    return projection(g_new.grid(), g_old.grid(),method);
}

///@copydoc projection(const RealGridX1d<real_type>&,const RealGridX1d<real_type>&,std::string)
template<class real_type>
dg::SparseMatrix< int, real_type, thrust::host_vector> projection( const aRealTopologyX3d<real_type>& g_new, const aRealTopologyX3d<real_type>& g_old,std::string method = "dg") {
    return projection(g_new.grid(), g_old.grid(),method);
}

///@}

}//namespace create
}//namespace dg
