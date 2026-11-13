#ifndef _DG_TYPEDEFS_CUH_
#define _DG_TYPEDEFS_CUH_
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "sparseblockmat.h"
#include "sparsematrix.h"

/*! @file
  @brief Useful typedefs of commonly used types.
  */

/*! @namespace dg::x
 * @brief Typedefs that resolve to MPI types if \c MPI_VERSION is defined
 * and to shared memory types else listed in the \ref dg_tensor_typedefs_table "tensor typedefs table", the \ref dg_grid_typedefs_table "grid typedefs table"
 * and the \ref dg_geometry_typedefs_table "geometry typedefs table".
 *
 * Code written with these typedefs compiles with and without MPI.
 */
/*!@addtogroup typedefs
 *
 * @anchor dg_tensor_typedefs_table
 * The typedefs in namespace \c dg::x expand to
 *
 * | Name | without MPI | with MPI |
 * |------|-------------|----------|
 * | `dg::x::HVec` | `dg::HVec` = `thrust::host_vector<double>` | `dg::MHVec` = `dg::MPI_Vector<thrust::host_vector<double>>` |
 * | `dg::x::cHVec` | `dg::cHVec` = `thrust::host_vector<thrust::complex<double>>` | `dg::cMHVec` = `dg::MPI_Vector<thrust::host_vector<thrust::complex<double>>>` |
 * | `dg::x::fHVec` | `dg::fHVec` = `thrust::host_vector<float>` | `dg::fMHVec` = `dg::MPI_Vector<thrust::host_vector<float>>` |
 * | `dg::x::iHVec` | `dg::iHVec` = `thrust::host_vector<int>` | `dg::iMHVec` = `dg::MPI_Vector<thrust::host_vector<int>>` |
 * | `dg::x::DVec` | `dg::DVec` = `thrust::device_vector<double>` | `dg::MDVec` = `dg::MPI_Vector<thrust::device_vector<double>>` |
 * | `dg::x::cDVec` | `dg::cDVec` = `thrust::device_vector<thrust::complex<double>>` | `dg::cMDVec` = `dg::MPI_Vector<thrust::device_vector<thrust::complex<double>>>` |
 * | `dg::x::fDVec` | `dg::fDVec` = `thrust::device_vector<float>` | `dg::fMDVec` = `dg::MPI_Vector<thrust::device_vector<float>>` |
 * | `dg::x::iDVec` | `dg::iDVec` = `thrust::device_vector<int>` | `dg::iMDVec` = `dg::MPI_Vector<thrust::device_vector<int>>` |
 * | `dg::x::HMatrix` | `dg::HMatrix` = `dg::EllSparseBlockMat<double, thrust::host_vector>` | `dg::MHMatrix` = `dg::MPISparseBlockMat<thrust::host_vector, dg::EllSparseBlockMat<double, thrust::host_vector>, dg::CooSparseBlockMat<double, thrust::host_vector>>` |
 * | `dg::x::fHMatrix` | `dg::fHMatrix` = `dg::EllSparseBlockMat<float, thrust::host_vector>` | `dg::fMHMatrix` = `dg::MPISparseBlockMat<thrust::host_vector, dg::EllSparseBlockMat<float, thrust::host_vector>, dg::CooSparseBlockMat<float, thrust::host_vector>>` |
 * | `dg::x::DMatrix` | `dg::DMatrix` = `dg::EllSparseBlockMat<double, thrust::device_vector>` | `dg::MDMatrix` = `dg::MPISparseBlockMat<thrust::device_vector, dg::EllSparseBlockMat<double, thrust::device_vector>, dg::CooSparseBlockMat<double, thrust::device_vector>>` |
 * | `dg::x::fDMatrix` | `dg::fDMatrix` = `dg::EllSparseBlockMat<float, thrust::device_vector>` | `dg::fMDMatrix` = `dg::MPISparseBlockMat<thrust::device_vector, dg::EllSparseBlockMat<float, thrust::device_vector>, dg::CooSparseBlockMat<float, thrust::device_vector>>` |
 * | `dg::x::IHMatrix` | `dg::IHMatrix` = `dg::SparseMatrix<int, double, thrust::host_vector>` | `dg::MIHMatrix` = `dg::MPIDistMat<thrust::host_vector, dg::SparseMatrix<int, double, thrust::host_vector>>` |
 * | `dg::x::IDMatrix` | `dg::IDMatrix` = `dg::SparseMatrix<int, double, thrust::device_vector>` | `dg::MIDMatrix` = `dg::MPIDistMat<thrust::device_vector, dg::SparseMatrix<int, double, thrust::device_vector>>` |
 */

namespace dg{

///@addtogroup typedefs
///@{
//vectors
///@name Vector typedefs
///@{
template<class T>
using HVec_t  = thrust::host_vector<T>; //!< Host Vector
using HVec  = thrust::host_vector<double>; //!< Host Vector
using cHVec = thrust::host_vector<thrust::complex<double>>; //!< complex Host Vector
using iHVec = thrust::host_vector<int>; //!< integer Host Vector
using fHVec = thrust::host_vector<float>; //!< Host Vector

using DVec  = thrust::device_vector<double>; //!< Device Vector. The device can be an OpenMP parallelized cpu or a gpu. This depends on the value of the macro THRUST_DEVICE_SYSTEM, which can be either THRUST_DEVICE_SYSTEM_OMP for openMP or THRUST_DEVICE_SYSTEM_CUDA for a gpu or THRUST_DEVICE_SYSTEM_CPP for a cpu.
using iDVec = thrust::device_vector<int>; //!< integer Device Vector
using cDVec = thrust::device_vector<thrust::complex<double>>; //!< complex Device Vector
using fDVec = thrust::device_vector<float>; //!< Device Vector. The device can be an OpenMP parallelized, a serial cpu or a gpu. This depends on the value of the macro THRUST_DEVICE_SYSTEM, which can be either THRUST_DEVICE_SYSTEM_OMP for openMP or THRUST_DEVICE_SYSTEM_CUDA for a gpu or THRUST_DEVICE_SYSTEM_CPP for a cpu.
///@}

///@name Derivative matrices
///@{
template<class T>
using HMatrix_t = EllSparseBlockMat<T, thrust::host_vector>; //!< Host Matrix for derivatives
using HMatrix = EllSparseBlockMat<double, thrust::host_vector>; //!< Host Matrix for derivatives
using fHMatrix = EllSparseBlockMat<float, thrust::host_vector>; //!< Host Matrix for derivatives
using DMatrix = EllSparseBlockMat<double, thrust::device_vector>; //!< Device Matrix for derivatives
using fDMatrix = EllSparseBlockMat<float, thrust::device_vector>; //!< Device Matrix for derivatives
///@}

///@name Interpolation matrices
///@{
template<class real_type>
using IHMatrix_t = dg::SparseMatrix<int, real_type, thrust::host_vector>; //!< Host Matrix for interpolation
template<class real_type>
using IDMatrix_t = dg::SparseMatrix<int, real_type, thrust::device_vector>; //!< Device Matrix for interpolation
using IHMatrix = IHMatrix_t<double>; //!< Host Matrix for interpolation
using IDMatrix = IDMatrix_t<double>; //!< Device Matrix for interpolation
///@}

///@}
}//namespace dg

#ifdef MPI_VERSION
#include "mpi_vector.h"
#include "mpi_matrix.h"

namespace dg{
///@addtogroup typedefs
///@{
///@name MPI vector typedefs
///@{
template<class T>
using MHVec_t   = dg::MPI_Vector<dg::HVec_t<T> >; //!< see the \ref dg_tensor_typedefs_table "tensor typedefs table"
using MHVec     = dg::MPI_Vector<dg::HVec >;
using cMHVec    = dg::MPI_Vector<dg::cHVec >;
using fMHVec    = dg::MPI_Vector<dg::fHVec >;
using iMHVec    = dg::MPI_Vector<dg::iHVec >;
using MDVec     = dg::MPI_Vector<dg::DVec >;
using cMDVec    = dg::MPI_Vector<dg::cDVec >;
using fMDVec    = dg::MPI_Vector<dg::fDVec >;
using iMDVec    = dg::MPI_Vector<dg::iDVec >;
///@}

///@cond
template< class T>
using CooMat_t  = dg::CooSparseBlockMat<T, thrust::host_vector>; //!< Host Matrix for the communication part of \c dg::MHMatrix_t
using CooMat    = dg::CooSparseBlockMat<double, thrust::host_vector>; //!< Host Matrix for the communication part of \c dg::MHMatrix
using fCooMat   = dg::CooSparseBlockMat<float, thrust::host_vector>; //!< Host Matrix for the communication part of \c dg::fMHMatrix
using DCooMat   = dg::CooSparseBlockMat<double, thrust::device_vector>; //!< Device Matrix for the communication part of \c dg::MDMatrix
using fDCooMat  = dg::CooSparseBlockMat<float, thrust::device_vector>; //!< Device Matrix for the communication part of \c dg::fMDMatrix
///@endcond

///@name MPI derivative matrix typedefs
///@{
template<class T>
using MHMatrix_t  = dg::MPISparseBlockMat<thrust::host_vector, dg::HMatrix_t<T>, dg::CooMat_t<T>>; //!< see the \ref dg_tensor_typedefs_table "tensor typedefs table"
using MHMatrix    = dg::MPISparseBlockMat<thrust::host_vector, dg::HMatrix, dg::CooMat>;
using fMHMatrix   = dg::MPISparseBlockMat<thrust::host_vector, dg::fHMatrix, dg::fCooMat>;
using MDMatrix    = dg::MPISparseBlockMat<thrust::device_vector, dg::DMatrix, dg::DCooMat>;
using fMDMatrix   = dg::MPISparseBlockMat<thrust::device_vector, dg::fDMatrix, dg::fDCooMat>;
///@}


///@name MPI interpolation matrix typedefs
///@{
template<class real_type>
using MIHMatrix_t = MPIDistMat< thrust::host_vector, IHMatrix_t<real_type> >; //!< see the \ref dg_tensor_typedefs_table "tensor typedefs table"
template<class real_type>
using MIDMatrix_t = MPIDistMat< thrust::device_vector, IDMatrix_t<real_type> >;
using MIHMatrix = MIHMatrix_t<double>;
using MIDMatrix = MIDMatrix_t<double>;
///@}


///@}
}//namespace dg
#endif //MPI_VERSION

//Definitions for platform independent code
namespace dg{
///@addtogroup typedefs
///@{
//vectors
namespace x{
///@name Typedefs for platform independent code
///@{
#ifdef MPI_VERSION
using HVec  = MHVec; ///< see the \ref dg_tensor_typedefs_table "tensor typedefs table"
using cHVec = cMHVec;
using fHVec = fMHVec;
using iHVec = iMHVec;

using DVec  = MDVec;
using cDVec = cMDVec;
using fDVec = fMDVec;
using iDVec = iMDVec;

//derivative matrices
using HMatrix = MHMatrix;
using fHMatrix = fMHMatrix;
using DMatrix = MDMatrix;
using fDMatrix = fMDMatrix;
//interpolation matrices
using IHMatrix = MIHMatrix;
using IDMatrix = MIDMatrix;
#else
using HVec  = HVec;
using cHVec = cHVec;
using fHVec = fHVec;
using iHVec = iHVec;

using DVec  = DVec;
using cDVec = cDVec;
using fDVec = fDVec;
using iDVec = iDVec;

//derivative matrices
using HMatrix = HMatrix;
using fHMatrix = fHMatrix;
using DMatrix = DMatrix;
using fDMatrix = fDMatrix;
//interpolation matrices
using IHMatrix = IHMatrix;
using IDMatrix = IDMatrix;
#endif //MPI_VERSION
///@}
}//namespace x
///@}
}//namespace dg

////CONVENIENCE MACRO////////
#ifdef MPI_VERSION
/**
 * @brief Expands to \c if(rank==0) if \c MPI_VERSION is defined else is empty
 *
 * Use to let only rank 0 write to \c std::cout, e.g. <tt>DG_RANK0 std::cout << "Hello\n";</tt>
 * @attention With MPI an \c int variable \c rank must be in scope
 * @ingroup macros
 */
#define DG_RANK0 if(rank==0)
#else //MPI_VERSION
#define DG_RANK0
#endif //MPI_VERSION

#endif//_DG_TYPEDEFS_CUH_
