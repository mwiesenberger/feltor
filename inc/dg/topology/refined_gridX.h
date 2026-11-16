#pragma once

#include "interpolation.h"
#include "evaluationX.h"
#include "weightsX.h"
#include "gridX.h"
#include "base_geometryX.h"
#include "refined_grid.h"


namespace dg
{
///@addtogroup generators
///@{
/**
* @brief Abstract base class for 2d grid refinement that increases the number of grid cells of a fixed basis grid
*/
template<class real_type>
struct aRealRefinementX2d
{
    /*! @brief Generate the grid transformation
     *
     * @param g_old The 1d grid to refine
     * @param weightsX A 2d vector of size nx_new()*ny_new(). These represent the Jacobian of the transformation \f[\frac{\partial \zeta}{\partial x} \f]. The new metric element has thus to be multiplied by weights^2 and the volume by 1/weights
     * @param weightsY A 2d vector of size nx_new()*ny_new(). These represent the Jacobian of the transformation \f[\frac{\partial \eta}{\partial y} \f]. The new metric element has thus to be multiplied by weights^2 and the volume by 1/weights
     * @param abscissasX A 2d vector of size nx_new()*ny_new(). These are the new abscissas \f$ x(\zeta) \f$ of the grid.
     * @param abscissasY A 2d vector of size nx_new()*ny_new(). These are the new abscissas \f$ y(\eta) \f$ of the grid.
    */
    void generate( const RealGridX2d<real_type>& g_old, thrust::host_vector<real_type>& weightsX, thrust::host_vector<real_type>& weightsY, thrust::host_vector<real_type>& abscissasX, thrust::host_vector<real_type>& abscissasY) const
    {
        thrust::host_vector<real_type> wx, ax, wy, ay;
        // The X-point is never on the periodic x boundary (fx < 1), so a
        // non-periodic bc refines only the left side of x1 if fx == 0
        RealGrid1d<real_type> gx( g_old.x0(), g_old.x1(), g_old.n(), g_old.Nx(), dg::DIR);
        do_generateX(gx,g_old.inner_Nx(), wx,ax);
        RealGridX1d<real_type> gy( g_old.y0(), g_old.y1(), g_old.fy(), g_old.n(), g_old.Ny(), g_old.bcy());
        do_generateY(gy,wy,ay);
        unsigned size=wx.size()*wy.size();
        weightsX.resize(size), weightsY.resize(size);
        abscissasX.resize(size), abscissasY.resize(size);
        //now make product space
        for( unsigned i=0; i<wy.size(); i++)
            for( unsigned j=0; j<wx.size(); j++)
            {
                weightsX[i*wx.size()+j] = wx[j];
                weightsY[i*wx.size()+j] = wy[i];
                abscissasX[i*wx.size()+j] = ax[j];
                abscissasY[i*wx.size()+j] = ay[i];
            }
    }

    /*! @brief the new number of cells
     * @param Nx_old the old number of cells
     * @param fx_old Factor to partition x-axis
     */
    unsigned nx_new( unsigned Nx_old, real_type fx_old) const
    {
        return do_Nx_new(Nx_old, fx_old);
    }
    /*! @brief the new number of cells in y
     * @param Ny_old the old number of cells
     * @param fy_old Factor to partition y-axis
     */
    unsigned ny_new( unsigned Ny_old, real_type fy_old) const
    {
        return do_Ny_new(Ny_old, fy_old);
    }
    /*! @brief the new factor to partition the x-axis
     * @param Nx_old the old number of cells
     * @param fx_old the old factor to partition x-axis
     */
    real_type fx_new( unsigned Nx_old, real_type fx_old) const
    {
        return do_fx_new(Nx_old, fx_old);
    }
    /*! @brief the new factor to partition the y-axis
     * @param Ny_old the old number of cells
     * @param fy_old the old factor to partition y-axis
     */
    real_type fy_new( unsigned Ny_old, real_type fy_old) const
    {
        return do_fy_new(Ny_old, fy_old);
    }
    ///@copydoc hide_clone
    virtual aRealRefinementX2d* clone()const=0;
    ///@brief Allow deletion through base class pointer
    virtual ~aRealRefinementX2d() = default;
    protected:
    aRealRefinementX2d() = default;
    ///@brief We do not allow object slicing so the copy is protected
    aRealRefinementX2d(const aRealRefinementX2d& src) = default;
    ///@brief We do not allow object slicing so the assignment is protected
    aRealRefinementX2d& operator=(const aRealRefinementX2d& src) = default;
    private:
    virtual void do_generateX( const RealGrid1d<real_type>& gx, unsigned nodeXX, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const =0;
    virtual void do_generateY( const RealGridX1d<real_type>& gy, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const =0;
    virtual unsigned do_Nx_new( unsigned Nx_old, real_type fx) const =0;
    virtual unsigned do_Ny_new( unsigned Ny_old, real_type fy) const =0;
    virtual real_type do_fx_new( unsigned Nx_old, real_type fx) const =0;
    virtual real_type do_fy_new( unsigned Ny_old, real_type fy) const =0;
};

/**
* @brief No refinement
*/
template<class real_type>
struct RealIdentityXRefinement : public aRealRefinementX2d<real_type>
{
    ///@copydoc hide_clone
    virtual RealIdentityXRefinement* clone()const override final{
        return new RealIdentityXRefinement(*this);
    }
    private:
    virtual void do_generateX( const RealGrid1d<real_type>& g, unsigned, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const override final{
        weights = thrust::host_vector<real_type>( g.size(), 1);
        abscissas = g.abscissas(0);
    }
    virtual void do_generateY( const RealGridX1d<real_type>& g, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const override final{
        weights = thrust::host_vector<real_type>( g.size(), 1);
        abscissas = g.abscissas(0);
    }
    virtual unsigned do_Nx_new( unsigned Nx_old, real_type) const override final{ return Nx_old; }
    virtual unsigned do_Ny_new( unsigned Ny_old, real_type) const override final{ return Ny_old; }
    virtual real_type do_fx_new( unsigned, real_type fx) const override final{ return fx; }
    virtual real_type do_fy_new( unsigned, real_type fy) const override final{ return fy; }
};

/**
 * @brief Equidistant cell refinement around the X-point
 */
template<class real_type>
struct RealEquidistXRefinement : public aRealRefinementX2d<real_type>
{
    /**
     * @brief Divide cells around the X-point into equidistant cells
     *
     * @param add_x number of additional cells in each refined cell in x
     * @param add_y number of additional cells in each refined cell in y
     * @param howmanyX number of cells on each side of the X-point to refine in x
     * @param howmanyY number of cells on each side of the X-point to refine in y
     * @sa RealEquidistRefinement
     */
    RealEquidistXRefinement( unsigned add_x, unsigned add_y, unsigned howmanyX=1, unsigned howmanyY=1): add_x_(add_x), howm_x_(howmanyX), add_y_(add_y), howm_y_(howmanyY){ }
    ///@copydoc hide_clone
    virtual RealEquidistXRefinement* clone()const override final{return new RealEquidistXRefinement(*this);}
    private:
    unsigned add_x_, howm_x_, add_y_, howm_y_;
    virtual void do_generateX( const RealGrid1d<real_type>& gx, unsigned nodeXX, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const override final
    {
        RealEquidistRefinement<real_type> equi(add_x_, nodeXX, howm_x_);
        equi.generate(gx,weights,abscissas);
    }
    virtual void do_generateY( const RealGridX1d<real_type>& g, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const override final
    {
        RealEquidistRefinement<real_type> equi0(add_y_,0,howm_y_);
        if( add_y_ == 0 || howm_y_ == 0) {
            equi0.generate(g.grid(), weights, abscissas);
            return;
        }
        if( g.f() == 0) {
            equi0.generate( RealGrid1d<real_type>( g.x0(), g.x1(), g.n(), g.N(), dg::PER), weights, abscissas);
            return;
        }
        thrust::host_vector<real_type> w1, w2, w3;
        thrust::host_vector<real_type> a1, a2, a3;
        unsigned node1 = g.outer_N();
        RealEquidistRefinement<real_type> equi1(add_y_, node1, howm_y_);
        equi1.generate( RealGrid1d<real_type>(g.x0()                  ,g.x0()+g.outer_N()*g.h(), g.n(), g.outer_N(), dg::DIR), w1,a1); //left side
        equi0.generate( RealGrid1d<real_type>(g.x0()+g.outer_N()*g.h(),g.x1()-g.outer_N()*g.h(), g.n(), g.inner_N(), dg::PER), w2,a2);//inner side
        equi0.generate( RealGrid1d<real_type>(g.x1()-g.outer_N()*g.h(),g.x1()                  , g.n(), g.outer_N(), dg::DIR), w3,a3);//right side
        //now combine and unnormalize weights
        thrust::host_vector<real_type> wtot( w1.size() + w2.size() + w3.size()), atot( wtot);
        for( unsigned i=0; i<w1.size() ; i++)
            wtot[i]                     = w1[i]*equi1.N_new( g.outer_N(), dg::DIR)/g.outer_N();
        for( unsigned i=0; i<w2.size(); i++)
            wtot[w1.size()+i]           = w2[i]*equi0.N_new( g.inner_N(), dg::PER)/g.inner_N();
        for( unsigned i=0; i<w3.size(); i++)
            wtot[w1.size()+w2.size()+i] = w3[i]*equi0.N_new( g.outer_N(), dg::DIR)/g.outer_N();
        weights = wtot;
        abscissas = detail::normalize_weights_and_compute_abscissas( g.grid(), weights);
    }
    virtual unsigned do_Ny_new( unsigned Ny, real_type fy) const override final{
        if( fy==0 ) return Ny + 2*add_y_*howm_y_;
        return Ny + 4*add_y_*howm_y_;
    }
    virtual unsigned do_Nx_new( unsigned Nx, real_type fx) const override final{
        if( fx==0 ) return Nx + add_x_*howm_x_;
        return Nx + 2*add_x_*howm_x_;
    }
    virtual real_type do_fx_new( unsigned Nx, real_type fx) const override final{
        if( fx==0 ) return 0;
        return (fx*(real_type)Nx + (real_type)(add_x_*howm_x_))/(real_type)(Nx+2.*add_x_*howm_x_);
    }
    virtual real_type do_fy_new( unsigned Ny, real_type fy) const override final{
        if( fy==0 ) return 0;
        return (fy*(real_type)Ny + (real_type)(add_y_*howm_y_))/(real_type)(Ny+4.*add_y_*howm_y_);
    }
};

/**
 * @brief The exponential refinement around the X-point
 *
 * The y-direction is refined exponentially, the x-direction equidistantly
 * as in \c RealEquidistXRefinement
 */
template<class real_type>
struct RealExponentialXRefinement : public aRealRefinementX2d<real_type>
{
    /**
     * @brief Refine cells around the X-point
     *
     * @param add_x number of additional cells in each refined cell in x
     * @param add_y number of additional cells in the cells next to the
     * X-point in y; every new cell halves the cell closest to the X-point
     * @param howmanyX number of cells on each side of the X-point to refine in x
     * @param howmanyY unused
     * @sa RealEquidistRefinement RealExponentialRefinement
     */
    RealExponentialXRefinement( unsigned add_x, unsigned add_y, unsigned howmanyX=1, unsigned howmanyY=1): add_x_(add_x), howm_x_(howmanyX), add_y_(add_y), howm_y_(howmanyY){ }
    ///@copydoc hide_clone
    virtual RealExponentialXRefinement* clone()const override final{return new RealExponentialXRefinement(*this);}
    private:
    unsigned add_x_, howm_x_, add_y_, howm_y_;
    virtual void do_generateX( const RealGrid1d<real_type>& gx, unsigned nodeXX, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const override final
    {
        RealEquidistRefinement<real_type> equi(add_x_, nodeXX, howm_x_);
        equi.generate(gx,weights,abscissas);
    }
    virtual void do_generateY( const RealGridX1d<real_type>& g, thrust::host_vector<real_type>& weights, thrust::host_vector<real_type>& abscissas) const override final
    {
        RealExponentialRefinement<real_type> expo0( add_y_, 0);
        if( add_y_ == 0) { return expo0.generate( g.grid(), weights, abscissas); }
        if( g.f()  == 0) { return expo0.generate( RealGrid1d<real_type>( g.x0(), g.x1(), g.n(), g.N(), dg::PER), weights, abscissas); }
        thrust::host_vector<real_type> w1, w2, w3;
        thrust::host_vector<real_type> a1, a2, a3;
        unsigned node1 = g.outer_N();
        RealExponentialRefinement<real_type> expo1( add_y_, node1);
        expo1.generate( RealGrid1d<real_type>(g.x0()                  ,g.x0()+g.outer_N()*g.h(), g.n(), g.outer_N(), dg::DIR), w1,a1); //left side
        expo0.generate( RealGrid1d<real_type>(g.x0()+g.outer_N()*g.h(),g.x1()-g.outer_N()*g.h(), g.n(), g.inner_N(), dg::PER), w2,a2); //inner side
        expo0.generate( RealGrid1d<real_type>(g.x1()-g.outer_N()*g.h(),g.x1()                  , g.n(), g.outer_N(), dg::DIR), w3,a3); //right side
        //now combine and unnormalize weights
        thrust::host_vector<real_type> wtot( w1.size() + w2.size() + w3.size());
        for( unsigned i=0; i<w1.size() ; i++)
            wtot[i]                     = w1[i]*expo1.N_new( g.outer_N(), dg::DIR)/g.outer_N();
        for( unsigned i=0; i<w2.size(); i++)
            wtot[w1.size()+i]           = w2[i]*expo0.N_new( g.inner_N(), dg::PER)/g.inner_N();
        for( unsigned i=0; i<w3.size(); i++)
            wtot[w1.size()+w2.size()+i] = w3[i]*expo0.N_new( g.outer_N(), dg::DIR)/g.outer_N();
        weights = wtot;

        abscissas = detail::normalize_weights_and_compute_abscissas( g.grid(), weights);
    }
    virtual unsigned do_Ny_new( unsigned Ny, real_type fy) const override final{
        if( fy==0 ) return Ny + 2*add_y_;
        return Ny + 4*add_y_;
    }
    virtual unsigned do_Nx_new( unsigned Nx, real_type fx) const override final{
        if( fx==0 ) return Nx + add_x_*howm_x_;
        return Nx + 2*add_x_*howm_x_;
    }
    virtual real_type do_fx_new( unsigned Nx, real_type fx) const override final{
        if( fx==0 ) return 0;
        return (fx*(real_type)Nx + (real_type)(add_x_*howm_x_))/(real_type)(Nx+2.*add_x_*howm_x_);
    }
    virtual real_type do_fy_new( unsigned Ny, real_type fy) const override final{
        if( fy==0 ) return 0;
        return (fy*(real_type)Ny + (real_type)add_y_)/(real_type)(Ny+4.*add_y_);
    }
};

///@name Double precision typedefs
///@{
using aRefinementX2d            = dg::aRealRefinementX2d<double>; ///< Double precision version of the class template with the \c Real prefix
using IdentityXRefinement       = dg::RealIdentityXRefinement<double>;
using EquidistXRefinement       = dg::RealEquidistXRefinement<double>;
using ExponentialXRefinement    = dg::RealExponentialXRefinement<double>;
///@}
///@}


/**
 * @brief Refined X-point grid
 * @ingroup basicXgeometry
 */
template<class real_type>
struct RealCartesianRefinedGridX2d : public dg::aRealGeometryX2d<real_type>
{
    ///@param ref the refinement (\c Nx, \c Ny, \c fx and \c fy are the values before refinement)
    ///@copydoc hide_gridX_parameters2d
    ///@copydoc hide_bc_parameters2d
    RealCartesianRefinedGridX2d( const aRealRefinementX2d<real_type>& ref,
            real_type x0, real_type x1, real_type y0, real_type y1,
            real_type fx, real_type fy,
            unsigned n, unsigned Nx, unsigned Ny,
            bc bcx = dg::PER, bc bcy = dg::PER) : dg::aGeometryX2d( x0, x1, y0, y1,
                ref.fx_new(Nx, fx), ref.fy_new(Ny, fy), n, ref.nx_new(Nx, fx), ref.ny_new(Ny, fy), bcx, bcy), w_(2), abs_(2)
    {
        RealGridX2d<real_type> g( x0,x1,y0,y1,fx,fy,n,Nx,Ny,bcx,bcy);
        ref.generate(g,w_[0],w_[1],abs_[0],abs_[1]);
    }

    ///@copydoc hide_clone
    virtual RealCartesianRefinedGridX2d* clone()const override final{return new RealCartesianRefinedGridX2d(*this);}
    private:
    std::vector<thrust::host_vector<real_type> > w_,abs_;
    virtual SparseTensor<thrust::host_vector<real_type> > do_compute_metric()const override final{
        SparseTensor<thrust::host_vector<real_type> > t(*this);
        t.values().push_back( w_[0]);
        t.values().push_back( w_[1]);
        dg::blas1::pointwiseDot( w_[0], w_[0], t.values()[2]);
        dg::blas1::pointwiseDot( w_[1], w_[1], t.values()[3]);
        t.idx(0,0)=2, t.idx(1,1)=3;
        return t;
    }
    virtual SparseTensor<thrust::host_vector<real_type> > do_compute_jacobian()const override final{
        SparseTensor<thrust::host_vector<real_type> > t(*this);
        t.values().push_back( w_[0]);
        t.values().push_back( w_[1]);
        t.idx(0,0)=2, t.idx(1,1)=3;
        return t;
    }
    virtual std::vector<thrust::host_vector<real_type> > do_compute_map()const override final{
        return abs_;
    }
};

/**
 * @brief Refined X-point grid
 * @ingroup basicXgeometry
 */
template<class real_type>
struct RealCartesianRefinedGridX3d : public dg::aRealGeometryX3d<real_type>
{
    ///@param ref the refinement (\c Nx, \c Ny, \c fx and \c fy are the values before refinement)
    ///@copydoc hide_gridX_parameters3d
    ///@copydoc hide_bc_parameters3d
    RealCartesianRefinedGridX3d( const aRealRefinementX2d<real_type>& ref,
            real_type x0, real_type x1, real_type y0, real_type y1, real_type z0, real_type z1,
            real_type fx, real_type fy,
            unsigned n, unsigned Nx, unsigned Ny, unsigned Nz,
            bc bcx = dg::PER, bc bcy = dg::PER, bc bcz = dg::PER) : dg::aGeometryX3d(
                x0, x1, y0, y1,z0,z1,
                ref.fx_new(Nx, fx), ref.fy_new(Ny, fy),
                n, ref.nx_new(Nx, fx), ref.ny_new(Ny, fy), Nz,
                bcx, bcy, bcz), w_(2), abs_(2)
    {
        RealGridX2d<real_type> g( x0,x1,y0,y1,fx,fy,n,Nx,Ny,bcx,bcy);
        ref.generate(g,w_[0],w_[1],abs_[0],abs_[1]);
        //lift to 3d
        w_[0].resize(this->size()), w_[1].resize(this->size()), abs_[0].resize(this->size()), abs_[1].resize(this->size());
        unsigned size2d=this->n()*this->n()*this->Nx()*this->Ny();
        for( unsigned i=1; i<Nz; i++)
            for(unsigned k=0; k<size2d; k++)
            {
                w_[0][i*size2d+k]=w_[0][k];
                w_[1][i*size2d+k]=w_[1][k];
                abs_[0][i*size2d+k]=abs_[0][k];
                abs_[1][i*size2d+k]=abs_[1][k];
            }
    }

    ///@copydoc hide_clone
    virtual RealCartesianRefinedGridX3d* clone()const{return new RealCartesianRefinedGridX3d(*this);}
    private:
    std::vector<thrust::host_vector<real_type> > w_,abs_;
    virtual SparseTensor<thrust::host_vector<real_type> > do_compute_metric()const override final {
        SparseTensor<thrust::host_vector<real_type> > t(*this);
        t.values().push_back( w_[0]);
        t.values().push_back( w_[1]);
        dg::blas1::pointwiseDot( w_[0], w_[0], t.values()[2]);
        dg::blas1::pointwiseDot( w_[1], w_[1], t.values()[3]);
        t.idx(0,0)=2, t.idx(1,1)=3;
        return t;
    }
    virtual SparseTensor<thrust::host_vector<real_type> > do_compute_jacobian()const override final {
        SparseTensor<thrust::host_vector<real_type> > t(*this);
        t.values().push_back( w_[0]);
        t.values().push_back( w_[1]);
        t.idx(0,0)=2, t.idx(1,1)=3;
        return t;
    }
    virtual std::vector<thrust::host_vector<real_type> > do_compute_map()const override final {
        return abs_;
    }
};

///@addtogroup gridtypes
///@{
///@name Double precision typedefs
///@{
using CartesianRefinedGridX2d = dg::RealCartesianRefinedGridX2d<double>; ///< Double precision version of the class template with the \c Real prefix
using CartesianRefinedGridX3d = dg::RealCartesianRefinedGridX3d<double>;
///@}
///@}


}//namespace dg
