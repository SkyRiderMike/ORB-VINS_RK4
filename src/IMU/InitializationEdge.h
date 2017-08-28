#pragma once

#include "Thirdparty/g2o/g2o/core/base_vertex.h"
#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "so3.h"
#include "NavState.h"
#include "IMUPreintegrator.h"
#include "Thirdparty/g2o/g2o/core/base_multi_edge.h"
#include "Thirdparty/g2o/g2o/core/base_binary_edge.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"

#include "mathfuncs.h"

// unused yet.
namespace g2o
{
using namespace Eigen;
using namespace Sophus;
using namespace ORB_SLAM2;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
/*
Vector :
s(scalar)
g(3x1)
*/
class VertexScaleAndGravity : public BaseVertex<3, Eigen::Matrix<double, 4, 1>>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VertexScaleAndGravity() : BaseVertex<3, Eigen::Matrix<double, 4, 1> >()
    {
        setToOriginImpl();
    }

    bool read(std::istream &is) { return true; }

    bool write(std::ostream &os) const { return true; }

    virtual void setToOriginImpl()
    {
        // scale should be initialized to 1.0
        // gravity should be 0,0,9.81
        _estimate = Eigen::Matrix<double, 4, 1>(1.0, 0, 0, 9.81);
    }

    virtual void oplusImpl(const double *update_);


};

class VertexVelocityAndBias : public BaseVertex<6, Eigen::Matrix<double, 6, 1>>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VertexVelocityAndBias() : BaseVertex<6, Eigen::Matrix<double, 6, 1> >()
    {
        setToOriginImpl();
    }

    bool read(std::istream &is) { return true; }
    bool write(std::ostream &os) const { return true; }

    virtual void setToOriginImpl() override
    {
        _estimate = Eigen::Matrix<double, 6, 1>::Zero();
    }

    virtual void oplusImpl(const double *update_) override
    {
        Eigen::Map<const Eigen::Matrix<double, 6, 1> > update(update_);
        _estimate += update;
    }
};

class EdgeNavStateInitialization : public BaseMultiEdge<9, IMUPreintegrator>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeNavStateInitialization() : BaseMultiEdge<9, IMUPreintegrator>()
    {
        resize(3);
    }

    bool read(std::istream &is) { return true; }

    bool write(std::ostream &os) const { return true; }

    virtual void computeError() override;

    virtual void linearizeOplus() override;

    void SetParams(const Eigen::Vector3d& Pi, // i is the last keyframe
                   const Eigen::Vector3d& Pj, // j is the current keyframe
                   const Eigen::Matrix3d& Ri, 
                   const Eigen::Matrix3d& Rj, 
                   const Eigen::Vector3d& Gbiasj)
    {
        // GravityVec = gw;
        _Pi = Pi;
        _Pj = Pj;
        _Ri = Ri;
        _Rj = Rj;
        _Gbiasj = Gbiasj;
    }


  protected:
    // Gravity vector in 'world' frame
    Vector3d _Pi;
    Vector3d _Pj;
    Matrix3d _Rj;
    Matrix3d _Ri;
    Vector3d _Gbiasj;

};
}