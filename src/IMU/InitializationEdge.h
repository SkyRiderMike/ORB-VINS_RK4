#pragma once

#include "Thirdparty/g2o/g2o/core/base_vertex.h"
#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "so3.h"
#include "NavState.h"
#include "IMUPreintegrator.h"
#include "Thirdparty/g2o/g2o/core/base_multi_edge.h"
#include "Thirdparty/g2o/g2o/core/base_binary_edge.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"

#define SMALL_EPS 1e-10

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

    VertexScaleAndGravity() : BaseVertex<3, Eigen::Matrix<double, 4, 1>>()
    {
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

    Eigen::Matrix<double, 3, 2> computeA(Eigen::Vector3d &n)
    {
        size_t minid;
        Eigen::Matrix<double, 3, 2> A;
        n.cwiseAbs().minCoeff(&minid);
        switch (minid)
        {
        case 0:
            A.col(0) = Eigen::Vector3d(0.0, -n(2), n(1));
            break;
        case 1:
            A.col(0) = Eigen::Vector3d(n(2), 0.0, n(0));
            break;
        case 2:
            A.col(0) = Eigen::Vector3d(n(1), -n(0), 0.0);
            break;
        }
        A.col(1) = n.cross(A.col(0));
        return A;
    }

    Eigen::Matrix3d Exp(const Eigen::Vector3d &dx)
    {
        double theta = dx.norm();
        if (theta < SMALL_EPS)
        {
            return Eigen::Matrix3d::Identity();
        }
        else
        {
            Eigen::Matrix3d hatdx = skew(dx / theta);
            return Eigen::Matrix3d::Identity() + sin(theta) * hatdx + (1 - cos(theta)) * hatdx * hatdx;
        }
    }
};

class VertexVelocityAndBias : public BaseVertex<6, Eigen::Matrix<double, 9, 1>>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VertexVelocityAndBias() : BaseVertex<6, Eigen::Matrix<double, 6, 1>>()
    {
    }

    bool read(std::istream &is) { return true; }
    bool write(std::ostream &os) { return true; }

    virtual void setToOriginImpl() override
    {
        _estimate = Eigen::Matrix<double, 6, 1>::Zero();
    }

    virtual void oplusImpl(const double *update_) override
    {
        Eigen::Map<const Eigen::Matrix<double, 6, 1>> update(update_);
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

    void SetParams(const Eigen::Vector3d& Pi,const Eigen::Vector3d& Pj, const Eigen::Matrix3d& Ri, const Eigen::Matrix3d& Rj, const Eigen::Vector3d& Gbias)
    {
        // GravityVec = gw;
        _Pi = Pi;
        _Pj = Pj;
        _Ri = Ri;
        _Rj = Rj;
        _Gbias = Gbias;
    }

    Eigen::Matrix<double, 3, 2> computeA(Eigen::Vector3d &n)
    {
        size_t minid;
        Eigen::Matrix<double, 3, 2> A;
        n.abs().minCoeff(&minid);
        switch (minid)
        {
        case 0:
            A.col(0) = Eigen::Vector3d(0.0, -n(2), n(1));
            break;
        case 1:
            A.col(0) = Eigen::Vector3d(n(2), 0.0, n(0));
            break;
        case 2:
            A.col(0) = Eigen::Vector3d(n(1), -n(0), 0.0);
            break;
        }
        A.col(1) = n.cross(A.col(0));
        return A;
    }

    Eigen::Matrix3d skew(const Eigen::Vector3d &v) {
        Eigen::Matrix3d mat;
        mat(0, 0) = 0.0;
        mat(0, 1) = -v(2);
        mat(0, 2) = v(1);
        mat(1, 0) = v(2);
        mat(1, 1) = 0.0;
        mat(1, 2) = -v(0);
        mat(2, 0) = -v(1);
        mat(2, 1) = v(0);
        mat(2, 2) = 0.0;

        return mat;
    }

  protected:
    // Gravity vector in 'world' frame
    Vector3d _Pi;
    Vector3d _Pj;
    Matrix3d _Rj;
    Matrix3d _Ri;
    Vector3d _Gbias;

};
}