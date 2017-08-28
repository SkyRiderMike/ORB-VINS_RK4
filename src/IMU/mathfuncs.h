#pragma once

#define SMALL_EPS 1e-10
#include <Eigen/Core>

inline Eigen::Matrix<double, 3, 2> computeA(Eigen::Vector3d& n)
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

inline Eigen::Matrix3d skew(const Eigen::Vector3d &v) {
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

inline Eigen::Matrix3d Exp(const Eigen::Vector3d &dx)
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
