#include "g2otypes.h"
#include "InitializationEdge.h"

namespace g2o
{


void VertexScaleAndGravity::oplusImpl(const double* update_)
{
    _estimate(0) += update_[0];
    Eigen::Map<const Eigen::Matrix<double,2,1> > update(update_ + 1);
    
    Eigen::Vector3d g = _estimate.segment<3>(1);
    Eigen::Vector3d _gravity_update = computeA(g) * update;
    _estimate.segment<3>(1) = Exp(_gravity_update) *g;
}

void EdgeNavStateInitialization::computeError() 
{
    // const VertexNavState* vNavState = static_cast<const VertexNavState*>(_vertices[0]);
    const VertexVelocityAndBias* vVRi = static_cast<const VertexVelocityAndBias*>(_vertices[0]);
    const VertexVelocityAndBias* vVRj = static_cast<const VertexVelocityAndBias*>(_vertices[1]);
    const VertexScaleAndGravity* vScaleGravity = static_cast<const VertexScaleAndGravity*>(_vertices[2]);

    const Eigen::Matrix<double,6,1> vrb_i  = vVRi->estimate();
    Eigen::Matrix<double,3,1> Vi = vrb_i.segment<3>(0);
    Eigen::Matrix<double,3,1> Bai = vrb_i.segment<3>(3);
    // Eigen::Matrix<double,3,1> Bwi = vrb_i.segment<3>(6);

    const Eigen::Matrix<double,6,1> vrb_j = vVRj->estimate();
    Eigen::Matrix<double,3,1> Vj = vrb_j.segment<3>(0);
    Eigen::Matrix<double,3,1> Baj = vrb_j.segment<3>(3);
    // Eigen::Matrix<double,3,1> Bwj = vrb_j.segment<3>(6);

    const Eigen::Matrix<double,4,1> vSv = vScaleGravity->estimate();
    double scale = vSv(0);
    Eigen::Matrix<double, 3,1> GravityVec = vSv.segment<3>(1);

    const IMUPreintegrator& M = _measurement;
    double dTij = M.getDeltaTime(); 
    double dT2 = dTij*dTij;
    Eigen::Vector3d dPij = M.getDeltaP();
    Eigen::Vector3d dVij = M.getDeltaV();
    
    // tmp variable, transpose of Ri
    Eigen::Matrix3d RiT = _Ri.transpose();

    // Eigen::Matrix<double, 12, 1> err;
    _error.setZero();

    _error.segment<3>(0) = RiT * (scale * _Pj - scale * _Pi - Vi * dTij - 0.5 * GravityVec * dT2) - (dPij + M.getJPBiasg() * _Gbiasj + M.getJPBiasa() * Baj);
    _error.segment<3>(3) = RiT * (Vj - Vi - GravityVec * dTij) - (dVij + M.getJVBiasg() * _Gbiasj + M.getJVBiasa() * Baj);
    _error.segment<3>(6) = Baj - Bai;
    // _error.segment<3>(9) = Bwj - Bwi;

}

void EdgeNavStateInitialization::linearizeOplus()
{
    // const VertexNavState* vNavStatex = static_cast<const VertexNavState*>(_vertices[0]);
    const VertexScaleAndGravity* vScaleGravity = static_cast<const VertexScaleAndGravity*>(_vertices[2]);

    const Eigen::Matrix<double,4,1> vSv = vScaleGravity->estimate();
    // double scale = vSv(0);
    Eigen::Matrix<double, 3,1> Gravity = vSv.segment<3>(1);

    const IMUPreintegrator& M = _measurement;
    double dTij = M.getDeltaTime(); 
    double dT2 = dTij*dTij;
    
    // tmp variable, transpose of Ri
    Eigen::Matrix<double, 9, 6> JVRi, JVRj;
    JVRi.setZero();
    JVRj.setZero();
    Eigen::Matrix<double, 9, 6> JscaleGravity;
    JscaleGravity.setZero();
    Eigen::Matrix3d RiT = _Ri.transpose(); 
    JVRi.block<3,3>(0,0) = -dTij * RiT;
    JVRi.block<3,3>(3,0) = -RiT;
    JVRi.block<3,3>(6,3) = -Eigen::Matrix3d::Identity();
    JVRi.block<3,3>(9,6) = -Eigen::Matrix3d::Identity();

    JVRj.block<3,3>(0,3) = -M.getJPBiasa();
    JVRj.block<3,3>(0,6) = -M.getJPBiasg();
    JVRj.block<3,3>(3,0) = RiT;
    JVRj.block<3,3>(3,3) = -M.getJVBiasa();
    JVRj.block<3,3>(3,6) = -M.getJVBiasg();
    JVRj.block<3,3>(6,3) = Eigen::Matrix3d::Identity();
    JVRj.block<3,3>(9,6) = Eigen::Matrix3d::Identity();

    Eigen::Matrix<double,3,2> skewgA  = skew(Gravity)*computeA(Gravity);

    JscaleGravity.block<3,1>(0,0) = RiT * (_Pj - _Pi);
    JscaleGravity.block<3,2>(0,1) = 0.5 * dT2 * RiT * skewgA;
    JscaleGravity.block<3,2>(3,1) = dTij * RiT * skewgA;

    // evaluate
    _jacobianOplus[0] = JVRi;
    _jacobianOplus[1] = JVRj;
    _jacobianOplus[2] = JscaleGravity;
    
}

}