#include <davinci-kinematics/davinci.h>
#include <cmath>
#include <iostream>
#include <pugixml.hpp>
#include <stdexcept>

namespace
{
  Eigen::Matrix4d xRotation(double theta)
  {
    Eigen::Matrix4d R = Eigen::Matrix4d::Identity();
    double c = cos(theta);
    double s = sin(theta);

    R(1,1) = c;
    R(2,1) = s;
  
    R(1,2) = -s;
    R(2,2) = c;

    return R;
  }

  Eigen::Matrix4d yRotation(double theta)
  {
    Eigen::Matrix4d R = Eigen::Matrix4d::Identity();
    double c = cos(theta);
    double s = sin(theta);

    R(0,0) = c;
    R(2,0) = -s;

    R(0,2) = s;
    R(2,2) = c;

    return R;
  }

  Eigen::Matrix4d zRotation(double theta)
  {
    Eigen::Matrix4d R = Eigen::Matrix4d::Identity();
    double c = cos(theta);
    double s = sin(theta);
    R(0,0) = c;
    R(1,0) = s;

    R(0,1) = -s;
    R(1,1) = c;

    return R;
  }
}

Eigen::Matrix4d getInverseTransform(const Eigen::Matrix4d& T)
{
  Eigen::Matrix3d R = T.topLeftCorner<3,3>();
  Eigen::Vector3d p = T.topRightCorner<3,1>();
  Eigen::Matrix4d T_inv;
  T_inv.topLeftCorner<3,3>() = R.transpose();
  T_inv.topRightCorner<3,1>() = -R.transpose()*p;

  T_inv(3,0) = T_inv(3,1) = T_inv(3,2) = 0.0;
  T_inv(3,3) = 1.0;

  return T_inv;
}

DavinciKinematics::DavinciKinematics(const std::string& inputFilename)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(inputFilename.c_str());

  if (!result)
    {
    std::cout << "Error in loading kinematics file: " << result.description() << std::endl;
    throw std::runtime_error("Kinematics file load error");
    }
  
  pugi::xml_node intra = doc.child("davinci_parameters").child("intracorporeal");
  params_.wristLength_ = intra.child("wrist_length").text().as_double();
  params_.gripperLength_ = intra.child("gripper_length").text().as_double();

  pugi::xml_node extra = doc.child("davinci_parameters").child("extracorporeal");
  params_.eL_ = extra.child("L").text().as_double();
  params_.el1_ = extra.child("l1").text().as_double();
  params_.el2_ = extra.child("l2").text().as_double();
  double l = extra.child("l").text().as_double();
  params_.eAlpha_ = asin(l / params_.eL_);
  params_.eRadius1_ = extra.child("radius1").text().as_double();
  params_.eRadius2_ = extra.child("radius2").text().as_double();

  pugi::xml_node passive = doc.child("davinci_parameters").child("passive");
  params_.pl2_ = passive.child("l2").text().as_double();
  params_.ph2_ = passive.child("h2").text().as_double();
  params_.pl3_ = passive.child("l3").text().as_double();
  params_.ph3_ = passive.child("h3").text().as_double();
  params_.pl4_ = passive.child("l4").text().as_double();
  params_.ph4_ = passive.child("h4").text().as_double();
  params_.pl5_ = passive.child("l5").text().as_double();
  
  pugi::xml_node rcmOffset = passive.child("rcm_offset");
  params_.pRCMOffset_(0) = rcmOffset.child("x").text().as_double();
  params_.pRCMOffset_(1) = rcmOffset.child("y").text().as_double();
  params_.pRCMOffset_(2) = rcmOffset.child("z").text().as_double();

  params_.pLinkRadius_ = passive.child("link_radius").text().as_double();
}

// \TODO reduce number of intermediate matrices for optimization
Eigen::Matrix4d DavinciKinematics::intraFK(const Eigen::Matrix4d& portFrame,
                                           const std::vector<double>& q) const
{
  assert(q.size() == 6);

  Eigen::Vector3d e1; e1(0) = 1.0; e1(1) = 0.0; e1(2) = 0.0;
  Eigen::Vector3d e2; e2(0) = 0.0; e2(1) = 1.0; e2(2) = 0.0;
  Eigen::Vector3d e3; e3(0) = 0.0; e3(1) = 0.0; e3(2) = 1.0;

  Eigen::Matrix4d T = portFrame;
  T *= zRotation(q[0]);

  Eigen::Matrix4d T_1_2 = Eigen::Matrix4d::Zero();
  T_1_2.block<3,1>(0,0) = e1;
  T_1_2.block<3,1>(0,1) = -e3;
  T_1_2.block<3,1>(0,2) = e2;
  T_1_2(3,3) = 1.0;
  T_1_2 *= zRotation(q[1]);
  T = T*T_1_2;

  Eigen::Matrix4d T_2_3 = Eigen::Matrix4d::Zero();
  T_2_3.block<3,1>(0,0) = -e2;
  T_2_3.block<3,1>(0,1) = e3;
  T_2_3.block<3,1>(0,2) = -e1;
  T_2_3(3,3) = 1.0;
  T_2_3 *= zRotation(q[2]);
  T = T*T_2_3;

  Eigen::Matrix4d T_3_4 = Eigen::Matrix4d::Zero();
  T_3_4.block<3,1>(0,0) = e1;
  T_3_4.block<3,1>(0,1) = e2;
  T_3_4.block<3,1>(0,2) = e3;
  T_3_4(2,3) = q[3];
  T_3_4(3,3) = 1.0;
  T = T*T_3_4;

  Eigen::Matrix4d T_4_5 = Eigen::Matrix4d::Zero();
  T_4_5.block<3,1>(0,0) = e3;
  T_4_5.block<3,1>(0,1) = e1;
  T_4_5.block<3,1>(0,2) = e2;
  T_4_5(3,3) = 1.0;
  T_4_5 *= zRotation(q[4]);
  T = T*T_4_5;

  Eigen::Matrix4d T_5_6 = Eigen::Matrix4d::Zero();
  T_5_6.block<3,1>(0,0) = e1;
  T_5_6.block<3,1>(0,1) = -e3;
  T_5_6.block<3,1>(0,2) = e2;
  T_5_6(0,3) = params_.wristLength_;
  T_5_6(3,3) = 1.0;
  T_5_6 *= zRotation(q[5]);
  T = T*T_5_6;

  // For gripper. Bring it in later when we have IK working to correspond to this
  // Eigen::Matrix4d T_6_e = Eigen::Matrix4d::Zero();
  // T_6_e.block<3,1>(0,0) = e3;
  // T_6_e.block<3,1>(0,1) = e2;
  // T_6_e.block<3,1>(0,2) = e1;
  // T_6_e(0,3) = gripperLength_;
  // T_6_e(3,3) = 1.0;
  // T = T*T_6_e;

  return T;
}

Eigen::Matrix4d DavinciKinematics::passiveFK(const Eigen::Matrix4d& baseFrame,
                                             const std::vector<double>& q) const
{
  std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d> > poses;
  this->passiveFK(baseFrame, q, &poses);

  return poses.back();
}

// TODO factor in length of gripper tool by just adjusting initial
// pose before doing anything
void DavinciKinematics::intraIK(const Eigen::Matrix4d& portFrame,
                                const Eigen::Matrix4d& worldPose,
                                std::vector<double>* q) const
{
  assert(q->size() == 6);

  // Express input pose in terms of arm's base frame
  Eigen::Matrix4d pose = getInverseTransform(portFrame) * worldPose;

  Eigen::Matrix3d R = pose.topLeftCorner<3,3>();
  Eigen::Vector3d p = pose.topRightCorner<3,1>();
  Eigen::Vector3d positionTrans = -R.transpose()*p;
  double u = positionTrans(0); 
  double v = positionTrans(1);
  double w = positionTrans(2);
 
  (*q)[5] = atan2(v, -u);
  double c5 = cos((*q)[5]);
  double s5 = sin((*q)[5]);
  (*q)[4] = atan2(w, -params_.wristLength_ - u*c5 + v*s5);
  double c4 = cos((*q)[4]);
  double s4 = sin((*q)[4]);
  (*q)[3] = (-params_.wristLength_ - u*c5 + v*s5)*c4 + w*s4;
  
  Eigen::Matrix3d m;
  m(0,0) = c5*s4;  m(0,1) = -s5;  m(0,2) = c4*c5;
  m(1,0) = -s4*s5; m(1,1) = -c5;  m(1,2) = -c4*s5;
  m(2,0) = c4;     m(2,1) = 0.0;  m(2,2) = -s4;
  Eigen::Matrix3d Rbar = R*m;

  (*q)[1] = atan2(Rbar(2,2), sqrt(Rbar(0,2)*Rbar(0,2) + Rbar(1,2)*Rbar(1,2)));
  double c1 = cos((*q)[1]);
  (*q)[2] = atan2(-Rbar(2,1) / c1, Rbar(2,0) / c1);
  (*q)[0] = atan2(-Rbar(1,2) / c1, -Rbar(0,2) / c1);
}

// \TODO allow for specification of q0 evaluation if computed
// beforehand
Eigen::Matrix<double, 3, 6> 
DavinciKinematics::passiveJacobian(const Eigen::Matrix4d& baseFrame,
                                   const std::vector<double>& q_in,
                                   double stepSize) const
{
  std::vector<double> q(q_in);

  Eigen::Matrix<double, 3, 6> jacobian;

  // Using forward difference
  Eigen::Vector3d p0 = this->passiveFK(baseFrame, q).topRightCorner<3,1>();

  for (std::size_t i = 0; i < q.size(); ++i)
    {
    double q0 = q[i];
    q[i] += stepSize;
    Eigen::Vector3d p1 = this->passiveFK(baseFrame, q).topRightCorner<3,1>();
    jacobian.col(i) = (p1 - p0) / stepSize;
    q[i] = q0;
    }

  return jacobian;
}

// Iteration step size calculation came from: 
// http://math.ucsd.edu.libproxy.lib.unc.edu/~sbuss/ResearchWeb/ikmethods/iksurvey.pdf
//
// TODO maybe try DLS from same paper
bool DavinciKinematics::passiveIK(const Eigen::Matrix4d& baseFrame,
                                  const Eigen::Vector3d& worldPosition,
                                  std::vector<double>* qOut,
                                  double stepSize) const
{
  std::vector<double> q(6, 0.0);
  Eigen::Matrix<double,6,1> q_vec;

  Eigen::Vector3d goalVec = 
    worldPosition - this->passiveFK(baseFrame, q).topRightCorner<3,1>();
  double error = goalVec.squaredNorm();
  unsigned iterations = 0;
  unsigned MAX_ITERATIONS = 1000;
  while (error > 0.001*0.001 && iterations < MAX_ITERATIONS) // TODO do tolerance better
    {
    Eigen::Matrix<double, 3, 6> jacobian = this->passiveJacobian(baseFrame,
                                                                 q, stepSize);
    for (unsigned i = 0; i < 6; ++i)
      q_vec(i) = q[i];
    
    Eigen::Vector3d JJTe = jacobian*jacobian.transpose()*goalVec;
    double step = goalVec.dot(JJTe) / JJTe.squaredNorm();
    q_vec += step*jacobian.transpose()*goalVec;

    for (unsigned i = 0; i < 6; ++i)
      q[i] = q_vec(i);

    goalVec = worldPosition - this->passiveFK(baseFrame, q).topRightCorner<3,1>();
    error = goalVec.squaredNorm();
    }

  qOut->swap(q);

  if (iterations >= MAX_ITERATIONS)
    {
    std::cout << "Warning: passive IK didn't converge!" << std::endl;
    return false;
    }
  else
    return true;  
}

void DavinciKinematics::getExtraCylispheres(const Eigen::Matrix4d& portFrame,
                                            const std::vector<double>& q,
                                            std::vector<Collisions::Cylisphere>* cylispheres) const
{
  double c1 = cos(q[0]);
  double c2 = cos(q[1]);
  double s1 = sin(q[0]);
  double s2 = sin(q[1]);
  double c2a = cos(q[1] + params_.eAlpha_);
  double s2a = sin(q[1] + params_.eAlpha_);

  // horizontal cylisphere
  Collisions::Cylisphere c_h;
  Eigen::Vector4d p;
  p(0) = params_.eL_*c2a*c1;
  p(1) = params_.eL_*c2a*s1;
  p(2) = -params_.eL_*s2a;
  p(3) = 1.0;
  c_h.p1.noalias() = (portFrame*p).head<3>();
  
  p(2) -= params_.el2_;
  c_h.p2.noalias() = (portFrame*p).head<3>();
  c_h.r = params_.eRadius2_;

  // vertical cylisphere 
  Collisions::Cylisphere c_v;
  c_v.p1 = c_h.p1;
  p(2) += params_.el2_; // undo z offset from p2h
  p(0) += params_.el1_*c2*c1;
  p(1) += params_.el1_*c2*s1;
  p(2) += -params_.el1_*s2;
  c_v.p2.noalias() = (portFrame*p).head<3>();
  c_v.r = params_.eRadius1_;

  cylispheres->push_back(c_h);
  cylispheres->push_back(c_v);
}

void DavinciKinematics::passiveFK(const Eigen::Matrix4d& baseFrame,
                                  const std::vector<double>& q,
                                  std::vector<Eigen::Matrix4d,
                                  Eigen::aligned_allocator<Eigen::Matrix4d> >* jointPoses) const
{
  assert(q.size() == 6);

  Eigen::Vector4d e1; e1(0) = 1.0; e1(1) = 0.0; e1(2) = 0.0; e1(3) = 0.0;
  Eigen::Vector4d e2; e2(0) = 0.0; e2(1) = 1.0; e2(2) = 0.0; e2(3) = 0.0;
  Eigen::Vector4d e3; e3(0) = 0.0; e3(1) = 0.0; e3(2) = 1.0; e3(3) = 0.0;

  Eigen::Matrix4d T = baseFrame;
  jointPoses->push_back(T);

  // \TODO add this joint pose and fix passive primitives accordingly
  Eigen::Matrix4d T_0_1 = Eigen::Matrix4d::Identity();
  T_0_1(2,3) = q[0];
  T = T*T_0_1;
  
  Eigen::Matrix4d T_1_2 = Eigen::Matrix4d::Identity();
  T_1_2 *= zRotation(q[1]);
  T = T*T_1_2;
  jointPoses->push_back(T);

  Eigen::Matrix4d T_2_3 = Eigen::Matrix4d::Identity();
  T_2_3(1,3) = params_.pl2_;
  T_2_3(2,3) = params_.ph2_;
  T_2_3 *= zRotation(q[2]);
  T = T*T_2_3;
  jointPoses->push_back(T);

  Eigen::Matrix4d T_3_4 = Eigen::Matrix4d::Identity();
  T_3_4(1,3) = params_.pl3_;
  T_3_4(2,3) = -params_.ph3_;
  T_3_4 *= zRotation(q[3]);
  T = T*T_3_4;
  jointPoses->push_back(T);

  Eigen::Matrix4d T_4_5 = Eigen::Matrix4d::Identity();
  T_4_5.col(0) = e1;
  T_4_5.col(1) = -e3;
  T_4_5.col(2) = e2;
  T_4_5(1,3) = params_.pl4_;
  T_4_5(2,3) = -params_.ph4_;
  T_4_5 *= zRotation(q[4]);
  T = T*T_4_5;
  jointPoses->push_back(T);

  Eigen::Matrix4d T_5_6 = Eigen::Matrix4d::Identity();
  T_5_6.col(0) = e1;
  T_5_6.col(1) = e3;
  T_5_6.col(2) = -e2;
  T_5_6(2,3) = params_.pl5_;
  T_5_6 *= zRotation(q[5]);
  T = T*T_5_6;
  jointPoses->push_back(T);

  Eigen::Matrix4d T_6_RCM = Eigen::Matrix4d::Identity();
  T_6_RCM.col(0) = e3;
  T_6_RCM.col(1) = e1;
  T_6_RCM.col(2) = e2;
  T_6_RCM.topRightCorner<3,1>() = params_.pRCMOffset_;
  T = T*T_6_RCM;
  jointPoses->push_back(T);
}

void DavinciKinematics::getPassivePrimitives(const Eigen::Matrix4d& baseFrame,
                                             const std::vector<double>& q,
                                             std::vector<Collisions::Cylisphere>* cylispheres,
                                             Collisions::Sphere* sphere) const
{
  std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d> > poses;
  this->passiveFK(baseFrame, q, &poses);

  Collisions::Cylisphere c;
  c.r = params_.pLinkRadius_;
  c.p1 = poses[0].topRightCorner<3,1>();
  c.p2 = poses[1].topRightCorner<3,1>();
  cylispheres->push_back(c);

  c.p1 = c.p2;
  c.p2.noalias() = c.p1 + poses[1].block<3,1>(0,1)*params_.pl2_;
  cylispheres->push_back(c);

  c.p1 = c.p2;
  c.p2 = poses[2].topRightCorner<3,1>();
  cylispheres->push_back(c);

  c.p1 = poses[2].topRightCorner<3,1>();
  c.p2.noalias() = c.p1 + poses[2].block<3,1>(0,1)*params_.pl3_;
  cylispheres->push_back(c);

  // skip next link because it's small

  c.p1 = poses[4].topRightCorner<3,1>();
  c.p2 = poses[5].topRightCorner<3,1>();
  cylispheres->push_back(c);

  // add the sphere at L units up in z-direction from the prior
  // cylisphere's p2
  sphere->r = params_.eL_;
  sphere->p.noalias() = poses[5].topRightCorner<3,1>() + poses[5].block<3,1>(0,2)*params_.eL_;
}

DavinciParameters DavinciKinematics::getParams() const
{
  return params_;
}

// \TODO only construct one triangle of covariance matrix (symmetric
// obviously)
void DavinciKinematics::unscentedIK(const Eigen::Matrix4d& portFrame,
                                    const Eigen::Matrix4d& meanTargetPose,
                                    const Eigen::Vector3d& positionVariance,
                                    const Eigen::Vector3d& orientVariance,
                                    std::vector<double>* mean_qOut,
                                    Eigen::Matrix<double,6,6>* covariance_qOut) const
{
  Eigen::Matrix4d targetPose = meanTargetPose;
  std::vector<std::vector<double> > sigmaPoints;
  std::vector<double> q(6);

  // Perturb translation elements
  for (unsigned i = 0; i < 3; ++i)
    {
    double delta = sqrt(2*positionVariance(i));

    targetPose(i,3) += delta;
    this->intraIK(portFrame, targetPose, &q);
    sigmaPoints.push_back(q);

    targetPose(i,3) -= 2*delta;
    this->intraIK(portFrame, targetPose, &q);
    sigmaPoints.push_back(q);
    }
  
  // Perturb orientation elements
  double delta = sqrt(2*orientVariance(0));    
  Eigen::Matrix4d R = xRotation(delta);
  targetPose.noalias() = R*meanTargetPose;
  this->intraIK(portFrame, targetPose, &q);
  sigmaPoints.push_back(q);

  targetPose.noalias() = R.transpose()*meanTargetPose;
  this->intraIK(portFrame, targetPose, &q);
  sigmaPoints.push_back(q);

  delta = sqrt(2*orientVariance(1));
  R = yRotation(delta);
  targetPose.noalias() = R*meanTargetPose;
  this->intraIK(portFrame, targetPose, &q);
  sigmaPoints.push_back(q);

  targetPose.noalias() = R.transpose()*meanTargetPose;
  this->intraIK(portFrame, targetPose, &q);
  sigmaPoints.push_back(q);

  delta = sqrt(2*orientVariance(2));
  R = zRotation(delta);
  targetPose.noalias() = R*meanTargetPose;
  this->intraIK(portFrame, targetPose, &q);
  sigmaPoints.push_back(q);

  targetPose.noalias() = R.transpose()*meanTargetPose;
  this->intraIK(portFrame, targetPose, &q);
  sigmaPoints.push_back(q);
  // end orientation perturbations
  
  // get moments of transformed sigma points
  std::vector<double> mean_q(6, 0.0);
  for (unsigned i = 0; i < 12; ++i)
  {
    for (unsigned j = 0; j < 6; ++j)
      mean_q[j] += sigmaPoints[i][j];
  }
  for (unsigned i = 0; i < 6; ++i)
    mean_q[i] /= 12.0;

  Eigen::Matrix<double,6,6> covariance_q = Eigen::Matrix<double,6,6>::Zero();
  for (unsigned i = 0; i < 6; ++i)
  {
    for (unsigned j = 0; j < 6; ++j)
    {
      for (unsigned k = 0; k < 12; ++k)
        covariance_q(i,j) += (sigmaPoints[k][i] - mean_q[i])*(sigmaPoints[k][j] - mean_q[j]);
      covariance_q(i,j) /= (12 - 1);
    }
  }

  // populate out-parameters
  mean_qOut->swap(mean_q);
  (*covariance_qOut) = covariance_q;
}

namespace
{
  double fromPoseToExtraClearance(const DavinciKinematics& kin,
                                  const Eigen::Matrix4d& portFrameL,
                                  const Eigen::Matrix4d& portFrameR,
                                  const Eigen::Matrix4d& targetPose)
  {
    std::vector<double> qL(6);
    std::vector<double> qR(6);
    std::vector<Collisions::Cylisphere> cylL, cylR;
    std::vector<Collisions::Sphere> sL, sR;

    kin.intraIK(portFrameL, targetPose, &qL);
    kin.intraIK(portFrameR, targetPose, &qR);
    kin.getExtraCylispheres(portFrameL, qL, &cylL);
    kin.getExtraCylispheres(portFrameR, qR, &cylR);
    return Collisions::distance(cylL, sL, cylR, sR);
  }
}

// \TODO for efficiency compute unscented IK and unscented clearance
// at the same time
void DavinciKinematics::unscentedClearance(const Eigen::Matrix4d& portFrameL,
                                           const Eigen::Matrix4d& portFrameR,
                                           const Eigen::Matrix4d& meanTargetPose,
                                           const Eigen::Vector3d& positionVariance,
                                           const Eigen::Vector3d& orientVariance,
                                           double* mean_clearanceOut,
                                           double* covariance_clearanceOut) const
{
  Eigen::Matrix4d targetPose = meanTargetPose;
  std::vector<double> sigmaPoints;
  std::vector<double> qL(6);
  std::vector<double> qR(6);
  std::vector<Collisions::Cylisphere> cylL, cylR;
  std::vector<Collisions::Sphere> sL, sR;

  // Perturb translation elements
  for (unsigned i = 0; i < 3; ++i)
    {
    double delta = sqrt(2*positionVariance(i));
    targetPose(i,3) += delta;    
    sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));

    targetPose(i,3) -= 2*delta;
    sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));
    }
  
  // Perturb orientation elements
  double delta = sqrt(2*orientVariance(0));    
  Eigen::Matrix4d R = xRotation(delta);
  targetPose.noalias() = R*meanTargetPose;
  sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));

  targetPose.noalias() = R.transpose()*meanTargetPose;
  sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));

  delta = sqrt(2*orientVariance(1));
  R = yRotation(delta);
  targetPose.noalias() = R*meanTargetPose;
  sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));

  targetPose.noalias() = R.transpose()*meanTargetPose;
  sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));

  delta = sqrt(2*orientVariance(2));
  R = zRotation(delta);
  targetPose.noalias() = R*meanTargetPose;
  sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));

  targetPose.noalias() = R.transpose()*meanTargetPose;
  sigmaPoints.push_back(fromPoseToExtraClearance(*this, portFrameL, portFrameR, targetPose));
  // end orientation perturbations
  
  // get moments of transformed sigma points
  double mean_clearance = 0.0;
  for (unsigned i = 0; i < 12; ++i)
    mean_clearance += sigmaPoints[i];
  mean_clearance /= 12.0;

  double covariance_clearance = 0.0;
  for (unsigned i = 0; i < 12; ++i)
    covariance_clearance += (sigmaPoints[i] - mean_clearance)*(sigmaPoints[i] - mean_clearance);
  covariance_clearance /= (12 - 1);

  // populate out-parameters
  *mean_clearanceOut = mean_clearance;
  *covariance_clearanceOut = covariance_clearance;
}
