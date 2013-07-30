#ifndef DAVINCI_KINEMATICS_H
#define DAVINCI_KINEMATICS_H

#include <Eigen/Dense>
#include <vector>

#include <collisions/collisions.h>

// Attributes prefixed with 'e' are extracorporeal mechanism
// parameters, and those prefixed with 'p' are passive arm parameters.
struct DavinciParameters
{ 
  double wristLength_;
  double gripperLength_;

  // dimensions of active mechanism (using terminology of Azimian's
  // thesis)
  double eL_;
  double el1_;
  double el2_;
  double eAlpha_;

  // radii for collision primitives
  double eRadius1_; // radius of primitive nearest the port
  double eRadius2_; // radius of primitive with double parallelogram

  double pl2_;
  double ph2_;
  double pl3_;
  double ph3_;
  double pl4_;
  double ph4_;
  double pl5_;
  Eigen::Vector3d pRCMOffset_;
  double pLinkRadius_;  
};

class DavinciKinematics
{
public:
  DavinciKinematics(const std::string& inputFilename);

  // Get frame at the end of the wrist (base of the gripper) given the
  // pose of the RCM and configuration q.
  //
  // Port frame is positioned at RCM of the instrument mechanism with
  // the pose as follows (at home position, i.e. all rotations zero):
  // +x : "up" the instrument from the RCM
  // +z : "forward" from mechanism's perspective, along double parallelogram
  Eigen::Matrix4d intraFK(const Eigen::Matrix4d& portFrame,
                          const std::vector<double>& q) const;

  // Get frame of the RCM at the end of a passive section with a given
  // baseFrame and configuration q.
  //
  // Where (0,0,0) is on the floor at the center of the entire da
  // Vinci's base, baseFrame's (x,y) is positioned so that the
  // z-component is centered at the first prismatic axis of motion
  // (therefore x-y will vary for each arm), and the z-component is 0.
  Eigen::Matrix4d passiveFK(const Eigen::Matrix4d& baseFrame,
                            const std::vector<double>& q) const;

  // Given that the intracorporeal mechanism is placed such that its
  // RCM is at portFrame, and given a goal pose worldPose for the
  // wrist, find the configuration q which gets us there.
  void intraIK(const Eigen::Matrix4d& portFrame,
               const Eigen::Matrix4d& worldPose,
               std::vector<double>* q) const;

  // Given the frame of the base of the robot and a configuration,
  // compute robot's positional Jacobian.
  //
  // Currently uses finite differences.
  Eigen::Matrix<double, 3, 6> passiveJacobian(const Eigen::Matrix4d& baseFrame,
                                              const std::vector<double>& q,
                                              double stepSize) const;

  bool passiveIK(const Eigen::Matrix4d& baseFrame,
                 const Eigen::Vector3d& worldPosition,
                 std::vector<double>* q,
                 double stepSize=0.001) const;

  // Get geometric primitives representing the active extracorporeal
  // mechanism
  void getExtraCylispheres(const Eigen::Matrix4d& portFrame,
                           const std::vector<double>& q,
                           std::vector<Collisions::Cylisphere>* cylispheres) const;

  void getPassivePrimitives(const Eigen::Matrix4d& baseFrame,
                            const std::vector<double>& q,
                            std::vector<Collisions::Cylisphere>* cylispheres,
                            Collisions::Sphere* sphere) const;

  DavinciParameters getParams() const;

private:
  DavinciParameters params_;

  void passiveFK(const Eigen::Matrix4d& baseFrame,
                 const std::vector<double>& q,
                 std::vector<Eigen::Matrix4d>* jointPoses) const;
};

#endif
