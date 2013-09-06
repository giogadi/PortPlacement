#ifndef PORT_PLACEMENT_AUTO_PORT_PLACEMENT_LOGIC_OPTIM_OPTIM_H
#define PORT_PLACEMENT_AUTO_PORT_PLACEMENT_LOGIC_OPTIM_OPTIM_H

#include <davinci-kinematics/davinci.h>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace Optim
{
  typedef std::vector<Eigen::Matrix4d,
    Eigen::aligned_allocator<Eigen::Matrix4d> > Matrix4dVec;

  bool findCollisionFreePassiveLR(const DavinciKinematics& kin,
                                  const Eigen::Matrix4d& baseFrameL,
                                  const Eigen::Matrix4d& baseFrameR,
                                  const Eigen::Vector3d& rcmL,
                                  const Eigen::Vector3d& rcmR,
                                  std::vector<double>* qL,
                                  std::vector<double>* qR);

  // Orientations of task frames are such that +x is the tip tangent,
  // +z is the last joint's axis
  bool findFeasiblePlan(const DavinciKinematics& kin,
                        const Eigen::Matrix4d& baseFrameL,
                        const Eigen::Matrix4d& baseFrameR,
                        const Matrix4dVec& taskFrames,
                        const Eigen::Vector3d& portCurvePoint1,
                        const Eigen::Vector3d& portCurvePoint2,
                        std::vector<double>* qL,
                        std::vector<double>* qR);
};

#endif
