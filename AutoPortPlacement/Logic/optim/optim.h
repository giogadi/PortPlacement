#ifndef PORT_PLACEMENT_AUTO_PORT_PLACEMENT_LOGIC_OPTIM_OPTIM_H
#define PORT_PLACEMENT_AUTO_PORT_PLACEMENT_LOGIC_OPTIM_OPTIM_H

#include <davinci-kinematics/davinci.h>
#include <Eigen/Dense>

namespace Optim
{
  bool findCollisionFreePassiveLR(const DavinciKinematics& kin,
                                  const Eigen::Matrix4d& baseFrameL,
                                  const Eigen::Matrix4d& baseFrameR,
                                  const Eigen::Vector3d& rcmL,
                                  const Eigen::Vector3d& rcmR,
                                  std::vector<double>* qL,
                                  std::vector<double>* qR);
};

#endif
