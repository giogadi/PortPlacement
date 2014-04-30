#include <optim/optim.h>

int main()
{
  DavinciKinematics kinematics;

  Eigen::Matrix4d baseFrameL = Eigen::Matrix4d::Identity();

  // Rotate about +z by 180 degrees
  Eigen::Matrix4d baseFrameR = Eigen::Matrix4d::Identity();
  baseFrameR(0,0) = -1.0;
  baseFrameR(1,1) = -1.0;

  Eigen::Vector3d rcm;
  rcm(0) = 0.0;
  rcm(1) = 0.8;
  rcm(2) = 0.7;

  std::vector<double> qL(6);
  std::vector<double> qR(6);
  if (Optim::findCollisionFreePassiveLR(kinematics, baseFrameL, baseFrameR, rcm, rcm, &qL, &qR))
    return 0;
  else
    return 1;
}
