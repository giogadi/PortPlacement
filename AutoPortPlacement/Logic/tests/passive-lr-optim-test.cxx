#include <optim/optim.h>

int main()
{
  DavinciKinematics kinematics;

  Eigen::Matrix4d baseFrameL = Eigen::Matrix4d::Identity();
  Eigen::Vector3d e1 = baseFrameL.block<3,1>(0,0);
  Eigen::Vector3d e2 = baseFrameL.block<3,1>(0,1);
  Eigen::Vector3d e3 = baseFrameL.block<3,1>(0,2);

  // Rotate about +z by 90 degrees, translate -0.5 in x
  baseFrameL.block<3,1>(0,0) = e2;
  baseFrameL.block<3,1>(0,1) = -e1;
  baseFrameL.block<3,1>(0,2) = e3;
  baseFrameL.topRightCorner<3,1>() = -0.5*e1;

  // Rotate about +z by -90 degrees, translate 0.5 in x
  Eigen::Matrix4d baseFrameR = Eigen::Matrix4d::Identity();
  baseFrameR.block<3,1>(0,0) = -e2;
  baseFrameR.block<3,1>(0,1) = e1;
  baseFrameR.block<3,1>(0,2) = e3;
  baseFrameR.topRightCorner<3,1>() = 0.5*e1;

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
};
