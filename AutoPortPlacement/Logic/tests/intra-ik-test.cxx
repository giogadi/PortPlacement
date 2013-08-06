#include <davinci-kinematics/davinci.h>
#include <iostream>

int main(int argc, char** argv)
{
  DavinciKinematics kinematics;

  Eigen::Matrix4d worldPose = Eigen::Matrix4d::Identity();
  worldPose(0,0) = -1.0;
  worldPose(1,1) = -1.0;
  worldPose(0,3) = -0.3;

  std::vector<double> q(6);
  kinematics.intraIK(Eigen::Matrix4d::Identity(),
                     worldPose,
                     &q);

  for (std::size_t i = 0; i < q.size(); ++i)
    std::cout << q[i] << std::endl;

  if (fabs(q[0]) <= 0.0000001 &&
      fabs(q[1]) <= 0.0000001 &&
      fabs(q[2]) <= 0.0000001 &&
      fabs(q[3] - 0.29) <= 0.0000001 &&
      fabs(q[4]) <= 0.0000001 &&
      fabs(q[5]) <= 0.0000001)
    return 0;
  else
    return 1;
}
