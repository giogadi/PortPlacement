#include <davinci-kinematics/davinci.h>
#include <iostream>

int main(int argc, char** argv)
{
  DavinciKinematics kinematics;

  Eigen::Matrix4d meanTargetPose = Eigen::Matrix4d::Identity();
  meanTargetPose(0,0) = -1.0;
  meanTargetPose(1,1) = -1.0;
  meanTargetPose(0,3) = -0.3;

  Eigen::Vector3d positionVariance = Eigen::Vector3d::Constant(0.0001);
  Eigen::Vector3d orientVariance = Eigen::Vector3d::Constant(0.007);
  std::vector<double> mean_q(6);
  Eigen::Matrix<double,6,6> cov_q;

  kinematics.unscentedIK(Eigen::Matrix4d::Identity(),
                         meanTargetPose,
                         positionVariance,
                         orientVariance,
                         &mean_q,
                         &cov_q);

  std::cout << "mu:" << std::endl;
  for (std::size_t i = 0; i < mean_q.size(); ++i)
    std::cout << mean_q[i] << std::endl;
  std::cout << std::endl;

  std::cout << "sigma:" << std::endl;
  std::cout << cov_q << std::endl;

  return 0;
}
