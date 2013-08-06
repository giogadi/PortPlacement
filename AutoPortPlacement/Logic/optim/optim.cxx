#include <optim/optim.h>

#include <nlopt.hpp>

#include <math.h>
#include <algorithm>
#include <iostream>

namespace
{
  struct PassiveLRProblem
  {
    // First 6 parameters of x are configuration variables for left
    // arm, second 6 are for right arm.
    double objective(const std::vector<double>& x,
                     std::vector<double>& grad);

    void constraint(std::vector<double>& c,
                    const std::vector<double>& x,
                    std::vector<double>& grad);

    const DavinciKinematics& kin_;
    const Eigen::Matrix4d& baseFrameL_;
    const Eigen::Matrix4d& baseFrameR_;
    const Eigen::Vector3d& rcmL_;
    const Eigen::Vector3d& rcmR_;

    static double wrapObj(const std::vector<double>& x,
                          std::vector<double>& grad,
                          void* data)
    {
      return reinterpret_cast<PassiveLRProblem*>(data)->objective(x, grad);
    }
    static void wrapConst(unsigned int m,
                          double* result,
                          unsigned int n,
                          const double* x,
                          double* grad,
                          void* data)
    {
      std::vector<double> c(m);
      std::vector<double> x_v(n);      
      for (unsigned i = 0; i < n; ++i)
        x_v[i] = x[i];
      std::vector<double> grad_v;

      reinterpret_cast<PassiveLRProblem*>(data)->constraint(c, x_v, grad_v);

      for (unsigned i = 0; i < m; ++i)
        result[i] = c[i];
    }
  };
}

bool Optim::findCollisionFreePassiveLR(const DavinciKinematics& kin,
                                       const Eigen::Matrix4d& baseFrameL,
                                       const Eigen::Matrix4d& baseFrameR,
                                       const Eigen::Vector3d& rcmL,
                                       const Eigen::Vector3d& rcmR,
                                       std::vector<double>* qLOut,
                                       std::vector<double>* qROut)
{
  PassiveLRProblem obj = { kin, baseFrameL, baseFrameR, rcmL, rcmR };

  // nlopt::opt opt(nlopt::LN_COBYLA, 12);
  nlopt::opt opt(nlopt::GN_ISRES, 12);

  std::vector<double> lb(12, -3.14);
  lb[0] = 0.1; lb[6] = 0.1;
  opt.set_lower_bounds(lb);

  std::vector<double> ub(12, 3.14);
  opt.set_upper_bounds(ub);

  opt.set_max_objective(&PassiveLRProblem::wrapObj, &obj);

  std::vector<double> tol(6, 1e-6);
  opt.add_equality_mconstraint(&PassiveLRProblem::wrapConst, &obj, tol); // 1mm^2

  opt.set_xtol_rel(-1e-8);
  opt.set_xtol_abs(-1);
  opt.set_maxtime(30.0);
  
  // Use IK to find initial passive configurations which reach the
  // required RCM's as nice initial guesses
  std::vector<double> qL(6, 0.0);
  qL[0] = 0.1;
  kin.passiveIK(baseFrameL, rcmL, &qL);

  std::cout << "Left IK done!" << std::endl; // debug
  std::cout << "left error: " 
            << (kin.passiveFK(baseFrameL, qL).topRightCorner<3,1>() - rcmL).norm()
            << std::endl;

  std::vector<double> qR(6, 0.0);
  qR[0] = 0.1;
  kin.passiveIK(baseFrameR, rcmR, &qR);

  std::cout << "Right IK done!" << std::endl; // debug
  std::cout << "right error: " 
            << (kin.passiveFK(baseFrameR, qR).topRightCorner<3,1>() - rcmR).norm()
            << std::endl;

  std::vector<double> x(12);
  std::copy(qL.begin(), qL.end(), x.begin());
  std::copy(qR.begin(), qR.end(), x.begin()+6);

  for (unsigned i = 0; i < 12; ++i)
    std::cout << x[i] << " ";
  std::cout << std::endl;

  std::vector<double> dummy;
  std::cout << "f(x0): " << obj.objective(x, dummy) << std::endl;
  // std::cout << "c1(x0): " << obj.constraintL(x, dummy) << std::endl;
  // std::cout << "c2(x0): " << obj.constraintR(x, dummy) << std::endl;

  double maxf;
  nlopt::result result = opt.optimize(x, maxf);

  std::cout << "result: " << " " << result << std::endl;

  std::cout << "maximum distance: " << maxf << std::endl; // debug

  std::vector<double> c(6);
  obj.constraint(c, x, dummy);
  std::cout << "c(x):";
  for (unsigned i = 0; i < 6; ++i)
    std::cout << " " << c[i];
  std::cout << std::endl;

  std::copy(x.begin(), x.begin()+6, qLOut->begin());
  std::copy(x.begin()+6, x.end(), qROut->begin());

  std::cout << "Left:" << std::endl;
  for (unsigned i = 0; i < 6; ++i)
    std::cout << " " << (*qLOut)[i];
  std::cout << std::endl;

  std::cout << "Right:" << std::endl;
  for (unsigned i = 0; i < 6; ++i)
    std::cout << " " << (*qROut)[i];
  std::cout << std::endl;

  return true;
}

double PassiveLRProblem::objective(const std::vector<double>& x,
                                   std::vector<double>& grad)
{
  // Get collision primitives for left arm
  std::vector<double> qL(x.begin(), x.begin()+6);
  std::vector<Collisions::Cylisphere> cL;
  std::vector<Collisions::Sphere> sL(1);
  kin_.getPassivePrimitives(baseFrameL_, qL, &cL, &(sL[0]));

  // Get collision primitives for right arm
  std::vector<double> qR(x.begin()+6, x.end());
  std::vector<Collisions::Cylisphere> cR;
  std::vector<Collisions::Sphere> sR(1);
  kin_.getPassivePrimitives(baseFrameR_, qR, &cR, &(sR[0]));
  
  return Collisions::distance(cL, sL, cR, sR);
}

void PassiveLRProblem::constraint(std::vector<double>& m,
                                  const std::vector<double>& x,
                                  std::vector<double>& grad)
{
  std::vector<double> qL(x.begin(), x.begin()+6);
  Eigen::Vector3d error = kin_.passiveFK(baseFrameL_, qL).topRightCorner<3,1>() - rcmL_;
  m[0] = error(0)*error(0);
  m[1] = error(1)*error(1);
  m[2] = error(2)*error(2);

  std::vector<double> qR(x.begin()+6, x.end());
  error = kin_.passiveFK(baseFrameR_, qR).topRightCorner<3,1>() - rcmR_;
  m[3] = error(0)*error(0);
  m[4] = error(1)*error(1);
  m[5] = error(2)*error(2);
}
