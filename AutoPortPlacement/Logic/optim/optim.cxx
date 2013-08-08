#include <optim/optim.h>

#include <nlopt.hpp>

#include <math.h>
#include <algorithm>
#include <iostream>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/constants/constants.hpp>

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
                          double* /*grad*/,
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
  lb[0] = 0.5; lb[6] = 0.5;
  opt.set_lower_bounds(lb);

  std::vector<double> ub(12, 3.14);
  ub[0] = 1.0;
  ub[6] = 1.0;
  opt.set_upper_bounds(ub);

  opt.set_max_objective(&PassiveLRProblem::wrapObj, &obj);

  std::vector<double> tol(6, 1e-6);
  opt.add_equality_mconstraint(&PassiveLRProblem::wrapConst, &obj, tol); // 1mm^2

  opt.set_xtol_rel(-1e-8);
  opt.set_xtol_abs(-1);
  opt.set_maxtime(10.0);
  
  // Use IK to find initial passive configurations which reach the
  // required RCM's as nice initial guesses
  std::vector<double> qL(6, 0.0);
  qL[0] = 0.5;
  kin.passiveIK(baseFrameL, rcmL, &qL);

  std::cout << "Left IK done!" << std::endl; // debug
  std::cout << "left error: " 
            << (kin.passiveFK(baseFrameL, qL).topRightCorner<3,1>() - rcmL).norm()
            << std::endl;

  std::vector<double> qR(6, 0.0);
  qR[0] = 0.5;
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
                                   std::vector<double>& /*grad*/)
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
                                  std::vector<double>& /*grad*/)
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

namespace
{
  // Represents a minimax optimization problem parameterized as in:
  // http://ab-initio.mit.edu/wiki/index.php/NLopt_Introduction#Equivalent_formulations_of_optimization_problems
  struct FeasiblePlanProblem
  {
    void portConstraint(const Eigen::Matrix4d& baseFrame,
                        const std::vector<double>& q,
                        double* c) const;

    void ikConstraint(const Eigen::Matrix4d& portFrame,
                      const Eigen::Matrix4d& taskFrame,
                      double spatialVariance,
                      double orientVariance,
                      std::vector<double>* c) const;

    void activeClearConstraint(const Eigen::Matrix4d& portFrameL,
                               const Eigen::Matrix4d& portFrameR,
                               const Eigen::Matrix4d& taskFrame,
                               double spatialVariance,
                               double orientVariance,
                               double* c) const;

    void passiveClearConstraint(const Eigen::Matrix4d& baseFrameL,
                                const Eigen::Matrix4d& baseFrameR,
                                const std::vector<double>& qL,
                                const std::vector<double>& qR,
                                double *c) const;

    const DavinciKinematics& kin;
    const Eigen::Matrix4d& baseFrameL;
    const Eigen::Matrix4d& baseFrameR;
    const Optim::Matrix4dVec& taskFrames;
    const Eigen::Vector3d& portCurvePoint1;
    const Eigen::Vector3d& portCurvePoint2;
    double chanceConstraint;

    static void wrapIneq(unsigned int m,
                         double* result,
                         unsigned int n,
                         const double* x,
                         double* grad,
                         void* data);

    static double feasibleMinimaxObj(const std::vector<double>& x,
                                     std::vector<double>& grad,
                                     void* data);
  };

  const double activeLowerBounds[] = {-boost::math::constants::pi<double>()/3,
                                      -boost::math::constants::pi<double>()/3,
                                      0.0,
                                      -boost::math::constants::pi<double>(),
                                      -boost::math::constants::pi<double>()/2,
                                      -boost::math::constants::pi<double>()/2};

  // For now I'm doing symmetric bounds, and making up a bound for the
  // last joint. Accoring to Azimian's implementation, q[1]'s limit is
  // actually pi/4, which is strange and should be checked up on.
  const double  activeUpperBounds[] = {boost::math::constants::pi<double>()/3,
                                       boost::math::constants::pi<double>()/3,
                                       0.2,
                                       boost::math::constants::pi<double>(),
                                       boost::math::constants::pi<double>()/2,
                                       boost::math::constants::pi<double>()/2};
}

bool Optim::findFeasiblePlan(const DavinciKinematics& kin,
                             const Eigen::Matrix4d& baseFrameL,
                             const Eigen::Matrix4d& baseFrameR,
                             const Matrix4dVec& taskFrames,
                             const Eigen::Vector3d& portCurvePoint1,
                             const Eigen::Vector3d& portCurvePoint2,
                             double chanceConstraint,
                             std::vector<double>* qL_out,
                             std::vector<double>* qR_out,
                             double* spatialVariance,
                             double* orientVariance)
{
  FeasiblePlanProblem problem = {kin, baseFrameL, baseFrameR, taskFrames, 
                                 portCurvePoint1, portCurvePoint2, chanceConstraint};

  // \TODO try another algorithm
  nlopt::opt opt(nlopt::GN_ISRES, 15);
  
  std::vector<double> lb(15);
  lb[0] = 0.5;
  double twoThirdsPi = (2.0/3.0)*boost::math::constants::pi<double>();
  std::fill(lb.begin()+1, lb.begin()+6, -twoThirdsPi);
  lb[6] = 0.5;
  std::fill(lb.begin()+7, lb.begin()+12, -twoThirdsPi);
  lb[12] = 0.0000001;
  lb[13] = 0.0000001;
  lb[14] = -100.0; // trying to make an order of magnitude or 2 greater than constraint values

  std::vector<double> ub(15);
  ub[0] = 1.0;
  std::fill(ub.begin()+1, ub.begin()+6, twoThirdsPi);
  ub[6] = 1.0;
  std::fill(ub.begin()+7, ub.begin()+12, twoThirdsPi);
  ub[12] = 0.0001;
  ub[13] = 0.0001;
  ub[14] = 100.0; // trying to make an order of magnitude or 2 greater than constraint values

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  opt.set_min_objective(&FeasiblePlanProblem::feasibleMinimaxObj, 0);

  std::size_t numIneqConstraints = 3 + 13*taskFrames.size();
  std::vector<double> ineqTol(numIneqConstraints, 0.0000001);
  opt.add_inequality_mconstraint(&FeasiblePlanProblem::wrapIneq, 
                                 (void*) &problem,
                                 ineqTol);

  // Stop the optimization once we've found a point that satisfies the constraints (t <= 0)
  opt.set_stopval(0.0);

  // Stop the optimization after 1000 seconds have passed no matter what
  opt.set_maxtime(100.0);

  // Use Jacobian IK to find a nice initial guess
  // Place RCM's at middle of port curve
  Eigen::Vector3d rcm = 0.5*(portCurvePoint1 + portCurvePoint2);
  std::vector<double> qL(6, 0.0);
  qL[0] = 0.5;
  kin.passiveIK(baseFrameL, rcm, &qL);

  std::cout << "Left IK done!" << std::endl; // debug
  std::cout << "left error: " 
            << (kin.passiveFK(baseFrameL, qL).topRightCorner<3,1>() - rcm).norm()
            << std::endl;

  std::vector<double> qR(6, 0.0);
  qR[0] = 0.5;
  kin.passiveIK(baseFrameR, rcm, &qR);

  std::cout << "Right IK done!" << std::endl; // debug
  std::cout << "right error: " 
            << (kin.passiveFK(baseFrameR, qR).topRightCorner<3,1>() - rcm).norm()
            << std::endl;

  std::vector<double> x(15);
  std::copy(qL.begin(), qL.end(), x.begin());
  std::copy(qR.begin(), qR.end(), x.begin()+6);
  x[12] = x[13] = 0.0000001;
  x[14] = 90.0;

  // Output initial objective and constraint costs
  std::vector<double> dummy;
  std::cout << "f(x0): " << FeasiblePlanProblem::feasibleMinimaxObj(x, dummy, (void*) &problem)  << std::endl;

  double* c_ineq = new double[numIneqConstraints];
  double x_array[15];
  for (unsigned i = 0; i < 15; ++i)
    x_array[i] = x[i];
  FeasiblePlanProblem::wrapIneq(numIneqConstraints, c_ineq, 15, x_array, 0, (void*) &problem);
  std::cout << "c_ineq(x0):";
  for (unsigned i = 0; i < numIneqConstraints; ++i)
    std::cout << " " << (c_ineq[i] + x[14]);
  std::cout << std::endl;

  double minf;
  nlopt::result result = opt.optimize(x, minf);
  std::cout << "result: " << result << std::endl;

  std::cout << "f(x): " << FeasiblePlanProblem::feasibleMinimaxObj(x, dummy, (void*) &problem) << std::endl;

  for (unsigned i = 0; i < 15; ++i)
    x_array[i] = x[i];
  FeasiblePlanProblem::wrapIneq(numIneqConstraints, c_ineq, 15, x_array, 0, (void*) &problem);
  std::cout << "c_ineq(x):";
  for (unsigned i = 0; i < numIneqConstraints; ++i)
    std::cout << " " << (c_ineq[i] + x[14]);
  std::cout << std::endl;

  std::cout << "x:";
  for (unsigned i = 0; i < x.size(); ++i)
    std::cout << " " << x[i];
  std::cout << std::endl;

  delete c_ineq;

  for (unsigned i = 0; i < 6; ++i)
    {
    (*qL_out)[i] = x[i];
    (*qR_out)[i] = x[i+6];
    }
  *spatialVariance = x[12];
  *orientVariance = x[13];

  return true;
}

void FeasiblePlanProblem::portConstraint(const Eigen::Matrix4d& baseFrame,
                                         const std::vector<double>& q,
                                         double* c) const
{
  Eigen::Matrix4d rcm = this->kin.passiveFK(baseFrame, q);
  double d = Collisions::distance(this->portCurvePoint1, this->portCurvePoint2,
                                  rcm.topRightCorner<3,1>());
  *c = d*d - 0.001*0.001;
}

void FeasiblePlanProblem::ikConstraint(const Eigen::Matrix4d& portFrame,
                                       const Eigen::Matrix4d& taskFrame,
                                       double spatialVariance,
                                       double orientVariance,
                                       std::vector<double>* c) const
{
  std::vector<double> mean_q;
  Eigen::Matrix<double,6,6> covariance_q;
  this->kin.unscentedIK(portFrame, taskFrame, 
                        Eigen::Vector3d::Constant(spatialVariance), 
                        Eigen::Vector3d::Constant(orientVariance),
                        &mean_q, &covariance_q);

  boost::math::normal n;
  double quantile = boost::math::quantile(n, 1 - this->chanceConstraint);
  for (unsigned i = 0; i < 6; ++i)
    {
    double midpt = 0.5*(activeUpperBounds[i] - activeLowerBounds[i]);
    (*c)[i] = fabs(mean_q[i] - midpt) + quantile*sqrt(covariance_q(i,i)) - midpt;
    }
}

void FeasiblePlanProblem::activeClearConstraint(const Eigen::Matrix4d& portFrameL,
                                                const Eigen::Matrix4d& portFrameR,
                                                const Eigen::Matrix4d& taskFrame,
                                                double spatialVariance,
                                                double orientVariance,
                                                double* c) const
{
  double mean_d, variance_d;
  this->kin.unscentedClearance(portFrameL, portFrameR, taskFrame,
                               Eigen::Vector3d::Constant(spatialVariance),
                               Eigen::Vector3d::Constant(orientVariance),
                               &mean_d, &variance_d);

  boost::math::normal n;
  double quantile = boost::math::quantile(n, 1 - this->chanceConstraint);
  *c = quantile*sqrt(orientVariance) - mean_d;
}

void FeasiblePlanProblem::passiveClearConstraint(const Eigen::Matrix4d& baseFrameL,
                                                 const Eigen::Matrix4d& baseFrameR,
                                                 const std::vector<double>& qL,
                                                 const std::vector<double>& qR,
                                                 double *c) const
{
  std::vector<Collisions::Cylisphere> cL, cR;
  std::vector<Collisions::Sphere> sL(1), sR(1);
  this->kin.getPassivePrimitives(baseFrameL, qL, &cL, &sL[0]);
  this->kin.getPassivePrimitives(baseFrameR, qR, &cR, &sR[0]);
  double d = Collisions::distance(cL, sL, cR, sR);
  *c = -d*d;
}

// inequality constraints:
// 1 passive clear constraint
// k*6*2 ik constraints
// k*1 active clear constraints
void FeasiblePlanProblem::wrapIneq(unsigned int m,
                                   double* result,
                                   unsigned int n,
                                   const double* x,
                                   double* /*grad*/,
                                   void* data)
{ 
  FeasiblePlanProblem* problem = reinterpret_cast<FeasiblePlanProblem*>(data);

  if (m != 3 + (problem->taskFrames.size())*13)
    throw std::runtime_error("wrapIneq Error: wrong number of constraints!");

  if (n != 15)
    throw std::runtime_error("wrapIneq Error: wrong number of variables!");
  
  std::vector<double> qL(&x[0], &x[6]);
  std::vector<double> qR(&x[6], &x[12]);

  // passive clear constraint
  problem->passiveClearConstraint(problem->baseFrameL, problem->baseFrameR,
                                  qL, qR, &result[0]); 

  // port constraints
  problem->portConstraint(problem->baseFrameL, qL, &result[1]);
  problem->portConstraint(problem->baseFrameR, qR, &result[2]);

  // ik constraints
  std::vector<double> c(6);
  Eigen::Matrix4d portFrameL = problem->kin.passiveFK(problem->baseFrameL, qL);
  Eigen::Matrix4d portFrameR = problem->kin.passiveFK(problem->baseFrameR, qR);
  double spatialVariance = x[12];
  double orientVariance = x[13];
  for (std::size_t k = 0; k < problem->taskFrames.size(); ++k)
    {
    problem->ikConstraint(portFrameL, 
                          problem->taskFrames[k], 
                          spatialVariance,
                          orientVariance,
                          &c);
    for (unsigned i = 0; i < 6; ++i)
      result[3+(k*12)+i] = c[i];

    problem->ikConstraint(portFrameR,
                          problem->taskFrames[k],
                          spatialVariance,
                          orientVariance,
                          &c);
    for (unsigned i = 0; i < 6; ++i)
      result[3+(k*12)+i+6] = c[i];
    }

  // active clearance constraints
  for (unsigned k = 0; k < problem->taskFrames.size(); ++k)
    {
    problem->activeClearConstraint(portFrameL, portFrameR, problem->taskFrames[k],
                                   spatialVariance, orientVariance,
                                   &result[3+(problem->taskFrames.size()*12)+k]);
    }

  // make sure to subtract away t from all constraints
  double t = x[14];
  for (unsigned i = 0; i < m; ++i)
    result[i] -= t;
}

double FeasiblePlanProblem::feasibleMinimaxObj(const std::vector<double>& x,
                                               std::vector<double>& /*grad*/,
                                               void* /*data*/)
{
  return x[14];
}
