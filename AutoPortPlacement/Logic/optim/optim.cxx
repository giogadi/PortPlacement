#include <optim/optim.h>

#include <nlopt.hpp>

#include <math.h>
#include <algorithm>
#include <iostream>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

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
                      std::vector<double>* c) const;

    void activeClearConstraint(const std::vector<double>& qL,
                               const std::vector<double>& qR,
                               const Eigen::Matrix4d& taskFrame,
                               std::vector<double>* c) const;

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

    static void wrapIneq(unsigned int m,
                         double* result,
                         unsigned int n,
                         const double* x,
                         double* grad,
                         void* data);

    static double feasibleMinimaxObj(const std::vector<double>& x,
                                     std::vector<double>& grad,
                                     void* data);

    void getBounds(std::vector<double>* lb, std::vector<double>* ub) const;

    void getInitialGuessIK(std::vector<double>* x0) const;
    void getInitialGuessExact(std::vector<double>* x0) const;

    void outputStateProperties(std::ostream& out, const std::vector<double>& x) const;

    unsigned getNumVariables() const;
    unsigned getNumConstraints() const;
  };

  const double activeLowerBounds[] = {-boost::math::constants::pi<double>()/3,
                                      -boost::math::constants::pi<double>()/3,
                                      -boost::math::constants::pi<double>(),
                                      0.0,
                                      -boost::math::constants::pi<double>()/2,
                                      -boost::math::constants::pi<double>()/2};

  // For now I'm doing symmetric bounds, and making up a bound for the
  // last joint. Accoring to Azimian's implementation, q[1]'s limit is
  // actually pi/4, which is strange and should be checked up on.
  const double  activeUpperBounds[] = {boost::math::constants::pi<double>()/3,
                                       boost::math::constants::pi<double>()/3,
                                       boost::math::constants::pi<double>(),
                                       0.2,
                                       boost::math::constants::pi<double>()/2,
                                       boost::math::constants::pi<double>()/2};

  const double passiveLowerBounds[] = {0.5,
                                       -(2.0/3.0)*boost::math::constants::pi<double>(),
                                       -(2.0/3.0)*boost::math::constants::pi<double>(),
                                       -(2.0/3.0)*boost::math::constants::pi<double>(),
                                       -(2.0/3.0)*boost::math::constants::pi<double>(),
                                       -(2.0/3.0)*boost::math::constants::pi<double>()};

  const double passiveUpperBounds[] = {1.0,
                                       (2.0/3.0)*boost::math::constants::pi<double>(),
                                       (2.0/3.0)*boost::math::constants::pi<double>(),
                                       (2.0/3.0)*boost::math::constants::pi<double>(),
                                       (2.0/3.0)*boost::math::constants::pi<double>(),
                                       (2.0/3.0)*boost::math::constants::pi<double>()};
}

bool Optim::findFeasiblePlan(const DavinciKinematics& kin,
                             const Eigen::Matrix4d& baseFrameL,
                             const Eigen::Matrix4d& baseFrameR,
                             const Matrix4dVec& taskFrames,
                             const Eigen::Vector3d& portCurvePoint1,
                             const Eigen::Vector3d& portCurvePoint2,
                             std::vector<double>* qL_out,
                             std::vector<double>* qR_out)
{
  FeasiblePlanProblem problem = {kin, baseFrameL, baseFrameR, taskFrames, 
                                 portCurvePoint1, portCurvePoint2};

  const unsigned numVariables = problem.getNumVariables();

  nlopt::opt opt(nlopt::LN_COBYLA, numVariables);
  
  std::vector<double> lb(numVariables);
  std::vector<double> ub(numVariables);
  problem.getBounds(&lb, &ub);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  opt.set_min_objective(&FeasiblePlanProblem::feasibleMinimaxObj, 0);

  std::size_t numIneqConstraints = problem.getNumConstraints();
  std::vector<double> ineqTol(numIneqConstraints, 0.0000001);
  opt.add_inequality_mconstraint(&FeasiblePlanProblem::wrapIneq, 
                                 (void*) &problem,
                                 ineqTol);

  // Stop the optimization once we've found a point that satisfies the constraints (t <= 0)
  // opt.set_stopval(0.0);

  // Stop the optimization after some seconds have passed no matter what
  opt.set_maxtime(100.0);

  std::vector<double> x(numVariables);
  problem.getInitialGuessIK(&x);

  // Output initial value, objective and constraint costs
  std::cout << "================== x0 =================" << std::endl;
  problem.outputStateProperties(std::cout, x);

  double minf;
  nlopt::result result = opt.optimize(x, minf);

  std::cout << "=================" << std::endl;
  std::cout << "Result: " << result << std::endl;
  std::cout << "=================" << std::endl;

  std::cout << "================== x_opt ==============" << std::endl;
  problem.outputStateProperties(std::cout, x);

  std::copy(x.begin(), x.begin()+6, qL_out->begin());
  std::copy(x.begin()+6, x.begin()+12, qR_out->begin());

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
                                       std::vector<double>* c) const
{
  std::vector<double> q(6);
  this->kin.intraIK(portFrame, taskFrame, &q);
  for (unsigned i = 0; i < 6; ++i)
    {
    (*c)[2*i] = activeLowerBounds[i] - q[i];
    (*c)[2*i+1] = q[i] - activeUpperBounds[i];
    }
}

void FeasiblePlanProblem::activeClearConstraint(const std::vector<double>& qL,
                                                const std::vector<double>& qR,
                                                const Eigen::Matrix4d& taskFrame,
                                                std::vector<double>* c) const
{
  std::vector<double> dists;
  this->kin.fullClearances(this->baseFrameL, this->baseFrameR, qL, qR, taskFrame, &dists);
  for (std::size_t i = 0; i < dists.size(); ++i)
    dists[i] = -dists[i]; // to make constraint in form of c(x) <= 0
  c->swap(dists);
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
  // *c = exp(-d) - 1;
  *c = -d;
}

// inequality constraints:
// 1 passive clear constraint
// 2 port constraints
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

  if (m != problem->getNumConstraints())
    throw std::runtime_error("wrapIneq Error: wrong number of constraints!");

  if (n != problem->getNumVariables())
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
  std::vector<double> c(6*2);
  Eigen::Matrix4d portFrameL = problem->kin.passiveFK(problem->baseFrameL, qL);
  Eigen::Matrix4d portFrameR = problem->kin.passiveFK(problem->baseFrameR, qR);
  for (std::size_t k = 0; k < problem->taskFrames.size(); ++k)
    {
    problem->ikConstraint(portFrameL, 
                          problem->taskFrames[k], 
                          &c);
    for (unsigned i = 0; i < 6*2; ++i)
      result[3+(k*2*6*2)+i] = c[i];

    problem->ikConstraint(portFrameR,
                          problem->taskFrames[k],
                          &c);
    for (unsigned i = 0; i < 6*2; ++i)
      result[3+(k*2*6*2)+i+6*2] = c[i];
    }

  // active clearance constraints
  for (unsigned k = 0; k < problem->taskFrames.size(); ++k)
    {
    std::vector<double> c(57);
    problem->activeClearConstraint(qL, qR, problem->taskFrames[k],
                                   &c);
    for (unsigned i = 0; i < 57; ++i)
      result[3+(problem->taskFrames.size()*2*6*2) + 57*k + i] = c[i];
    }

  // make sure to subtract away t from all constraints
  double t = x[12];
  for (unsigned i = 0; i < m; ++i)
    result[i] -= t;
}

double FeasiblePlanProblem::feasibleMinimaxObj(const std::vector<double>& x,
                                               std::vector<double>& /*grad*/,
                                               void* /*data*/)
{
  return x[12];
}

void FeasiblePlanProblem::getBounds(std::vector<double>* lb, 
                                    std::vector<double>* ub) const
{
  for (unsigned i = 0; i < 6; ++i)
    {
    (*lb)[i] = (*lb)[6+i] = passiveLowerBounds[i];
    (*ub)[i] = (*ub)[6+i] = passiveUpperBounds[i];
    }
  (*lb)[12] = -100.0; // trying to make an order of magnitude or 2 greater than constraint values
  (*ub)[12] = 100.0; // trying to make an order of magnitude or 2 greater than constraint values
}

void FeasiblePlanProblem::getInitialGuessIK(std::vector<double>* x) const
{
  boost::mt19937 rng(2);
  boost::uniform_01<> uniDist;

  // Use Jacobian IK to find a nice initial guess
  // Place RCM's at some points on port curve
  double curveParamL = 0.5;
  double curveParamR = 0.5;
  Eigen::Vector3d rcmL = 
    this->portCurvePoint1 + curveParamL*(this->portCurvePoint2 - this->portCurvePoint1);
  Eigen::Vector3d rcmR = 
    this->portCurvePoint1 + curveParamR*(this->portCurvePoint2 - this->portCurvePoint1);
  std::vector<double> qL(6), qR(6);
  for (unsigned i = 0; i < 6; ++i)
    {
    qL[i] = passiveLowerBounds[i] + uniDist(rng)*(passiveUpperBounds[i] - passiveLowerBounds[i]);
    qR[i] = passiveLowerBounds[i] + uniDist(rng)*(passiveUpperBounds[i] - passiveLowerBounds[i]);
    }
  kin.passiveIK(this->baseFrameL, rcmL, &qL);
  std::cout << "Left IK done!" << std::endl; // debug
  std::cout << "left error: " 
            << (kin.passiveFK(baseFrameL, qL).topRightCorner<3,1>() - rcmL).norm()
            << std::endl;

  kin.passiveIK(this->baseFrameR, rcmR, &qR);
  std::cout << "Right IK done!" << std::endl; // debug
  std::cout << "right error: " 
            << (kin.passiveFK(baseFrameR, qR).topRightCorner<3,1>() - rcmR).norm()
            << std::endl;


  std::copy(qL.begin(), qL.end(), x->begin());
  std::copy(qR.begin(), qR.end(), x->begin()+6);
  (*x)[12] = 90.0;
}

void FeasiblePlanProblem::getInitialGuessExact(std::vector<double>* x) const
{
  (*x)[0] = 0.988914;
  (*x)[1] = -0.058769;
  (*x)[2] = -1.71274;
  (*x)[3] = -0.270708;
  (*x)[4] = -0.54294;
  (*x)[5] = 1.15854;
  (*x)[6] = 0.630289;
  (*x)[7] = 1.91681;
  (*x)[8] = 1.15106;
  (*x)[9] = -1.00725;
  (*x)[10] = 0.290322;
  (*x)[11] = 2.01327;
  (*x)[12] = -1.08085e-06;         
}

void FeasiblePlanProblem::outputStateProperties(std::ostream& out, 
                                                const std::vector<double>& x) const
{
  out << "x:";
  for (std::size_t i = 0; i < x.size(); ++i)
    out << " " << x[i];
  out << std::endl;

  std::vector<double> dummy;
  out << "f(x): " << FeasiblePlanProblem::feasibleMinimaxObj(x, dummy, (void*) this)  << std::endl;

  double* c_ineq = new double[this->getNumConstraints()];
  double x_array[this->getNumVariables()];
  for (unsigned i = 0; i < this->getNumVariables(); ++i)
    x_array[i] = x[i];
  FeasiblePlanProblem::wrapIneq(this->getNumConstraints(), c_ineq, this->getNumVariables(), x_array, 0, (void*) this);
  out << "c_ineq(x):";
  for (unsigned i = 0; i < this->getNumConstraints(); ++i)
    out << " " << (c_ineq[i] + x[12]);
  out << std::endl;

  std::vector<double> qpL(x.begin(), x.begin()+6);
  std::vector<double> qpR(x.begin()+6, x.begin()+12);
  std::vector<double> qaL(6);
  std::vector<double> qaR(6);
  this->kin.intraIK(this->kin.passiveFK(this->baseFrameL, qpL), this->taskFrames[0], &qaL);
  this->kin.intraIK(this->kin.passiveFK(this->baseFrameR, qpR), this->taskFrames[0], &qaR);
  out << "aL:";
  for (unsigned i = 0; i < 6; ++i)
    out << " " << qaL[i];
  out << std::endl << "aR:";
  for (unsigned i = 0; i < 6; ++i)
    out << " " << qaR[i];
  out << std::endl;

  delete c_ineq;
}

unsigned FeasiblePlanProblem::getNumVariables() const
{
  return 13;
}

unsigned FeasiblePlanProblem::getNumConstraints() const
{
  return 3 + (24+57)*taskFrames.size();
}
