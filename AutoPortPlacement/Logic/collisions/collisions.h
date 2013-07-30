#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <Eigen/Dense>
#include <vector>

namespace Collisions
{
  struct Cylisphere
  {
    Eigen::Vector3d p1;
    Eigen::Vector3d p2;
    double r;
  };

  struct Sphere
  {
    Eigen::Vector3d p;
    double r;
  };

  double distance(const Cylisphere& c, const Sphere& s);
  double distance(const Cylisphere& c1, const Cylisphere& c2);
  double distance(const Sphere& s1, const Sphere& s2);

  // Compute the minimum distance between one set of
  // cylispheres+sphere and another set
  double distance(const std::vector<Cylisphere>& c1,
                  const std::vector<Sphere>& s1,
                  const std::vector<Cylisphere>& c2,
                  const std::vector<Sphere>& s2);
}

#endif
