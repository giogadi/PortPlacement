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

  // Compute distance between a segment and a point
  double distance(const Eigen::Vector3d& segmentPoint1, 
                  const Eigen::Vector3d& segmentPoint2,
                  const Eigen::Vector3d& point);

  double distance(const Cylisphere& c, const Sphere& s);
  double distance(const Cylisphere& c1, const Cylisphere& c2);
  double distance(const Sphere& s1, const Sphere& s2);

  // Compute the minimum distance between one set of
  // cylispheres+spheres and another set
  double distance(const std::vector<Cylisphere>& c1,
                  const std::vector<Sphere>& s1,
                  const std::vector<Cylisphere>& c2,
                  const std::vector<Sphere>& s2);

  // Compute all the distances between one set of cylispheres+spheres
  // to another set
  double distances(const std::vector<Cylisphere>& c1,
                   const std::vector<Sphere>& s1,
                   const std::vector<Cylisphere>& c2,
                   const std::vector<Sphere>& s2,
                   std::vector<double>* d);
}

#endif
