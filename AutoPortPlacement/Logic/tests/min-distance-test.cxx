#include <collisions/collisions.h>

#include <iostream>

int main(int, char* [])
{
  std::vector<Collisions::Cylisphere> c1, c2;
  std::vector<Collisions::Sphere> s1, s2;

  Collisions::Cylisphere c;
  Collisions::Sphere s;

  // cylisphere (0,0,0) - (1,0,0); r = 1.0
  c.p1 = Eigen::Vector3d::Zero();
  c.p2(0) = 1.0; c.p2(1) = 0.0; c.p2(2) = 0.0;
  c.r = 1.0;
  c1.push_back(c);
  
  // sphere (2,0,0); r = 1.0
  s.p(0) = 2.0; s.p(1) = 0.0; s.p(2) = 0.0;
  s.r = 1.0;
  s1.push_back(s);

  // cylisphere (0,3,0) - (1,3,0); r = 1.0
  c.p1(0) = 0.0; c.p1(1) = 3.0; c.p1(2) = 0.0;
  c.p2(0) = 1.0; c.p2(1) = 3.0; c.p2(2) = 0.0;
  c.r = 1.0;
  c2.push_back(c);

  // sphere (2,3,0); r = 1.0
  s.p(0) = 2.0; s.p(1) = 3.0; s.p(2) = 0.0;
  s.r = 1.0;
  s2.push_back(s);

  std::vector<double> dists;
  double dist = Collisions::distances(c1,s1,c2,s2, &dists);
  std::cout << dist << std::endl;
  std::cout << "num dists: " << dists.size() << std::endl;

  if (fabs(dist - 1.0) < 0.0000001)
    return 0;
  else
    return 1;
}
