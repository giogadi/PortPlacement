#include <collisions/collisions.h>
#include <iostream>

// Segment-segment distance computation adapted from:
// http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
double Collisions::distance(const Cylisphere& c1, const Cylisphere& c2)
{
  // \TODO do this smarter. maybe use some multiple of
  // numeric_limits<double>::epsilon()
  const double SMALL_NUM = 0.00000001;

  Eigen::Vector3d u = c1.p2 - c1.p1;
  Eigen::Vector3d v = c2.p2 - c2.p1;
  Eigen::Vector3d w = c1.p1 - c2.p1;
  double a = u.dot(u);         // always >= 0
  double b = u.dot(v);
  double c = v.dot(v);         // always >= 0
  double d = u.dot(w);
  double e = v.dot(w);
  double D = a*c - b*b;        // always >= 0
  double sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
  double tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

  // compute the line parameters of the two closest points
  if (D < SMALL_NUM)  // the lines are almost parallel
    {
    sN = 0.0;         // force using point P0 on segment S1
    sD = 1.0;         // to prevent possible division by 0.0 later
    tN = e;
    tD = c;
    }
  else // get the closest points on the infinite lines
    {
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0.0) // sc < 0 => the s=0 edge is visible
      {
      sN = 0.0;
      tN = e;
      tD = c;
      }
    else if (sN > sD) // sc > 1  => the s=1 edge is visible
      {
      sN = sD;
      tN = e + b;
      tD = c;
      }
    }

  if (tN < 0.0) // tc < 0 => the t=0 edge is visible
    {
    tN = 0.0;
    // recompute sc for this edge
    if (-d < 0.0)
      sN = 0.0;
    else if (-d > a)
      sN = sD;
    else
      {
      sN = -d;
      sD = a;
      }
    }
  else if (tN > tD) // tc > 1  => the t=1 edge is visible
    {
    tN = tD;
    // recompute sc for this edge
    if ((-d + b) < 0.0)
      sN = 0;
    else if ((-d + b) > a)
      sN = sD;
    else
      {
      sN = (-d +  b);
      sD = a;
      }
    }

  // finally do the division to get sc and tc
  sc = (fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
  tc = (fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);

  // get the difference of the two closest points
  Eigen::Vector3d dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)

  return dP.norm() - (c1.r + c2.r);   // return the closest distance
}

// Adapted from point-to-segment distance calculation from:
// http://geomalgorithms.com/a02-_lines.html#Distance-to-Ray-or-Segment
double Collisions::distance(const Eigen::Vector3d& segmentPoint1, 
                            const Eigen::Vector3d& segmentPoint2,
                            const Eigen::Vector3d& point)
{
  Eigen::Vector3d v = segmentPoint2 - segmentPoint1;
  Eigen::Vector3d w = point - segmentPoint1;

  double c1 = w.dot(v);
  if ( c1 <= 0 )
    return (segmentPoint1 - point).norm();

  double c2 = v.dot(v);
  if ( c2 <= c1 )
    return (segmentPoint2 - point).norm();

  double b = c1 / c2;
  Eigen::Vector3d Pb = segmentPoint1 + b*v;
  return (Pb - point).norm();
}

double Collisions::distance(const Cylisphere& c, const Sphere& s)
{
  return distance(c.p1, c.p2, s.p) - (c.r + s.r);
}

double Collisions::distance(const Sphere& s1, const Sphere& s2)
{
  return (s2.p - s1.p).norm() - (s1.r + s2.r);
}

double Collisions::distance(const std::vector<Cylisphere>& c1,
                            const std::vector<Sphere>& s1,
                            const std::vector<Cylisphere>& c2,
                            const std::vector<Sphere>& s2)
{
  return distances(c1, s1, c2, s2, 0);
}

double Collisions::distances(const std::vector<Cylisphere>& c1,
                             const std::vector<Sphere>& s1,
                             const std::vector<Cylisphere>& c2,
                             const std::vector<Sphere>& s2,
                             std::vector<double>* d)
{
 double minDistance = std::numeric_limits<double>::infinity();
  
  // Compute distances of all cylispheres in c1 against all
  // cylispheres in c1 and spheres in c2
  for (std::vector<Cylisphere>::const_iterator c1iter = c1.begin();
       c1iter != c1.end();
       ++c1iter)
    {
    for (std::vector<Cylisphere>::const_iterator c2iter = c2.begin();
         c2iter != c2.end();
         ++c2iter)
      {
      double dist = distance(*c1iter, *c2iter);
      minDistance = std::min(minDistance, dist);
      if (d)
        d->push_back(dist);
      }
    for (std::vector<Sphere>::const_iterator s2iter = s2.begin();
         s2iter != s2.end();
         ++s2iter)
      {
      double dist = distance(*c1iter, *s2iter);
      minDistance = std::min(minDistance, dist);
      if (d)
        d->push_back(dist);
      }
    }

  for (std::vector<Sphere>::const_iterator s1iter = s1.begin();
       s1iter != s1.end();
       ++s1iter)
    {
    for (std::vector<Cylisphere>::const_iterator c2iter = c2.begin();
         c2iter != c2.end();
         ++c2iter)
      {
      double dist = distance(*c2iter, *s1iter);
      minDistance = std::min(minDistance, dist);
      if (d)
        d->push_back(dist);
      }
    for (std::vector<Sphere>::const_iterator s2iter = s2.begin();
         s2iter != s2.end();
         ++s2iter)
      {
      double dist = distance(*s1iter, *s2iter);
      minDistance = std::min(minDistance, dist);
      if (d)
        d->push_back(dist);
      }
    }

  return minDistance;
}
