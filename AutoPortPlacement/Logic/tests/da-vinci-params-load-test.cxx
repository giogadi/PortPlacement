#include <davinci-kinematics/davinci.h>
#include <cmath>

int main(int argc, char** argv)
{
  DavinciKinematics kinematics;

  DavinciParameters p = kinematics.getParams();

  if (p.wristLength == 0.008 &&
      p.gripperLength == 0.01 &&
      p.el1 == 0.126 &&
      p.el2 == 0.04 &&
      p.el3 == 0.445 &&
      p.el4 == 0.043 &&
      p.el5 == 0.52 &&
      p.er1 == 0.005 &&
      p.er2 == 0.05 &&
      p.pl2 == 0.419 &&
      p.ph2 == 0.145 &&
      p.pl3 == 0.425 &&
      p.ph3 == 0.129 &&
      p.pl4 == 0.07 &&
      p.ph4 == 0.07 &&
      p.pl5 == 0.3 &&
      p.pRCMOffset(0) == 0.0 &&
      p.pRCMOffset(1) == 0.32 &&
      p.pRCMOffset(2) == 0.15 &&
      p.pLinkRadius == 0.07)
    {
    return 0;
    }
  else
    return 1;
}
