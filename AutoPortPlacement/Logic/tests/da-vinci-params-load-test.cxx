#include <davinci-kinematics/davinci.h>
#include <cmath>

int main(int argc, char** argv)
{
  DavinciKinematics kinematics(argv[1]);

  DavinciParameters p = kinematics.getParams();

  if (p.wristLength_ == 0.01 &&
      p.gripperLength_ == 0.01 &&
      p.eL_ == 0.15 &&
      p.eAlpha_ == asin(0.01 / 0.15) &&
      p.el1_ == 0.3 &&
      p.el2_ == 0.32 &&
      p.eRadius1_ == 0.05 &&
      p.eRadius2_ == 0.05 &&
      p.pl2_ == 0.3 &&
      p.ph2_ == 0.15 &&
      p.pl3_ == 0.3 &&
      p.ph3_ == 0.07 &&
      p.pl4_ == 0.07 &&
      p.ph4_ == 0.07 &&
      p.pl5_ == 0.3 &&
      p.pRCMOffset_(0) == 0.0 &&
      p.pRCMOffset_(1) == 0.32 &&
      p.pRCMOffset_(2) == 0.15 &&
      p.pLinkRadius_ == 0.07)
    {
    return 0;
    }
  else
    return 1;
}
