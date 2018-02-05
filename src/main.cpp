#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"
#include "spline.h"

//int main() {

//}

//#define PYTHONPLOTTING
//#define QWTPLOTTING
const bool DEBUGTRANSFORMS=false;
const bool DEBUGPATHSPLINE=false;
const bool DEBUGVEHICLEPATH=false;
const bool DEBUGLANESELECTION=false;
const bool DEBUGSENSORREADINGS=false;
const bool DEBUGSACCELERATION=false;

using namespace std;
#ifdef PYTHONPLOTTING
namespace plt = matplotlibcpp;
#endif

#ifdef QWTPLOTTING
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_point_data.h>
#endif

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

string printVector(const vector<double>& theVector, const int precision=-1) {
  std::ostringstream oss;
  if (precision >= 0) {
    //cout << "precision:" << precision << std::endl;
    oss << std::setprecision(precision);
  }
  for (int i = 0; i < theVector.size(); i++) {
    oss << theVector[i];
    if (i <  theVector.size()-1)
      oss << ", ";
  }
  return oss.str();
}

string printVector(const vector<Eigen::Vector2d> theVector, const int precision=-1) {
  std::ostringstream oss;
  if (precision >= 0) {
    //cout << "precision:" << precision << std::endl;
    oss << std::setprecision(precision);
  }
  for (int i = 0; i < theVector.size(); i++) {
    oss << i << ":\n";
    Eigen::Vector2d vector = theVector[i];
    oss << vector.matrix();
    if (i <  theVector.size()-1)
      oss << "\n";
  }
  return oss.str();
}

string printVector(const vector<vector<double>> theVector, const int precision=-1) {
  std::ostringstream oss;
  if (precision >= 0) {
    //cout << "precision:" << precision << std::endl;
    oss << std::setprecision(precision);
  }
  for (int i = 0; i < theVector.size(); i++) {
    oss << i << ":\n";
    oss << printVector(theVector[i]);
    if (i <  theVector.size()-1)
      oss << "\n";
  }
  return oss.str();
}

string printVector(const vector<vector<vector<double>>> theVector, const int precision=-1) {
  std::ostringstream oss;
  if (precision >= 0) {
    //cout << "precision:" << precision << std::endl;
    oss << std::setprecision(precision);
  }
  for (int i = 0; i < theVector.size(); i++) {
    oss << i << ":\n";
    oss << printVector(theVector[i]);
    if (i <  theVector.size()-1)
      oss << "\n";
  }
  return oss.str();
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{
  
  double closestLen = 100000; //large number
  int closestWaypoint = 0;
  
  for(int i = 0; i < maps_x.size(); i++)
  {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen)
    {
      closestLen = dist;
      closestWaypoint = i;
    }
    
  }
  
  return closestWaypoint;
  
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);
  
  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];
  
  double heading = atan2((map_y-y),(map_x-x));
  
  double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);
  
  if(angle > pi()/4)
  {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size())
    {
      closestWaypoint = 0;
    }
  }
  
  return closestWaypoint;
}

vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);
  
  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size()-1;
  }
  
  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];
  
  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;
  
  double frenet_d = distance(x_x,x_y,proj_x,proj_y);
  
  //see if d value is positive or negative by comparing it to a center point
  
  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);
  
  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }
  
  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }
  
  frenet_s += distance(0,0,proj_x,proj_y);
  
  return {frenet_s,frenet_d};
  
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  //cout << "maps_s:" << printVector(maps_s);
  int prev_wp = -1;
  
  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }
  
  int wp2 = (prev_wp+1)%maps_x.size();
  
  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);
  
  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);
  
  double perp_heading = heading-pi()/2;
  
  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);
  
  return {x,y};
  
}

const int LANESIZE=4;//meters
const int RIGHTLANE=0;
const int MIDDLELANE=1;
const int LEFTLANE=2;

int determineLane(double theCarD) {
  /*
   The d vector has a magnitude of 1 and points perpendicular to the road in the direction of the right-hand side of the road. The d
   vector can be used to calculate lane positions. For example, if you want to be in the left lane at some waypoint just add the
   waypoint's (x,y) coordinates with the d vector multiplied by 2. Since the lane is 4 m wide, the middle of the left lane (the lane
   closest to the double-yellow diving line) is 2 m from the waypoint.
   */
  return int(theCarD)/LANESIZE;
}

double laneCenter(int theLane) {
  /*
   If you would like to be in the middle lane, add the waypoint's coordinates to the d vector multiplied by 6 = (2+4), since the
   center of the middle lane is 4 m from the center of the left lane, which is itself 2 m from the double-yellow diving line and the
   waypoints.
   */
  return double(theLane*LANESIZE)+(LANESIZE/2.);
}

const double HOURSPERMINUTE = 1./60.;
const double MINUTESPERSECOND = 1./60.;
const double HOURSPERSECOND = HOURSPERMINUTE*MINUTESPERSECOND;

const double FEETPERMILE= 5280.;
const double METERSPERFOOT = .3048;
const double MILESPERHOURTOMETERSPERSECOND = HOURSPERSECOND*FEETPERMILE*METERSPERFOOT;

const static double milesPerHourToMetersPerSecond(const double theMilesPerHour) {
  return theMilesPerHour*MILESPERHOURTOMETERSPERSECOND;
}

const double MAXIMUMMILESPERHOUR=50.;// miles per hour
const int CARMOVESPERSECOND=50;// car moves at a well known fixed rate of 50 times per second
const int LENGTHOFCARPATHTRAJECTORY=CARMOVESPERSECOND*1/*second*/;
const double PATHLOOKAHEADXDISTANCE=40.;// meters
const double PATHLOOKAHEADSDISTANCE=40.;// meters

const double MAXIMUMMETERSPERSECOND=milesPerHourToMetersPerSecond(MAXIMUMMILESPERHOUR);
// In order for the passenger to have an enjoyable ride both the jerk and the total acceleration should not exceed 10 m/s^2.
const double MAXIMUMACCERLERATION=10.;// m/s^2

const bool isInitialPath(const vector<double> thePath) {
  return thePath.size()==0;
}

// 2 vectors holding 2 x points & 2 y points — for a path size >= 2
const vector<vector<double>>& last2PathPoints(const vector<double> thePathX, const vector<double> thePathY) {
  vector<double> xPoints;
  vector<double> yPoints;
  const int pathSize=thePathX.size();
  if (pathSize >= 2) {
    // previous first
    const double secondX=thePathX[pathSize-2];
    const double secondY=thePathY[pathSize-2];
    xPoints.push_back(secondX);
    yPoints.push_back(secondY);
    // then the current
    const double firstX=thePathX[pathSize-1];
    const double firstY=thePathY[pathSize-1];
    xPoints.push_back(firstX);
    yPoints.push_back(firstY);
  }
  vector<vector<double>>* xyPoints= new vector<vector<double>>;
  (*xyPoints).push_back(xPoints);
  (*xyPoints).push_back(yPoints);
  return *xyPoints;
}

// 2 vectors holding 2 x points & 2 y points  — for the initial path
const vector<vector<double>>& last2PathPoints(const double theCurrentX, const double theCurrentY, const double theCurrentYaw) {
  vector<double> xPoints;
  vector<double> yPoints;
  // previous first
  // interpolate a prior xy point to the current one based on yaw
  xPoints.push_back(theCurrentX-cos(theCurrentYaw));
  yPoints.push_back(theCurrentY-sin(theCurrentYaw));
  // then the current
  xPoints.push_back(theCurrentX);
  yPoints.push_back(theCurrentY);
  
  vector<vector<double>>* xyPoints= new vector<vector<double>>;
  (*xyPoints).push_back(xPoints);
  (*xyPoints).push_back(yPoints);
  return *xyPoints;
}

const double yawFromXY(const vector<vector<double>>& theXYs) {
  vector<double> xPoints=theXYs[0];
  vector<double> yPoints=theXYs[1];
  double yaw = atan2((yPoints[1]-yPoints[0]), (xPoints[1]-xPoints[0]));
  return yaw;
}

const vector<double> sFromXY(vector<vector<double>>& theXYs, const double theTheta, const vector<double> &theMapsX, const vector<double> &theMapsY) {
  const vector<double> xPoints=theXYs[0];
  const vector<double> yPoints=theXYs[1];
  vector<double> sValues;
  for (int xy=0; xy<theXYs.size(); xy++) {
    // 0:s & 1:d
    const vector<double> sdValue = getFrenet(xPoints[xy], yPoints[xy], theTheta, theMapsX, theMapsY);
    sValues.push_back(sdValue[0]/* s value*/);
  }
  cout << "sFromXY-theXYs:" << printVector(theXYs) << ", sValues:" << printVector(sValues) << std::endl;
  return sValues;
}

const vector<double> velocityFromXY(const double theDeltaT, vector<vector<double>> theXYs, const double theTheta, const vector<double> &theMapsX, const vector<double> &theMapsY) {
  //const vector<double> sValues=sFromXY(theXYs, theTheta, theMapsX, theMapsY);
  const vector<double> xPoints=theXYs[0];
  const vector<double> yPoints=theXYs[1];
  vector<double> vValues;
  double x0 = xPoints[0];
  double y0 = yPoints[0];
  for (int s=1; s<theXYs.size(); s++) {
    const double x = xPoints[s];
    const double y = yPoints[s];
    const double deltaX = x-x0;
    const double deltaY = y-y0;
    const double deltaD = sqrt(deltaX*deltaX+deltaY*deltaY);
    const double velocity=deltaD/theDeltaT;
    if (DEBUGVEHICLEPATH) {
      cout << "velocityFromXY-velocity:" << velocity << ", deltaD:" << deltaD << ", deltaX:" << deltaX << ", deltaY:" << deltaY << ", theXYs[" << s << "]:" << printVector(theXYs[s]) << std::endl;
    }
    vValues.push_back(velocity);
    x0 = xPoints[s];
    y0 = yPoints[s];
  }
  if (DEBUGVEHICLEPATH) {
    cout << "velocityFromXY-theXYs:" << printVector(theXYs) << ", vValues:" << printVector(vValues) << std::endl;
  }
  return vValues;
}

const bool static areSame(const double a, const double b, const double epsilon) {
  return fabs(a - b) < epsilon;
}

const bool static isZero(const double a, const double epsilon) {
  return fabs(a) < epsilon;
}

const bool static isNotZero(const double a, const double epsilon) {
  return !isZero(a, epsilon);
}

const bool static areSame(const  Eigen::VectorXd a, const  Eigen::VectorXd b, const double epsilon) {
  for (int r=0; r<a.rows(); r++) {
    if (!areSame(a(r),b(r), epsilon)) {
      std::cout << std::endl
      << "a(" << r << "):" << a(r) << " != "
      << "b(" << r << "):" << b(r) << " == "
      << "fabs(a-b):" << fabs(a(r) - b(r))
      << std::endl;
      return false;
    }
  }
  return true;
}

const bool static areSame(const  Eigen::MatrixXd a, const  Eigen::MatrixXd b, const double epsilon) {
  for (int r=0; r<a.rows(); r++) {
    for (int c=0; c<a.cols(); c++) {
      if (!areSame(a(r,c),b(r,c), epsilon)) {
        std::cout << std::endl
        << "a(" << r << "," << c << "):" << a(r,c) << " != "
        << "b(" << r << "," << c << "):" << b(r,c) << " == "
        << "fabs(a-b):" << fabs(a(r,c) - b(r,c))
        << std::endl;
        return false;
      }
    }
  }
  return true;
}

const bool static areSame(const Eigen::Matrix3d a, const  Eigen::Affine2d b, const double epsilon) {
  return areSame(a, (Eigen::MatrixXd) b.matrix(), epsilon);
}

static const double EPSILON=1.e-5;

const bool static areSame(const Eigen::Matrix3d a, const  Eigen::Affine2d b) {
  return areSame(a,b,EPSILON);
}

const bool static areSame(const  Eigen::VectorXd a, const  Eigen::VectorXd b) {
  return areSame(a, b, EPSILON);
}

const bool static areSame( Eigen::MatrixXd a,  Eigen::MatrixXd b) {
  return areSame(a, b, EPSILON);
}

const bool static areSame(double a, double b) {
  return areSame(a, b, EPSILON);
}

const Eigen::Affine2d static rtAffineTransform(const double theX, const double theY, const double theTheta) {
  // translation followed by rotation == RT
  if (DEBUGTRANSFORMS)
    cout << "xyAffineTransform-theX:" << std::setprecision(10) << theX << ", theY:" << theY << ", theTheta:" << theTheta << std::endl;
  const Eigen::Affine2d translateToOrigin(Eigen::Translation<double,2>(theX,theY));
  const Eigen::Affine2d rotation((Eigen::Rotation2D<double>(theTheta)));
  if (DEBUGTRANSFORMS) {
    cout  << "xyAffineTransform-translateToOrigin:" << std::setprecision(10) << std::endl << translateToOrigin.matrix() << " ," << std::endl
    << "rotation:" << std::setprecision(10) <<std::endl << rotation.matrix() << std::endl << " ," << std::endl
    << "rotation*translateToOrigin:" << std::endl << (rotation*translateToOrigin).matrix()
    << std::endl;
  }
  return rotation*translateToOrigin;
}

const Eigen::Affine2d static trAffineTransform(const double theX, const double theY, const double theTheta) {
  // rotation followed by translation
  if (DEBUGTRANSFORMS)
    cout << "xyAffineTransform-theX:" << std::setprecision(10) << theX << ", theY:" << theY << ", theTheta:" << theTheta << std::endl;
  const Eigen::Affine2d translateToOrigin(Eigen::Translation<double,2>(theX,theY));
  const Eigen::Affine2d rotation((Eigen::Rotation2D<double>(theTheta)));
  if (DEBUGTRANSFORMS) {
    cout  << "xyAffineTransform-translateToOrigin:" << std::setprecision(10) << std::endl << translateToOrigin.matrix() << " ," << std::endl
    << "rotation:" << std::setprecision(10) <<std::endl << rotation.matrix() << std::endl << " ," << std::endl
    << "rotation*translateToOrigin:" << std::endl << (rotation*translateToOrigin).matrix()
    << std::endl;
  }
  return translateToOrigin*rotation;
}

const Eigen::Affine2d static rtInverseTransform(const double theX, const double theY, const double theTheta) {
  // rotation followed by translation == (RT)-1
  // (RT)-1 == (T-1)(R-1) & T-1(x,y)=T(-x, -y), R-1(theta)= R(-theta) => (RT)-1=T(-x,-y)R(-theta)
  if (DEBUGTRANSFORMS)
    cout << "rtInverseTransform-theX:" << std::setprecision(10) << theX << ", theY:" << theY << ", theTheta:" << theTheta << std::endl;
  return trAffineTransform(-theX, -theY, -theTheta);
}

const Eigen::Affine2d static trInverseTransform(const double theX, const double theY, const double theTheta) {
  // translation followed by rotation
  // (TR)-1 == (R-1)(T-1) & T-1(x,y)=T(-x, -y), R-1(theta)= R(-theta) => (TR)-1=R(-theta)T(-x,-y)
  if (DEBUGTRANSFORMS)
    cout << "xyAffineTransform-theX:" << std::setprecision(10) << theX << ", theY:" << theY << ", theTheta:" << theTheta << std::endl;
  return rtAffineTransform(-theX, -theY, -theTheta);
}

const Eigen::Matrix3d static rtMatrixTransform(const double theX, const double theY, const double theTheta) {
  if (DEBUGTRANSFORMS) cout << std::setprecision(10) << "xyMatrixTransform-theX:" << theX << ", theY:" << theY << ", theTheta:" << theTheta << std::endl;
  const double sinMapRotation=sin(theTheta);
  const double cosMapRotation=cos(theTheta);
  const double sinRotationXtranslation=sinMapRotation*(theX);
  const double sinRotationYtranslation=sinMapRotation*(theY);
  const double cosRotationXtranslation=cosMapRotation*(theX);
  const double cosRotationYtranslation=cosMapRotation*(theY);
  Eigen::Matrix3d rotationAndTranslation;
  rotationAndTranslation <<
  cosMapRotation, -sinMapRotation,  +cosRotationXtranslation-sinRotationYtranslation,
  sinMapRotation, cosMapRotation,   +sinRotationXtranslation+cosRotationYtranslation,
  0.,             0.,               1;
  if (DEBUGTRANSFORMS) cout << std::setprecision(10) << "xyMatrixTransform-translationFollowedByRotation:" << std::endl << rotationAndTranslation.matrix() << std::endl;
  return rotationAndTranslation;
}

const Eigen::Affine2d static rtAffineTransform(const Eigen::Vector2d theXY, const double theTheta) {
  return rtAffineTransform(theXY[0], theXY[1], theTheta);
}

const Eigen::Affine2d static rtInverseTransform(const Eigen::Vector2d theXY, const double theTheta) {
  return rtInverseTransform(theXY[0], theXY[1], theTheta);
}

const Eigen::Affine2d static trAffineTransform(const Eigen::Vector2d theXY, const double theTheta) {
  return trAffineTransform(theXY[0], theXY[1], theTheta);
}

const Eigen::Affine2d static trInverseTransform(const Eigen::Vector2d theXY, const double theTheta) {
  return trInverseTransform(theXY[0], theXY[1], theTheta);
}

const Eigen::Matrix3d static rtMatrixTransformation(const Eigen::Vector2d theXY, const double theTheta) {
  return rtMatrixTransform(theXY[0], theXY[1], theTheta);
}

const Eigen::Vector3d static xyAffineTransform(const Eigen::Affine2d theAffineTransform, const Eigen::Vector3d theHomogeneousXY) {
  if (DEBUGTRANSFORMS)
    cout << "xyAffineTransform-theAffineTransform:\n" << std::setprecision(10) << theAffineTransform.matrix()
    << "\n, theHomogeneousXY:\n" << theHomogeneousXY.matrix() << std::endl;
  return theAffineTransform*theHomogeneousXY;
}

const Eigen::Vector3d static xyMatrixTransform(const Eigen::Matrix3d theMatrixTransform, const Eigen::Vector3d theHomogeneousXY) {
  return theMatrixTransform*theHomogeneousXY;
}

const Eigen::Vector2d static xyAffineTransform(const Eigen::Affine2d theAffineTransform, const Eigen::Vector2d theXY) {
  const Eigen::Vector3d homogeneousPosition(theXY[0], theXY[1], 1.);
  const Eigen::Vector3d transformedPosition=xyAffineTransform(theAffineTransform,homogeneousPosition);
  const Eigen::Vector2d xyPosition(transformedPosition[0], transformedPosition[1]);
  return xyPosition;
}

const void testXY() {
  const bool DEBUGTRANSFORMS=false;
  const double PI=M_PI;
  Eigen::Vector2d xy = Eigen::Vector2d(1., -1.);
  
  Eigen::Affine2d rtAffine=rtAffineTransform(xy, PI);
  
  Eigen::Vector3d xyTransform = xyAffineTransform(rtAffine, Eigen::Vector3d(-xy[0], -xy[1], 1.));
  if (DEBUGTRANSFORMS)
    cout << "testXY-rtAffine:\n" << std::setprecision(10) << rtAffine.matrix()
    << "\n, xyTransform:\n" << xyTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xyTransform, Eigen::Vector3d(0., 0., 1.)));
  
  Eigen::Affine2d rtInverse=rtInverseTransform(xy, PI);
  xyTransform = xyAffineTransform(rtInverse, Eigen::Vector3d(-xy[0], -xy[1], 1.));
  if (DEBUGTRANSFORMS)
    cout << "testXY-rtInverse:\n" << std::setprecision(10) << rtInverse.matrix()
    << "\n, xyTransform:\n" << xyTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xyTransform, Eigen::Vector3d(0., 0., 1.)));
  
  Eigen::Affine2d trAffine=trAffineTransform(xy, PI);
  
  xyTransform = xyAffineTransform(trAffine, Eigen::Vector3d(-xy[1], -xy[0], 1.));
  if (DEBUGTRANSFORMS)
    cout << "testXY-trAffine:\n" << std::setprecision(10) << trAffine.matrix()
    << "\n, xyTransform:\n" << xyTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xyTransform, Eigen::Vector3d(0., 0., 1.)));
  
  Eigen::Affine2d trInverse=rtInverseTransform(xy, PI);
  xyTransform = xyAffineTransform(trInverse, Eigen::Vector3d(-xy[0], -xy[1], 1.));
  if (DEBUGTRANSFORMS)
    cout << "testXY-trInverse:\n" << std::setprecision(10) << trInverse.matrix()
    << "\n, xyTransform:\n" << xyTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xyTransform, Eigen::Vector3d(0., 0., 1.)));
  
  Eigen::Vector2d xy2dTransform = xyAffineTransform(rtAffine, Eigen::Vector2d(-xy[0], -xy[1]));
  if (DEBUGTRANSFORMS)
    cout << "testXY-rtAffine:\n" << std::setprecision(10) << rtAffine.matrix()
    << "\n, xy2dTransform:\n" << xy2dTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xy2dTransform, Eigen::Vector2d(0., 0.)));
  
  xy2dTransform = xyAffineTransform(rtInverse, Eigen::Vector2d(-xy[0], -xy[1]));
  if (DEBUGTRANSFORMS)
    cout << "testXY-rtInverse:\n" << std::setprecision(10) << rtInverse.matrix()
    << "\n, xy2dTransform:\n" << xy2dTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xy2dTransform, Eigen::Vector2d(0., 0.)));
  
  xy2dTransform = xyAffineTransform(trAffine, Eigen::Vector2d(-xy[1], -xy[0]));
  if (DEBUGTRANSFORMS)
    cout << "testXY-trAffine:\n" << std::setprecision(10) << trAffine.matrix()
    << "\n, xy2dTransform:\n" << xy2dTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xy2dTransform, Eigen::Vector2d(0., 0.)));
  
  xy2dTransform = xyAffineTransform(trInverse, Eigen::Vector2d(-xy[0], -xy[1]));
  if (DEBUGTRANSFORMS)
    cout << "testXY-trInverse:\n" << std::setprecision(10) << trInverse.matrix()
    << "\n, xy2dTransform:\n" << xy2dTransform.matrix() << std::endl;
  assert(areSame((Eigen::VectorXd)xy2dTransform, Eigen::Vector2d(0., 0.)));
  
}

const Eigen::Vector2d static xyMatrixTransform(const Eigen::Matrix3d theMatrixTransform, const Eigen::Vector2d theXY) {
  const Eigen::Vector3d homogeneousPosition(theXY[0], theXY[1], 1.);
  const Eigen::Vector3d transformedPosition=xyMatrixTransform(theMatrixTransform,homogeneousPosition);
  const Eigen::Vector2d xyPosition(transformedPosition[0], transformedPosition[1]);
  return xyPosition;
}

const vector<Eigen::Vector2d> static xyAffineTransform(const Eigen::Affine2d theAffineTransform, const vector<Eigen::Vector2d> theXYs) {
  vector<Eigen::Vector2d> transformedXYs;
  for (int p=0; p<theXYs.size(); p++) {
    Eigen::Vector2d transformedXY = xyAffineTransform(theAffineTransform, theXYs[p]);
    if (DEBUGTRANSFORMS) cout << "xyAffineTransform-theMapPositions[" << p << "]:" << theXYs[p].matrix() << std::setprecision(10) << ", transformedXY:" << transformedXY.matrix() << std::endl;
    transformedXYs.push_back(transformedXY);
  }
  return transformedXYs;
}

const vector<Eigen::Vector2d> static rtAffineTransform(const double theX, const double theY, const double theTheta, const vector<double> theXs, const vector<double> theYs) {
  assert(theXs.size() == theYs.size());
  vector<Eigen::Vector2d> xy;
  for (int p=0; p<theXs.size(); p++) {
    xy.push_back(Eigen::Vector2d(theXs[p], theYs[p]));
  }
  const Eigen::Affine2d rtAffineTransformation=rtAffineTransform(theX, theY, theTheta);
  const Eigen::Matrix3d rtMatrixTransformation=rtMatrixTransform(theX, theY, theTheta);
  assert(areSame(rtMatrixTransformation, rtAffineTransformation));
  return xyAffineTransform(rtAffineTransformation, xy);
}

const vector<Eigen::Vector2d> static rtInverseTransform(const double theX, const double theY, const double theTheta, const vector<double> theXs, const vector<double> theYs) {
  assert(theXs.size() == theYs.size());
  vector<Eigen::Vector2d> xy;
  for (int p=0; p<theXs.size(); p++) {
    xy.push_back(Eigen::Vector2d(theXs[p], theYs[p]));
  }
  const Eigen::Affine2d rtInverseTransformation=rtInverseTransform(theX, theY, theTheta);
  //const Eigen::Matrix3d rtMatrixTransformation=rtMatrixTransform(theX, theY, theTheta);
  //assert(areSame(rtMatrixTransformation, rtAffineTransformation));
  return xyAffineTransform(rtInverseTransformation, xy);
}

const static vector<vector<double>> unpackV2D(const vector<Eigen::Vector2d> theXYEigenVector) {
  vector<vector<double>> vectorCopy;
  vector<double> xValues;
  vector<double> yValues;
  for (int v=0; v<theXYEigenVector.size(); v++) {
    const Eigen::Vector2d eigenVector=theXYEigenVector[v];
    xValues.push_back(eigenVector[0]);
    yValues.push_back(eigenVector[1]);
  }
  vectorCopy.push_back(xValues);
  vectorCopy.push_back(yValues);
  return vectorCopy;
}

const void testRT() {
  const bool DEBUGTRANSFORMS=false;
  const double PI=M_PI;
  Eigen::Affine2d rtAffine=rtAffineTransform(1., -1., PI);
  Eigen::Affine2d rtInverse=rtInverseTransform(1., -1., PI);
  Eigen::Affine2d rtAndRtInverse=rtAffine*rtInverse;
  if (DEBUGTRANSFORMS)
    cout << "rtAndRtInverse:\n" << rtAndRtInverse.matrix() << ", rows:" << rtAndRtInverse.matrix().rows()
    << ", columns:" << rtAndRtInverse.matrix().cols() << std::endl;
  Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(rtAndRtInverse.matrix().rows(),rtAndRtInverse.matrix().cols());
  if (DEBUGTRANSFORMS)
    cout << "testRT-rtAffine:\n" << std::setprecision(10) << rtAffine.matrix() << "\n, rtInverse:\n" << rtInverse.matrix()
    << ", rtAndRtInverse:\n" << rtAndRtInverse.matrix() << "\n, identity:\n" << identity.matrix() << std::endl;
  assert(areSame(identity, rtAndRtInverse));
  
  
  rtAffine=rtAffineTransform(Eigen::Vector2d(1., -1.), PI);
  rtInverse=rtInverseTransform(Eigen::Vector2d(1., -1.), PI);
  rtAndRtInverse=rtAffine*rtInverse;
  if (DEBUGTRANSFORMS)
    cout << "rtAndRtInverse:\n" << rtAndRtInverse.matrix() << ", rows:" << rtAndRtInverse.matrix().rows()
    << ", columns:" << rtAndRtInverse.matrix().cols() << std::endl;
  identity=Eigen::MatrixXd::Identity(rtAndRtInverse.matrix().rows(),rtAndRtInverse.matrix().cols());
  if (DEBUGTRANSFORMS)
    cout << "testRT-rtAffine:\n" << std::setprecision(10) << rtAffine.matrix() << "\n, rtInverse:\n" << rtInverse.matrix()
    << ", rtAndRtInverse:\n" << rtAndRtInverse.matrix() << "\n, identity:\n" << identity.matrix() << std::endl;
  assert(areSame(identity, rtAndRtInverse));
  
  vector<double> theXs;
  theXs.push_back(1.);
  theXs.push_back(-1.);
  vector<double> theYs;
  theYs.push_back(-1.);
  theYs.push_back(1.);
  
  vector<Eigen::Vector2d> rtAffineXYs = rtAffineTransform(1., -1., PI, theXs, theYs);
  const vector<vector<double>> xyVector = unpackV2D(rtAffineXYs);
  const vector<double> xS = xyVector[0];
  const vector<double> yS = xyVector[1];
  vector<Eigen::Vector2d> rtInverseXYs = rtInverseTransform(1., -1., PI, xS, yS);
  if (DEBUGTRANSFORMS)
    cout << "testRT-rtAffineXYs:\n" << std::setprecision(10) << printVector(rtAffineXYs)
    << "\rtInverseXYs:\n" << printVector(rtInverseXYs) << std::endl;
  for (int xy=0; xy<rtAffineXYs.size(); xy++) {
    const Eigen::Vector2d inverseXY=rtInverseXYs[xy];
    assert(areSame(theXs[xy], inverseXY[0]));
    assert(areSame(theYs[xy], inverseXY[1]));
  }
}

const vector<Eigen::Vector2d> static trAffineTransform(const double theX, const double theY, const double theTheta, const vector<double> theXs, const vector<double> theYs) {
  assert(theXs.size() == theYs.size());
  vector<Eigen::Vector2d> xy;
  for (int p=0; p<theXs.size(); p++) {
    xy.push_back(Eigen::Vector2d(theXs[p], theYs[p]));
  }
  const Eigen::Affine2d trAffineTransformation=trAffineTransform(theX, theY, theTheta);
  return xyAffineTransform(trAffineTransformation, xy);
}

const vector<Eigen::Vector2d> static trInverseTransform(const double theX, const double theY, const double theTheta, const vector<double> theXs, const vector<double> theYs) {
  assert(theXs.size() == theYs.size());
  vector<Eigen::Vector2d> xy;
  for (int p=0; p<theXs.size(); p++) {
    xy.push_back(Eigen::Vector2d(theXs[p], theYs[p]));
  }
  const Eigen::Affine2d trInverse=trInverseTransform(theX, theY, theTheta);
  return xyAffineTransform(trInverse, xy);
}


const void testTR() {
  const bool DEBUGTRANSFORMS=false;
  const double PI=M_PI;
  Eigen::Affine2d trAffine=trAffineTransform(1., -1., PI);
  Eigen::Affine2d trInverse=trInverseTransform(1., -1., PI);
  Eigen::Affine2d trAndTrInverse=trAffine*trInverse;
  if (DEBUGTRANSFORMS)
    cout << "trAndTrInverse:\n" << trAndTrInverse.matrix() << ", rows:" << trAndTrInverse.matrix().rows()
    << ", columns:" << trAndTrInverse.matrix().cols() << std::endl;
  Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(trAndTrInverse.matrix().rows(),trAndTrInverse.matrix().cols());
  if (DEBUGTRANSFORMS)
    cout << "testTR-rtAffine:\n" << std::setprecision(10) << trAffine.matrix() << "\n, trInverse:\n" << trInverse.matrix()
    << ", trAndTrInverse:\n" << trAndTrInverse.matrix() << "\n, identity:\n" << identity.matrix() << std::endl;
  assert(areSame(identity, trAndTrInverse));
  
  trAffine=trAffineTransform(Eigen::Vector2d(1., -1.), PI);
  trInverse=trInverseTransform(Eigen::Vector2d(1., -1.), PI);
  trAndTrInverse=trAffine*trInverse;
  if (DEBUGTRANSFORMS)
    cout << "trAndTrInverse:\n" << trAndTrInverse.matrix() << ", rows:" << trAndTrInverse.matrix().rows()
    << ", columns:" << trAndTrInverse.matrix().cols() << std::endl;
  identity=Eigen::MatrixXd::Identity(trAndTrInverse.matrix().rows(),trAndTrInverse.matrix().cols());
  if (DEBUGTRANSFORMS)
    cout << "testTR-trAffine:\n" << std::setprecision(10) << trAffine.matrix() << "\n, trInverse:\n" << trInverse.matrix()
    << ", trAndTrInverse:\n" << trAndTrInverse.matrix() << "\n, identity:\n" << identity.matrix() << std::endl;
  assert(areSame(identity, trAndTrInverse));
  
  vector<double> theXs;
  theXs.push_back(1.);
  theXs.push_back(-1.);
  vector<double> theYs;
  theYs.push_back(-1.);
  theYs.push_back(1.);
  
  vector<Eigen::Vector2d> trAffineXYs = trAffineTransform(1., -1., PI, theXs, theYs);
  const vector<vector<double>> xyVector = unpackV2D(trAffineXYs);
  const vector<double> xS = xyVector[0];
  const vector<double> yS = xyVector[1];
  vector<Eigen::Vector2d> trInverseXYs = trInverseTransform(1., -1., PI, xS, yS);
  if (DEBUGTRANSFORMS)
    cout << "testTR-trAffineXYs:\n" << std::setprecision(10) << printVector(trAffineXYs)
    << "\ntrInverseXYs:\n" << printVector(trInverseXYs) << std::endl;
  for (int xy=0; xy<trAffineXYs.size(); xy++) {
    const Eigen::Vector2d inverseXY=trInverseXYs[xy];
    assert(areSame(theXs[xy], inverseXY[0]));
    assert(areSame(theYs[xy], inverseXY[1]));
  }
  
}

/*
const static Eigen::Vector2d trAffineTransformXY(const double theX, const double theY, const double theTheta) {
  Eigen::Vector2d *xy = new Eigen::Vector2d;
  (*xy) << theX, theY;
  const Eigen::Affine2d trAffineTransformation=rtAffineTransform(theX, theY, theTheta);
  return xyAffineTransform(trAffineTransformation, *xy);
}
*/
const vector<Eigen::Vector2d> static transformFromMapToCar(const double theCarX, const double theCarY, const double theCarTheta, const vector<double> theMapXs, const vector<double> theMapYs) {
  return trInverseTransform(theCarX, theCarY, theCarTheta, theMapXs, theMapYs);
}

const static vector<Eigen::Vector2d> transformFromCarToMap(const double theMapX, const double theMapY, const double theMapTheta, const vector<double> theCarXs, const vector<double> theCarYs) {
  return trAffineTransform(theMapX, theMapY, theMapTheta, theCarXs, theCarYs);
}

const void testCarMap() {
  const bool DEBUGTRANSFORMS=false;
  const double PI=M_PI;
  
  vector<double> theCarXs;
  theCarXs.push_back(1.);
  theCarXs.push_back(-1.);
  vector<double> theCarYs;
  theCarYs.push_back(-1.);
  theCarYs.push_back(1.);
  
  vector<Eigen::Vector2d> mapXYs = transformFromCarToMap(1., -1., PI, theCarXs, theCarYs);
  const vector<vector<double>> xyVector = unpackV2D(mapXYs);
  const vector<double> mapXs = xyVector[0];
  const vector<double> mapYs = xyVector[1];
  vector<Eigen::Vector2d> carXYs = transformFromMapToCar(1., -1., PI, mapXs, mapYs);
  if (DEBUGTRANSFORMS)
    cout << "testCarMap-mapXYs:\n" << std::setprecision(10) << printVector(mapXYs)
    << "\ncarXYs:\n" << printVector(carXYs) << std::endl;
  for (int xy=0; xy<theCarXs.size(); xy++) {
    const Eigen::Vector2d carXY=carXYs[xy];
    assert(areSame(theCarXs[xy], carXY[0]));
    assert(areSame(theCarYs[xy], carXY[1]));
  }
  
}
const static vector<vector<double>>& rtClassTransform(const double refX, const double refY, const double refYaw, const vector<double>& ptsX, const vector<double>& ptsY) {
  // translation folloed by a rotation => RT
  if (DEBUGTRANSFORMS)
    cout << std::setprecision(10) << "trClassTransform-refX:" << refX << ", refY:" << refY
    << ", refYaw:" << refYaw << ", cos(-refYaw):" << cos(0-refYaw) << ", sin(-refYaw):" << sin(0-refYaw) <<  std::endl;
  vector<double> xPoints;
  vector<double> yPoints;
  
  for (int i=0; i<ptsX.size(); i++) {
    const double shiftX=ptsX[i]-refX;
    const double shiftY=ptsY[i]-refY;
    
    xPoints.push_back(shiftX*cos(0-refYaw)-shiftY*sin(0-refYaw));
    yPoints.push_back(shiftX*sin(0-refYaw)+shiftY*cos(0-refYaw));
    if (DEBUGTRANSFORMS)
      cout << "trClassTransform-i:" << i << std::setprecision(10) << ", shiftX:" << shiftX << ", shiftY:" << shiftY << std::endl;
  }
  
  vector<vector<double>>* xyPoints= new vector<vector<double>>;
  (*xyPoints).push_back(xPoints);
  (*xyPoints).push_back(yPoints);
  return *xyPoints;
}

const static vector<double>& trClassTransform(const double refX, const double refY, const double refYaw, const double theX, const double theY) {
  // rotation followed by translation => TR
  if (DEBUGTRANSFORMS)
    cout << std::setprecision(10) << "rtClassTransform-refX:" << refX << ", refY:" << refY
    << ", refYaw:" << refYaw << ", cos(refYaw):" << cos(refYaw) << ", sin(refYaw):" << sin(refYaw) <<  std::endl;
  
  const double xPoint=theX*cos(refYaw)-theY*sin(refYaw);
  const double yPoint=theX*sin(refYaw)+theY*cos(refYaw);
  
  vector<double> *xy = new vector<double>;
  (*xy).push_back(xPoint+refX);
  (*xy).push_back(yPoint+refY);
  
  if (DEBUGTRANSFORMS)
    cout << "rtClassTransform:" << std::setprecision(10) << ", xPoint:" << xPoint << ", yPoint:" << yPoint << ", xy:" <<  printVector(*xy, 10) << std::endl;
  
  return *xy;
}

const static vector<vector<double>>& trClassTransform(const double refX, const double refY, const double refYaw, const vector<double>& ptsX, const vector<double>& ptsY) {
  // rotation followed by translation => TR
  vector<double> *xPoints=new vector<double>;
  vector<double> *yPoints=new vector<double>;
  for (int i=0; i<ptsX.size(); i++) {
    vector<double> rtXY=trClassTransform(refX, refY, refYaw, ptsX[i], ptsY[i]);
    (*xPoints).push_back(rtXY[0]);
    (*yPoints).push_back(rtXY[1]);
  }
  
  vector<vector<double>> *xyPoints=new vector<vector<double>>;
  (*xyPoints).push_back(*xPoints);
  (*xyPoints).push_back(*yPoints);
  return *xyPoints;
}

// The data format for each car is: [ id, x, y, vx, vy, s, d].
// The id is a unique identifier for that car. The x, y values are in global map coordinates,
// and the vx, vy values are the velocity components, also in reference to the global map.
// Finally s and d are the Frenet coordinates for that car.
const static double getId(const vector<double> theSensorFusion) {
  return theSensorFusion[0];
}

const static int XINDEX = 1;
const static double getX(const vector<double> theSensorFusion) {
  return theSensorFusion[XINDEX];
}

const static int YINDEX = 2;
const static double getY(const vector<double> theSensorFusion) {
  return theSensorFusion[YINDEX];
}

const static int VXINDEX = 3;
const static double getVx(const vector<double> theSensorFusion) {
  return theSensorFusion[VXINDEX];
}

const static int VYINDEX = 4;
const static double getVy(const vector<double> theSensorFusion) {
  return theSensorFusion[VYINDEX];
}

const static int SINDEX = 5;
const static double getS(const vector<double> theSensorFusion) {
  return theSensorFusion[SINDEX];
}

const static int DINDEX = 6;
const static double getD(const vector<double> theSensorFusion) {
  return theSensorFusion[DINDEX];
}

const static double getV(const double theVx, const double theVy) {
  return sqrt(theVx*theVx+theVy*theVy);
}

const static double getV(const vector<double> theSensorFusion) {
  const double vX = getVx(theSensorFusion);
  const double vY = getVy(theSensorFusion);
  return getV(vX, vY);
}

const static int TINDEX = 7;
const static double getT(const vector<double> theSensorFusion) {
  return theSensorFusion[TINDEX];
}

const void setT(vector<double>& theSensorFusion, const double theTime) {
  if (theSensorFusion.size() < (TINDEX+1)) {
    theSensorFusion.resize(TINDEX+1);
  }
  theSensorFusion[TINDEX]=theTime;
}


const static std::string toString(const vector<double> theSensorFusion) {
  std::ostringstream oss;
  if (theSensorFusion.size()>0) {
    oss << "id:" << getId(theSensorFusion) << ", x:"  << getX(theSensorFusion) << ", y:"  << getY(theSensorFusion) << ", vX:"  << getVx(theSensorFusion) << ", vY:"  << getVy(theSensorFusion) << ", s:"  << getS(theSensorFusion) << ", d:"  << getD(theSensorFusion);
    if (theSensorFusion.size() > TINDEX) {
      oss << ", t:" << getT(theSensorFusion);
    }
  }
  return oss.str();
}

const static std::string toString(const vector<vector<double>> theSensorFusionVector) {
  std::ostringstream oss;
  for (int vector=0; vector<theSensorFusionVector.size(); vector++) {
    oss << "#:" << vector << "\n" << toString(theSensorFusionVector[vector]) << std::endl;
  }
  return oss.str();
}

const static string trueFalse(const bool theBoolean) {
  return theBoolean==0?"false":"true";
}

const static std::string toString(const vector<bool> theBooleanVector) {
  std::ostringstream oss;
  for (int boolean=0; boolean<theBooleanVector.size(); boolean++) {
    if (boolean>0) {
      oss << ", ";
    }
    oss << boolean << "? " << trueFalse(theBooleanVector[boolean]);
  }
  return oss.str();
}

const void printSensors(const vector<vector<double>> theSensors) {
  for (int s=0; s<theSensors.size(); s++) {
    cout << "#:" << s << " " << toString(theSensors[s]) << std::endl;
  }
}

const static vector<vector<double>>& sensedInSameLane(const int theTargetLane, const vector<vector<double>> theSensed) {
  vector<vector<double>> *inSameLane = new vector<vector<double>>;
  for (int s=0; s<theSensed.size(); s++) {
    int sensedLane=determineLane(getD(theSensed[s]));
    if (sensedLane == theTargetLane) {
      if (DEBUGVEHICLEPATH) {
        cout << "id:" << getId(theSensed[s]) << ", sensedLane:" << sensedLane << ", theTargetLane:" << theTargetLane << std::endl;
      }
      (*inSameLane).push_back(theSensed[s]);
    }
  }
  return *inSameLane;
}

const static bool isInFrontOfS(const double theS, const vector<double> theSensed) {
  const double sensedS=getS(theSensed);
  return (sensedS>theS);// may be way far away
}

const static bool isBehindS(const double theS, const vector<double> theSensed) {
  const double sensedS=getS(theSensed);
  return (sensedS<theS);// may be way far away
}

const static vector<double> closestSensedInFront(const int theLane, const double theS, const vector<vector<double>> theSensed) {
  const vector<vector<double>> inSameLane = sensedInSameLane(theLane, theSensed);
  vector<double> closestSensed;
  if (DEBUGVEHICLEPATH) {
    cout << "closestInFrontInSameLane-closestSensed.size:" << closestSensed.size() << std::endl;
  }
  for (int s=0; s<inSameLane.size(); s++) {
    const bool inFront=isInFrontOfS(theS, inSameLane[s]);
    if (inFront) {
      if (closestSensed.size()==0) {
        closestSensed=inSameLane[s];
      } else {
        const double sensedS = getS(inSameLane[s]);
        const double closestS = getS(closestSensed);
        if (DEBUGVEHICLEPATH) {
          cout << "closestInFrontInSameLane-s:" << s << ",sensedS:" << sensedS << ", closestS:" << closestS << std::endl;
        }
        if (sensedS < closestS) {
          closestSensed=inSameLane[s];
        }
      }
    }
  }
  if (DEBUGVEHICLEPATH) {
    cout << "closestInFrontInSameLane-closestSensed:" << toString(closestSensed) << std::endl;
  }
  return closestSensed;
}

const static vector<double> closestSensedInFront(const double theD, const double theS, const vector<vector<double>> theSensed) {
  return closestSensedInFront(determineLane(theD), theS, theSensed);
}

const static vector<vector<double>> closestInFrontInLane(const int theNumberOfLanes, const double theS, const vector<vector<double>> theSensed) {
  vector<vector<double>> closestInLane;
  closestInLane.resize(theNumberOfLanes);
  for (int lane=0; lane<theNumberOfLanes; lane++) {
    closestInLane[lane] = closestSensedInFront(lane, theS, theSensed);
  }
  cout << "closestInFrontInLane-closestInLane:" << toString(closestInLane) << std::endl;
  return closestInLane;
  
}

const static vector<double> slowestSensedInFront(const int theTargetLane, const double theS,
                                                 const vector<vector<double>> theSensed) {
  const vector<vector<double>> inSameLane = sensedInSameLane(theTargetLane, theSensed);
  vector<double> slowestSensed;
  if (DEBUGVEHICLEPATH) {
    cout << "slowestSensedInFront-closestSensed.size:" << slowestSensed.size() << std::endl;
  }
  for (int s=0; s<inSameLane.size(); s++) {
    const bool inFront=isInFrontOfS(theS, inSameLane[s]);
    if (inFront) {
      if (slowestSensed.size()==0) {
        slowestSensed=inSameLane[s];
      } else {
        const double sensedV = getV(inSameLane[s]);
        const double slowestV = getV(slowestSensed);
        if (DEBUGVEHICLEPATH) {
          cout << "slowestSensedInFront-theTargetLane:" << theTargetLane << ", s:" << s << ",sensedV:" << sensedV << ", slowestV:" << slowestV << std::endl;
        }
        if (sensedV < slowestV) {
          slowestSensed=inSameLane[s];
        }
      }
    }
  }
  if (DEBUGVEHICLEPATH) {
    cout << "slowestSensedInFront-closestSensed:" << toString(slowestSensed) << std::endl;
  }
  return slowestSensed;
}

const static vector<vector<double>> slowestInFrontInAnyLane(const int theNumberOfLanes, const double theS, const vector<vector<double>> theSensed) {
  vector<vector<double>> slowestInLane;
  slowestInLane.resize(theNumberOfLanes);
  for (int lane=0; lane<theNumberOfLanes; lane++) {
    slowestInLane[lane] = slowestSensedInFront(lane, theS, theSensed);
  }
  if (DEBUGVEHICLEPATH) {
    cout << "slowestInFrontInAnyLane-slowestInLane:" << toString(slowestInLane) << std::endl;
  }
  return slowestInLane;
}

const static vector<double> closestBehindInLane(const int theTargetLane, const double theS, const vector<vector<double>> theSensed) {
  const vector<vector<double>> inTargetLane = sensedInSameLane(theTargetLane, theSensed);
  vector<double> closestSensed;
  cout << "closestBehindInLane-closestSensed.size:" << closestSensed.size() << std::endl;
  for (int s=0; s<inTargetLane.size(); s++) {
    const bool isBehind=isBehindS(theS, inTargetLane[s]);
    if (isBehind) {
      if (closestSensed.size()==0) {
        closestSensed=inTargetLane[s];
      } else {
        const double sensedS = getS(inTargetLane[s]);
        const double closestS = getS(closestSensed);
        cout << "closestBehindInLane-theTargetLane:" << theTargetLane << ", s:" << s << ",sensedS:" << sensedS << ", closestS:" << closestS << std::endl;
        if (sensedS > closestS) {
          closestSensed=inTargetLane[s];
        }
      }
    }
  }
  cout << "closestBehindInLane-closestSensed:" << toString(closestSensed) << std::endl;
  return closestSensed;
}

const static vector<vector<double>> closestBehindInAnyLane(const int theNumberOfLanes, const double theS, const vector<vector<double>> theSensed) {
  vector<vector<double>> closestInLane;
  closestInLane.resize(theNumberOfLanes);
  for (int lane=0; lane<theNumberOfLanes; lane++) {
    closestInLane[lane] = closestBehindInLane(lane, theS, theSensed);
  }
  cout << "closestBehindInAnyLane-closestInLane:" << toString(closestInLane) << std::endl;
  return closestInLane;
}

const vector<bool> findFasterLane(const int theCurrentLane, const vector<vector<double>> theSensedSlowest) {
  const double SPEEDFACTORTOBEFASTER=1.2;
  const double benchmarkSpeed = getV(theSensedSlowest[theCurrentLane]);
  vector<bool> fasterLane(theSensedSlowest.size(), false);
  fasterLane.resize(theSensedSlowest.size());
  for (int lane=0; lane<theSensedSlowest.size(); lane++) {
    const bool isCurrentLane = (lane==theCurrentLane);
    if (!isCurrentLane) {
      const bool isEmptyLane = theSensedSlowest[lane].size()==0;
      if (isEmptyLane) {
        fasterLane[lane] = true;
      } else {// if it is not empty, is it better than the benchmark (.i.e. the current lane)
        const bool isBetterThanBenchmark = getV(theSensedSlowest[lane]) > (SPEEDFACTORTOBEFASTER*benchmarkSpeed);
        if (isBetterThanBenchmark) {
          fasterLane[lane]=true;
        }
      }
    }
  }
  if (DEBUGVEHICLEPATH) {
    cout << "findFasterLane-fasterLane:" << toString(fasterLane) << std::endl;
  }
  return fasterLane;
}

const vector<bool> findSafeLane(const int theCurrentLane, const double theFutureLocation, const double theSafeTimeGap, const vector<vector<double>> theSensedBehind, const double theCurrentLocation, const double theCurrentSpeed, const vector<vector<double>> theSensedInFront) {
  cout << "findSafeLane-theSensedBehind?" << printVector(theSensedBehind) << std::endl;
  cout << "findSafeLane-theSensedInFront?" << printVector(theSensedInFront) << std::endl;
  vector<bool> safeLane(theSensedBehind.size(), false);
  for (int lane=0; lane<theSensedBehind.size(); lane++) {
    const bool isCurrentLane = (lane==theCurrentLane);
    if (isCurrentLane) {// the current lane is always safe
      safeLane[lane]=true;
    } else {
      const bool isEmptyLane = theSensedBehind[lane].size()==0;
      if (isEmptyLane) {
        safeLane[lane] = true;
      } else {// if it is not empty
        const double sBehind = getS(theSensedBehind[lane]);
        const double vBehind = getV(theSensedBehind[lane]); // car can accelerate which is not being tracked
        const double deltaS = theSafeTimeGap*vBehind; // dS=v*dT
        const double sBehindAtT = sBehind+deltaS;// location of this car at the end of theSafeTimeGap
        const bool inFrontOfCarBehind = (sBehindAtT<theFutureLocation);
        const double sInFront = getS(theSensedInFront[lane])+(theCurrentSpeed*1./*second*/);
        const bool isBehindCarInFront = sInFront>theCurrentLocation;
        cout << "findSafeLane-inFrontOfCarBehind?" << trueFalse(inFrontOfCarBehind) << ", isBehindCarInFront? "
          << trueFalse(isBehindCarInFront) << ", sBehind:" << sBehind << ", vBehind:" << vBehind << ", sBehindAtT:"
          << sBehindAtT  << std::endl;
        if (inFrontOfCarBehind && isBehindCarInFront) {
          safeLane[lane]=true;
        }
      }
    }
  }
  cout << "findSafeLane-theCurrentLane:" << theCurrentLane << ", theFutureLocation:" << theFutureLocation <<  ", theSafeTimeGap:" << theSafeTimeGap << std::endl;
  cout << "findSafeLane-safeLane:" << toString(safeLane) << std::endl;
  return safeLane;
}

const int findFasterSafeLane(const int theCurrentLane, const vector<bool>& theSafeLanes, const vector<bool>& theFasterLanes) {
  for (int lane=0; lane<theSafeLanes.size(); lane++) {
    const bool isCurrentLane = (lane==theCurrentLane);
    const bool isOnlyOneLaneOver = abs(lane-theCurrentLane) == 1;
    if (!isCurrentLane && theSafeLanes[lane] && isOnlyOneLaneOver && theFasterLanes[lane]) {
      return lane;
    }
  }
  return -1;
}

static vector<vector<vector<double>>> sensedVehicles; // by car id, by rolling buffer size, by sensor id
static int SENSEDHISTORYLENGTH = 10; // approx seconds

const static void recordSensorReading(const vector<double> theSensed) {
  const int sensorId = getId(theSensed);
  vector<vector<double>>& sensedReadings=sensedVehicles[sensorId];
  //cout << "recordSensorReading-before-sensedReadings[" << sensorId << "]:" << printVector(sensedReadings) << std::endl;
  if (sensedReadings.size() >= SENSEDHISTORYLENGTH) {// sensed buffer is full, remove oldest;
    sensedReadings.erase(sensedReadings.begin());
  }
  sensedReadings.push_back(theSensed);
  //cout << "recordSensorReading-after-sensedVehicles[" << sensorId << "]:" << toString(sensedVehicles[sensorId]) << std::endl;
}

const static void recordSensorReading(vector<double> theSensed, const double theTime) {
  setT(theSensed, theTime);
  const int sensorId = getId(theSensed);
  const bool newlySensedVehicle = sensedVehicles.size()<(sensorId+1);// sensorId is index, so size must be at least id+1
  //cout << "recordSensorReading-newlySensedVehicle:" << newlySensedVehicle << ", sensorId:" << sensorId << ", sensedVehicles.size:" << sensedVehicles.size() << std::endl;
  if (newlySensedVehicle) {
    sensedVehicles.resize(sensedVehicles.size()+20); // increase number of sensedVehicles by 100
    //cout << "recordSensorReading-sensedVehicles.(re)size:" << sensedVehicles.size() << std::endl;
  }
  setT(theSensed, theTime); // record sensor time
  recordSensorReading(theSensed);
  //cout << "recordSensorReading-sensedVehicles" << "[" << sensorId << "]:" << toString(sensedVehicles[sensorId]) << std::endl;
  //cout << "recordSensorReading-sensedVehicles.size:" << sensedVehicles.size() << std::endl;
}

const static void recordSensorReadings(const vector<vector<double>> theSensed, const double theTime) {
  //cout << "recordSensorReadings-before-sensedVehicles:" << printVector(sensedVehicles) << std::endl;
  for (int reading=0; reading<theSensed.size(); reading++) {
    recordSensorReading(theSensed[reading], theTime);
  }
  if (DEBUGSENSORREADINGS) {
    cout << "recordSensorReadings-after-sensedVehicles:" << printVector(sensedVehicles) << std::endl;
  }
}

const void printVehiclePath(const tk::spline& thePath, const double theStartingX, const double theEndingX, const double theDelta) {
  for (double xValue=theStartingX; xValue<=theEndingX; xValue+=theDelta) {
    const double yValue=thePath(xValue);
    cout << "printVehiclePath-recordSensorReadings-xValue:" << xValue << " -> " << yValue << std::endl;
  }
}

const tk::spline createVehiclePath(const vector<double>& theSensedVehicle) {
  const int vehicleId = getId(theSensedVehicle);
  const vector<vector<double>> vehicleHistory = sensedVehicles[vehicleId];
  vector<double> tPoints;// x
  vector<double> sPoints;// y
  for (int history=0; history<vehicleHistory.size(); history++) {
    tPoints.push_back(getT(vehicleHistory[history]));
    sPoints.push_back(getS(vehicleHistory[history]));
  }
  tk::spline vehiclePath;
  if (DEBUGLANESELECTION) {
    cout << "createVehiclePath-tPoints:" << printVector(tPoints) << "\nsPoints:" << printVector(sPoints) << std::endl;
  }
  vehiclePath.set_points(tPoints, sPoints);
  const double startingX=getT(vehicleHistory[0]);
  const double endingX=getT(vehicleHistory[vehicleHistory.size()-1]);
  const double deltaX=1.;
  if (DEBUGLANESELECTION) {
    printVehiclePath(vehiclePath, startingX, endingX, deltaX);
  }
  return vehiclePath;
}

const tk::spline createVehiclePath(const double sStart, const double sVelocity, double theStartTime, const double theEndTime, const double theDeltaTime) {
  vector<double> tPoints;// x
  vector<double> sPoints;// y
  double sVehicle=sStart;
  for (int time=theStartTime; time<theEndTime; time+=theDeltaTime) {
    tPoints.push_back(time);
    sPoints.push_back(sVehicle);
    sVehicle+=sVelocity*theDeltaTime;
  }
  tk::spline vehiclePath;
  if (DEBUGLANESELECTION) {
    cout << "createVehiclePath\ntPoints:" << printVector(tPoints) << "\nsPoints:" << printVector(sPoints) << ", sStart:" << sStart << ", sVelocity:" << sVelocity << ", theDeltaTime:" << theDeltaTime << std::endl;
  }
  vehiclePath.set_points(tPoints, sPoints);
  if (DEBUGLANESELECTION) {
    printVehiclePath(vehiclePath, theStartTime, theEndTime, theDeltaTime);
  }
  return vehiclePath;
}

const bool pathsAreSeparated(const tk::spline& theFirstPath, const tk::spline& theSecondPath, const double thePathSeparation, const double theStartTime, const double theEndTime, const double theDeltaTime) {
  for (double time=theStartTime; time<theEndTime; time+=theDeltaTime ) {
    const double s1 = theFirstPath(time);
    const double s2 = theSecondPath(time);
    const bool isSeparated = abs(s2-s1)>abs(thePathSeparation);
    if (DEBUGLANESELECTION) {
      cout << "pathsAreSeparated-s1[" << time << "]:" << s1 << ", s2[" << time << "]:" << s2 << ", isSeparated? " << isSeparated << ", thePathSeparation:" << thePathSeparation << std::endl;
    }
    if (!isSeparated) {
      return false;
    }
  }
  return true;
}

const bool isLaneSafe(const int theLane, const tk::spline& theVehiclePath, const double theVehicleSize, const double theStartTime, const double theEndTime, const double theDeltaTime, const vector<vector<double>> theSensed) {
  vector<vector<double>> inTheLane=sensedInSameLane(theLane, theSensed); // all sensed vehicles in this lane
  for (int vehicle=0; vehicle<inTheLane.size(); vehicle++) {// for all vehicles in lane
    const vector<double> otherVehicle=inTheLane[vehicle];
    const tk::spline theOtherVehiclePath = createVehiclePath(otherVehicle);
    const bool separated=pathsAreSeparated(theVehiclePath, theOtherVehiclePath, theVehicleSize, theStartTime, theEndTime, theDeltaTime);
    const int vehicleId = getId(otherVehicle);
    if (DEBUGLANESELECTION) {
      cout << "isLaneSafe-vehicleId:" << vehicleId << ", separated? " << trueFalse(separated) << ", theLane:" << theLane << std::endl;
    }
    if (!separated) {
      return false;
    }
  }
  return true;

}

const vector<bool> isLaneSafe(const double theVehicleS, const double theVehicleSVelocity, const double theVehicleSize, const double theStartTime, const double theEndTime, const double theDeltaTime, const vector<vector<double>> theSensed) {
  vector<bool> isSafe(3,false);
  if (DEBUGLANESELECTION) {
    cout << "isLaneSafe-theVehicleS:" << theVehicleS << ", theVehicleSVelocity: " << theVehicleSVelocity << ", theVehicleSize:" << theVehicleSize << std::endl;
  }
  const tk::spline vehiclePath = createVehiclePath(theVehicleS, theVehicleSVelocity, theStartTime, theEndTime, theDeltaTime);
  for (int lane=0; lane<2; lane++) {
    isSafe[lane]=isLaneSafe(lane, vehiclePath, theVehicleSize, theStartTime, theEndTime, theDeltaTime, theSensed);
  }
  return isSafe;
}
  
static bool CANCHANGELANES=true;
static bool isChangingLanes=false;
static int targetLane=-1;
static double lastRecordingT = -1;
static double currentT = 0;
static double v0;
static double v1;
static bool stateTransition=true;

int main(int argc, char* argv[]) {
  
  testXY();
  testRT();
  testTR();
  testCarMap();
  
  const int NUMBEROFLANES=3;

  uWS::Hub h;
  
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  /*
   Each waypoint has an (x,y) global map position, and a Frenet s value and Frenet d unit normal vector (split up into the x
   component, and the y component). The s value is the distance along the direction of the road. The first waypoint has an s value of
   0 because it is the starting point.
   */
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;
  
  // Waypoint map to read from
  std::string argv_str(argv[0]);
  std::string base = argv_str.substr(0, argv_str.find_last_of("/"));
  cout << "argv_str:" << argv_str << std::endl;
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  // The track is 6945.554 meters around (about 4.32 miles)
  // If the car averages near 50 MPH, then it should take a little more than 5 minutes
  double max_s = 6945.554;
  
  ifstream in_map_(map_file_.c_str(), ifstream::in);
  
  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }
  
  //cout << "map_waypoints_s:" << printVector(map_waypoints_s) << std::endl;
  
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                                                                                                       uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      
      auto s = hasData(data);
      
      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];// in degrees, use deg2rad(car_yaw) to convert
          double car_speed = j[1]["speed"];// miles per hour
          
          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];
          
          if (car_speed > MAXIMUMMILESPERHOUR) {
            cout << std::setprecision(6) << "\n\n-------->car_speed:" << car_speed << ", car_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s
            << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << std::endl;
          }
          
          if (DEBUGSACCELERATION) {
            cout << std::setprecision(6) << "\n\ncar_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << std::endl;
          }
          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
          if (DEBUGLANESELECTION) {
            printSensors(sensor_fusion);
          }
          const vector<double> closestInLane = closestSensedInFront(car_d, car_s, sensor_fusion);
          const bool isCarInFrontInLane = closestInLane.size()>0;
          if (DEBUGVEHICLEPATH) {
            cout << "isCarInFrontInLane?" << isCarInFrontInLane << ", closestInLane:" << toString(closestInLane) << std::endl;
          }

          json msgJson;
          
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          
          const double PATHPOINTDELTAT=1./double(LENGTHOFCARPATHTRAJECTORY);
          if (previous_path_x.size() > 0) {// if it is not the initial message, update the current time
            const int movesCompleted=LENGTHOFCARPATHTRAJECTORY-previous_path_x.size();
            const double timePassed = PATHPOINTDELTAT*movesCompleted;
            currentT += timePassed;
            if (DEBUGSENSORREADINGS) {
              cout << std::setprecision(6) << "\ncurrentT:" << currentT <<  ", timePassed:" << timePassed << ", movesCompleted:" << movesCompleted << std::endl;
            }
            if ((currentT-lastRecordingT >= 1.)) {// record every second
              recordSensorReadings(sensor_fusion, currentT);
              lastRecordingT=currentT;
            }
          }

          // (re)use what is left in previous path
          for (int i=0; i<previous_path_x.size(); i++) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }
          if (DEBUGSENSORREADINGS) {
            cout << "\nrecovered:" << next_x_vals.size() << " values from previous_path" << ", currentT:" << currentT << std::endl;
            if (previous_path_x.size() > 2) {
              cout << "next_x_vals:" << next_x_vals[0] << "," << next_x_vals[1] << "..." <<  next_x_vals[next_x_vals.size()-2] << ", "
              <<  next_x_vals[next_x_vals.size()-1]
              << "\nnext_y_vals:" << next_y_vals[0] << "," << next_y_vals[1] << "..." <<  next_y_vals[next_y_vals.size()-2]  <<  ", "
              << next_y_vals[next_y_vals.size()-1]  << std::endl;
            }
          }


          /*
           Start by simply trying to move the car forward in a straight line at a constant 50 MPH velocity. Use the car's (x, y) localization information and its heading direction to create a simple, straight path that is drawn directly in front of the car.
           
           In main.cpp, instead of setting the speed directly, we pass next_x_vals, and next_y_vals to the simulator. We will set the points 0.5 m apart. Since the car moves 50 times a second, a distance of 0.5m per move will create a velocity of 25 m/s. 25 m/s is close to 50 MPH.
           
           (.5 meters/move)(50 moves/second)=25 meters/second
           (25 m/sec)(3.2804 ft/m)(1/5280 mi/ft)(60 sec/min)(60 min/hr)=55.92 mi/hr
           */
          targetLane=(targetLane==-1)?determineLane(car_d):targetLane;//set targetLane on start-up
          const vector<vector<double>> lastXyPathPoints =
          previous_path_x.size()<2?last2PathPoints(car_x, car_y, deg2rad(car_yaw)):last2PathPoints(previous_path_x, previous_path_y);
          if (DEBUGPATHSPLINE) {
            cout << std::setprecision(6) << "car_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << ", car_speed:" << car_speed << " mph" << std::endl;
          }
          vector<double> xMapPoints=lastXyPathPoints[0];
          vector<double> yMapPoints=lastXyPathPoints[1];
          const double pathInterpolatedYaw=yawFromXY(lastXyPathPoints);
          if (DEBUGPATHSPLINE) {
            cout << "xMapPoints:" << printVector(xMapPoints, 10) << std::endl;
            cout << "yMapPoints:" << printVector(yMapPoints, 10) << std::endl;
            cout << std::setprecision(6) << "pathInterpolatedYaw:" << pathInterpolatedYaw << ", car_yaw:" << deg2rad(car_yaw) << std::endl;
          }
          // add additional waypoints in front of car in map coordinates, the spline needs more than 2 points
          //const int currentLane=determineLane(car_d);
          const int NUMBEROFSPLINEPOINTS=3;
          //const double deltaS = (3.*PATHLOOKAHEADSDISTANCE)/double(NUMBEROFSPLINEPOINTS);
          const double pathLookaheadDistance=isChangingLanes?PATHLOOKAHEADSDISTANCE:PATHLOOKAHEADSDISTANCE;
          for (int p=0; p<NUMBEROFSPLINEPOINTS; p++) {
            double lookAheadS = car_s+(p+1)*pathLookaheadDistance;
            //double lookAheadS = car_s+(p+1)*deltaS;
            vector<double> mapXY = getXY(lookAheadS, laneCenter(targetLane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            xMapPoints.push_back(mapXY[0]);
            yMapPoints.push_back(mapXY[1]);
          }
          // xPoints[1], yPoints[1] is last projected location for the car (usually) using the previous path
          const vector<vector<double>>& classTransformedXY=rtClassTransform(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xMapPoints, yMapPoints);// car_yaw for 1st 2 points is around 0 & should stay 0 for last 3 added for spline to be flat
          const vector<double> xTransformedPoints=classTransformedXY[0];
          const vector<double> yTransformedPoints=classTransformedXY[1];
          if (DEBUGPATHSPLINE) {
            cout << "class-xTransformedPoints:" << printVector(xTransformedPoints, 10) << std::endl;
            cout << "class-yTransformedPoints:" << printVector(yTransformedPoints, 10) << std::endl;
            cout << "affine-xPoints:" << printVector(xMapPoints, 10) << std::endl;
            cout << "affine-yPoints:" << printVector(yMapPoints, 10) << std::endl;
            cout << std::setprecision(6) << "pathInterpolatedYaw:" << pathInterpolatedYaw << ", car_yaw:" << deg2rad(car_yaw) << std::endl;
          }
          const vector<Eigen::Vector2d> carXYFromMapXY=transformFromMapToCar(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xMapPoints, yMapPoints);
          //cout << "affine-carXYFromMapXY:\n" << printVector(carXYFromMapXY, 10) << std::endl;
          
          // from v1 to an eignen vector
          // float* ptr_data = &v1[0];
          // Eigen::VectorXf v2 = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(v1.data(), v1.size());
          // from the eigen vector to the std vector
          // std::vector<float> v3(&v2[0], v2.data()+v2.cols()*v2.rows());
          // to check
          // for(int i = 0; i < v1.size() ; i++){
          //   std::cout << std::to_string(v1[i]) << " | " << std::to_string(v2[i]) << " | " << std::to_string(v3[i]) << std::endl;
          // }
          //from the eigen vector to the std vector
          //cout << "carPointsFromMap:\n" << printVector(carPointsFromMap, 10) << std::endl;
          vector<vector<double>> carXY=unpackV2D(carXYFromMapXY);
          if (DEBUGPATHSPLINE) {
            cout << "carXY[0]:\n" << printVector(carXY[0], 10) << std::endl;
            cout << "carXY[1]:\n" << printVector(carXY[1], 10) << std::endl;
          }
          
          const double patLookaheadPlanningDistance=PATHLOOKAHEADSDISTANCE;// mixing s & x==PATHLOOKAHEADSDISTANCE
          //const double targetVelocity=isCarInFrontInLane? min(MAXIMUMMETERSPERSECOND, getV(closestInLane)): MAXIMUMMETERSPERSECOND;
          double targetVelocity=.90*MAXIMUMMETERSPERSECOND;
          double targetAcceleration=0.75*MAXIMUMACCERLERATION;
          
          double projectedPathVelocity=previous_path_x.size()<2?
          0. : velocityFromXY(PATHPOINTDELTAT/*theDeltaT*/, lastXyPathPoints, pathInterpolatedYaw, map_waypoints_x, map_waypoints_y)[0];
          if (DEBUGVEHICLEPATH) {
            cout << "projectedPathVelocity:" << projectedPathVelocity << ", pathInterpolatedYaw:" << pathInterpolatedYaw << std::endl;
            cout << std::setprecision(6) << "car_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << ", car_speed:" << car_speed << " mph, closestInLane:" << (closestInLane.size()>0?getV(closestInLane):-1.) << std::endl;
          }
          
          if (projectedPathVelocity > MAXIMUMMETERSPERSECOND) {
            cout << "\n\ncar_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s
            << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << ", projectedPathVelocity:" << projectedPathVelocity << std::endl;
          }

          
          if (isCarInFrontInLane) {
            const double closestInLaneDeltaS = getS(closestInLane)-car_s;
            const bool closeToPlanningHorizon=closestInLaneDeltaS<1.35*patLookaheadPlanningDistance;
            if (DEBUGVEHICLEPATH) {
            cout << "closeToPlanningHorizon:" << closeToPlanningHorizon << ", closestInLaneDeltaS:" << closestInLaneDeltaS
            << ", getS:" << getS(closestInLane) << ", (2x)patLookaheadPlanningDistance:" << 2.*patLookaheadPlanningDistance
            << ", closestInLaneDeltaS:" << closestInLaneDeltaS << std::endl;
            }
            if (closeToPlanningHorizon) {
              if (DEBUGVEHICLEPATH) {
              cout << "before-isChangingLanes:" << isChangingLanes << ", in targetLane? " << (determineLane(car_d)==targetLane) << ", targetVelocity:" << targetVelocity <<  ", targetLane:" << targetLane << std::endl;
              }
              if (CANCHANGELANES && isChangingLanes) {// if i'm changing lanes
                if (determineLane(car_d)==targetLane) {// have already changed?
                  isChangingLanes=false;// stop lane changing
                  stateTransition=false;
                } else {// changing, maintain velocity
                  //targetVelocity=currentCarPassingVelocity; // maintain car passing velocity to avoid acceleration errors
                }
              } else {// if i'm not changing & i'm inside the planning horizon, can i change?
                //const vector<vector<double>> closestSensedBehind = closestBehindInAnyLane(NUMBEROFLANES, car_s, sensor_fusion);
                const vector<vector<double>> slowestSensedInFront = slowestInFrontInAnyLane(NUMBEROFLANES, car_s, sensor_fusion);
                //const vector<vector<double>> closestInFront = closestInFrontInLane(NUMBEROFLANES, car_s, sensor_fusion);
                const vector<bool> fasterLanes = findFasterLane(determineLane(car_d), slowestSensedInFront);
                //const vector<bool> safeLanes = findSafeLane(determineLane(car_d), end_path_s, PATHLOOKAHEADSDISTANCE/milesPerHourToMetersPerSecond(car_speed), closestSensedBehind, car_s, milesPerHourToMetersPerSecond(car_speed), closestInFront);
                const vector<bool> safeLanes=isLaneSafe(car_s, targetVelocity, 25./*size of vehicle*/, currentT, currentT+PATHLOOKAHEADSDISTANCE/PATHPOINTDELTAT, 60.*(50.*PATHPOINTDELTAT/* 1 second */) /* 1 minute */ , sensor_fusion);
                vector<bool> dummy(3,false);
                const int safeLane=findFasterSafeLane(determineLane(car_d), safeLanes, fasterLanes);
                if (CANCHANGELANES && safeLane>=0) {
                  targetLane=safeLane;
                  isChangingLanes=true;
                  stateTransition=true;
                  if (DEBUGVEHICLEPATH) {
                    cout << "changing-isChangingLanes:" << isChangingLanes << ", targetVelocity:" << targetVelocity << std::endl;
                  }
                } else {// inside the planning & not changing lanes, track vehicle in front
                  targetVelocity=getV(closestInLane);// track vehicle in front, unless i'm changing lanes
                  if (DEBUGVEHICLEPATH) {
                    cout << "not changing-isChangingLanes:" << isChangingLanes << ", targetVelocity:" << targetVelocity << std::endl;
                  }
                }
              }
              if (DEBUGVEHICLEPATH) {
                cout << "after-isChangingLanes:" << isChangingLanes << ", targetVelocity:" << targetVelocity << std::endl;
              }
              const double endDeltaT=PATHPOINTDELTAT*next_x_vals.size()/*number of previous points*/;// time in previous path = dt*N
              const double endDeltaS = endDeltaT*getV(closestInLane);// dT*V=dS, distance traveled in dT
              const double closestEndS = getS(closestInLane)+endDeltaS;// s+dS, position of closest car t the end of the previous path
              const double closestDeltaS = closestEndS-end_path_s;// s to closest car at the end of the previous path (current planning cycle)
              //const double deltaD = deltaS-patLookaheadPlanningDistance;// compare distance to closest & planning horizon
              const bool insidePlanningHorizon = closestDeltaS<patLookaheadPlanningDistance;// compare distance to closest & planning horizon
              if (DEBUGVEHICLEPATH) {
                cout << "insidePlanningHorizon:" << insidePlanningHorizon /*<< ", deltaD:" << deltaD*/ << ", closestDeltaS:" << closestDeltaS
                  << ", endDeltaT:" << endDeltaT << ", end_path_s:" << end_path_s << ", patLookaheadPlanningDistance:"
                  << patLookaheadPlanningDistance << ", endDeltaS:" << endDeltaS << ", getV:" << getV(closestInLane)
                  << ", getS:" << getS(closestInLane) << ", closestEndS:" << closestEndS << std::endl;
              }
              if (insidePlanningHorizon) {// adjust velocity
                const double deltaS = patLookaheadPlanningDistance-endDeltaS;// shortfall
                const double percent = deltaS/patLookaheadPlanningDistance;
                const double velocityFactor = 1.-max(0.,percent/* 0 to 1*/);// 1 to 0
                //const double deltaV = getV(closestInLane)-projectedPathVelocity;// speed difference
                targetVelocity = velocityFactor*getV(closestInLane);
                //const double tAtEndOfPreviousPath = next_x_vals.size()*PATHPOINTDELTAT;
                //const double deltaT = LENGTHOFCARPATHTRAJECTORY*PATHPOINTDELTAT-tAtEndOfPreviousPath;
                //const double accelerationToMatchV = deltaV/deltaT;
                const bool isSlowingDown = getV(closestInLane)<projectedPathVelocity;
                const bool isCloseToCarInFront = closestDeltaS < .5*patLookaheadPlanningDistance;
                if ( isSlowingDown && isCloseToCarInFront) {// getting desperate?
                  targetAcceleration=MAXIMUMACCERLERATION;
                  const double timeToIntercept = closestInLaneDeltaS/milesPerHourToMetersPerSecond(car_speed);// miles per hour!!
                  cout << "!!------------>targetVelocity:" << targetVelocity << ", timeToIntercept:" << timeToIntercept << std::endl;
                }
                if (DEBUGVEHICLEPATH) {
                  cout << "insidePlanningHorizon-targetVelocity:" << targetVelocity << ", percent:" << percent
                    << ", velocityFactor:" << velocityFactor
                    << ". targetAcceleration:" << targetAcceleration << ", deltaS:" << deltaS << ", isSlowingDown:"
                    << isSlowingDown << ", isCloseToCarInFront" << isCloseToCarInFront << std::endl;
                }
              }
              //const int addtionalPathPoints=LENGTHOFCARPATHTRAJECTORY-next_x_vals.size();
              //const double gapDeltaT=PATHPOINTDELTAT*addtionalPathPoints;// time
              //const double numberOfPointsToMatch = deltaD/deltaT;
              //cout << "velocityToMatchS:" << velocityToMatchS << ", deltaD:" << deltaD << ", deltaS:" << deltaS << ", end_path_s:" << end_path_s << ", endDeltaS:" << endDeltaS << ", endDeltaT" << endDeltaT << std::endl;
            }
           
            // if distanceTo<PATHLOOKAHEADSDISTANCE is close then use min, is less then reduce
            //car_targetVelocity=(isCarInFrontInLane && insideThePlanningHorizon)?min(MAXIMUMMETERSPERSECOND, getV(closestInLane)):targetVelocity;
            if (DEBUGVEHICLEPATH) {
              cout << "insidePlanningHorizon-deltaD:" << closestInLaneDeltaS << targetVelocity << std::endl;
            }
          }

          
          const double TARGETMILESPERHOUR=MAXIMUMMILESPERHOUR;// miles per hour
          const double TARGETMETERSPERSECOND=milesPerHourToMetersPerSecond(TARGETMILESPERHOUR);
          const double PATHDELTADISTANCEMAX=PATHPOINTDELTAT*MAXIMUMMETERSPERSECOND;
          if (DEBUGVEHICLEPATH) {
            std::cout << "MAXIMUMMILESPERHOUR:" << MAXIMUMMILESPERHOUR << ", MAXIMUMMETERSPERSECOND:" << MAXIMUMMETERSPERSECOND <<std::endl;
            std::cout << "LENGTHOFCARPATHTRAJECTORY:" << LENGTHOFCARPATHTRAJECTORY << ", PATHPOINTDELTAT:" << PATHPOINTDELTAT << ", PATHDELTADISTANCEMAX:" << PATHDELTADISTANCEMAX << std::endl;
          }
        
          tk::spline sPathFromCar;
          if (DEBUGVEHICLEPATH) {
            cout << "carXY" << printVector(carXY) << std::endl;
          }
          sPathFromCar.set_points(carXY[0], carXY[1]);
          const double sAtPATHLOOKAHEADPLANNINGDISTANCE=sPathFromCar(PATHLOOKAHEADXDISTANCE);
          const double lengthOfLinearFitToSpline=sqrt(PATHLOOKAHEADXDISTANCE*PATHLOOKAHEADXDISTANCE+sAtPATHLOOKAHEADPLANNINGDISTANCE*sAtPATHLOOKAHEADPLANNINGDISTANCE);
          // spline points in car coordinates that will fill in next_?_vals after those copied over from previous_path
          int addtionalCarPoints=LENGTHOFCARPATHTRAJECTORY-next_x_vals.size();
          vector<double> xAdditionalCarPoints;
          vector<double> yAdditionalCarPoints;
          double carPointX=0.;
          const int ACCELERATIONSTEPS=.25*LENGTHOFCARPATHTRAJECTORY;
          double pathPointVelocity=projectedPathVelocity;// initial value (reducing acceleration discontinuity)
          for (int additionalPoint=0; additionalPoint<addtionalCarPoints; additionalPoint++) {
            
            double deltaV=targetVelocity-pathPointVelocity; // m/s
            double signOfDeltaV = deltaV<0.?-1.:1.;
            double transitionAcceleration = targetAcceleration;
            
            if (stateTransition) {
              transitionAcceleration=isChangingLanes?transitionAcceleration:.5*transitionAcceleration;
              if (DEBUGVEHICLEPATH) {
                cout << "stateTransition:" << stateTransition << ", transitionAcceleration:" << transitionAcceleration << ", targetAcceleration:" << targetAcceleration << std::endl;
              }
            }
            
            pathPointVelocity+=(signOfDeltaV*transitionAcceleration*PATHPOINTDELTAT); // v1=v0+a*dt
            
            if (pathPointVelocity > MAXIMUMMETERSPERSECOND) {
              cout << std::setprecision(6) << "\n\n-------->pathPointVelocity:" << pathPointVelocity << ", car_x:" << car_x
                << ", car_y:" << car_y << ", car_s:" << car_s
                << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << std::endl;
            }
            if (DEBUGVEHICLEPATH) {
              cout << "pathPointVelocity:" << pathPointVelocity << ", targetVelocity:" << targetVelocity << ", transitionAcceleration:" << transitionAcceleration << ", targetAcceleration:" << targetAcceleration << ", deltaV:" << deltaV << ", projectedPathVelocity:" << projectedPathVelocity << std::endl;
            }
            //pathPointVelocity=targetVelocity;
            //const int numberOfSamples=lengthOfLinearFitToSpline/(PATHPOINTDELTAT*TARGETMETERSPERSECOND);// meters/(seconds*metes/second)
            const int numberOfSamples=lengthOfLinearFitToSpline/(PATHPOINTDELTAT*abs(pathPointVelocity));// meters/(seconds*metes/second)
            if (DEBUGVEHICLEPATH) {
              //const int classN = lengthOfLinearFitToSpline/(.02*49.5/2.24);
              const int classN = lengthOfLinearFitToSpline/(.02*abs(pathPointVelocity)/2.24);
              cout << "numberOfSamples:" << numberOfSamples << ", next_x_vals.size:" << next_x_vals.size()
              << ", targetVelocity:" << targetVelocity << ", pathPointVelocity:" << pathPointVelocity << ", classN:" << classN
                <<std::endl;
            }
            const double deltaX=PATHLOOKAHEADXDISTANCE/double(numberOfSamples);
            
            // additionalPoint+1: so deltaX is never 0, so next point added will not match last point in previous path
            carPointX+=deltaX;
            xAdditionalCarPoints.push_back(carPointX);
            const double carPointY=sPathFromCar(carPointX);
            yAdditionalCarPoints.push_back(carPointY);
            
            v0=v1;
            v1=pathPointVelocity;
            const double acceleration=(v1-v0)/PATHPOINTDELTAT;// a=dv/dt
            if (DEBUGSACCELERATION) {
              cout << "acceleration:" << acceleration << ", v1:" << v1 << ", v0:" << v0 << "(v1-v0):" << (v1-v0) << ", deltaV:" << deltaV <<std::endl;
            
              if (abs(acceleration) > abs(MAXIMUMACCERLERATION)) {
                cout << "------------------->acceleration:" << acceleration << ", pathPointVelocity:" << pathPointVelocity << ", v1:" << v1 << ", v0:" << v0 << ", deltaV:" << deltaV << ", targetAcceleration:" << targetAcceleration << ", PATHPOINTDELTAT:" << PATHPOINTDELTAT << "(v1-v0):" << (v1-v0) << std::endl;
              }
            }
          }
          if (projectedPathVelocity > .9*targetVelocity) {
            stateTransition=false;
          }
          // transform the additional car xy point back to the map
          const vector<vector<double>>& classMapPointsXY=trClassTransform(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xAdditionalCarPoints, yAdditionalCarPoints);
          const vector<double> xClassMapPoints=classMapPointsXY[0];
          const vector<double> yClassMapPoints=classMapPointsXY[1];
          if (DEBUGPATHSPLINE) {
            cout << "class-classMapPointsXY:" << printVector(classMapPointsXY, 10) << std::endl;
            cout << "class-xClassMapPoints:" << printVector(xClassMapPoints, 10) << std::endl;
            cout << "class-yClassMapPoints:" << printVector(yClassMapPoints, 10) << std::endl;
          }
          const vector<Eigen::Vector2d> additionalMapXYPoints=transformFromCarToMap(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xAdditionalCarPoints, yAdditionalCarPoints);
          //cout << "affine-additionalMapXYPoints:\n" << printVector(additionalMapXYPoints, 10) << std::endl;
          
          // add new xy points to previous path from car
          vector<vector<double>> mapXYPoints=unpackV2D(additionalMapXYPoints);
          for (int additionalPoint=0; additionalPoint<addtionalCarPoints; additionalPoint++) {
            next_x_vals.push_back(additionalMapXYPoints[additionalPoint][0]);
            next_y_vals.push_back(additionalMapXYPoints[additionalPoint][1]);
          }
          if (DEBUGPATHSPLINE) {
            cout << "mapXYPoints[0]:\n" << printVector(mapXYPoints[0], 10) << std::endl;
            cout << "mapXYPoints[1]:\n" << printVector(mapXYPoints[1], 10) << std::endl;
            cout << "next_x_vals.size:" << next_x_vals.size() << ", next_y_vals.size:" << next_y_vals.size() << std::endl;
          }
#ifdef xx
          //double dist_inc = 0.5;// meters, speed is set by
          double dist_inc = PATHDELTADISTANCEMAX;
          
          double same_d=laneCenter(currentLane);
          cout << "currentLane:" << currentLane << ", same_d:" << same_d << std::endl;
          for(int i = 0; i < LENGTHOFCARPATHTRAJECTORY; i++)
          {
            double next_s=car_s+(dist_inc*(i+1));
            vector<double> next_xy=getXY(next_s, same_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            //cout << "currentLane:" << currentLane << ", same_d:" << same_d << ", same_d:" << same_d << ", next_xy:" << printVector(next_xy) << std::endl;
            next_x_vals.push_back(next_xy[0]);
            next_y_vals.push_back(next_xy[1]);
            //next_x_vals.push_back(car_x/*current car x position*/+(dist_inc*i)*cos(deg2rad(car_yaw)));
            //next_y_vals.push_back(car_y/*current car y position*/+(dist_inc*i)*sin(deg2rad(car_yaw)));
          }
#endif
          if (DEBUGVEHICLEPATH) {
            cout << std::setprecision(6) << "car_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << ", car_speed:" << car_speed << std::endl;
            cout << "next_x_vals:" << next_x_vals[0] << "," << next_x_vals[1] << "..." <<  next_x_vals[next_x_vals.size()-2] << ", "
            <<  next_x_vals[next_x_vals.size()-1]
            << "\nnext_y_vals:" << next_y_vals[0] << "," << next_y_vals[1] << "..." <<  next_y_vals[next_y_vals.size()-2]  <<  ", "
            << next_y_vals[next_y_vals.size()-1]  << std::endl;
          }
          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;
          
          auto msg = "42[\"control\","+ msgJson.dump()+"]";
          
          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });
  
  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });
  
  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });
  
  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });
  
  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
#ifdef PYTHONPLOTTING
  //cout << "glfwGetVersionString: " << glfwGetVersionString() << std::endl;
  plt::plot({1,3,2,4});
  plt::show();
#endif
#ifdef QWTPLOTTING
  QwtPlot plot(QwtText("QwtExample1"));
  plot.setGeometry(0,0,640,400);
  plot.setAxisScale(QwtPlot::xBottom, 0.0,2.0 * M_PI);
  plot.setAxisScale(QwtPlot::yLeft,-1.0,1.0);
  
  QwtPlotCurve curve("Sine");
  std::vector<double> xs;
  std::vector<double> ys;
  for (double x = 0; x < 2.0 * M_PI; x+=(M_PI / 100.0))
  {
    xs.push_back(x);
    ys.push_back(std::sin(x) * std::cos(x));
  }
  QwtPointArrayData * const data = new QwtPointArrayData(&xs[0],&ys[0],xs.size());
  curve.setData(data);
  curve.attach(&plot);
  
  plot.show();
#endif
  h.run();
}
