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

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
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

const double MAXIMUMMILESPERHOUR=49.5;// miles per hour
const double TARGETMILESPERHOUR=MAXIMUMMILESPERHOUR;// miles per hour
const int CARMOVESPERSECOND=50;// car moves at a well known fixed rate of 50 times per second
const int LENGTHOFCARPATHTRAJECTORY=CARMOVESPERSECOND*1/*second*/;
const double PATHLOOKAHEADXDISTANCE=30.;// meters

const double MAXIMUMMETERSPERSECOND=milesPerHourToMetersPerSecond(MAXIMUMMILESPERHOUR);
const double TARGETMETERSPERSECOND=milesPerHourToMetersPerSecond(TARGETMILESPERHOUR);
const double PATHPOINTDELTAT=1./double(LENGTHOFCARPATHTRAJECTORY);
const double PATHDELTADISTANCEMAX=PATHPOINTDELTAT*MAXIMUMMETERSPERSECOND;

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

const double yawFromXY(vector<vector<double>> theXYs) {
  vector<double> xPoints=theXYs[0];
  vector<double> yPoints=theXYs[1];
  double yaw = atan2(yPoints[1]-yPoints[0], xPoints[1]-xPoints[0]);
  return yaw;
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
    << "\carXYs:\n" << printVector(carXYs) << std::endl;
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
  if (DEBUGTRANSFORMS || true)
    cout << std::setprecision(10) << "rtClassTransform-refX:" << refX << ", refY:" << refY
    << ", refYaw:" << refYaw << ", cos(refYaw):" << cos(refYaw) << ", sin(refYaw):" << sin(refYaw) <<  std::endl;
  
  const double xPoint=theX*cos(refYaw)-theY*sin(refYaw);
  const double yPoint=theX*sin(refYaw)+theY*cos(refYaw);
  
  vector<double> *xy = new vector<double>;
  (*xy).push_back(xPoint+refX);
  (*xy).push_back(yPoint+refY);
  
  if (DEBUGTRANSFORMS || true)
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

int main(int argc, char* argv[]) {
  
  testXY();
  testRT();
  testTR();
  testCarMap();

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
  
  std::cout << "MAXIMUMMILESPERHOUR:" << MAXIMUMMILESPERHOUR << ", MAXIMUMMETERSPERSECOND:" << MAXIMUMMETERSPERSECOND <<std::endl;
  std::cout << "LENGTHOFCARPATHTRAJECTORY:" << LENGTHOFCARPATHTRAJECTORY << ", PATHPOINTDELTAT:" << PATHPOINTDELTAT << ", PATHDELTADISTANCEMAX:" << PATHDELTADISTANCEMAX << std::endl;
  
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
          double car_speed = j[1]["speed"];
          
          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];
          
          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
          
          json msgJson;
          
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          
          // (re)use what is left in previous path
          for (int i=0; i<previous_path_x.size(); i++) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }
          cout << "\n\nrecovered:" << next_x_vals.size() << " values from previous_path" << std::endl;
          if (previous_path_x.size() > 2) {
            cout << "next_x_vals:" << next_x_vals[0] << "," << next_x_vals[1] << "..." <<  next_x_vals[next_x_vals.size()-2] << ", "
            <<  next_x_vals[next_x_vals.size()-1]
            << "\nnext_y_vals:" << next_y_vals[0] << "," << next_y_vals[1] << "..." <<  next_y_vals[next_y_vals.size()-2]  <<  ", "
            << next_y_vals[next_y_vals.size()-1]  << std::endl;
          }

          /*
           Start by simply trying to move the car forward in a straight line at a constant 50 MPH velocity. Use the car's (x, y) localization information and its heading direction to create a simple, straight path that is drawn directly in front of the car.
           
           In main.cpp, instead of setting the speed directly, we pass next_x_vals, and next_y_vals to the simulator. We will set the points 0.5 m apart. Since the car moves 50 times a second, a distance of 0.5m per move will create a velocity of 25 m/s. 25 m/s is close to 50 MPH.
           
           (.5 meters/move)(50 moves/second)=25 meters/second
           (25 m/sec)(3.2804 ft/m)(1/5280 mi/ft)(60 sec/min)(60 min/hr)=55.92 mi/hr
           */
          const vector<vector<double>>& lastXyPathPoints =
          previous_path_x.size()<2?last2PathPoints(car_x, car_y, deg2rad(car_yaw)):last2PathPoints(previous_path_x, previous_path_y);
          cout << std::setprecision(6) << "car_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << std::endl;
          vector<double> xMapPoints=lastXyPathPoints[0];
          cout << "xMapPoints:" << printVector(xMapPoints, 10) << std::endl;
          vector<double> yMapPoints=lastXyPathPoints[1];
          cout << "yMapPoints:" << printVector(yMapPoints, 10) << std::endl;
          const double pathInterpolatedYaw=yawFromXY(lastXyPathPoints);
          cout << std::setprecision(6) << "pathInterpolatedYaw:" << pathInterpolatedYaw << ", car_yaw:" << deg2rad(car_yaw) << std::endl;
          
          // add additional waypoints in front of car in map coordinates, the spline needs more than 2 points
          const int currentLane=determineLane(car_d);
          for (int p=0; p<3; p++) {
            double lookAheadS = car_s+(p+1)*PATHLOOKAHEADXDISTANCE;
            vector<double> mapXY = getXY(lookAheadS, laneCenter(currentLane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            xMapPoints.push_back(mapXY[0]);
            yMapPoints.push_back(mapXY[1]);
          }
          // xPoints[1], yPoints[1] is last projected location for the car (usually) using the previous path
          const vector<vector<double>>& classTransformedXY=rtClassTransform(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xMapPoints, yMapPoints);
          const vector<double> xTransformedPoints=classTransformedXY[0];
          cout << "class-xTransformedPoints:" << printVector(xTransformedPoints, 10) << std::endl;
          const vector<double> yTransformedPoints=classTransformedXY[1];
          cout << "class-yTransformedPoints:" << printVector(yTransformedPoints, 10) << std::endl;
          
          cout << "affine-xPoints:" << printVector(xMapPoints, 10) << std::endl;
          cout << "affine-yPoints:" << printVector(yMapPoints, 10) << std::endl;
          cout << std::setprecision(6) << "pathInterpolatedYaw:" << pathInterpolatedYaw << ", car_yaw:" << deg2rad(car_yaw) << std::endl;
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
          cout << "carXY[0]:\n" << printVector(carXY[0], 10) << std::endl;
          cout << "carXY[1]:\n" << printVector(carXY[1], 10) << std::endl;

          tk::spline s;
          s.set_points(carXY[0], carXY[1]);
          const double yAtLookAheadDistance=s(PATHLOOKAHEADXDISTANCE);
          const double lengthOfLinearFitToSpline=sqrt(PATHLOOKAHEADXDISTANCE*PATHLOOKAHEADXDISTANCE+yAtLookAheadDistance*yAtLookAheadDistance);
          const int numberOfSamples=lengthOfLinearFitToSpline/(PATHPOINTDELTAT*TARGETMETERSPERSECOND);// meters/(seconds*metes/second)
          const int classN = lengthOfLinearFitToSpline/(.02*49.5/2.24);
          cout << "numberOfSamples:" << numberOfSamples << ", classN:" << classN << ", next_x_vals.size:" << next_x_vals.size() << std::endl;
          
          const double deltaX=PATHLOOKAHEADXDISTANCE/double(numberOfSamples);
          // spline points in car coordinates that will fill in next_?_vals after those copied over from previous_path
          int addtionalCarPoints=LENGTHOFCARPATHTRAJECTORY-next_x_vals.size();
          vector<double> xAdditionalCarPoints;
          vector<double> yAdditionalCarPoints;
          for (int additionalPoint=0; additionalPoint<addtionalCarPoints; additionalPoint++) {
            // additionalPoint+1: so deltaX is never 0, so next point added will not match last point in previous path
            const double carPointX=(additionalPoint+1)*deltaX;
            xAdditionalCarPoints.push_back(carPointX);
            const double carPointY=s(carPointX);
            yAdditionalCarPoints.push_back(carPointY);
          }
          
          // transform the additional car xy point back to the map
          const vector<vector<double>>& classMapPointsXY=trClassTransform(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xAdditionalCarPoints, yAdditionalCarPoints);
          cout << "class-classMapPointsXY:" << printVector(classMapPointsXY, 10) << std::endl;
          const vector<double> xClassMapPoints=classMapPointsXY[0];
          cout << "class-xClassMapPoints:" << printVector(xClassMapPoints, 10) << std::endl;
          const vector<double> yClassMapPoints=classMapPointsXY[1];
          cout << "class-yClassMapPoints:" << printVector(yClassMapPoints, 10) << std::endl;
          
          const vector<Eigen::Vector2d> additionalMapXYPoints=transformFromCarToMap(xMapPoints[1], yMapPoints[1], deg2rad(car_yaw), xAdditionalCarPoints, yAdditionalCarPoints);
          //cout << "affine-additionalMapXYPoints:\n" << printVector(additionalMapXYPoints, 10) << std::endl;
          
          // add new xy points to previous path from car
          vector<vector<double>> mapXYPoints=unpackV2D(additionalMapXYPoints);
          cout << "mapXYPoints[0]:\n" << printVector(mapXYPoints[0], 10) << std::endl;
          cout << "mapXYPoints[1]:\n" << printVector(mapXYPoints[1], 10) << std::endl;
          for (int additionalPoint=0; additionalPoint<addtionalCarPoints; additionalPoint++) {
            next_x_vals.push_back(additionalMapXYPoints[additionalPoint][0]);
            next_y_vals.push_back(additionalMapXYPoints[additionalPoint][1]);
          }
          cout << "next_x_vals.size:" << next_x_vals.size() << ", next_y_vals.size:" << next_y_vals.size() << std::endl;

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
          cout << std::setprecision(6) << "car_x:" << car_x << ", car_y:" << car_y << ", car_s:" << car_s << ", car_d:" << car_d << ", car_yaw:"<< deg2rad(car_yaw) << std::endl;
          cout << "next_x_vals:" << next_x_vals[0] << "," << next_x_vals[1] << "..." <<  next_x_vals[next_x_vals.size()-2] << ", "
          <<  next_x_vals[next_x_vals.size()-1]
          << "\nnext_y_vals:" << next_y_vals[0] << "," << next_y_vals[1] << "..." <<  next_y_vals[next_y_vals.size()-2]  <<  ", "
          << next_y_vals[next_y_vals.size()-1]  << std::endl;
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
