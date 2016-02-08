#include <cstdlib>
#include <math.h>
#include <iostream>
#include <string>

using namespace std;

class Generator
{
 public:
  // new type for four-momenta, and define +, -, * on them. The latter contracts indices using the (+---) metric.
  struct momentum {double E, px, py, pz;};
  friend momentum operator+( const momentum& left, const momentum& right ) {
    momentum result = {left.E + right.E, left.px + right.px, left.py + right.py, left.pz + right.pz};
    return result;
  }
  friend momentum operator-( const momentum& left, const momentum& right ) {
    momentum result = {left.E - right.E, left.px - right.px, left.py - right.py, left.pz - right.pz};
    return result;
  }
  friend double operator*( const momentum& left, const momentum& right ) {
    return left.E * right.E - left.px * right.px - left.py * right.py - left.pz * right.pz;
  }
  friend std::ostream& operator<<(std::ostream &out, const momentum& right)
  {
    return out << "(" << right.E << ", " << right.px << ", " << right.py << ", " << right.pz << ")";
  }

  // Constructor and destructor
  Generator(int NEventsIn, double sqrtSMinIn = 10., string filenameIn = "events.dat");
  ~Generator();

  // Calculate new phase-space point described by four random numbers in [0,1]
  bool NewEvent(double z1, double z2, double z3, double z4);

  // Get total XSec (once all events are generated)
  double GetXSec();

  
 private:
  // Settings for event generation
  int NEvents;
  double sqrtSMin;
  string filename;
  double xsec;
  int eventCounter;

  // Internal functions
  double CalculateMSquared(double s, double t, double u);
  double GetPDFValue(int pdgid, double x, double q);
  bool SaveEvent(double weight, momentum p3, momentum p4);
};
