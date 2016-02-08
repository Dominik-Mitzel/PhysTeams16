#include "Generator.h"



//================================================================================//
//                                                                                //
// Constructor: set number of events, sqrt(s) cut, and filename for event file    //
//                                                                                //
//================================================================================//

Generator::Generator(int NEventsIn, double sqrtSMinIn, string filenameIn)
{

  // Settings: total number of events, sqrt(S) cut, filename
  if (NEventsIn > 0) {
    NEvents = NEventsIn;
  } else {
    cout << "Warning: Number of events given as " << NEventsIn << ", use default 1 instead!";
    NEvents = 1;
  }

  if (sqrtSMinIn > 0.)
    sqrtSMin = sqrtSMinIn;
  else {
    cout << "Warning: sqrtSMin given as " << sqrtSMinIn << ", use default 10 instead!";
    sqrtSMin = 10.;
  }

  filename = filenameIn;

  // Initialize other variables
  eventCounter = 0;
  xsec = 0.;
}



//================================================================================//
//                                                                                //
// Destructor                                                                     //
//                                                                                //
//================================================================================//

Generator::~Generator()
{
  // Not sure if anything has to be done at all?    
}






//================================================================================//
//                                                                                //
// NewEvent: given four random numbers in [0,1], calculates the corresponding     //
//           phase-space point, the Jacobian, matrix element squared and pdfs.    //
//           The resulting weight is saved together with the momenta in the       //
//           event file. Returns true if succesful (hopefully always!).           //
//                                                                                //
//================================================================================//

bool Generator::NewEvent(double z1, double z2, double z3, double z4)
{
  // Momenta in the proton-proton center-of-mass frame, all in GeV.
  // p1 and p2 are incoming momenta of u and ubar.
  // p3 and p4 are outgoing momenta of d and dbar.
  momentum p1, p2, p3, p4;

  // azimuthal angle of outgoing d quark
  double phi3;

  // Mandelstam variables (in parton-level process)
  double s, t, u;

  // Momentum fractions carried by initial quarks
  double x1, x2;
  double protonE = 6500.;
  double protonS = pow(2.*protonE,2);

  // Check input variables are in (0,1) range
  if (z1 < 0 || z2 < 0 || z3 < 0 || z4 < 0 || z1 > 1 || z2 > 1 || z3 > 1 || z4 > 1) {
    cout << "Error: input numbers (" << z1 << ", " << z2 << ", " << z3 << ", " << z4
	 << ") not all in [0,1] range! Ignoring event!" << endl;
    return false;
  }
  
  // Transform z1 to z4 into kinematics
  // for now, just use a very-not-optimal scheme:
  // z1 and z2 directly give x1 and x2, thus setting s
  // z3 sets t, normalized such that z3 = 0 -> t = -s and z3 = 1 -> t = 0
  // z4 sets phi of p3, normalised such that z4 = 1 -> phi = 2pi
  x1 = z1;
  x2 = z2;

  // avoid division by zero later
  if (x1 < 0.000000001) x1 = 0.000000001;
  if (x2 < 0.000000001) x2 = 0.000000001;
  
  s = x1*x2*protonS;
  t = -s + z3*s;
  phi3 = 2.*M_PI*z4;

  // Debug output
  cout << endl;
  cout << "Input kinematics:" << endl;
  cout << "  x1 = " << x1 << endl;
  cout << "  x2 = " << x2 << endl;
  cout << "  s = " << s << endl;
  cout << "  t = " << t << endl;
  cout << "  phi3 = " << phi3 << endl;
  cout << endl;

  // Calculate four-momenta based on x1, x2, s, t, phi3 (to do)
  double pt = sqrt(-t*(t+s)) / (sqrt(s));
  
  p1 = (momentum){x1*protonE,
		  0.,
		  0.,
		  x1*protonE};
  p2 = (momentum){x2*protonE,
		  0.,
		  0.,
		  -x2*protonE};
  p3 = (momentum){(s*x1+t*x1-t*x2)*protonE/s,
		  pt*cos(phi3),
		  pt*sin(phi3),
		  (s*x1+t*x1+t*x2)*protonE/s,};
  p4 = (momentum){(s*x2-t*x1+t*x2)*protonE/s,
		  -pt*cos(phi3),
		  -pt*sin(phi3),
		  (-s*x2-t*x1-t*x2)*protonE/s};

  u = (p4-p1)*(p4-p1);

  // debug output
  cout << endl;
  cout << "Momentum conservation checks (should be zero):" << endl;
  cout << "  (p3+p4) - (p1+p2) = " << p3 + p4 - p1 - p2 << endl;
  cout << "Mass-shell checks (should be zero):" << endl;
  cout << "  p1^2 = " << p1*p1 << endl;
  cout << "  p2^2 = " << p2*p2 << endl;
  cout << "  p3^2 = " << p3*p3 << endl;
  cout << "  p4^2 = " << p4*p4 << endl;
  cout << "Mandelstam checks (should be zero):" << endl;
  cout << "  s - (p1+p2)^2 = " << s - (p1+p2)*(p1+p2) << endl;
  cout << "  s - (p3+p4)^2 = " << s - (p3+p4)*(p3+p4) << endl;
  cout << "  t - (p3-p1)^2 = " << t - (p3-p1)*(p3-p1) << endl;
  cout << "  t - (p4-p2)^2 = " << t - (p4-p2)*(p4-p2) << endl;
  cout << "  u - (p4-p1)^2 = " << u - (p4-p1)*(p4-p1) << endl;
  cout << "  u - (p3-p2)^2 = " << u - (p3-p2)*(p3-p2) << endl;
  cout << endl;

  // Calculate Jacobian (to do)
  double JacobiDeterminant = 1.;

  // Overall factor for hard matrix element (to do)
  double overallConstant = 1./s;

  // Calculate hard matrix elements (to do)
  double MSquared = CalculateMSquared(s, t, u);

  // Pdfs
  double q = sqrt(s);
  int pdgid1 = 1; // up quark
  int pdgid2 = -1; // anti-up quark
  double pdfFactor1 = GetPDFValue(pdgid1, x1, q);
  double pdfFactor2 = GetPDFValue(pdgid2, x2, q);

  // Plug everything together
  double weight = 1 / (double) NEvents * overallConstant * JacobiDeterminant * MSquared * pdfFactor1 * pdfFactor2;

  // Save result
  return SaveEvent(weight, p3, p4);
}



//================================================================================//
//                                                                                //
// CalculateMSquared: Returns matrix element squared for u ubar -> d dbar,        //
//                    for the phase-space point given by the Mandelstam variables //
//                    s, t, u                                                     //
//                                                                                //
//================================================================================//

double Generator::CalculateMSquared(double s, double t, double u)
{
  return 1./pow(t,2); // to do
}



//================================================================================//
//                                                                                //
// GetPDFValue: Returns the value of the pdf                                      //
//                                                                                //
//================================================================================//

double Generator::GetPDFValue(int pdgid, double x, double q)
{
  return 1.;
}



//================================================================================//
//                                                                                //
// SaveEvent: Adds a line to the output file with the weight and momenta given.   //
//            Returns true if succesful.                                          //
//                                                                                //
//================================================================================//

bool Generator::SaveEvent(double weight, momentum p3, momentum p4)
{
  // Add weight to total xsec so far and increment event counter
  xsec += weight;
  eventCounter++;
  
  // for now: write on screen
  cout << "Event " << eventCounter << ": p3 = (" << p3.E << ", " << p3.px << ", " << p3.py << ", " << p3.pz
       << "); p4 = (" << p4.E << ", " << p4.px << ", " << p4.py << ", " << p4.pz
       << ") -> weight = " << weight << endl;
  return true;
}




//================================================================================//
//                                                                                //
// GetXSec: Returns the sum of all event weights calculated so far                //
//          (corrected for missing events)                                        //
//                                                                                //
//================================================================================//

double Generator::GetXSec()
{
  if (eventCounter > 0)
    return xsec * (double) NEvents / (double) eventCounter;
  else // no events generated yet
    return -1.;
}
