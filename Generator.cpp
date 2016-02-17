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
  filenameKin = "debug.dat";

  // Initialize other variables
  eventCounter = 0;
  xsec = 0.;

  // delete old output files
  if( remove( filename.c_str() ) != 0 )
    cout << "Error deleting file old output File.  Maybe it didn't exist?\n";
  else
    cout << "Old output file successfully deleted \n";
  if( remove( filenameKin.c_str() ) != 0 )
    cout << "Error deleting file old output File.  Maybe it didn't exist?\n";
  else
    cout << "Old output file successfully deleted \n";
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
  // debug output
  cerr << endl;
  cerr << "------------------------------------------------------------" << endl;

  
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
  if (x1 < 0.000000001) x1 = 0.000000001;
  x2 = z2;
  if (x2 < 0.000000001) x2 = 0.000000001;
  
  s = x1*x2*protonS;

  // Check if cut is satisfied
  if (s < sqrtSMin*sqrtSMin) {
    return SaveEvent(0., p3, p4);
  }

  
  t = -s + z3*s;
  phi3 = 2.*M_PI*z4;

  // Debug output
  cerr << endl;
  cerr << "Input kinematics:" << endl;
  cerr << "  x1 = " << x1 << endl;
  cerr << "  x2 = " << x2 << endl;
  cerr << "  s = " << s << endl;
  cerr << "  t = " << t << endl;
  cerr << "  phi3 = " << phi3 << endl;
  cerr << endl;

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
  cerr << endl;
  cerr << "Four-momenta:" << endl;
  cerr << "  p1 = " << p1 << endl;
  cerr << "  p2 = " << p2 << endl;
  cerr << "  p3 = " << p3 << endl;
  cerr << "  p4 = " << p4 << endl;
  cerr << "Mandelstams:" << endl;
  cerr << "  s = " << s << endl;
  cerr << "  t = " << t << endl;
  cerr << "  u = " << u << endl;

  // debug output
  cerr << endl;
  cerr << "Momentum conservation checks (should be zero):" << endl;
  cerr << "  (p3+p4) - (p1+p2) = " << p3 + p4 - p1 - p2 << endl;
  cerr << "Mass-shell checks (should be zero):" << endl;
  cerr << "  p1^2 = " << p1*p1 << endl;
  cerr << "  p2^2 = " << p2*p2 << endl;
  cerr << "  p3^2 = " << p3*p3 << endl;
  cerr << "  p4^2 = " << p4*p4 << endl;
  cerr << "Mandelstam checks (should be zero):" << endl;
  cerr << "  s - (p1+p2)^2 = " << s - (p1+p2)*(p1+p2) << endl;
  cerr << "  s - (p3+p4)^2 = " << s - (p3+p4)*(p3+p4) << endl;
  cerr << "  t - (p3-p1)^2 = " << t - (p3-p1)*(p3-p1) << endl;
  cerr << "  t - (p4-p2)^2 = " << t - (p4-p2)*(p4-p2) << endl;
  cerr << "  u - (p4-p1)^2 = " << u - (p4-p1)*(p4-p1) << endl;
  cerr << "  u - (p3-p2)^2 = " << u - (p3-p2)*(p3-p2) << endl;
  
  // Calculate Jacobian
  double JacobiDeterminant = s;

  // Calculate hard matrix elements
  double dsigmadt = CalculateDSigmaDt(s, t, u);

  // Pdfs
  double q = sqrt(s);
  int pdgid1 = 1; // up quark
  int pdgid2 = -1; // anti-up quark
  double pdfFactor1 = GetPDFValue(pdgid1, x1, q);
  double pdfFactor2 = GetPDFValue(pdgid2, x2, q);

  // Conversion from GeV^{-2}  to pb
  double conversionFactor = 389370000.;

  // two assignments of quarks to protons
  double multiplicityFactor = 2.;

  // Plug everything together
  double weight = multiplicityFactor * conversionFactor / (double) NEvents * JacobiDeterminant * dsigmadt * pdfFactor1 * pdfFactor2;

  // Save data to be plotted
  SaveKinematics(weight, z1, z2, z3, z4, sqrt(s), sqrt(-t), sqrt(p3.px*p3.px + p3.py * p3.py));

  // Save result
  return SaveEvent(weight, p3, p4);
}



//================================================================================//
//                                                                                //
// CalculateDSigmaDt: Returns d sigma / d t squared for u ubar -> d dbar,         //
//                    for the phase-space point given by the Mandelstam variables //
//                    s, t, u                                                     //
//                                                                                //
//================================================================================//

double Generator::CalculateDSigmaDt(double s, double t, double u)
{
  const double alphaS = 0.1;
  if (t != 0. && s != 0. && u != 0.)
    return 4.*M_PI*alphaS*alphaS/(9.*s*s) * (t*t + u*u)/(s*s);
  return 0.;
      
}



//================================================================================//
//                                                                                //
// GetPDFValue: Returns the value of the pdf, following some weird paper          //
//                                                                                //
//================================================================================//

double Generator::GetPDFValue(int pdgid, double x, double q)
{
  double a = 0;
  double b = 0;
  double n = 1;

  if (x < 0.05) {
    switch (pdgid) {
    case 1:
      a = 2.85;
      b = 1.02;
      n = 1.24;
    case -1:
      a = 6.26;
      b = 1.50;
      n = 19.46;
    }
  } else {
    switch (pdgid) {
    case 1:
      a = 4.11;
      b = 0.52;
      n = 0.38;
    case -1:
      a = 6.69;
      b = 1.48;
      n = 17.0;
    }
  }

  
  return 1./n * pow(x,-b) * pow((1-x),a);
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

  if (weight > 0.) {
    ofstream outFile;
    outFile.open(filename, ios::app);
    
    outFile << weight << " " << p3.E << " " << p3.px << " " << p3.py << " " << p3.pz << " " << p4.E << " " << p4.px << " " << p4.py << " " << p4.pz << "\n";
    outFile.close();
    
    return true;
  } else {
    return true;
  }
}



//================================================================================//
//                                                                                //
// SaveEvent: Adds a line to the output file with the weight and momenta given.   //
//            Returns true if succesful.                                          //
//                                                                                //
//================================================================================//

bool Generator::SaveKinematics(double weight, double y1, double y2, double y3, double y4, double y5, double y6, double y7)
{
  string s = " ";

  if (weight > 0.) {
    
    ofstream outFile;
    outFile.open(filenameKin, ios::app);
    
    outFile << weight << "\t" << y1
	 << "\t"  << y2
	 << "\t"  << y3
	 << "\t"  << y4
	 << "\t"  << y5
	 << "\t"  << y6
	 << "\t"  << y7 << "\n";
    outFile.close();
  } 
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
