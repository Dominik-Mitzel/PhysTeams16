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
  momentum p1 = {1,0,0,1};
  momentum p2 = {1,0,0,-1};
  momentum p3 = {1,0,1,0};
  momentum p4 = {1,0,-1,0};

  // Transform z1 to z4 into four-momenta, taking into account the sqrtSMin cut
  // and the fact that the quarks should not have more energy than 6500 GeV
  // ...

  // Mandelstam variables
  double s = (p1+p2)*(p1+p2);
  double t = (p3-p1)*(p3-p1);
  double u = (p4-p1)*(p4-p1);

  // Calculate Jacobian
  double JacobiDeterminant = 1.;

  // Overall factor for hard matrix element
  double overallConstant = 1.;

  // Calculate hard matrix elements
  double MSquared = CalculateMSquared(s, t, u);

  // Momentum fractions carried by quark 1 and 2
  double protonE = 6500.; // in GeV
  double x1 = p1.E / protonE;
  double x2 = p2.E / protonE;
  // Energy scale of hard process, important for  
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
  return 1.;
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
