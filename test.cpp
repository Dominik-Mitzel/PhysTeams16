#include <cstdlib>
#include "Generator.h"
#include <fstream>

using namespace std;



//================================================================================//
//                                                                                //
//   main: test Generator class                                                   //
//                                                                                //
//================================================================================//

int main(int argc, char* argv[])
{

  
  // for testing, let's just generate events with the C++ random number generator -- the inputFile does not contain enough events
  int nEvents = 1000000;
  
  Generator* gen = new Generator(nEvents);

  cout << "Starting event generation!" << endl;
  
  for (int i = 0; i < nEvents; i++) {
    gen->NewEvent( (double) rand() / ((double) RAND_MAX + 1.),
     		   (double) rand() / ((double) RAND_MAX + 1.),
     		   (double) rand() / ((double) RAND_MAX + 1.),
 		   (double) rand() / ((double) RAND_MAX + 1.));
    
    if (i % ((int) nEvents / 10) == 1) cout << "  " << i - 1 << " phase-space points done..." << endl;
  }
  

  // // alternative: use input file
  // std::ifstream infile("inputFile.dat");
  // int nEvents = 0;

  // double in1, in2, in3, in4;
  // while (infile >> in1 >> in2 >> in3 >> in4)
  // {
  //    nEvents++;
  // }

  // infile.clear();
  // infile.seekg(0, ios::beg);
  
  // cout << "nEvents = " << nEvents << "\n";

  // Generator* gen = new Generator(nEvents);
  // while (infile >> in1 >> in2 >> in3 >> in4)
  // {
  //    cout << "in1 = " << in1 << ", in2 = " << in2 << ", in3 = " << in3 << ", in4 = " << in4 << "\n";
  //    gen->NewEvent(in1, in2, in3, in4);
  // }

  cout << endl;
  cout << "--------------------------------------------------------------------------------" << endl;
  cout << "  Generated " << gen->GetNEvents() << " events passing the sqrt(s) cut." << endl;
  cout << "  Total xsec = " << gen->GetXSec() << " pb." << endl;
  cout << "--------------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << endl;
  cout << endl;
  
  return 0;
}
