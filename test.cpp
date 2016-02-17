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
  std::ifstream infile("inputFile.dat");
  int nEvents = 0;

  double in1, in2, in3, in4;
  while (infile >> in1 >> in2 >> in3 >> in4)
  {
     nEvents++;
  }

  infile.clear();
  infile.seekg(0, ios::beg);
  
  cout << "nEvents = " << nEvents << "\n";

  Generator* gen = new Generator(nEvents);

  while (infile >> in1 >> in2 >> in3 >> in4)
  {
     cout << "in1 = " << in1 << ", in2 = " << in2 << ", in3 = " << in3 << ", in4 = " << in4 << "\n";
     gen->NewEvent(in1, in2, in3, in4);
     
  }

  cout << "------------------------------------------------------------" << endl;
  cout << " Total xsec = " << gen->GetXSec() << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << endl;
  
  return 0;
}
