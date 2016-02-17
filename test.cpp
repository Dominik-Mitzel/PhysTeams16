#include <cstdlib>
#include "Generator.h"

using namespace std;



//================================================================================//
//                                                                                //
//   main: test Generator class                                                   //
//                                                                                //
//================================================================================//

int main(int argc, char* argv[])
{
  const int NEvents = 1000;

  Generator* gen = new Generator(NEvents,100.);
  
  cout << endl;
  cout << "------------------------------------------------------------" << endl;

  for (int i = 0; i < NEvents; i++) {
    gen->NewEvent( (double) rand() / ((double) RAND_MAX + 1.),
    		   (double) rand() / ((double) RAND_MAX + 1.),
    		   (double) rand() / ((double) RAND_MAX + 1.),
		   (double) rand() / ((double) RAND_MAX + 1.));
  }

  cout << "------------------------------------------------------------" << endl;
  cout << " Total xsec = " << gen->GetXSec() << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << endl;
  
  return 0;
}
