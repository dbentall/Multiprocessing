#include <omp.h>
#include <iostream>
#include <sstream>
#include <time.h>

#include "wiener.h"

#define FILE_STR "images/hubble-20g-"
#define EXT ".png"
#define K 1.0
#define SIGMA 5.0

using namespace std;

void println(string out){
  cout << out;
}

int main(int argc, char* argv[]){
  time_t startTime = time(NULL);
  int numImages = 8;
  int fileNo = 0;
  #pragma omp parallel for shared(fileNo)
  for(int i = 0; i < numImages; i++){
    while(1){
      int privateNo;
      #pragma omp critical
      {
        privateNo = fileNo++;
      }
      if(privateNo >= numImages){
        break;
      }
      ostringstream inFile, outFile;
      inFile << FILE_STR << privateNo << EXT;
      outFile << FILE_STR << privateNo << "out-" << EXT;
      cout << inFile << endl << outFile << endl;
      wFilter(inFile.str().c_str(), outFile.str().c_str(), SIGMA, K);
      #pragma omp critical
      {
        cout << FILE_STR << privateNo << EXT << " thread: "
                << omp_get_thread_num() << endl;
      }
    }
  }
  return 0;
}
