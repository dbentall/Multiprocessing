#include <mpi.h>
using namespace std;

int main(int argc, char* argv[]){
  cout << "s\n";
  int rank;
  MPI::Init(argc, argv);
  
  rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  cout << "rank: " << rank << " size: " << size << endl;
  
  MPI::Finalize();
  cout << "e\n";
  return 0;
}
