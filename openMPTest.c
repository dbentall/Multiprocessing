#include <omp.h>

int main(int argc, char* argv[]){

  int thID;
  int i;
#pragma omp parallel private(thID)
{
  thID = omp_get_thread_num();
  for(i = 0; i < 10; i++){
    printf("i: %d, id: %d\n", i, thID);
  }
}
  return 0;
}
