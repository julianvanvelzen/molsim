// #include "../system.h"
#include <stdio.h>
#define GRIDSIZE 3

int main(int argc, char** argv) {


  double x,y;
  int cell;
  x = 2.2;
  y = 1.5;

  printf("x %d y %d cell %d\n", (int)x  % GRIDSIZE, (int)y * GRIDSIZE, (int)x  % GRIDSIZE + (int)y * GRIDSIZE );


}