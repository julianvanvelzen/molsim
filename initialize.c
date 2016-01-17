#include "system.h"

void Initialize (){

  int i;
  int c;
  double x;
  double y;

  for (i=0;i<NUMBER_OF_PARTICLES;i++){
    c = RandomNumber() * NUMBER_OF_PROCESSORS;
    y = RandomNumber() + c / GRIDSIZE;
    x = RandomNumber() + c % GRIDSIZE;

    (particlelist + i)->position.x = x;
    (particlelist + i)->position.y = y;

    (particlelist + i)->velocity.x = RandomVelocity(TEMPERATURE);
    (particlelist + i)->velocity.y = RandomVelocity(TEMPERATURE);
    
    (particlelist + i)->cellnumber = c ;
    (particlelist + i)->force.x = 0.0;
    (particlelist + i)->force.y = 0.0;
  }

}