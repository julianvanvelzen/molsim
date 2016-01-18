#include "system.h"

void Initialize (){

  int i;

  for (i=0;i<NUMBER_OF_PARTICLES;i++){

    (particlelist + i)->position.x = 0.95 + i/30.0;
    (particlelist + i)->position.y = 0.95 + i/30.0;

    // (particlelist + i)->position.x = RandomNumber() * GRIDSIZE;
    // (particlelist + i)->position.y = RandomNumber() * GRIDSIZE;

    (particlelist + i)->velocity.x = /*RandomVelocity(TEMPERATURE)*/ 0;
    (particlelist + i)->velocity.y = /*RandomVelocity(TEMPERATURE)*/ 0;
    
    AssignCellnumber(i);

    (particlelist + i)->force.x = 0.0;
    (particlelist + i)->force.y = 0.0;
  }

}