#include "system.h"

void Initialize (){

  int i,j;
  // Particle particle;

  // particle.position.x = RandomNumber();
  // particle.position.y = RandomNumber();
  // particle.velocity.x = RandomVelocity(Temperature);
  // particle.velocity.y = RandomVelocity(Temperature);


  for (i=0;i<NUMBER_OF_PARTICLES;i++){
    (particlelist + i)->position.x = RandomNumber()-0.5 + (i/5)%GRIDSIZE;
    (particlelist + i)->position.y = RandomNumber()-0.5 + (i/5)/GRIDSIZE;
    (particlelist + i)->cellnumber = i/5;
    (particlelist + i)->force.x = 0.0;
    (particlelist + i)->force.y = 0.0;
  }

}