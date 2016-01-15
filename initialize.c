#include "system.h"

void Initialize (){

  int i,j;
  // Particle particle;

  // particle.position.x = RandomNumber();
  // particle.position.y = RandomNumber();
  // particle.velocity.x = RandomVelocity(Temperature);
  // particle.velocity.y = RandomVelocity(Temperature);


  for (i=0;i<NUMBER_OF_PARTICLES;i++){
    (particlelist + i)->position.x = 0.5 + i%GRIDSIZE;
    (particlelist + i)->position.y = 0.5 + i/GRIDSIZE;
    (particlelist + i)->cellnumber = 1;
  }

}