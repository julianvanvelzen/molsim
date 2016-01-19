#include "system.h"

void Initialize (){
  int i, j;
  double dE;
  Vector forceVector;

  for (i=0;i<NUMBER_OF_PARTICLES;i++)
  {
    (particlelist + i)->position.x = RandomNumber() * GRIDSIZE;
    (particlelist + i)->position.y = RandomNumber() * GRIDSIZE;

    (particlelist + i)->velocity.x = /*RandomVelocity(TEMPERATURE)*/ 0;
    (particlelist + i)->velocity.y = /*RandomVelocity(TEMPERATURE)*/ 0;
    
    AssignCellnumber(i);
  }

  // a single double loop is necessary in initialisation to determine the forces during the first timestep
  for (i = 0; i < NUMBER_OF_PARTICLES-1; i++)
  {
    for(j = i + 1; j < NUMBER_OF_PARTICLES; j++)
    {
      ForceEnergy((particlelist + i)->position, (particlelist + j)->position, &forceVector, &dE);
      (particlelist + i)->force[0] = VectorAddition((particlelist + i)->force[0], forceVector);     
      (particlelist + j)->force[0] = VectorAddition((particlelist + j)->force[0], VectorMultiplication(forceVector, -1));
    }
    (particlelist + i)->force[1].x = 0.0;
    (particlelist + i)->force[1].y = 0.0;
  }
}