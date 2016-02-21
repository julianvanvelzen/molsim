#include "system.h"

void Initialize (){
  int i, j;
  Vector momentum, temp;

  for (i=0;i<NUMBER_OF_PARTICLES;i++)
  {
    (particlelist + i)->position.x = RandomNumber() * GRIDSIZE;
    (particlelist + i)->position.y = RandomNumber() * GRIDSIZE;

    (particlelist + i)->velocity = VectorScalar(RanUnit(), RandomVelocity(TEMPERATURE));// 0;
    //(particlelist + i)->velocity.y = RandomVelocity(TEMPERATURE);// 0;
    
    momentum = VectorAddition(momentum, (particlelist + i)->velocity);

    AssignCellnumber(i);
  }
  
  momentum = VectorScalar(momentum, 1.0/NUMBER_OF_PARTICLES);

  for(i=0; i<NUMBER_OF_PARTICLES; i++){
   (particlelist + i)->velocity = VectorAddition((particlelist +i)->velocity, VectorScalar(momentum, -1));
   temp = VectorAddition(temp, (particlelist + i)->velocity);
  }

  printf("momentum: %lf %lf\n", temp.x, temp.y);

  // one double loop over all particles is necessary in initialisation 
  // to determine the forces during the first timestep
  for (i = 0; i < NUMBER_OF_PARTICLES-1; i++)
  {
    for(j = i + 1; j < NUMBER_OF_PARTICLES; j++)
    {
      ForceEnergy((particlelist + i), (particlelist + j));
    }
  }
  for (i = 0; i < NUMBER_OF_PARTICLES; i++){
    (particlelist + i)->force[0] = (particlelist + i)->force[1];  
    (particlelist + i)->force[1].y = 0.0;
    (particlelist + i)->force[1].x = 0.0;
    for (j = 0; j < NUMBER_OF_BINS; j++){
      (particlelist + i)->radial_distribution[j] = 0;
    }
  }
}