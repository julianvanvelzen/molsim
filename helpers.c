#include "system.h"

void CheckInputErrors(){
  int Error = 0;
  if(NUMBER_OF_PROCESSORS != SQR(GRIDSIZE)){
    printf("Number of processors must be equal to Gridsize^2\n");
    Error = 1;
  }
  if(TEMPERATURE <= 0){
    printf("Temperature has to be positive and non-zero\n");
    Error = 1;
  }
  if(NUMBER_OF_PARTICLES <= 0){
    printf("There must be at least one particle in the box\n");
    Error = 1;
  }
  if(NUMBER_OF_CYCLES < 300){
    printf("Recommended to simulate at least 300 cycles\n");
    Error = 1;
  }
  if(RCUT > 1 || RCUT <=0){
    printf("Rcut must be higher than 0 and may not exceed 1\n");
    Error = 1;
  }
  if(DELTAT <= 0){
    printf("DeltaT must be higher than 0\n");
    Error = 1;
  }
  if(INITIALISATION_STEPS < 0){
    printf("Number of initialisation steps may not be negative\n");
    Error = 1;
  }
  if(NUMBER_OF_CYCLES <= INITIALISATION_STEPS){
    printf("Simulation must be longer than number of initialisation steps (%d)\n", INITIALISATION_STEPS);
    Error = 1;
  }
  if(Error == 1){
    printf("\n");
    exit(0);
  }
}

void ForceEnergy(Particle *p1, Particle *p2, int pbc_x, int pbc_y){
  if(pbc_x != 0) p1->position.x -= GRIDSIZE;
  if(pbc_y != 0) p1->position.y -= GRIDSIZE;
  
  double distance = VectorDistance(p1->position, p2->position);
  if (distance >= RCUT) {
    if(pbc_x != 0) p1->position.x += GRIDSIZE;
    if(pbc_y != 0) p1->position.y += GRIDSIZE;
    return;
  }
  double force = (2.0*REPULSIVE_CST*sqrt(SQR(distance-RCUT))/SQR(RCUT));
  
  Vector relative_position;
  relative_position.x = (p1->position.x - p2->position.x);
  relative_position.y = (p1->position.y - p2->position.y);

  Vector forceVector;
  forceVector.x = relative_position.x * force/distance;
  forceVector.y = relative_position.y * force/distance;
  
  p1->force[1].x += forceVector.x;
  p1->force[1].y += forceVector.y;

  p2->force[1].x -= forceVector.x;
  p2->force[1].y -= forceVector.y;

  if (forceVector.x != - VectorScalar(forceVector, -1.0).x)
    printf("f1 != f2\n");

  // Potential and pressure are summed over all particles later,
  // there is no need to explicitly save p2->potential or p2->pressure
  p1->potential += REPULSIVE_CST*SQR(distance-RCUT)/SQR(RCUT);
  p1->pressure_contribution += 2*(forceVector.x*relative_position.x + forceVector.y*relative_position.y) / (2 * SQR(GRIDSIZE));
  
  int rdf_bin_index = NUMBER_OF_BINS*distance/RCUT;
  p1->radial_distribution[rdf_bin_index] += 2;

  if(pbc_x != 0) p1->position.x += GRIDSIZE;
  if(pbc_y != 0) p1->position.y += GRIDSIZE;
}

Vector VectorAddition(Vector v1, Vector v2){
  Vector v3;
  v3.x = v1.x + v2.x;
  v3.y = v1.y + v2.y;
  return v3;
}

Vector VectorScalar(Vector vector, double scalar){
  vector.x *= scalar;
  vector.y *= scalar;
  return vector;
}

double VectorDistance(Vector v1, Vector v2){
  double distance = sqrt(SQR(v1.x - v2.x) + SQR(v1.y - v2.y));
  return distance;
}

Vector RanUnit(void){
  double random = RandomNumber() * 2 * M_PI;
  Vector unit;
  unit.x = cos(random);
  unit.y = sin(random);
  return unit;
}

void getNearbyCoordinates(Cell *cells, int currentPosition){
  (cells+currentPosition)->neighbouringcells[0] = currentPosition + GRIDSIZE;
  (cells+currentPosition)->neighbouringcells[1] = currentPosition + GRIDSIZE + 1;
  (cells+currentPosition)->neighbouringcells[2] = currentPosition + 1;
  (cells+currentPosition)->neighbouringcells[3] = currentPosition - GRIDSIZE + 1;
  (cells+currentPosition)->neighbouringcells[4] = currentPosition - GRIDSIZE;
  (cells+currentPosition)->neighbouringcells[5] = currentPosition - GRIDSIZE - 1;
  (cells+currentPosition)->neighbouringcells[6] = currentPosition - 1;
  (cells+currentPosition)->neighbouringcells[7] = currentPosition + GRIDSIZE - 1;

  // top
  if (currentPosition+GRIDSIZE > SQR(GRIDSIZE)-1){
    (cells+currentPosition)->neighbouringcells[7] = (currentPosition + GRIDSIZE-1)%GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[0] = (currentPosition + GRIDSIZE)%GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[1] = (currentPosition + GRIDSIZE+1)%GRIDSIZE;
  }
  // right 
  if ((currentPosition+1)%GRIDSIZE == 0 ){
    (cells+currentPosition)->neighbouringcells[1] = (currentPosition + 1)%(SQR(GRIDSIZE));
    (cells+currentPosition)->neighbouringcells[2] -= GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[3] -= GRIDSIZE;
  }
  // bottom
  if (currentPosition-GRIDSIZE < 0){
    (cells+currentPosition)->neighbouringcells[3] += SQR(GRIDSIZE);
    (cells+currentPosition)->neighbouringcells[4] += SQR(GRIDSIZE);
    (cells+currentPosition)->neighbouringcells[5] += SQR(GRIDSIZE);
  }
  // left
  if(currentPosition%GRIDSIZE==0){
    (cells+currentPosition)->neighbouringcells[5] += GRIDSIZE;   
    (cells+currentPosition)->neighbouringcells[6] += GRIDSIZE; 
    (cells+currentPosition)->neighbouringcells[7] += GRIDSIZE;
  }
}

int cmpfunc (const void * a, const void * b){
  Particle *A = (Particle *)a;
  Particle *B = (Particle *)b;
  return ( A->cellnumber - B->cellnumber );
} 

void setindeces(Particle *particlelist, Cell *cells){
  int i = 0, j;
  int currentcell = 0;
  
  // finds which particles belong to which cells,
  // assumes the particles are sorted by cellnumber
  while (i<NUMBER_OF_PARTICLES) {
    while ((particlelist+i)->cellnumber > currentcell ) {
      (cells+currentcell)->start = i;
      (cells+currentcell)->end = i;
      currentcell++;
    }
  
    (cells+currentcell)->start = i;
    while((particlelist+i)->cellnumber == currentcell) i++;
    (cells+currentcell)->end = i;
    
    if(i >= NUMBER_OF_PARTICLES){
      (cells+currentcell)->end = NUMBER_OF_PARTICLES;
      for(j=currentcell+1;j<NUMBER_OF_PROCESSORS;j++){
        (cells+j)->start = NUMBER_OF_PARTICLES;
        (cells+j)->end   = NUMBER_OF_PARTICLES;
      }
    }
    currentcell++;
  }
}

void loopforces(Cell *cells, int world_rank){
  int i,k,l,m;
  int mstart, mend;
  int pbc_x = 0;
  int pbc_y = 0;

  // loop over all particles in your own cell
  int pstart = (cells + world_rank)->start;
  int pend = (cells + world_rank)->end;
  for (i = pstart; i < pend; i++){
    
    // loop over all other particles in your own cell
    for (k = i + 1; k < pend; k++)
      ForceEnergy((particlelist + i), (particlelist + k), 0, 0);
    
    // loop over all neighbouring cells
    for (l=0; l<4; l++) {
      // test for periodic boundary conditions        
      if ( (l >= 1) && (cells+world_rank)->neighbouringcells[l] % GRIDSIZE == 0)      pbc_x = 1;
      if ( (l == 0 || l == 1) && (cells+world_rank)->neighbouringcells[l] < GRIDSIZE) pbc_y = 1;

      // loop over all particles in neighbouring cells
      mstart = (cells + (cells+world_rank)->neighbouringcells[l])->start;
      mend   = (cells + (cells+world_rank)->neighbouringcells[l])->end;
      for (m  = mstart; m < mend; m++)
        ForceEnergy((particlelist + i), (particlelist + m), pbc_x, pbc_y);

      pbc_x = 0;
      pbc_y = 0;
    }
  }
}

void sum_apply_contributions(Cell *cells, Particle *gather, int cycle){
  int i,j,k, neighbour_offset, current_box_offset;
  double Ek = 0;
  double Ev = 0;
  double scale;
  double current_pressure = 0;
  double current_temperature = 0;

  // loop contributions of every particle
  for (i = 0; i < NUMBER_OF_PARTICLES; i++)
  {
    current_box_offset = (particlelist+i)->cellnumber;
    // Apply contributions of neighboring cells
    // Only the contributions of bottom and left neighbours need to be considered
    for (j = 4; j < 8; j++) {   
      neighbour_offset = (cells+current_box_offset)->neighbouringcells[j];
      if(neighbour_offset == 0) continue; // (particlelist + i) already includes the contribution of box 0
      (particlelist+i)->force[1].x += (gather + (NUMBER_OF_PARTICLES*neighbour_offset) + i)->force[1].x;
      (particlelist+i)->force[1].y += (gather + (NUMBER_OF_PARTICLES*neighbour_offset) + i)->force[1].y;
    }
    // Apply contributions of own cells
    // (particlelist + i) already includes the contribution of box 0
    if(current_box_offset != 0){    
      (particlelist +i)->force[1].x += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->force[1].x;
      (particlelist +i)->force[1].y += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->force[1].y;
      (particlelist +i)->pressure_contribution += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->pressure_contribution;
    }
  
    // Calculate stuff
    (particlelist + i)->velocity.x += (particlelist + i)->force[1].x*0.5*DELTAT;
    (particlelist + i)->velocity.y += (particlelist + i)->force[1].y*0.5*DELTAT;
    (particlelist + i)->force[0].x = (particlelist + i)->force[1].x;
    (particlelist + i)->force[0].y = (particlelist + i)->force[1].y;
    (particlelist + i)->force[1].x = 0;
    (particlelist + i)->force[1].y = 0;

    Ev += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->potential;
    Ek += (SQR((particlelist + i)->velocity.x) + SQR((particlelist + i)->velocity.y)); // divided by 2 later for efficiency
    current_pressure += (particlelist + i)->pressure_contribution;
    
    if(cycle > INITIALISATION_STEPS){
      for(k = 0; k < NUMBER_OF_BINS; k++){
        rdf_total[k] += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->radial_distribution[k]; 
        (particlelist + i)->radial_distribution[k] = 0;
      }
    }
  }
  
  Ek *= 0.5;

  // scale velocities to temperature
  if(cycle < INITIALISATION_STEPS){
    scale = sqrt(TEMPERATURE*(2.0*NUMBER_OF_PARTICLES-2.0)/(2.0*Ek));
    Ek = 0;
    if (cycle%50 ==0 && cycle > 51 ) printf("Temperature scale %lf", scale);
    for(i = 0; i < NUMBER_OF_PARTICLES; i++){
      (particlelist + i)->velocity = VectorScalar((particlelist + i)->velocity, scale);
      Ek += (SQR((particlelist + i)->velocity.x) + SQR((particlelist + i)->velocity.y));
    }
    Ek *= 0.5;
  }

  current_pressure += 2.0*Ek*NUMBER_OF_PARTICLES/(SQR(GRIDSIZE)*(2.0*NUMBER_OF_PARTICLES-2));
  current_temperature = Ek/(NUMBER_OF_PARTICLES-1.0);

  kinetic_energy_array[cycle]   = Ek;
  potential_energy_array[cycle] = Ev;
  pressure_array[cycle]         = current_pressure;
  temperature_array[cycle]      = current_temperature;

  if(cycle == INITIALISATION_STEPS) Energy_Reference = Ek + Ev;
  if(cycle > INITIALISATION_STEPS){
    averages[0] += Ek;
    averages[1] += Ev;
    averages[2] += fabs(Energy_Reference - (Ek + Ev))/Energy_Reference;
    averages[3] += current_pressure;
    averages[4] += current_temperature;
  }
}


void gnuprint(FILE *gp){
  int i;

  // fack c
  char options[200] = "unset autoscale\nset xrange [0:";
  char a[] = "]\nset yrange [0:";
  char b[] = "]\nplot '-'\n";
  char c[3];
  sprintf(c, "%d", GRIDSIZE);
  strcat(options, c);
  strcat(options, a);
  strcat(options, c);
  strcat(options, b);
   
  fprintf(gp, options);

  for (i=0; i<NUMBER_OF_PARTICLES; i++) fprintf(gp, "%g %g\n", (particlelist + i)->position.x , (particlelist + i)->position.y );

  fflush(gp);
  fprintf(gp, "e\n");
}

void LiveLinePrint(FILE *gp, int i){
  int j;

  char options[400] = "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 4 ps 1.5 \nset style line 13 lc rgb '#09ad00' lt 1 lw 1.5\nset style line 8  lc rgb '#00ad88' lt 1 lw 1.5\nplot '-' using 1:2 with linespoints title 'Potential Energy', '-' using 1:3 with linespoints title 'Kinetic Energy', '-' using 1:4 with linespoints title 'Total Energy'\n ";

  // printf("%d\n", i);
  fprintf(gp, options);

  for (j=0; j<i; j+=10){
    fprintf(gp, "%d %g %g %g\n", j , potential_energy_array[j],kinetic_energy_array[j],potential_energy_array[j]+kinetic_energy_array[j]  );
    // printf("%d %lf %lf %lf\n",j , potential_energy_array[j],kinetic_energy_array[j],potential_energy_array[j]+kinetic_energy_array[j]  );
  } 

  fflush(gp);
  fprintf(gp, "e\n");
}

void WriteToFile(FILE *gp, int i){
  fprintf(gp, "%d %g %g %g\n", i , potential_energy_array[i],kinetic_energy_array[i],potential_energy_array[i]+kinetic_energy_array[i]  );
  fflush(gp);
}

void AssignCellnumber(int ParticleIndex){
  int poep = (particlelist + ParticleIndex)->position.x;
  int poeper = (particlelist + ParticleIndex)->position.y;
  (particlelist + ParticleIndex)->cellnumber =  poep + poeper * GRIDSIZE;
}

void displace_particles(){
  int i;
  for (i = 0; i < NUMBER_OF_PARTICLES; i++) {    
    // Velocity Verlet. force[0] = force at previous timestep. force[1] = force at current timestep.
    (particlelist + i)->velocity = VectorAddition((particlelist + i)->velocity, VectorScalar((particlelist + i)->force[0], (0.5*DELTAT)));
    (particlelist + i)->position = VectorAddition((particlelist + i)->position, VectorScalar((particlelist + i)->velocity, DELTAT));
    (particlelist + i)->potential = 0;
    (particlelist + i)->pressure_contribution = 0;
    
    // apply pbc and assign to correct cell
    if( (particlelist + i)->position.y > GRIDSIZE ){
      (particlelist + i)->position.y -= GRIDSIZE; 
    } else if( (particlelist + i)->position.y < 0 ){
      (particlelist + i)->position.y += GRIDSIZE;
    }
    if( (particlelist + i)->position.x > GRIDSIZE ){
      (particlelist + i)->position.x -= GRIDSIZE;
    } else if( (particlelist + i)->position.x < 0 ){
      (particlelist + i)->position.x += GRIDSIZE;
    }
    AssignCellnumber(i);
  }
}

void clean_exit_on_sig(int sig_num){
  // printf ("\n Signal %d received",sig_num);
}