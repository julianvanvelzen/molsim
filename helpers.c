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

void ForceEnergy(Particle *p1, Particle *p2){
  double distance = VectorDistance(p1->position, p2->position);
  if (distance > RCUT) return;

  double force = (2.0*REPULSIVE_CST*sqrt(SQR(distance-RCUT))/SQR(RCUT));
  
  Vector relative_position;
  relative_position.x = (p1->position.x - p2->position.x);
  relative_position.y = (p1->position.y - p2->position.y);

  Vector forceVector;
  forceVector = VectorScalar(relative_position, force/distance);
  
  p1->force[1] = VectorAddition(p1->force[1], forceVector);     
  p2->force[1] = VectorAddition(p2->force[1], VectorScalar(forceVector, -1)); 

  double potential = REPULSIVE_CST*SQR(distance-RCUT)/SQR(RCUT);
  p1->potential += potential;
  p2->potential += potential;

  double pressure = (forceVector.x*relative_position.x + forceVector.y*relative_position.y) / (2 * SQR(GRIDSIZE));
  p1->pressure_contribution += pressure;
  p2->pressure_contribution += pressure;

  int rdf_bin_index = NUMBER_OF_BINS*distance/RCUT;
  p1->radial_distribution[rdf_bin_index] += 2;
}

Vector VectorAddition(Vector v1, Vector v2){
  Vector v3;
  v3.x = v1.x + v2.x;
  v3.y = v1.y + v2.y;
  return v3;
}

Vector VectorScalar(Vector vector, double factor){
  vector.x *= factor;
  vector.y *= factor;
  return vector;
}

double VectorDistance(Vector v1, Vector v2){
  double distance;
  distance = sqrt(SQR(v1.x - v2.x) + SQR(v1.y - v2.y));
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
  (cells+currentPosition)->totalcount = 0;

  (cells+currentPosition)->neighbouringcells[0] = currentPosition + GRIDSIZE;
  (cells+currentPosition)->neighbouringcells[1] = currentPosition + GRIDSIZE + 1;
  (cells+currentPosition)->neighbouringcells[2] = currentPosition + 1;
  (cells+currentPosition)->neighbouringcells[3] = currentPosition - GRIDSIZE + 1;
  (cells+currentPosition)->neighbouringcells[4] = currentPosition - GRIDSIZE;
  (cells+currentPosition)->neighbouringcells[5] = currentPosition - GRIDSIZE - 1;
  (cells+currentPosition)->neighbouringcells[6] = currentPosition - 1;
  (cells+currentPosition)->neighbouringcells[7] = currentPosition + GRIDSIZE - 1;

  // boven
  if (currentPosition+GRIDSIZE > SQR(GRIDSIZE)-1){
    (cells+currentPosition)->neighbouringcells[7] = (currentPosition + GRIDSIZE-1)%GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[0] = (currentPosition + GRIDSIZE)%GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[1] = (currentPosition + GRIDSIZE+1)%GRIDSIZE;
  }
  // rechts 
  if ((currentPosition+1)%GRIDSIZE == 0 ){
    (cells+currentPosition)->neighbouringcells[1] = (currentPosition + 1)%(SQR(GRIDSIZE));
    (cells+currentPosition)->neighbouringcells[2] -= GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[3] = currentPosition - 2*GRIDSIZE+1;
  }
  // onder
  if (currentPosition-GRIDSIZE < 0){
    (cells+currentPosition)->neighbouringcells[3] += SQR(GRIDSIZE);
    (cells+currentPosition)->neighbouringcells[4] += SQR(GRIDSIZE);
    (cells+currentPosition)->neighbouringcells[5] += SQR(GRIDSIZE);
  }
  // links
  if(currentPosition%GRIDSIZE==0){
    (cells+currentPosition)->neighbouringcells[5] += GRIDSIZE;   
    (cells+currentPosition)->neighbouringcells[6] += GRIDSIZE; 
    (cells+currentPosition)->neighbouringcells[7] += GRIDSIZE;
  }
}

int cmpfunc (const void * a, const void * b){
  Particle *A = (Particle *)a;
  Particle *B = (Particle *)b;
  // int c = A->cellnumber - B->cellnumber;
  // free(A);
  // free(B);
  return ( A->cellnumber - B->cellnumber );
} 

void setindeces(Particle *particlelist, Cell *cells){
  int i = 0, j;
  int currentcell = 0;
  
  while (i<NUMBER_OF_PARTICLES)
  {
    while ((particlelist+i)->cellnumber > currentcell )
    {
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

  // loop over all particles in your own cell
  i = (cells + world_rank)->start-1;
  while(1) 
  {
    i++;
    if(i == (cells + world_rank)->end) break;

    // loop over all other particles in your own cell
    k = i;
    while(1)
    {
      k++;
      if(k >= (cells + world_rank)->end) break;
      ForceEnergy((particlelist + i), (particlelist + k));
    }

    // loop over all neighbouring cells
    for (l=0; l<4; l++)
    {
      // loop over all particles in neighbouring cells
      m = (cells + (cells+world_rank)->neighbouringcells[l])->start-1;
      while (1)
      {
        m++;
        if (m >= (cells + (cells+world_rank)->neighbouringcells[l])->end) break;

        // apply periodic boundary conditions
        if ( (l == 0 || l == 1) && (cells+world_rank)->neighbouringcells[l] < GRIDSIZE) (particlelist + i)->position.y -= GRIDSIZE;
        if ( (l >= 1) && (cells+world_rank)->neighbouringcells[l] % GRIDSIZE == 0)      (particlelist + i)->position.x -= GRIDSIZE;

        ForceEnergy((particlelist + i), (particlelist + m));

        if ( (l == 0 || l == 1) && (cells+world_rank)->neighbouringcells[l] < GRIDSIZE) (particlelist + i)->position.y += GRIDSIZE;              
        if ( (l >= 1) && (cells+world_rank)->neighbouringcells[l] %GRIDSIZE == 0)       (particlelist + i)->position.x += GRIDSIZE;
      }
    }
  }
}

void sum_apply_contributions(Cell *cells, Particle *gather, int cycle){
  int i,j,k, neighbour_offset, current_box_offset;
  double Ek = 0;
  double Ev = 0;
  double current_pressure = 0;

  for (i = 0; i < NUMBER_OF_PARTICLES; i++)
  {
    current_box_offset = (particlelist+i)->cellnumber;
    for (j = 7; j >= 4 ; j--) {   // Only the contribution of bottom and left neighbours needs to be considered
      neighbour_offset = (cells+current_box_offset)->neighbouringcells[j];
      if(neighbour_offset != 0){  // (particlelist + i) already includes the contribution of box 0
        (particlelist+i)->force[1] = VectorAddition((particlelist+i)->force[1], (gather + (NUMBER_OF_PARTICLES*neighbour_offset) + i)->force[1]);
        (particlelist+i)->potential += (gather + (NUMBER_OF_PARTICLES*neighbour_offset) + i)->potential;
        (particlelist+i)->pressure_contribution += (gather + (NUMBER_OF_PARTICLES*neighbour_offset) + i)->pressure_contribution;
      } 
    }
    if(current_box_offset != 0){  // (particlelist + i) already includes the contribution of box 0
      (particlelist+i)->force[1] = VectorAddition((particlelist+i)->force[1], (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->force[1]);
      (particlelist+i)->potential += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->potential;
      (particlelist+i)->pressure_contribution += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->pressure_contribution;
    }
  
    (particlelist + i)->velocity = VectorAddition((particlelist + i)->velocity, VectorScalar((particlelist + i)->force[1], (0.5*DELTAT)));
    (particlelist + i)->force[0] = (particlelist + i)->force[1];
    (particlelist + i)->force[1].x = 0;
    (particlelist + i)->force[1].y = 0;

    if(cycle >= INITIALISATION_STEPS){
      for(k=0;k<NUMBER_OF_BINS;k++) rdf_total[k] += (gather + (NUMBER_OF_PARTICLES*current_box_offset) + i)->radial_distribution[k]; 
      Ek += ( SQR((particlelist + i)->velocity.x) + SQR((particlelist + i)->velocity.y) ) / 2.0;
      Ev += (particlelist + i)->potential;
      current_pressure += (particlelist + i)->pressure_contribution;
    }

    for(k=0;k<NUMBER_OF_BINS;k++) (particlelist + i)->radial_distribution[k] = 0;
  }

  *(kinetic_energy_array + cycle) = Ek;
  *(potential_energy_array + cycle) = Ev;
  current_pressure += 2.0*Ek*NUMBER_OF_PARTICLES/(SQR(GRIDSIZE)*(2.0*NUMBER_OF_PARTICLES-2));

  if(cycle == INITIALISATION_STEPS) Energy_Reference = Ek + Ev;
  if(cycle >= INITIALISATION_STEPS){
    averages[0] += Ek;
    averages[1] += Ev;
    averages[2] += sqrt(SQR(Energy_Reference - (Ek + Ev)))/Energy_Reference;
    averages[3] += current_pressure;
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

void HistPrint(FILE *gp, int i){
  int j;

  char options[400] = "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 4 ps 1.5 \nset style line 13 lc rgb '#09ad00' lt 1 lw 1.5\nset style line 8  lc rgb '#00ad88' lt 1 lw 1.5\nplot '-' using 1:2 with linespoints, '-' using 1:3 with linespoints, '-' using 1:4 with linespoints\n ";

  // printf("%d\n", i);
  fprintf(gp, options);

  for (j=INITIALISATION_STEPS; j<i; j+=20) fprintf(gp, "%d %g %g %g\n", j , potential_energy_array[j],kinetic_energy_array[j],potential_energy_array[j]+kinetic_energy_array[j]  );

  fflush(gp);
  fprintf(gp, "e\n");
}

void AssignCellnumber(int ParticleIndex){
  (particlelist + ParticleIndex)->cellnumber = (int)(particlelist + ParticleIndex)->position.x + GRIDSIZE*(int)(particlelist + ParticleIndex)->position.y;
}

void displace_particles(){
  int i;
  for (i = 0; i < NUMBER_OF_PARTICLES; i++)
  {    
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

/*
char* VECTOR_DUMP(Vector d){ 
  char str[100];
  sprintf(str, "x: %lf y: %lf", d.x, d.y );
  return &str;
}

char* PARTICLE_DUMP(Particle d){ 
  char str[800];
  sprintf(str, "\
              position: %s\n \
              Velocity: %s\n \
              Force t0: %s\n \
              Force t1: %s\n \
              cellnumber: %d\n \
              potential: %lf\n \
              radial_distribution: %lf", \
              VECTOR_DUMP(d.position), \
              VECTOR_DUMP(d.velocity), \
              VECTOR_DUMP(d.force[0]), \
              VECTOR_DUMP(d.force[1]), \
              d.cellnumber, \ 
              d.potential, \ 
              d.radial_distribution );

  return &str;
}

char* INT_ARRAY_DUMP(int length, int data[]  ){
  char str[length];
  char buf[10];
  int i;

  for (i = 0; i < length; i++) {
    sprintf(buf, "%d", data[i]); 
    strcat(str, buf);
  }  
  return &str;
}

char* CELL_DUMP(Cell d){ 
  char str[800];
  sprintf(str, "\n\
               start              %d\n \
              end                %d\n \
              totalcount         %d\n \
              neighbouringcells  %s\n\n", \              
              d.start, \
              d.end, \
              d.totalcount, \
              INT_ARRAY_DUMP(8, d.neighbouringcells) );
  fflush;
  return &str;
}*/
