#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ran_uniform.h"
#include <unistd.h>
#include <signal.h>

typedef struct {
  double x;
  double y;
} Vector;

typedef struct {
  Vector position;  
  Vector velocity;

  // force[0] = force at previous timestep. force[1] = force at current timestep
  Vector force[2];

  int cellnumber;
  double potential;
  int radial_distribution[21]; // radial_distribution[20] = outside Rcut
  double pressure_contribution;
} Particle;

typedef struct {
  int start;
  int end;
  int totalcount;
  int neighbouringcells[8];
} Cell;

extern Particle *particlelist;
Particle *particlelist;

int main(int argc, char** argv) {

  int i,j,k;
  int world_rank, world_size;
  int NUMBER_OF_PARTICLES = 10;
  double size;
  particlelist = malloc(sizeof(Particle)*NUMBER_OF_PARTICLES);
  Particle *gather;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  for (i=0;i<NUMBER_OF_PARTICLES;i++){
    for (j=0;j<21;j++){
      (particlelist+i)->radial_distribution[j] = 6;
    }
  }

  if(world_rank == 0) 
    gather = calloc(NUMBER_OF_PARTICLES * 2, sizeof(Particle));

  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  for (i=0;i<1000; i++){
    for (k=0;k<NUMBER_OF_PARTICLES;k++){
      for (j=0;j<21;j++){
        (particlelist+k)->radial_distribution[j] = 6;
        printf("%d\n",(particlelist+k)->radial_distribution[j]  );
      }
    }
    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);
    
  }
  if(world_rank == 0)
    free(gather);
  free(particlelist);

}
