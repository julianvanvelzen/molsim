#include "system.h"

void Mdloop(world_rank){
  int i,j,k,l,m;
  double size;
  Vector dF;
  double dE;

  Cell cells[NUMBER_OF_PROCESSORS];
 
  for (i = 0; i < NUMBER_OF_PROCESSORS; i++) getNearbyCoordinates(&cells, i);

  Particle *gather;
  FILE * gp = popen ("gnuplot -persist", "w");
  gather = malloc(sizeof(Particle) * NUMBER_OF_PARTICLES * NUMBER_OF_PROCESSORS );

  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  for(i = 0; i < NUMBER_OF_CYCLES; i++){

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0){
        displace_particles();
        qsort(particlelist, NUMBER_OF_PARTICLES, sizeof(Particle), cmpfunc);
    }

    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);

    setindeces(particlelist, &cells);

    loopforces(&cells, world_rank);

    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);

    if (world_rank == 0){
        sum_contributions(&cells, gather);
        ApplyNewForces();
    }
    if (world_rank == 1) gnuprint(gp);

  }
  free(gather);  
}
