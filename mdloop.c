#include "system.h"

void Mdloop(world_rank){
  int i,j,k,l,m;
  double size;
  FILE * gp = popen ("gnuplot -persistent", "w");

  Cell cells[NUMBER_OF_PROCESSORS];
  Particle *gather;
  gather = malloc(sizeof(Particle) * NUMBER_OF_PARTICLES * NUMBER_OF_PROCESSORS );
 
  for (i = 0; i < NUMBER_OF_PROCESSORS; i++)
    getNearbyCoordinates(&cells, i);



  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  for(i = 0; i < NUMBER_OF_CYCLES; i++){

    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);

    setindeces(particlelist, &cells);

    loopforces(&cells, world_rank);

    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);

    if (world_rank == 0){
        sum_contributions(&cells, gather);
        gnuprint(gp);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);

    displace_particles();

    if (world_rank == 0)
        qsort(particlelist, NUMBER_OF_PARTICLES, sizeof(Particle), cmpfunc);

  }
  free(gather);
}

