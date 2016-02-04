#include "system.h"

void Mdloop(world_rank){
  int i,j,k,l,m;
  double size;
  Vector dF;
  double dE;
  double startwtime, endwtime;

  Cell cells[NUMBER_OF_PROCESSORS];
 
  startwtime = MPI_Wtime();
  for (i = 0; i < NUMBER_OF_PROCESSORS; i++) getNearbyCoordinates(&cells, i);
  endwtime = MPI_Wtime();
  // printf("Time elapsed getNearbyCoordinates: %fms\n", (endwtime-startwtime)*1000);

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

    startwtime = MPI_Wtime();
    setindeces(particlelist, &cells);
    endwtime = MPI_Wtime();
    // printf("Time elapsed setindeces: %fms\n", (endwtime-startwtime)*1000);

    startwtime = MPI_Wtime();
    loopforces(&cells, world_rank);
    endwtime = MPI_Wtime();
    // printf("Time elapsed loopforces: %fms\n", (endwtime-startwtime)*1000);

    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);

    
    if (world_rank == 0){
        // printf("####################\nCycle %d\n####################\n", i);
        // for (j = 0; j< NUMBER_OF_PARTICLES; j++)
        //   printf("Particle %d\n %s\n", j,  PARTICLE_DUMP( *(particlelist+j)) );
        // for (j = 0; j< 8; j++)
        //   printf("cell %d\n %s\n",j, CELL_DUMP(cells[j]) );


        startwtime = MPI_Wtime();
        sum_contributions(&cells, gather);
        ApplyNewForces();
        endwtime = MPI_Wtime();
        // printf("Time elapsed sum_contributions & ApplyNewForces: %fms\n", (endwtime-startwtime)*1000);
    }
    // if (world_rank == 1) gnuprint(gp);
  }
  free(gather);  
}
