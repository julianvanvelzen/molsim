#include "system.h"

double averages[4] = {0, 0, 0, 0}; // [0] = kinetic, [1] = potential, [2] = sum, [3] = pressure  
double *kinetic_energy_array;
double *potential_energy_array;
double rdf_total[NUMBER_OF_BINS] = { 0 };


void Mdloop(world_rank){
  int i, j,k;
  double size, normalisation;
  double startwtime, endwtime;
  double time1, time2, time3, time4, time5;

  Cell cells[NUMBER_OF_PROCESSORS];
 
  startwtime = MPI_Wtime();
  for (i = 0; i < NUMBER_OF_PROCESSORS; i++) getNearbyCoordinates(cells, i);
  endwtime = MPI_Wtime();
  time1 =  (endwtime-startwtime)*1000;

  FILE * gp = popen ("gnuplot -persist", "w");
  FILE * gphist = popen ("gnuplot -persist", "w");

  Particle *gather;

  // allocate memory
  if(world_rank == 0) {
    gather = calloc(NUMBER_OF_PARTICLES * NUMBER_OF_PROCESSORS, sizeof(Particle));
    kinetic_energy_array = malloc(sizeof(double) * NUMBER_OF_CYCLES);
    potential_energy_array = malloc(sizeof(double) * NUMBER_OF_CYCLES);

    if(gather == NULL || kinetic_energy_array == NULL || potential_energy_array == NULL){
      printf("Error during memory allocation");
      exit(0);
    }
  }



  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);
  


  // main loop
  for(i = 0; i < NUMBER_OF_CYCLES; i++){

    // for (j = 0; j < NUMBER_OF_PARTICLES; j++){
    //   for(k=0;k<NUMBER_OF_BINS;k++){
    //     (particlelist+i)->radial_distribution[k] += 1;
    //     printf("moree %d\n", (particlelist+i)->radial_distribution[k]);
    //   }

    // }    
  

    if (world_rank == 0){
      displace_particles();
      qsort(particlelist, NUMBER_OF_PARTICLES, sizeof(Particle), cmpfunc);
    }

    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);


    startwtime = MPI_Wtime();
    setindeces(particlelist, cells);
    endwtime = MPI_Wtime();
    time2 +=  (endwtime-startwtime)*1000;

    startwtime = MPI_Wtime();
    loopforces(cells, world_rank);
    endwtime = MPI_Wtime();
    time3 +=  (endwtime-startwtime)*1000;

    startwtime = MPI_Wtime();    
    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);
    endwtime = MPI_Wtime();
    time5 +=  (endwtime-startwtime)*1000;

    if (world_rank == 0){
      startwtime = MPI_Wtime();
      sum_apply_contributions(cells, gather, i, world_rank);
      endwtime = MPI_Wtime();
      time4 =  (endwtime-startwtime)*1000;
      // if (i%20 == 0 && i>INITIALISATION_STEPS) HistPrint(gphist, i);
      // gnuprint(gp);
    }

  }

  if(world_rank == 0){
    normalisation = 1.0/(NUMBER_OF_CYCLES-INITIALISATION_STEPS);
    for (i=0; i<4;  i++) averages[i] *= normalisation;
    printf("Averages:\
            \nKinetic:      %lf\
            \nPotential:    %lf\
            \nEnergy drift: %lf\
            \nPressure:     %lf\n", averages[0], averages[1], averages[2], averages[3]);
    for (i=0; i<NUMBER_OF_BINS; i++) {
      printf("\nrdf total %d\n", rdf_total[i]);
      rdf_total[i] *= ((1.0/(M_PI * (SQR((i+1)*RCUT/NUMBER_OF_BINS)-SQR(i*RCUT/NUMBER_OF_BINS)))) * (normalisation/NUMBER_OF_PARTICLES));
      printf("rcut20 %lf sqr-sqr%lf norm %lf rdf %lf\n", (i+1)*RCUT/NUMBER_OF_BINS, (SQR((i+1)*RCUT/NUMBER_OF_BINS)-SQR(i*RCUT/NUMBER_OF_BINS)), normalisation/NUMBER_OF_PARTICLES, rdf_total[i]);
    }
  }
  // printf("\n\n world rank%d\ngetNearbyCoordinates %lf\n setindeces %lf\n loopforces %lf\n ApplyNewForces %lf\n gather %lf\n sum %lf\n",world_rank, time1, time2, time3, time4, time5, time1 +  time2 +  time3 +  time4 );

  // clear dynamic allocated memory
  if (world_rank == 0){
    free(gather);
    free(kinetic_energy_array);
    free(potential_energy_array);
  }
}
