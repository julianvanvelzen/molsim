#include "system.h"

double averages[5] = {0, 0, 0, 0, 0}; // [0] = kinetic, [1] = potential, [2] = sum, [3] = pressure, [4] = temperature
double *kinetic_energy_array;
double *potential_energy_array;
double *pressure_array;
double *temperature_array;
double rdf_total[NUMBER_OF_BINS] = { 0 };

void Mdloop(world_rank){
  int i, j, k;
  double size, normalisation, shell_area;
  double startwtime, endwtime;
  double time1, time2, time3, time4, time5;

  Cell cells[NUMBER_OF_PROCESSORS];
  for (i = 0; i < NUMBER_OF_PROCESSORS; i++) getNearbyCoordinates(cells, i);

  char filepathEnergy[250]; 
  char filepathRadial[250]; 
  char filepathMomentum[250]; 
  char filepathForce[250]; 

  Particle *gather;
  FILE * gp = popen ("gnuplot -persist", "w");
  FILE * gphist = popen ("gnuplot -persist", "w");
  FILE * fpRadial;
  FILE * fpMomentum;
  FILE * fpEnergy;
  FILE * fpForce;  

  // allocate memory
  if(world_rank == 0) {

    sprintf(filepathMomentum, "data/MomentumForT%lfRpCst%lfRCUT%lfDT%lfNC%dNPar%dNpro%d.dat", TEMPERATURE, REPULSIVE_CST, RCUT, DELTAT, NUMBER_OF_CYCLES, NUMBER_OF_PARTICLES, NUMBER_OF_PROCESSORS);
    sprintf(filepathEnergy, "data/EnergyForT%lfRpCst%lfRCUT%lfDT%lfNC%dNPar%dNpro%d.dat", TEMPERATURE, REPULSIVE_CST, RCUT, DELTAT, NUMBER_OF_CYCLES, NUMBER_OF_PARTICLES, NUMBER_OF_PROCESSORS);
    sprintf(filepathRadial, "data/RadialForT%lfRpCst%lfRCUT%lfDT%lfNC%dNPar%dNpro%d.dat", TEMPERATURE, REPULSIVE_CST, RCUT, DELTAT, NUMBER_OF_CYCLES, NUMBER_OF_PARTICLES, NUMBER_OF_PROCESSORS);
    sprintf(filepathForce, "data/ForceForT%lfRpCst%lfRCUT%lfDT%lfNC%dNPar%dNpro%d.dat", TEMPERATURE, REPULSIVE_CST, RCUT, DELTAT, NUMBER_OF_CYCLES, NUMBER_OF_PARTICLES, NUMBER_OF_PROCESSORS);

    fpMomentum = fopen(filepathMomentum, "w+");
    fpRadial = fopen(filepathRadial, "w+");
    fpEnergy = fopen(filepathEnergy, "w+");
    fpForce = fopen(filepathForce, "w+");
    
    gather = calloc(NUMBER_OF_PARTICLES * NUMBER_OF_PROCESSORS, sizeof(Particle));
    kinetic_energy_array   = malloc(sizeof(double) * NUMBER_OF_CYCLES);
    potential_energy_array = malloc(sizeof(double) * NUMBER_OF_CYCLES);
    pressure_array         = malloc(sizeof(double) * NUMBER_OF_CYCLES);
    temperature_array      = malloc(sizeof(double) * NUMBER_OF_CYCLES);

    if(gather == NULL || kinetic_energy_array == NULL || potential_energy_array == NULL || pressure_array == NULL || temperature_array == NULL){
      printf("Error during memory allocation");
      exit(0);
    }
  }

  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);

  // main loop
  for(i = 0; i < NUMBER_OF_CYCLES; i++){

    if (world_rank == 0){
      startwtime = MPI_Wtime();
      displace_particles();
      endwtime   = MPI_Wtime(); time1 += (endwtime-startwtime)*1000;

      startwtime = MPI_Wtime();
      qsort(particlelist, NUMBER_OF_PARTICLES, sizeof(Particle), cmpfunc);
      endwtime   = MPI_Wtime(); time2 += (endwtime-startwtime)*1000;

    }

    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);

    setindeces(particlelist, cells);
    
    startwtime = MPI_Wtime();
    loopforces(cells, world_rank);
    endwtime = MPI_Wtime(); time3 += (endwtime-startwtime)*1000;

    startwtime = MPI_Wtime();    
    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);
    endwtime = MPI_Wtime(); time4 += (endwtime-startwtime)*1000;

    if (world_rank == 0){
      startwtime = MPI_Wtime();
      sum_apply_contributions(cells, gather, i);
      endwtime = MPI_Wtime(); time5 += (endwtime-startwtime)*1000;

      if (i%10 == 0 && i > 0 ) LiveLinePrint(gphist, i);     
      if (i%10 == 0 && i > 0 ) WriteToFile(fpEnergy, i);
      if (i%50 == 0 && i > 0 ) printf("\nKinetic Energy: %lf \t Potential Energy: %lf \t Sum: %lf \t Pressure: %lf \t Temperature: %lf \t ", \
        kinetic_energy_array[i], potential_energy_array[i], kinetic_energy_array[i] + potential_energy_array[i], pressure_array[i], temperature_array[i]);
    }
    // if (world_rank == 1) gnuprint(gp);
  }

  if(world_rank == 0){
    normalisation = 1.0/(NUMBER_OF_CYCLES-INITIALISATION_STEPS);
            
    for (i=0; i<4;  i++) averages[i] *= normalisation;
    printf("\n\nAverages:\
            \nKinetic energy:      %lf\
            \nPotential energy:    %lf\
            \nEnergy drift:        %lf\
            \nPressure:            %lf\n\n", averages[0], averages[1], averages[2], averages[3]);
    for (i=0; i<NUMBER_OF_BINS; i++) {
      rdf_total[i] *= normalisation/(M_PI * (SQR((i+1)*RCUT/NUMBER_OF_BINS)-SQR(i*RCUT/NUMBER_OF_BINS))*SQR(NUMBER_OF_PARTICLES)/SQR(GRIDSIZE));
      fprintf(fpRadial, "%d %g\n", i , rdf_total[i]  );
      printf("RDF bin %d:\t%lf\n", i, rdf_total[i]);
    }
    printf("\n\nTime spent:\ndisplace_particles %lf\nqsort %lf\nloopforces %lf\nMPI_gather %lf\nsum_apply_contributions %lf\n", time1, time2, time3, time4, time5);
  }
  

  // clear dynamic allocated memory
  if (world_rank == 0){
    free(gather);
    free(kinetic_energy_array);
    free(potential_energy_array);
    free(pressure_array);
    free(temperature_array);
    fclose( fpEnergy );
    fclose( fpRadial );
    fclose( fpMomentum );
    fclose( fpForce );
  }
}
