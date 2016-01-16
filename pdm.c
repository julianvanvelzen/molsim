#include "system.h"

int main(int argc, char** argv) {

  // initialisation
  InitializeRandomNumberGenerator(time(0l));
  int i,j,k;
  int world_rank, world_size;
  double startwtime, endwtime;
  if(argc>3){
    NUMBER_OF_PROCESSORS = atoi(argv[1]);
    TEMPERATURE          = atoi(argv[2]);
    NUMBER_OF_CYCLES     = atoi(argv[3]);
    NUMBER_OF_PARTICLES  = atoi(argv[4]);
    RCUT                 = atof(argv[5]);
    REPULSIVE_CST        = atoi(argv[6]);
    GRIDSIZE             = atoi(argv[7]);
  } else {
    printf("Not enough parameters\n");
    exit(1);
  }

  particlelist = malloc(sizeof(Particle)*NUMBER_OF_PARTICLES);

  // Initialize the MPI environment
  MPI_Init(&argc, &argv);

  // Find out rank, size
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // initialize particles
  if(world_rank==0){
    Initialize();
    qsort(particlelist, NUMBER_OF_PARTICLES, sizeof(Particle), cmpfunc);
    startwtime = MPI_Wtime();

    printf("argc:                 %d\
          \nNUMBER_OF_PROCESSORS: %d\
          \nTEMPERATURE:          %d\
          \nNUMBER_OF_CYCLES:     %d\
          \nNUMBER_OF_PARTICLES:  %d\
          \nRCUT:                 %lf\
          \nREPULSIVE_CST         %d\
          \nGRIDSIZE:             %d\ 
          \n\n\n", argc, NUMBER_OF_PROCESSORS, TEMPERATURE, NUMBER_OF_CYCLES, NUMBER_OF_PARTICLES, RCUT, REPULSIVE_CST, GRIDSIZE );
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  Mdloop(world_rank);

  if(world_rank==0) {
    endwtime = MPI_Wtime();
  


    // for(i=0; i< NUMBER_OF_PARTICLES; i++)
    //   printf("Particle %d Cellnumber: %d\n",i, (particlelist+i)->cellnumber );
    

    printf("Time elapsed: %fms\n", (endwtime-startwtime)*1000);
  }

  MPI_Finalize();

}



