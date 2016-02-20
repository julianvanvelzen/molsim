#include "system.h"

double TEMPERATURE;
double RCUT;
double REPULSIVE_CST;
double DELTAT;
double initialisation_sum;
int NUMBER_OF_CYCLES;
int NUMBER_OF_PROCESSORS;
int NUMBER_OF_PARTICLES;
int GRIDSIZE;
Particle *particlelist;

int main(int argc, char** argv) {

  // initialisation
  InitializeRandomNumberGenerator(time(0l));
  int world_rank, world_size;
  double startwtime, endwtime;
  if(argc>3){
    NUMBER_OF_PROCESSORS = atoi(argv[1]);
    TEMPERATURE          = atoi(argv[2]);
    NUMBER_OF_CYCLES     = atoi(argv[3]);
    NUMBER_OF_PARTICLES  = atoi(argv[4]);
    RCUT                 = atof(argv[5]);
    REPULSIVE_CST        = atof(argv[6]);
    GRIDSIZE             = atoi(argv[7]);
    DELTAT		 = atof(argv[8]);
  } else {
    printf("Not enough parameters\n");
    exit(1);
  }

  CheckInputErrors();

  particlelist = malloc(sizeof(Particle)*NUMBER_OF_PARTICLES);

  // Initialize the MPI environment
  MPI_Init(&argc, &argv);

  // Find out rank, size
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);


  // initialize particles
  if(world_rank==0){
    startwtime = MPI_Wtime();
    Initialize();
    qsort(particlelist, NUMBER_OF_PARTICLES, sizeof(Particle), cmpfunc);

    printf("argc:                 %d\
          \nNUMBER_OF_PROCESSORS: %d\
          \nTEMPERATURE:          %lf\
          \nNUMBER_OF_CYCLES:     %d\
          \nNUMBER_OF_PARTICLES:  %d\
          \nRCUT:                 %lf\
          \nREPULSIVE_CST:        %lf\
          \nGRIDSIZE:             %d\
          \nDELTAT:		  %lf\
          \n\n\n", argc, NUMBER_OF_PROCESSORS, TEMPERATURE, NUMBER_OF_CYCLES, NUMBER_OF_PARTICLES, RCUT, REPULSIVE_CST, GRIDSIZE, DELTAT);
  }
  
  Mdloop(world_rank);

  if(world_rank==0) {
    endwtime = MPI_Wtime();
    printf("\nTime elapsed: %fms\n",(endwtime-startwtime)*1000);
  }

  // signal(SIGSEGV, clean_exit_on_sig); 
  // printf("Processor %d is klaar!\n", world_rank );
  MPI_Finalize();
  return 0;
}



