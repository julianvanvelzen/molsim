#include "system.h"

void Mdloop(world_rank){
  int i,j,k,l,m;
  int dummy;
  double size;
  Vector dF;
  double dE;
  double energy;
  double distance; 
  Cell cell[NUMBER_OF_PROCESSORS];
  Particle *gather;
  gather = malloc(sizeof(Particle) * NUMBER_OF_PARTICLES * NUMBER_OF_PROCESSORS );
 
  for (i = 0; i < NUMBER_OF_PROCESSORS; i++)
    getNearbyCoordinates(&cell, i);



  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  for(i = 0; i < NUMBER_OF_CYCLES; i++){
    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);
    setindeces(particlelist, &cell);

    // loop over all particles in your own cell
    for(i = cell[world_rank].start; i < cell[world_rank].end; i++){

        // loop over all other particles in your own cell
        for(k = cell[world_rank].start; k < cell[world_rank].end; k++){
            if(k <= i)
                continue;
            dF.x = 0;
            dF.y = 0;
            ForceEnergy((particlelist + i)->position, (particlelist + k)->position,&dF,&dE);
            (particlelist + i)->force = VectorAddition((particlelist + i)->force, dF);
            (particlelist + k)->force = VectorAddition((particlelist + k)->force, VectorFlip(dF));
        }

        // loop over all particles in neighbouring cells
        for (l=0; l<5; l++){    
            for (m = cell[cell[world_rank].neighbouringcells[l]].start; m < cell[cell[world_rank].neighbouringcells[l]].end; m++){
                dF.x = 0;
                dF.y = 0;
                ForceEnergy(particlelist[i].position, particlelist[m].position,&dF,&dE);
                particlelist[i].force = VectorAddition(particlelist[i].force, dF);
                particlelist[m].force = VectorAddition(particlelist[m].force, VectorFlip(dF));
            }
        }

    }


    MPI_Gather(particlelist, size, MPI_BYTE, gather, size, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (world_rank == 0){
        for (k = 0; k < NUMBER_OF_PARTICLES; k++){
            for (j = 0; j < 8; j++){
                (particlelist + k)->force = VectorAddition( (particlelist + k)->force, (gather+NUMBER_OF_PARTICLES*cell[(particlelist+k)->cellnumber].neighbouringcells[j] + k)->force );  
            }
        }
        for (j =0; j< NUMBER_OF_PARTICLES; j++){
            printf("particle %d: force: %lf\n", j, (particlelist + j)->force.x);
        }
    }  


    MPI_Barrier(MPI_COMM_WORLD);
  }
  free(gather);
}

