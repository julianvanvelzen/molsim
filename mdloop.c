#include "system.h"

void Mdloop(world_rank){
  int i,j,k;
  int dummy;
  double size;
  Vector dF;
  Vector force;
  double dE;
  double energy;
  double distance; 
  int neighbouringcells[4];
  int indices[2*NUMBER_OF_PROCESSORS];

  neighbouringcells[0] = world_rank + GRIDSIZE;
  neighbouringcells[1] = world_rank + GRIDSIZE + 1;
  neighbouringcells[2] = world_rank + 1;
  neighbouringcells[3] = world_rank - GRIDSIZE + 1;

  // rechts 
  if ((world_rank+1)%GRIDSIZE == 0 ){
    neighbouringcells[1] = (world_rank + 1)%(SQR(GRIDSIZE));
    neighbouringcells[2] -= GRIDSIZE;
    neighbouringcells[3] = world_rank - 2*GRIDSIZE+1;
  }
  // boven
  if (world_rank+GRIDSIZE > SQR(GRIDSIZE)-1){
    neighbouringcells[0] = (world_rank + GRIDSIZE)%GRIDSIZE;
    neighbouringcells[1] = neighbouringcells[1] % SQR(GRIDSIZE);
  }
  // onder
  if (neighbouringcells[3] < 0)
    neighbouringcells[3] += SQR(GRIDSIZE);

  size = NUMBER_OF_PARTICLES * sizeof(Particle);
  for(i = 0; i < NUMBER_OF_CYCLES; i++){
    MPI_Bcast(particlelist, size , MPI_BYTE, 0, MPI_COMM_WORLD);

    if (world_rank == 0){
        printf("indices \n");
        getindeces(particlelist, &indices);
        for (j =0; j<(2*NUMBER_OF_PROCESSORS); j++){
            printf("%d\n", indices[j]);
        }
    }

    for( j=0; j< NUMBER_OF_PARTICLES; j++){
        if ( (particlelist + j)->cellnumber == world_rank){
            for(k=0;k<NUMBER_OF_PARTICLES;k++){
                ForceEnergy((particlelist + j)->position, (particlelist + j)->position,&dF,&dE);
                distance += VectorDistance((particlelist + j)->position, (particlelist + k)->position);
                force = VectorAddition(force, dF);
                // printf("world = %d distance = %lf %d %d\n", world_rank, , j,k);
            }

            (particlelist + j)->force = force;
            dE   = 0;
            dF.x = 0;
            dF.y = 0;

        }

    }

    //pak je eigen particles
    //bereken vector distance met alle andere
    //print
    
    // printf("Data: %d, Cycle: %d Process:%d size: %d\n", Map[2][2].particles[2].position.x, i, world_rank, size );    
    

    // MPI_Reduce(&i, &dummy, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
    MPI_Barrier(MPI_COMM_WORLD);

  }
}



    // 1: gehele map gaat naar helft van processors
    // 2.1: cell reset forces
    // 2: berekent forces voor eigen cell en cell +- 1 
    // 3: alle copien van map gaat terug naar master
    // 4: voegt alle copien samen
    // 5: stuurt terug naar nu alle processors
    // 6: update posities
    // 7: slave sturen alles terug
    // 8: master doet particles in juiste cell (als ze versprongen)

    // to do:
    // 1: voor map[x][y] (afhankelijk van world rank) doe je berekeningen
    //      -bewegingsstap
    //      -radial distr
    //      -
    // 2: als deeltje cell uit gaat, communiceer je naar andere cell dat er eentje bij komt
    // 3: geef met mpi_reduce je nieuwe cell terug
    // 4: master collecteerd ze allemaal en maakt er nieuwe Map van

    // if(world_rank==0) printf("sum of Processes: %d\n", dummy );