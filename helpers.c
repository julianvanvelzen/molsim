#include "system.h"

void ForceEnergy(Vector v1, Vector v2,Vector *dF, double *dE){
	double distance = VectorDistance(v1, v2);
	Vector v3;
	if (distance > RCUT) return;
	v3.x = (v1.x - v2.x)/distance;
	v3.y = (v1.y - v2.y)/distance;
	*dE  = REPULSION_CST*SQR(distance-RCUT)/SQR(RCUT);
	dF->x = (2*REPULSION_CST*(distance-RCUT)/SQR(RCUT))*v3.x;
	dF->y = (2*REPULSION_CST*(distance-RCUT)/SQR(RCUT))*v3.y;
}

Vector VectorAddition(Vector v1, Vector v2){
  Vector v3;
  v3.x = v1.x+v2.x;
  v3.y = v1.y+v2.y;
  return v3;
}

Vector VectorFlip(Vector vector){
  vector.x *= -1;
  vector.x *= -1;
  return vector;
}

double VectorDistance(Vector v1, Vector v2){
  double distance;
  distance =  sqrt(SQR(v1.x - v2.x) + SQR(v1.y - v2.y));
  return distance;
}

Vector RanUnit(void){
	double random = RandomNumber() * 2 * M_PI;
	Vector unit;
	unit.x = cos(random);
	unit.y = sin(random);
	return unit;

}

void getNearbyCoordinates(Cell *cells, int currentPosition){
  (cells+currentPosition)->neighbouringcells[0] = currentPosition + GRIDSIZE;
  (cells+currentPosition)->neighbouringcells[1] = currentPosition + GRIDSIZE + 1;
  (cells+currentPosition)->neighbouringcells[2] = currentPosition + 1;
  (cells+currentPosition)->neighbouringcells[3] = currentPosition - GRIDSIZE + 1;
  (cells+currentPosition)->neighbouringcells[4] = currentPosition - GRIDSIZE;
  (cells+currentPosition)->neighbouringcells[5] = currentPosition - GRIDSIZE - 1;
  (cells+currentPosition)->neighbouringcells[6] = currentPosition - 1;
  (cells+currentPosition)->neighbouringcells[7] = currentPosition + GRIDSIZE - 1;

  // links
  if(currentPosition%GRIDSIZE==0){
    (cells+currentPosition)->neighbouringcells[5] += GRIDSIZE;   
    (cells+currentPosition)->neighbouringcells[6] += GRIDSIZE; 
    (cells+currentPosition)->neighbouringcells[7] += GRIDSIZE;
  }

  // rechts 
  if ((currentPosition+1)%GRIDSIZE == 0 ){
    (cells+currentPosition)->neighbouringcells[1] = (currentPosition + 1)%(SQR(GRIDSIZE));
    (cells+currentPosition)->neighbouringcells[2] -= GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[3] = currentPosition - 2*GRIDSIZE+1;
  }
  // boven
  if (currentPosition+GRIDSIZE > SQR(GRIDSIZE)-1){
    (cells+currentPosition)->neighbouringcells[0] = (currentPosition + GRIDSIZE)%GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[1] = (currentPosition + GRIDSIZE+1)%GRIDSIZE;
    (cells+currentPosition)->neighbouringcells[7] = (currentPosition + GRIDSIZE-1)%GRIDSIZE;
  }
  // onder
  if (currentPosition-GRIDSIZE < 0){
    (cells+currentPosition)->neighbouringcells[3] += SQR(GRIDSIZE);
    (cells+currentPosition)->neighbouringcells[4] += SQR(GRIDSIZE);
    (cells+currentPosition)->neighbouringcells[5] += SQR(GRIDSIZE);
  }
}

int cmpfunc (const void * a, const void * b){
  Particle *A = (Particle *)a;
  Particle *B = (Particle *)b;
  // int c = A->cellnumber - B->cellnumber;
  // free(A);
  // free(B);
  return ( A->cellnumber - B->cellnumber );
} 

void setindeces(Particle *particlelist, Cell *cells){
  int i, j;
  int currentcell = 0;


  for(i=0;i<NUMBER_OF_PROCESSORS;i++){
    (cells+i)->start = 0;
    (cells+i)->end = 0;
  }
  
  i = 0;
  
  while (i<NUMBER_OF_PARTICLES){
    while ((particlelist+i)->cellnumber > currentcell ){
      (cells+currentcell)->start = i;
      (cells+currentcell)->end = i;
      currentcell++;
    }
  
    (cells+currentcell)->start = i;
  
    while((particlelist+i)->cellnumber == currentcell){
      i++;
    }
    (cells+currentcell)->end = i;
    currentcell++;


    if(i >= NUMBER_OF_PARTICLES){
      for(j=currentcell;j<NUMBER_OF_PROCESSORS;j++){
        (cells+j)->start = NUMBER_OF_PARTICLES;
        (cells+j)->end   = NUMBER_OF_PARTICLES;
      }
    }
  }
}

void loopforces(Cell *cells, int world_rank){
  int i,k,l,m;
  Vector dF;
  double dE;

  // loop over all particles in your own cell
  for(i = (cells + world_rank)->start; i < (cells + world_rank)->end; i++){

      // loop over all other particles in your own cell
      for(k = (cells + world_rank)->start; k < (cells + world_rank)->end; k++){
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
          for (m = (cells + (cells + world_rank)->neighbouringcells[l])->start; m < (cells + (cells + world_rank)->neighbouringcells[l])->end; m++){
              dF.x = 0;
              dF.y = 0;
              ForceEnergy(particlelist[i].position, particlelist[m].position,&dF,&dE);
              particlelist[i].force = VectorAddition(particlelist[i].force, dF);
              particlelist[m].force = VectorAddition(particlelist[m].force, VectorFlip(dF));
          }
      }
  }
}

void sum_contributions(Cell *cells, Particle *gather){
  int k,j;

  for (k = 0; k < NUMBER_OF_PARTICLES; k++){
      for (j = 0; j < 8; j++){
          // printf("neighbouringcells %d %d\n", (cells + (particlelist+k)->cellnumber)->neighbouringcells[j], (particlelist+k)->cellnumber );
          (particlelist + k)->force = VectorAddition( (particlelist + k)->force, (gather+NUMBER_OF_PARTICLES*cells[(particlelist+k)->cellnumber].neighbouringcells[j] + k)->force );  
      }
  }
}

void gnuprint(FILE *gp){

  int i;
  fprintf(gp, "plot '-'\n");

  for (i=0; i<NUMBER_OF_PARTICLES; i++){
    fprintf(gp, "%g %g\n", (particlelist + i)->position.x , (particlelist + i)->position.y );
  }

  fflush(gp);
  fprintf(gp, "e\n");
}

void displace_particles(){
  int i;
  double x,y;
  int c;

  for (i = 0; i < NUMBER_OF_PARTICLES; i++){    
    x = (particlelist + i)->position.x;
    y = (particlelist + i)->position.y;

    // update positions
    x += (particlelist + i)->velocity.x * DELTAT;
    y += (particlelist + i)->velocity.y * DELTAT;

    // check for pbc
    if ( y > GRIDSIZE ){
      printf("before %lf \n ", y);
      y = y - GRIDSIZE; 
      printf("after %lf \n", y);
    }
    if ( y < 0 ){
      y = GRIDSIZE - y;
    }
    if ( x > GRIDSIZE ){
      x -= GRIDSIZE; 
    }
    if ( x < 0 ){
      x = GRIDSIZE - x;
    }

    c = (int)x  % GRIDSIZE + (int)y * GRIDSIZE;

    if (c > 8 || c < 0){
      printf("y %lf true %d \n", y, y > (double)GRIDSIZE );
      printf("cellnumber %d xc %d yc %d x %lf y %lf\n", c, (int)x  % GRIDSIZE ,  (int)y * GRIDSIZE, x, y);
      
    }

    (particlelist + i)->position.x = x;
    (particlelist + i)->position.y = y;
    (particlelist + i)->cellnumber = c;

  }
}