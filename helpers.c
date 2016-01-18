#include "system.h"

void ForceEnergy(Vector v1, Vector v2,Vector *dF, double *dE){
	double distance = VectorDistance(v1, v2);
	Vector v3;
	if (distance > RCUT){ 
    return;
  }
	v3.x = (v1.x - v2.x)/distance;
	v3.y = (v1.y - v2.y)/distance;
	*dE  = REPULSION_CST*SQR(distance-RCUT)/SQR(RCUT);
	dF->x = (2.0*REPULSION_CST*(distance-RCUT)/SQR(RCUT))*v3.x;
	dF->y = (2.0*REPULSION_CST*(distance-RCUT)/SQR(RCUT))*v3.y;
}

Vector VectorAddition(Vector v1, Vector v2){
  Vector v3;
  v3.x = v1.x + v2.x;
  v3.y = v1.y + v2.y;
  return v3;
}

Vector VectorFlip(Vector vector){
  vector.x = -vector.x;
  vector.y = -vector.y;
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
    
    if(i >= NUMBER_OF_PARTICLES){

      (cells+currentcell)->end = NUMBER_OF_PARTICLES+1;

      for(j=currentcell+1;j<NUMBER_OF_PROCESSORS;j++){
        (cells+j)->start = NUMBER_OF_PARTICLES+1;
        (cells+j)->end   = NUMBER_OF_PARTICLES+1;
      }

    }

    currentcell++;
  }
}

void loopforces(Cell *cells, int world_rank){
  int i,k,l,m;
  Vector dF;
  double dE;

  // loop over all particles in your own cell
  for(i = (cells + world_rank)->start; i < (cells + world_rank)->end-1; i++){

      // loop over all other particles in your own cell
      k = (cells + world_rank)->start;
      do {
          k++;
          if(k <= i)
              continue;

          dF.x = 0;
          dF.y = 0;
          ForceEnergy((particlelist + i)->position, (particlelist + k)->position,&dF,&dE);
          (particlelist + i)->force = VectorAddition((particlelist + i)->force, VectorFlip(dF));     
          (particlelist + k)->force = VectorAddition((particlelist + k)->force, dF); 

          printf("i %d, k %d if %lf kf %lf \n", i, k, (particlelist + i)->force.x, (particlelist + k)->force.x);

      } while (k < (cells + world_rank)->end-2);

      // loop over all particles in neighbouring cells
      // for (l=0; l<=4; l++){
      //     for (m = (cells + (cells+world_rank)->neighbouringcells[l])->start; m < (cells + (cells + world_rank)->neighbouringcells[l])->end-1; m++){
      //         dF.x = 0;
      //         dF.y = 0;
      //         ForceEnergy((particlelist + i)->position, (particlelist + m)->position,&dF,&dE);
      //         (particlelist + i)->force = VectorAddition((particlelist + i)->force, VectorFlip(dF));
      //         (particlelist + m)->force = VectorAddition((particlelist + i)->force, dF);
      //     }
      // } 
  }
}

void sum_contributions(Cell *cells, Particle *gather){
  int k,j,offset, box;
  Particle sum;

  for (k = 0; k < NUMBER_OF_PARTICLES; k++){
      sum.force.x = 0;
      sum.force.y = 0;

      for (j = 0; j < 8 ; j++){
        offset = (cells+((particlelist+k)->cellnumber))->neighbouringcells[j];
        sum.force.x += (gather + (NUMBER_OF_PARTICLES*offset) + k)->force.x;
        
        // sum.force.y += (gather+(NUMBER_OF_PARTICLES*(cells+((particlelist+k)->cellnumber))->neighbouringcells[j]) + k)->force.y;
      }
      box = (particlelist+k)->cellnumber;
      (particlelist+k)->force.x = (gather + (NUMBER_OF_PARTICLES*box) + k)->force.x + sum.force.x;
      (particlelist+k)->force.y = (gather + (NUMBER_OF_PARTICLES*box) + k)->force.y + sum.force.y;
  }
}

void gnuprint(FILE *gp){

  int i;

  // fack c
  char options[100] = "unset autoscale\nset xrange [0:";
  char a[] = "]\nset yrange [0:";
  char b[] = "]\nplot '-'\n";
  char c[3];
  sprintf(c, "%d", GRIDSIZE);
  strcat(options, c);
  strcat(options, a);
  strcat(options, c);
  strcat(options, b);
   

  fprintf(gp, options);

  for (i=0; i<NUMBER_OF_PARTICLES; i++){
    fprintf(gp, "%g %g\n", (particlelist + i)->position.x , (particlelist + i)->position.y );
  }

  fflush(gp);
  fprintf(gp, "e\n");
}

void AssignCellnumber(int Particlenumber){
  (particlelist + Particlenumber)->cellnumber = (int)(particlelist + Particlenumber)->position.x + GRIDSIZE*(int)(particlelist + Particlenumber)->position.y;
}

void displace_particles(){
  int i;
  double x,y;

  for (i = 0; i < NUMBER_OF_PARTICLES; i++){    
    x = (particlelist + i)->position.x;
    y = (particlelist + i)->position.y;


    // update positions
    x += (particlelist + i)->velocity.x * DELTAT + 0.5*(particlelist + i)->force.x * SQR(DELTAT);
    y += (particlelist + i)->velocity.y * DELTAT + 0.5*(particlelist + i)->force.y * SQR(DELTAT);

    (particlelist + i)->velocity.x += (particlelist + i)->force.x * DELTAT;
    (particlelist + i)->velocity.y += (particlelist + i)->force.y * DELTAT;

    // if(i == 1){
    //   printf("box: %d v: %lf f: %lf\n", (particlelist+1)->cellnumber, (particlelist)->velocity.x, (particlelist)->force.x );
    //   sleep(1);
    // }

    // check for pbc
    if( y > GRIDSIZE )
      y -= GRIDSIZE; 

    if( x > GRIDSIZE )
      x -= GRIDSIZE; 

    if( y < 0        )
      y += GRIDSIZE;

    if( x < 0        )
      x += GRIDSIZE;

    (particlelist + i)->position.x = x;
    (particlelist + i)->position.y = y;

    AssignCellnumber(i);

  }

  for(i = 0; i<NUMBER_OF_PARTICLES; i++){
    (particlelist + i)->force.x = 0;
    (particlelist + i)->force.y = 0;
  }

}