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

void getNearbyCoordinates(Cell *cell, int currentPosition){
  (cell+currentPosition)->neighbouringcells[0] = currentPosition + GRIDSIZE;
  (cell+currentPosition)->neighbouringcells[1] = currentPosition + GRIDSIZE + 1;
  (cell+currentPosition)->neighbouringcells[2] = currentPosition + 1;
  (cell+currentPosition)->neighbouringcells[3] = currentPosition - GRIDSIZE + 1;
  (cell+currentPosition)->neighbouringcells[4] = currentPosition - GRIDSIZE;
  (cell+currentPosition)->neighbouringcells[5] = currentPosition - GRIDSIZE - 1;
  (cell+currentPosition)->neighbouringcells[6] = currentPosition - 1;
  (cell+currentPosition)->neighbouringcells[7] = currentPosition + GRIDSIZE - 1;

  // links
  if(currentPosition%GRIDSIZE==0){
    (cell+currentPosition)->neighbouringcells[5] += GRIDSIZE;   
    (cell+currentPosition)->neighbouringcells[6] += GRIDSIZE; 
    (cell+currentPosition)->neighbouringcells[7] += GRIDSIZE;
  }

  // rechts 
  if ((currentPosition+1)%GRIDSIZE == 0 ){
    (cell+currentPosition)->neighbouringcells[1] = (currentPosition + 1)%(SQR(GRIDSIZE));
    (cell+currentPosition)->neighbouringcells[2] -= GRIDSIZE;
    (cell+currentPosition)->neighbouringcells[3] = currentPosition - 2*GRIDSIZE+1;
  }
  // boven
  if (currentPosition+GRIDSIZE > SQR(GRIDSIZE)-1){
    (cell+currentPosition)->neighbouringcells[0] = (currentPosition + GRIDSIZE)%GRIDSIZE;
    (cell+currentPosition)->neighbouringcells[1] = (currentPosition + GRIDSIZE+1)%GRIDSIZE;
    (cell+currentPosition)->neighbouringcells[7] = (currentPosition + GRIDSIZE-1)%GRIDSIZE;
  }
  // onder
  if (currentPosition-GRIDSIZE < 0){
    (cell+currentPosition)->neighbouringcells[3] += SQR(GRIDSIZE);
    (cell+currentPosition)->neighbouringcells[4] += SQR(GRIDSIZE);
    (cell+currentPosition)->neighbouringcells[5] += SQR(GRIDSIZE);
  }

  // printf("w:%d  %d %d %d %d %d %d %d %d\n", currentPosition, cell[currentPosition].neighbouringcells[0], cell[currentPosition].neighbouringcells[1], cell[currentPosition].neighbouringcells[2], cell[currentPosition].neighbouringcells[3], cell[currentPosition].neighbouringcells[4], cell[currentPosition].neighbouringcells[5], cell[currentPosition].neighbouringcells[6], cell[currentPosition].neighbouringcells[7]);

}

int cmpfunc (const void * a, const void * b){
  Particle *A = (Particle *)a;
  Particle *B = (Particle *)b;
  return ( A->cellnumber - B->cellnumber );
} 

void setindeces(Particle *particlelist, Cell *cell){
  int i, j;
  int currentcell = 0;


  for(i=0;i<NUMBER_OF_PROCESSORS;i++){
    (cell+i)->start = 0;
    (cell+i)->end = 0;
  }
  
  i = 0;
  
  while (i<NUMBER_OF_PARTICLES){
    while ((particlelist+i)->cellnumber > currentcell ){
      (cell+currentcell)->start = i;
      (cell+currentcell)->end = i;
      currentcell++;
    }
  
    (cell+currentcell)->start = i;
  
    while((particlelist+i)->cellnumber == currentcell){
      i++;
    }
    (cell+currentcell)->end = i;
    currentcell++;


    if(i >= NUMBER_OF_PARTICLES){
      for(j=currentcell;j<NUMBER_OF_PROCESSORS;j++){
        (cell+j)->start = NUMBER_OF_PARTICLES;
        (cell+j)->end   = NUMBER_OF_PARTICLES;
      }
    }
  }

  // for(i=0;i<NUMBER_OF_PROCESSORS;i++){
  //   printf("start %d \tend %d \n", (cell+i)->start, (cell+i)->end);
  // }

  //   // als particle in lijst nieuwe waarde heeft tov cell number waar je mee bezig bent, 
  //   // dan is huidig punt het einde van vorige cell reeks, en (i+1) begin van nieuwe
  //   while ((particlelist+i)->cellnumber > currentcell+1 ){
  //     currentcell++;
  //   }
  //   printf("%d %d %d\n", i, (particlelist+i)->cellnumber, currentcell );
  //   if ((particlelist+i)->cellnumber > currentcell ){
  //     (cell + currentcell )->start = i;
  //     (cell + currentcell )->end   = i;
  //     currentcell ++;
  //   }
  // }
  // for(i=currentcell;i<NUMBER_OF_PROCESSORS;i++){
  //   (cell+i)->end   = NUMBER_OF_PROCESSORS-1;
  //   (cell+i)->start = NUMBER_OF_PROCESSORS-1;
  // }
  // for (j =0; j<NUMBER_OF_PROCESSORS; j++){
  //   printf("%d %d\n", (cell+j)->start, (cell+j)->end  );
  // }

  // *(cell) = 0;
  // int currentcell = 1;
  // for (i=0;i<NUMBER_OF_PARTICLES;i++){
  //   if ( (particlelist+i)->cellnumber > processor ){
  //     *(cell+currentcell) = i-1;
  //     currentcell++;
  //     *(cell+currentcell) = i;
  //     currentcell++;  
  //     processor++;
  //   }
  // }
  // for(i=currentcell;i<2*NUMBER_OF_PROCESSORS;i++){
  //   *(cell+i)=NUMBER_OF_PROCESSORS;
  // }
}