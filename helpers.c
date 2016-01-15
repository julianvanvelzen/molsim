#include "system.h"

void ForceEnergy(Vector v1, Vector v2,Vector *dF, double *dE){
	double distance = VectorDistance(v1, v2);
	// printf("%lf\n", distance );
	Vector v3;
	if (distance > RCUT) return;
	v3.x = (v1.x - v2.x)/distance;
	v3.y = (v1.y - v2.y)/distance;
	*dE  = REPULSION_CST*SQR(distance-RCUT)/SQR(RCUT);
	dF->x = 2*REPULSION_CST*(distance-RCUT)/SQR(RCUT)*v3.x;
	dF->y = 2*REPULSION_CST*(distance-RCUT)/SQR(RCUT)*v3.y;
}

Vector VectorAddition(Vector v1, Vector v2){
  Vector v3;
  v3.x = v1.x+v2.x;
  v3.y = v1.y+v2.y;
  return v3;
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

int cmpfunc (const void * a, const void * b){
  Particle *A = (Particle *)a;
  Particle *B = (Particle *)b;
  return ( A->cellnumber - B->cellnumber );
} 

void getindeces(Particle *particlelist, int *indices){
  int i;
  int processor=0;
  *(indices) = 0;
  for (i=0;i<NUMBER_OF_PARTICLES;i++){
    if ( (particlelist+i)->cellnumber > processor ){
      *(indices+(processor)+1) = i-1;
      processor++;
      *(indices+processor+1) = i;
      processor++;
    }
  }
  for(i=processor;i<2*NUMBER_OF_PROCESSORS;i++){
    *(indices+i)=NUMBER_OF_PROCESSORS;
  }
}