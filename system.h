#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ran_uniform.h"


// constants
#define MAX_NUMBER_OF_PARTICLES 100
#define MAX_ROWS 3
#define MAX_COLLUMNS 3
#define REPULSION_CST 1
#define SQR(x) ((x)*(x))

// Structs
typedef struct {
	double x;
	double y;
} Vector;

typedef struct {
	Vector position;
	Vector velocity;
	Vector force;
	int cellnumber;
	double  radial_distribution;
} Particle;


// typedef struct {
// 	int particles[MAX_NUMBER_OF_PARTICLES];
// 	int number_of_particles;
// 	double pressure;
// } Cell;

// prototypes
void Initialize(void);
Vector VectorAddition(Vector v1, Vector v2);
Vector RanUnit(void);
void ForceEnergy(Vector v1, Vector v2,Vector *dF, double *dE);
double VectorDistance(Vector v1, Vector v2);
void Mdloop(int world_rank);
int cmpfunc (const void * a, const void * b);
void getindeces(Particle *particlelist, int *indices);


// globals
extern int TEMPERATURE;
extern int RCUT;
extern int REPULSIVE_CST;
extern int NUMBER_OF_CYCLES;
extern int NUMBER_OF_PROCESSORS;
extern int NUMBER_OF_PARTICLES;
extern int GRIDSIZE;
extern Particle *particlelist;
// extern Cell Map[MAX_COLLUMNS][MAX_ROWS];

int NUMBER_OF_PARTICLES;
int NUMBER_OF_CYCLES;
int NUMBER_OF_PROCESSORS;
int GRIDSIZE;
int RCUT;
int TEMPERATURE;
int REPULSIVE_CST;
Particle *particlelist;
// Cell Map[MAX_COLLUMNS][MAX_ROWS];