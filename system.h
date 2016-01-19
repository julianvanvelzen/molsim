#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ran_uniform.h"
#include <unistd.h>
#include <signal.h>

// constants
#define MAX_NUMBER_OF_PARTICLES 100
#define MAX_ROWS 3
#define MAX_COLLUMNS 3
#define SQR(x) ((x)*(x))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DELTAT 0.001

// Structs
typedef struct {
	double x;
	double y;
} Vector;

typedef struct {
	Vector position;	
	Vector velocity;

	// force[0] = force at previous timestep. force[1] = force at current timestep
	Vector force[2];

	int cellnumber;
	double potential;
	double radial_distribution;
} Particle;

typedef struct {
  int start;
  int end;
  int totalcount;
  int neighbouringcells[8];
  int world_rank;
} Cell;

// prototypes
void CheckInputErrors();
void Initialize(void);
Vector VectorAddition(Vector v1, Vector v2);
Vector VectorMultiplication(Vector vector, double factor);
Vector RanUnit(void);
void ForceEnergy(Vector v1, Vector v2,Vector *forceVector, double *dE);
double VectorDistance(Vector v1, Vector v2);
void getNearbyCoordinates(Cell *cell, int currentPosition);
void Mdloop(int world_rank);
int cmpfunc (const void * a, const void * b);
void setindeces(Particle *particlelist, Cell *indices);
void AssignCellnumber(int Particlenumber);
void loopforces(Cell *cells, int world_rank);
void sum_contributions(Cell *cells, Particle *gather);
void gnuprint(FILE *gp);
void ApplyNewForces();
void displace_particles();
void clean_exit_on_sig(int sig_num);

// globals
extern int TEMPERATURE;
extern double RCUT;
extern double REPULSIVE_CST;
extern int NUMBER_OF_CYCLES;
extern int NUMBER_OF_PROCESSORS;
extern int NUMBER_OF_PARTICLES;
extern int GRIDSIZE;
extern Particle *particlelist;

int NUMBER_OF_PARTICLES;
int NUMBER_OF_CYCLES;
int NUMBER_OF_PROCESSORS;
int GRIDSIZE;
double RCUT;
double REPULSIVE_CST;
int TEMPERATURE;
Particle *particlelist;
