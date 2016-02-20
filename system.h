#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ran_uniform.h"
#include <unistd.h>
#include <signal.h>

// constants
#define INITIALISATION_STEPS 10
#define SQR(x) ((x)*(x))
#define NUMBER_OF_BINS 10

// Structs
typedef struct {
  double x;
  double y;
} Vector;

typedef struct {
  Vector position;	
  Vector velocity;
  Vector force[2]; // force[0] = force at previous timestep. force[1] = force at current timestep
  int cellnumber;
  double potential;
  int radial_distribution[NUMBER_OF_BINS]; // radial_distribution[NUMBER_OF_BINS] = outside Rcut
  double pressure_contribution;
} Particle;

typedef struct {
  int start;
  int end;
  int totalcount;
  int neighbouringcells[8];
} Cell;

// prototypes
void CheckInputErrors();
void Initialize(void);
Vector VectorAddition(Vector v1, Vector v2);
Vector VectorScalar(Vector vector, double factor);
Vector RanUnit(void);
void ForceEnergy(Particle *p1, Particle *p2);
double VectorDistance(Vector v1, Vector v2);
void getNearbyCoordinates(Cell *cell, int currentPosition);
void Mdloop(int world_rank);
int cmpfunc (const void * a, const void * b);
void setindeces(Particle *particlelist, Cell *indices);
void AssignCellnumber(int Particlenumber);
void loopforces(Cell *cells, int world_rank);
void sum_apply_contributions(Cell *cells, Particle *gather, int cycle, int worl);
void gnuprint(FILE *gp);
void displace_particles();
void clean_exit_on_sig(int sig_num);
char* VECTOR_DUMP(Vector d);
char* PARTICLE_DUMP(Particle d);
char* CELL_DUMP(Cell d);
char* INT_ARRAY_DUMP(int length, int data[]  );
void HistPrint(FILE *gp, int i);

// globals
extern double TEMPERATURE;
extern double RCUT;
extern double REPULSIVE_CST;
extern double DELTAT;
extern double *kinetic_energy_array;
extern double *potential_energy_array; 
extern double averages[4]; 
extern double pressure;
extern double initialisation_sum;
extern double rdf_total[NUMBER_OF_BINS];
extern int NUMBER_OF_CYCLES;
extern int NUMBER_OF_PROCESSORS;
extern int NUMBER_OF_PARTICLES;
extern int GRIDSIZE;
extern Particle *particlelist;
