#include <stdio.h>
#include <stdlib.h>


typedef struct {
	float x;
	float *y;
} Vector;

typedef struct {
	Vector position;
	Vector *positionptr;
} World;


int main(int argc, char** argv) {
	float number = 30;
	Vector *posptr;
	Vector pos;
	World wrld;
	World *wrldptr;
	World *lptr;

  printf("%d %d %d %d %d\n", sizeof(int), sizeof(double), sizeof(float), sizeof(char), sizeof(_Bool) );

	// arrays
	lptr = malloc(10*sizeof(World));
	(lptr+2)->positionptr = malloc(sizeof(Vector));
	(lptr+2)->positionptr->y = malloc(sizeof(float));
	(lptr+2)->positionptr->x = 2;
	*((lptr+2)->positionptr->y) = 3;
	printf("\npointers to arrays met pointers\n");
	printf("struct pointer %lf\n", (*(lptr+2)).positionptr->x );
	printf("struct pointer %lf\n", *((lptr+2)->positionptr->y) );


	// allocate memory when assigning a pointer
	wrldptr = malloc(sizeof(World));
	posptr  = malloc(sizeof(Vector));

	// Nested structs (pointers)
	printf("\nnested structs\n");
	wrldptr->position.x = 1;
	printf("wrldptr x: %lf\n", wrldptr->position.x );
	wrldptr->position.y = malloc(sizeof(float));
	*(wrldptr->position.y) = 2;
	printf("wrldptr y: %lf\n", *(wrldptr->position.y) );
	wrldptr->positionptr = malloc(sizeof(Vector));
	wrldptr->positionptr->y = malloc(sizeof(float));
	*(wrldptr->positionptr->y) = 3;
	printf("wrldptr y: %lf\n", *(wrldptr->positionptr->y) );


	// Pointer structs and Points in structs
	printf("\nPointer in structs\n");
	posptr->x = 10;
	printf("ptr x: %lf\n", posptr->x );
	posptr->y = malloc(sizeof(float));
	*(posptr->y) = 20;
	printf("ptr y: %lf\n", *(posptr->y) );
	posptr->y = &number;
	printf("ptr y: %lf\n", *(posptr->y) );
	pos.x = 100;
	printf("pos x: %lf\n", pos.x );
	pos.y = malloc(sizeof(float));
	*(pos.y) = 200;
	printf("pos x: %lf\n", *(pos.y) );



}
