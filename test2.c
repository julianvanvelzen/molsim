#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

	int particlelist[9] =  {0,0,2,2,2,4,4,4,4};
	// int particlelist[9] =  {0,1,2,3,4,5,6,7,8};
	// int particlelist[9] =  {0,0,0,0,0,0,0,0,0};
	// int particlelist[9] =  {1,1,1,1,1,1,1,1,1};
	// int particlelist[9] =  {5,5,5,5,5,5,5,5,5};



	int processors = 9;
	int numberofparticles = 9;
	int i;
	int starts[processors];
	int ends[processors];
	int current = -1;
	int count = 0;

	for (i=0; i<numberofparticles;i++){
		starts[i] = processors;
		ends[i]   = processors;
	}

	for (i = 0; i < numberofparticles; i++){
		if( particlelist[i] ==  current )
			continue;
		if( particlelist[i] > count){
			printf("true %d\n", count );
			starts[count] = 0;
			count ++;
			continue;
		}

		current = particlelist[i];
		starts[count]  = i;
		count ++;
	}
	for (i=0; i < processors; i++){
		ends[i] = starts[i+1];
		if (starts[i+1] <= 0)
			ends[i] = numberofparticles;
		if (starts[i+1] == 0)
			ends[i] = 0;
		if (i == processors -1)
			ends[i] = processors;
	}

	printf("\n");

	for (i=0; i< processors; i++)
		printf("%d %d\n", starts[i], ends[i]);


}

