#####   Potentiele verbeteringen  ######

forces optellen per particle/core -> checks in vector functie als iets t zelfde is
sorting parrallel
displacement parralel
get indices

###### VRAGEN #############
Hoe cell lists te implementeren?? two dimensional array of linked list (met rma) --> array
ghost cells / pbc?
welke MPI? openmpi?
hoe computers te verbinden? ssh? science park computers root?

#######  Run & compile #######
compile: mpicc [file] -o [output]
run: mpirun [options] [file]


#### Install: #########
sudo apt-get install libcr-dev mpich2 mpich2-doc









// create MPI object
// MPI_Aint Cell_size;
// MPI_Datatype Cell;
// {
//   int blocklen[] = {1, 1};
//   MPI_Aint addr[3];
//   MPI_Address(*Map, &addr[0]);
//   MPI_Address(Map.test, &addr[1]);
//   MPI_Address(Map.number_of_particles, &addr[2]);   
//   MPI_Aint disp[] = { addr[1] - addr[0], addr[2] - addr[0] };        
//   MPI_Datatype types[]  = {MPI_INT, MPI_INT};
//   MPI_Type_create_struct(3, blocklen, disp, types, &Cell);
//   MPI_Type_commit(&Cell);
// }
// MPI_Type_extent(Cell, &Cell_size);
// Map = malloc(Cell_size);    