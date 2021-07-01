#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
   int rank, size, i, ierror;
   
   ierror = MPI_Init(&argc, &argv);
   ierror = MPI_Comm_size(MPI_COMM_WORLD, &size);
   ierror = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   printf("Hello World, I am rank %d of %d procs\n", rank, size);
   ierror = MPI_Finalize();
   return 0;
}
// ## Variables
// # rank   = machine number
// # size   = number of machines
// # i      = index of current machine
// # ierror = integer indicating success or failure returned by
// #          MPI_Init, MPI_Comm_size, MPI_Comm_rank

// ## Functions
// # MPI_Init      = Initializes MPI, taking number of machines from
// #                 command line
// # MPI_Comm_size = Puts communication size in the variable `size`
// # MPI_Comm_rank = Puts rank of current machine in variable `rank`
// # MPI_Finalize  = Terminates MPI execution environment, all processes must call
// #                 call this routing before exiting.                  

// ## Explanation
/* # Everytime we run an MPI code, it runs the same version of a.out on all
the different machines. Each version knows its rank, so it can print it
on screen. That is all this code is doing. In the background, upon,
upon initialization by MPI_Init(), MPI sets up the so-called 
'world communicator', which is called MPI_COMM_WORLD. A communicator defines
a group of mpi processes that can be referred to by a communicator handle.
The MPI_COMM_WORLD handle describes all processes that have been started as 
part of the parallel program. It is already required as an argument in 
MPI_COMM_size, MPI_COMM_rank and will continue to be required in most functions
*/ 