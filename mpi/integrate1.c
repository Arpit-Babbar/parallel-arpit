// \int_a^b cos(x) dx by splitting domain over multiple processors

#include <stdio.h>
#include <mpi.h>
#include <math.h>

# define M_PI		3.14159265358979323846	/* pi */

// function integrating cos from a to b
// TODO - Replace with quadrature
double integrate(double a, double b)
{
   return sin(b) - sin(a);
}

int main(int argc, char** argv) {
   int rank, size, ierror;
   MPI_Status status; // CLUELESS. Source at bottom
   
   ierror = MPI_Init(&argc, &argv);
   ierror = MPI_Comm_size(MPI_COMM_WORLD, &size);
   ierror = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // integration limits
   double a = 0.0, b = 0.5 * M_PI;

   // Each machine will integrate over a particular part, as specified by the
   // rank of the machine. The limits for current machine(called "me")
   double mya = a + rank*(b-a)/size;
   double myb = mya + (b-a)/size;
   // integrate f(x) over my own chunk - actual work
   double psum = integrate(mya, myb);

   // rank 0 collects partial results
   if (rank==0)
   {
      printf("Number of procs = %d\n", size);
      double res = psum;// Total result stored here
      double tmp;       // Temporary variable to recieve results from other 
                        // processors
      for (int i = 1; i < size; i++)
      {
         MPI_Recv(&tmp,           // variable to recieve buffer in
                  1,              // maximum number of elements in receive buffer
                  MPI_DOUBLE,     // MPI double(recieve type)
                  i,              // rank of source
                  0,              // message tag(integer). Unused here
                  MPI_COMM_WORLD, // MPI Communicator. Visit link at bottom
                  &status         // CLUELESS. Source at bottom
                  );
         res = res + tmp;
      }
       printf("Result = %f\n",res);
   }
   else
   {
      MPI_Send(&psum,         // variable of sent buffer
               1,             // maximum number of elements in send buffer
               MPI_DOUBLE,    // data type of send buffer
               0,             // rank of destination processor
               0,             // message tag(integer). Unused here
               MPI_COMM_WORLD // MPI Communicator. Visit link at bottom
               );
   }
   ierror = MPI_Finalize();
   return 0;
}
/*
Variables
rank   = processor number
size   = number of processors
i      = index of current processor
ierror = integer indicating success or failure returned by
         MPI_Init, MPI_Comm_size, MPI_Comm_rank
*/
/*Functions
MPI_Init      = Initializes MPI, taking number of processors from
                command line
MPI_Comm_size = Puts communication size in the variable `size`
MPI_Comm_rank = Puts rank of current processor in variable `rank`
MPI_Finalize  = Terminates MPI execution environment, all processes must call
                call this routing before exiting.                  
*/

/* Structure
All processors compute a part of the integral and then send their result to 
rank 0 processor, which adds them all to give the final value
*/

/*
Misc
Info about communicators
https://www.codingame.com/playgrounds/47058/have-fun-with-mpi-in-c/mpi-communicators
MPI_Status source - 
https://www.geeksforgeeks.org/sum-of-an-array-using-mpi/
MPI Status info
https://docs.microsoft.com/en-us/message-passing-interface/mpi-status-structure
*/
