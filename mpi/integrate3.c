// Integrate cos(x) from a to b by splitting domain over multiple processors
// In this code, we do a very simplified version using MPI_Reduce

#include <stdio.h>
#include <mpi.h>
#include <math.h>

// function integrating cos from a to b
double integrate(double a, double b)
{
   return sin(b) - sin(a);
} // TODO - Replace with quadrature!!

int main(int argc, char** argv) {
   int rank, size, ierror;
   
   ierror = MPI_Init(&argc, &argv);
   ierror = MPI_Comm_size(MPI_COMM_WORLD, &size);
   ierror = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // integration limits
   double a = 0.0, b = 0.5 * M_PI;

   // Each machine will integrate over a particular part, as specified by the
   // rank of the machine. Here, we specify limits for current machine(called "me")
   double mya = a + rank*(b-a)/size;
   double myb = mya + (b-a)/size;

   // integrate f(x) over my own chunk - actual work
   double psum = integrate(mya, myb);
   double res = 0.0;
   
   // MPI_Reduce will do the operation MPI_Op between send and receive buffer
   // and store the output in receive buffer. The receive buffer is only in 
   // root processor
   ierror = MPI_Reduce(&psum,       // Send buffer
              &res,        // Receive buffer
              1,           // Buffer size
              MPI_DOUBLE,  // MPI_Datatype
              MPI_SUM,     // MPI_Op
              0,           // Rank of root process
              MPI_COMM_WORLD
              );

   // rank 0 prints the final result
   if (rank==0)
   {
      printf("Result = %f\n", res);
   }
}
// Variables
// rank   = processor number
// size   = number of processors
// i      = index of current processor
// ierror = integer indicating success or failure returned by
//          MPI_Init, MPI_Comm_size, MPI_Comm_rank

// Functions
// MPI_Init      = Initializes MPI, taking number of processors from
//                 command line
// MPI_Comm_size = Puts communication size in the variable `size`
// MPI_Comm_rank = Puts rank of current processor in variable `rank`
// MPI_Finalize  = Terminates MPI execution environment, all processes must call
//                 call this routing before exiting.                  

// Structure
// In this code, processor 0 doesn't iteratively wait for each processor to 
// send the result. Instead, processor 0 waits for ALL the processors to
// give their results. When that is available, it adds them to get
// the final result

// Misc
// Info about communicators
// https://www.codingame.com/playgrounds/47058/have-fun-with-mpi-in-c/mpi-communicators
// MPI_Status source - 
// https://www.geeksforgeeks.org/sum-of-an-array-using-mpi/