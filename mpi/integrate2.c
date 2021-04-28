// Integrate cos(x) from a to b by splitting domain over multiple processors
// In this code, we do a non-blocking recieve. See this link for info
// https://stackoverflow.com/questions/10017301/mpi-blocking-vs-non-blocking#:~:text=Similarly%2C%20MPI_Recv()%20returns%20when,MPI_Isend()%20and%20MPI_Irecv()%20.&text=Non%2Dblocking%20communication%20is%20used,computations%2C%20then%20do%20MPI_Wait()%20.

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

   // Set rank 0 to do non-blocking recieve
   MPI_Status  statuses[size];  // CLUELESS. Source at bottom
   MPI_Request requests[size];
   double      tmp[size];       // Temporary variable to recieve 
                                 // results from other processors
  // 0th element will not be used.
   if (rank==0)
   {
      printf("Number of procs = %d\n", size);
      for (int i = 1; i < size; i++)
      {
         // non-blocking receive
         ierror = MPI_Irecv(&tmp[i],        // variable to recieve buffer in
                           1,              // maximum number of elements in receive buffer
                           MPI_DOUBLE,     // MPI double
                           i,              // rank of source
                           0,              // message tag(integer). Unused here
                           MPI_COMM_WORLD, // MPI Communicator. Visit link at bottom
                           &requests[i]
                           // &status         // CLUELESS. Source at bottom
                           );
         // int sender = status.MPI_SOURCE; // CLUELESS. Source at bottom
      } 
   }
   // Each machine will integrate over a particular part, as specified by the
   // rank of the machine. Here, we specify limits for current machine(called "me")
   double mya = a + rank*(b-a)/size;
   double myb = mya + (b-a)/size;

   // integrate f(x) over my own chunk - actual work
   double psum = integrate(mya, myb);
   // rank 0 collects partial results
   if (rank==0)
   {
      double res = psum;   // Total result stored here
      ierror = MPI_Waitall(size-1, requests, statuses);
      for (int i=0; i< size-1; i++)
         res = res + tmp[i+1];
      printf("Result = %f\n", res);
   }
   // other ranks send their results to rank 0
   else
   {
      ierror = MPI_Send(&psum,         // variable of sent buffer
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