// ut + a ux = 0 using Lax-Wendroff

#include <stdio.h>
#include <stdlib.h> // for exit(0)
#include <string.h>
#include <mpi.h>
#include <math.h>

double init_condn(double x)
{
  double value = sin(2.0 * M_PI * x);
  return value;
}
void output_vectors_to_file(double grid[], double u[], double u_exact[],
                            int myN, int rank)
{
    FILE *fptr;
    char file_name[] = "Solution_";
    char index[30];
    sprintf(index, "%d", rank);
    strcat(file_name, index);
    strcat(file_name,".txt");
    fptr = fopen(file_name, "w");
    if (fptr == NULL)
    {
        printf("Could not open file, exiting...\n");
        exit(0);
    }
    for (int j = 1; j <= myN; j++)
        fprintf(fptr, "%f %f %f\n", grid[j], u[j], u_exact[j]);
    fclose(fptr);
    // Source - https://www.programiz.com/c-programming/c-file-input-output
    // Better source - https://www.geeksforgeeks.org/fprintf-in-c/
}

// send u[source_index] to some variable in processor dest, non-blocking
void send_value(int dest, int source_index,
                MPI_Request *request,
                double u[]
                /*double &source_buff, double &rec_buff*/)
{
  int ierror;
  ierror = MPI_Isend(&u[source_index],
                     1,
                     MPI_DOUBLE,
                     dest,
                     0,
                     MPI_COMM_WORLD,
                     request
                    );
  return;
}

// recieve some variable in u[dest_index], non-blocking
void recv_value(int source, int dest_index,
                   MPI_Request *request,
                   double u[])
{
  int ierror;
  ierror = MPI_Irecv(&u[dest_index],
                     1,
                     MPI_DOUBLE,
                     source,
                     0,
                     MPI_COMM_WORLD,
                     request);
  return;
} 

void update_ghost(double u0[], int myN, // Maybe find size yourself?
                                         // That'll fix the uneven size of last array
                  int rank, int size,
                  MPI_Request recv_req[2], // left, right
                  MPI_Request send_req[2]  // left, right
                  )// Do we really need 4 request arrays? Must check
{
  if (rank==0)
  {
    recv_value(size-1, 0    , &recv_req[0], u0); // fill my u[ 0] by size-1
    recv_value(1     , myN+1, &recv_req[1], u0); // fill my u[-1] by rank+1

    send_value(1     , myN, &send_req[1], u0); // send my u[-2] to rank+1
    send_value(size-1, 1  , &send_req[0], u0); // send my u[ 1] to size-1
  }
  else if (rank < size-1)
  {
    recv_value(rank-1, 0    , &recv_req[0], u0); // fill my u[ 0] by rank-1
    recv_value(rank+1, myN+1, &recv_req[1], u0); // fill my u[-1] by rank+1

    send_value(rank+1, myN, &send_req[1], u0); // send my u[-2] to rank+1
    send_value(rank-1, 1  , &send_req[0], u0); // send my u[ 1] to rank-1
  }
  else
  {
    recv_value(rank-1, 0    , &recv_req[0], u0); // fill my u[ 0] by rank-1
    recv_value(0     , myN+1, &recv_req[1], u0); // fill my u[-1] by rank=0

    send_value(0     , myN, &send_req[1], u0); // send my u[-2] to rank=0
    send_value(rank-1, 1  , &send_req[0], u0); // send my u[ 1] to rank-1
  }
}

void update_soln(double u0[], double u[], double myN, double cfl,
                 double dt, 
                 int rank, int size)
{
  MPI_Request send_req[2], recv_req[2];
  MPI_Status statuses[4];
  int ierror;
  for (int i = 1; i<= myN; i++)
    u0[i] = u[i];
  update_ghost(u0, myN, rank, size, recv_req, send_req);
  ierror = MPI_Waitall(2, recv_req, statuses);
  for (int i = 1; i <= myN; i++)
    u[i] = u0[i] - 0.5 * cfl * (u0[i+1]  - u0[i-1])
            + 0.5 * cfl * cfl * (u0[i-1] - 2.0 * u0[i] + u0[i+1]);
}

int main(int argc, char** argv)
{
  int ierror;
  int rank, size;
  ierror = MPI_Init(&argc, &argv);
  ierror = MPI_Comm_size(MPI_COMM_WORLD, &size);
  ierror = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // PDE parameters
  const double xmin = 0.0, xmax = 1.0;
  const double coeff = 1.0;
  // Scheme parameters
  const int N = 50;                                // number of points
  const double h = (xmax-xmin) / N;                          // grid size                 
  double cfl = 0.95;
  double dt = cfl * h / abs(coeff);
  double Tf = 1.0; 

  // Points per processor
  int myN = (int) floor(N/size);                     
  // Account for 1 potential uncounted point
  if (rank==size-1)
    myN = floor(N/size) + (N-floor(N/size)*size); 
  double myh = (xmax-xmin)/myN;
  double u[myN+2], u0[myN+2], u_exact[myN+2];              // Soln with 2 ghosts/proc
  printf("My rank is %d, size is %d\n", rank, size);
  double mya = xmin + rank * (xmax-xmin)/size;      
  printf("mya = %f\n", mya);

  double grid[myN+1];                                    // CHECK

  for (int i = 1; i <= myN+1; i++)
  {
    grid[i] = mya + (i-1) * h;
    u[i]  = init_condn(grid[i]); // updating non-ghost values
    u_exact[i] = init_condn(grid[i] - coeff * Tf);
  }

  int it = 0;
  double t = 0;
  while (t<Tf)
  {
    update_soln(u0, u, myN, cfl, dt, rank, size);
    t  += dt;
    it += 1;
  }
  output_vectors_to_file(grid, u, u_exact, myN, rank);
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

// Periodic boundary. Right end point not plotted

// Maybe write in such a way that all the parallel computing stuff is put
// outside as functions

// TODO -  Replace MPI_Irecv with MPI_Reduce?
// TODO -  Read parameters from a header file.