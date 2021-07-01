// Tutorial 1

// A 2D tutorial of MPI_Cart_create() and illustration of its functionality
// using MPI_Cart_coords() and MPI_Cart_rank().

// The user specifies a non-prime number N of processes and we order those
// processes into a 2D array of size N = k*l with largest possible k.

// What's lacking is a way to visualize the periodicity and open-ness.
// Hager and Wellein has a good picture, but there must be a way to print it
// on screen.

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

int check_prime(int n);
void factorize2d(int n, int dims[2]); // split n = k*l with largest k

int main(int argc, char** argv)
{
  int rank, size, ierr;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  assert (check_prime(size)==1 && "Enter a non-prime np\n");
  int ndims  = 2; // Actual dimension of the cartesian mesh
  int dims[2];    // Number of processors in each dimension
  factorize2d(size, dims);
  int periodicity[] = {0,1}; // Open x, periodic along y

  if (rank == 0)
    printf("dims = (%d X %d)\n", dims[0], dims[1]);

  MPI_Comm comm_cart;

  ierr = MPI_Cart_create(MPI_COMM_WORLD, // input communicator
                         2,              // number of dimensions
                         dims,           // dims[0] x dims[1] grid
                         periodicity,    // open in x, periodic in y
                         0,              // Don't reorder processes/ranks
                         &comm_cart      // output communicator
                         );

  int coords[2];
  // rank -> cartesian coordinates
  MPI_Cart_coords(comm_cart,  // needed to tell coordinates structure
                  rank,       // current rank
                  2,          // two dimensions
                  coords      // output cartesian coordinates of current rank
                  );
  printf("Co-ordinates of rank %d are %d x %d\n", rank, coords[0], coords[1]);
  // cartesian coordinates -> rank
  if (rank == size-1)
    for (int i = 0; i < dims[0]; i++)
      for (int j = 0; j < dims[1]; j++)
      {
        int corr_rank;
        coords[0] = i, coords[1] = j;
        ierr = MPI_Cart_rank(comm_cart, coords, &corr_rank);
        printf("Rank corresponding to %d x %d is %d\n",i, j, corr_rank);
      }
  ierr = MPI_Finalize();
  return 0;
}

int check_prime(int n)
{
  int flag = 0; // 0 indicates prime number
  for (int i=2; i<= n/2; i++)
  {
    if (n % i == 0)
    {
      flag = 1;
      break;
    }
  }
  return flag;
}

void factorize2d(int n, int dims[2])
{
  dims[0] = (int) ceil(n/2);
  for (int i = dims[0]; i>1; i--)
  {
    if (n % i == 0)
    {
      dims[0] = i;
      dims[1] = (int) n / i;
      break;
    }
  }
  return;
}

/*

// Variable names(bad looking names!)

ndims       = Actual dimension of the cartesian mesh
dims[ndims] = Specifying the number of processors in each dimension

// Summary

In a 2D problem, it can be hard to keep track of which process is
the neighbouring process of which, info that is naturally needed.
That is why we renumber the processes in cartesian coordinates.
0,1,2,3 can be renumbered as {(0,0),(0,1),(1,0),(1,1)}
MPI_Cart_create() does it for us.

// Communicator

Communicator is the set of cpu allocated to the problem.
MPI_COMM_WORLD is the default communicator which contains all cpu allocated to our mpirun.
I guess a communicator is an object containing ranks and their ordering

MPI_Cart_create generates a new "cartesian" communicator which can be used to
refer to the topology. Its syntax from mpich.org is

int MPI_Cart_create(
                    MPI_Comm comm_old,   // Input communicator (handle. ?)
                    int ndims,           // Number of dimensions of cartesian grid
                    const int dims[],    // Array of size ndims specifying
                                         // no. of processes in each dimension
                                         // If ndims = 2, it'd be something like
                                         // (k,l) representing  k X l
                    const int periods[], // logical array of size ndims telling
                                         // where the grid is periodic(true)
                                         // or not (false) in each dimension
                    int reorder,         // ranking may be reordered (true) or
                                         // not(false) (logical)
                                         // i.e., rank of process in comm_cart
                                         // may be different from comm_old


                    MPI_Comm *comm_cart  // the communicator output with new
                                         // cartesian topology (handle. ?)
                    )

// Notice that this does not mention the actual problem size. Apparently,
'it's the user's job to take care of data distribution. All MPI can do is
keep track of data distribution'.


// This function determines process coordinates of the particular rank
int MPI_Cart_coords(
                    MPI_Comm comm, // Communicator of the cartesian topology
                    int rank,      // Rank of process within cartesian topology
                    int maxdims,   // length of vector 'coords'
                                   // (in the calling program. ?)
                    int coords[]   // Vector to store the coordinates in
                    )
NOTE - If rank reordering is allowed, a process should first obtain its rank
from MPI_Comm_rank(MPI_COMM_WORLD, &rank).


// This simple function gives rank corresponding to a particular set of coordinates

int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);

bool is not in C by default, so I guess they are using int?
*/
