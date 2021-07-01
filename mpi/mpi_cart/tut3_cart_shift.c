// Tutorial 2

// In this tutorial, we see the neighbours of a particular point in grid
// along a specified axis towards specifed positive/negative side with
// specified distance.

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

void set2zero(int arr[], int size);

int main(int argc, char** argv)
{
  int rank, size, ierr;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int ndims = 2;
  int dims[ndims]; set2zero(dims,ndims);
  int periodicity[] = {0,0};

  MPI_Comm comm_cart;

  ierr = MPI_Dims_create(size,  // number of nodes in grid
                         ndims, // Actual dimension of grid
                         dims   // Array to store processors in each dimension
                         );

  if (rank == 0)
    printf("dims = (%d X %d)\n", dims[0], dims[1]);

  ierr = MPI_Cart_create(MPI_COMM_WORLD, // input communicator
                         2,              // number of dimensions
                         dims,            // k x l grid
                         periodicity,    // open in x, periodic in y
                         0,              // Don't reorder processes
                         &comm_cart      // output communicator
                        );
  int direction = 0;
  int disp = -1;
  int rank_source, rank_dest;

  ierr = MPI_Cart_shift(comm_cart,
                        direction,
                        disp,
                        &rank_source, &rank_dest);

  printf("My rank is %d.\n", rank);
  if (rank_source != MPI_PROC_NULL)
    printf("Along direction %d and displacement %d, my neighbour's rank is %d.\n",
            direction, disp, rank_source);
  else
    printf("Along direction %d and displacement %d, I don't have a neighbour.\n",
            direction, disp);
  if (rank_dest != MPI_PROC_NULL)
    printf("Going opposite, my neighbour has rank %d.\n",
            rank_dest);
  else
    printf("Going opposite, I don't have a neighbour.\n");

  int coords[2];
  MPI_Cart_coords(comm_cart,
                  rank,
                  2,
                  coords
                  );
  printf("Co-ordinates of my rank are %d x %d\n", coords[0], coords[1]);


  ierr = MPI_Finalize();
  return 0;
}

void set2zero(int arr[], int size)
{
  for(int i=0;i<size;i++)
    arr[i]=0;
}

/*

int MPI_Dims_create(int nnodes, // Number of nodes in grid
                    int ndims,  // Actual dimension of grid
                    int dims[]  // Array with #(processes) in each dimension
                    )

Interesting how this has nothing to do with MPI. It'll just always work. Should
we only be doing this with one processor btw?

NOTE - If we leave an entry as non-zero, MPI_Dims_create() would take it to
mean that we want the number of processors in that dimension to be that
non-zero entry. With the non-zero entries fixed, the function fills the zero
entries in the optimum way.

// MPI_Cart_shift

Returns the shifted source and destination ranks, given a shift direction and
amount.

MPI_Cart_shift(MPI_Comm comm,
              int direction,    // The axis along which we'd move
              int disp,         // The sign and amount of displacement
                                // along the axis specified in direction
              int *rank_source, // My rank, i.e., current rank
              int *rank_dest)   // Rank of my neighbour along specified
                                // direction and displacement.

The book says incorrect wrong will be returned as MPI_PROC_NULL and no
communication will take place if those are given as arguments. But, I am
getting negative values. Isn't that different? Hm...
*/