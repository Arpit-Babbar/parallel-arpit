// Tutorial 2

// In the previous code, we chose a specific way to split processes. However,
// there's no reason for that choice to be optimum. To get the optimum choice,
// we use the MPI_Dims_create() function.

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
  int dims[ndims]; set2zero(dims,ndims); // Initialize array and set it to zero
  int periodicity[] = {0,1};             // open, periodic

  MPI_Comm comm_cart;

  ierr = MPI_Dims_create(size,  // number of nodes in grid, i.e., processes
                         ndims, // Input actual dimension of grid(2,3,etc.)
                         dims   // Output array(Ex:- (2,2) reflecting 2 X 2)
                         // dims is also input, see 'NOTE' at bottom.
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

  int coords[2];
  MPI_Cart_coords(comm_cart,
                  rank,
                  2,
                  coords
                  );
  printf("Co-ordinates of rank %d are %d x %d\n", rank, coords[0], coords[1]);
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

*/