#include <stdio.h>
#include <stdlib.h>

/* We'll be using MPI routines, definitions, etc. */
#include <mpi.h>

/* Calculate local integral  */
double Trap(double left_endpt, double right_endpt, int trap_count,
            double base_len);

/* Function we're integrating */
double f(double x);

int main(int argc, char **argv) {

  int my_rank, comm_sz, n = 100000000, local_n;
  double a = 0.0, b = 3.0, h, local_a, local_b;
  double local_int = 0.0, total_int = 0.0;
  int source;

  /* Let the system do what it needs to start up MPI */
  MPI_Init(&argc, &argv);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  h = (b - a) / n;       /* h is the same for all processes */
  local_n = n / comm_sz; /* So is the number of trapezoids  */

  /* Length of each process' interval of
   * integration = local_n * h.  So my interval
   * starts at: */
  local_a = a + my_rank * local_n * h;
  local_b = local_a + local_n * h;
  local_int = Trap(local_a, local_b, local_n, h);

  /* Add up the integrals calculated by each process */
  // 如果用点对点阻塞
#if M == 1
  if (my_rank != 0) {
    MPI_Send(&local_int, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
  } else {
    for (int i = 1; i < comm_sz; i++) {
      MPI_Recv(&local_int, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      total_int += local_int;
    }
  }
#elif M == 2
  if (my_rank != 0) {
		MPI_Request request = MPI_REQUEST_NULL;
    MPI_Isend(&local_int, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD,
              &request);
  } else {
  	MPI_Request *requests = (MPI_Request *)malloc(sizeof(MPI_REQUEST_NULL) * (comm_sz - 1));
    double *all_int = (double *)malloc(sizeof(double) * comm_sz);
    all_int[0] = local_int;
    for (int i = 1; i < comm_sz; i++) {
      MPI_Irecv(&all_int[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD,
                &requests[i - 1]);
    }
    MPI_Waitall((comm_sz - 1), requests, MPI_STATUSES_IGNORE);
    for (int i = 0; i < comm_sz; i++) {
      total_int += all_int[i];
    }
    free(all_int);
    free(requests);
  }
#elif M == 3
  MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#elif M == 4
	double* all_int;
  if (my_rank == 0) {
    all_int = (double *)malloc(sizeof(double) * comm_sz);
	}
	MPI_Gather(&local_int, 1, MPI_DOUBLE, all_int, 1, MPI_DOUBLE, 0,
						 MPI_COMM_WORLD);
	if (my_rank == 0) {
    for (int i = 0; i < comm_sz; i++) {
      total_int += all_int[i];
    }
    free(all_int);
  }

#endif

  /* Print the result */
  if (my_rank == 0) {
    printf("With n = %d trapezoids, our estimate\n", n);
    printf("of the integral from %f to %f = %.15e\n", a, b, total_int);
  }

  /* Shut down MPI */
  MPI_Finalize();

  return 0;
} /*  main  */

// 调用的例子比如: Trap(local_a, local_b, local_n, h); 但h真的不是多余的么?
// 我感觉h=(b-a)/n 推导可得, 是h*(a+b+f(中间这些x))
double Trap(double left_endpt /* in */, double right_endpt /* in */,
            int trap_count /* in */, double base_len /* in */) {

  double estimate, x;
  int i;

  estimate = (f(left_endpt) + f(right_endpt)) / 2.0;
  for (i = 1; i <= trap_count - 1; i++) {
    x = left_endpt + i * base_len;
    estimate += f(x);
  }
  estimate = estimate * base_len;

  return estimate;

} /*  Trap  */

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) { return x * x; } /* f */
