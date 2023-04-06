#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

int schedule_thr;
int numthr;
double get_wtime(void){ 
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}

void compute_max_density(double *rho_, int N)
{
	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (rho_[i*N + j] > max_rho)
			{
				max_rho = rho_[i*N + j];
				max_i = i;
				max_j = j;
			}
		}

	printf("=====================================\n");
	printf("Output of compute_max_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);
}
void compute_max_density_omp(double *rho_, int N)
{
	// TODO 2: parallelize this function with OpenMP. (The provided code is identical to that of compute_max_density()).

	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;
	#pragma omp parallel for num_threads(numthr)
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (rho_[i*N + j] > max_rho)
			{
				max_rho = rho_[i*N + j];
				max_i = i;
				max_j = j;
			}
		}

	printf("=====================================\n");
	printf("Output of compute_max_omp_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);
}

void compute_max_density_omp_wss(double *rho_, int N)
{
	// TODO 2: parallelize this function with OpenMP. (The provided code is identical to that of compute_max_density()).

	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;
	#pragma omp parallel for schedule(static,schedule_thr) num_threads(numthr)
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (rho_[i*N + j] > max_rho)
			{
				max_rho = rho_[i*N + j];
				max_i = i;
				max_j = j;
			}
		}

	printf("=====================================\n");
	printf("Output of compute_max_omp_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);
}


void compute_max_density_omp_wsd(double *rho_, int N)
{
	// TODO 2: parallelize this function with OpenMP. (The provided code is identical to that of compute_max_density()).

	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;
	#pragma omp parallel for schedule(dynamic,schedule_thr) num_threads(numthr)
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (rho_[i*N + j] > max_rho)
			{
				max_rho = rho_[i*N + j];
				max_i = i;
				max_j = j;
			}
		}
	printf("=====================================\n");
	printf("Output of compute_max_omp_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);

}
void compute_max_density_omp_wsg(double *rho_, int N)
{
	// TODO 2: parallelize this function with OpenMP. (The provided code is identical to that of compute_max_density()).

	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;
	#pragma omp parallel for schedule(guided,schedule_thr) num_threads(16)
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (rho_[i*N + j] > max_rho)
			{
				max_rho = rho_[i*N + j];
				max_i = i;
				max_j = j;
			}
		}
	printf("=====================================\n");
	printf("Output of compute_max_omp_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);
}


int main(int argc, char *argv[])
{
	int N = 4096;
	double t1,t0;

	double *rho = (double *)malloc(N*N*sizeof(double));	// NxN 'density' matrix
	
	printf("\nGive number of threads __: ");
	scanf("%d", &numthr);

	printf("\nGive number for openmp scheduling __: ");
	scanf("%d", &schedule_thr);
	
	printf("\n");
	fflush(stdout);
	
	// matrix initialization
	srand48(1);
	for (int i = 0 ; i < N; i++) 
	for (int j = 0 ; j < N; j++) 
		rho[i*N+j] = drand48();			// these will be density values when we study diffusion

	t0=get_wtime();
	compute_max_density(rho, N);			// sequential version 
	t1=get_wtime();
	printf("\nSerial search time taken %.4f\n\n", t1-t0);
	fflush(stdout);
	
	t0=get_wtime();
	compute_max_density_omp(rho, N);		// your parallel version: must produce the same results
	t1=get_wtime();
	printf("\nParallel search time taken without schedule %.4f\n\n", t1-t0);
	fflush(stdout);
	
	t0=get_wtime();
	compute_max_density_omp_wss(rho, N);		// your parallel version: must produce the same results
	t1=get_wtime();
	printf("\nParallel search time taken with static schedule %.4f\n\n", t1-t0);
	fflush(stdout);
	
	t0=get_wtime();
	compute_max_density_omp_wsd(rho, N);		// your parallel version: must produce the same results
	t1=get_wtime();
	printf("\nParallel search time taken with dynamic schedule %.4f\n\n", t1-t0);

	t0=get_wtime();
	compute_max_density_omp_wsg(rho, N);		// your parallel version: must produce the same results
	t1=get_wtime();
	printf("\nParallel search time taken with guided schedule %.4f\n\n", t1-t0);



	free(rho);

	return 0;
}

